//
// Created by danielbord on 12/12/17.
//

#ifndef AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H
#define AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/math/special_functions/erf.hpp>
#include <ceres/ceres.h>
#include <unsupported/Eigen/Polynomials>

namespace Eigen {
    namespace internal {

        template<class T, int N, typename NewType>
        struct cast_impl<ceres::Jet<T, N>, NewType> {
            EIGEN_DEVICE_FUNC
            static inline NewType run(ceres::Jet<T, N> const &x) {
                return static_cast<NewType>(x.a);
            }
        };

    }  // namespace internal
}


namespace undistortion_utils {

    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> Points;
    using StdVector = std::vector<double, Eigen::aligned_allocator<double>>;

    template<typename T>
    void normalizePoints(auto &points, T w, T h, T r) {
        points.row(0).array() -= w / 2.0;
        points.row(1).array() -= h / 2.0;
        points /= r;
    }

    template<typename T>
    T undistortionDenominator(const T &r_distorted, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) {
        T denominator(1.0);
        T r_distorted2 = r_distorted * r_distorted;
        T r_distorted2_pow = r_distorted2;
        for (int i = 0; i < lambdas.rows(); ++i) {
            denominator += lambdas[i] * r_distorted2_pow;
            r_distorted2_pow *= r_distorted2;
        }

        return denominator;
    }

    template<typename T>
    T undistortionDenominator(const T &r_distorted, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas1,
                              const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas2) {
        auto denominator = undistortionDenominator<T>(r_distorted, lambdas1);
        auto nominator = undistortionDenominator<T>(r_distorted, lambdas2);

        return denominator / nominator;
    }

    template<typename T>
    std::pair<T, T>
    computeEpipolarLineDistanceError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas,
                                     const Eigen::Matrix<T, 3, 3> &hyp_F,
                                     const Eigen::Matrix<T, 2, 1> &u1d,
                                     const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h,
                                     Eigen::Matrix<T, 2, 1> &u1,
                                     Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {
        T wT(w), hT(h);
        T d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T alpha = undistortionDenominator(dr, hyp_lambdas);


        T r1d, r2d;
        r1d = u1d.norm();
        r2d = u2d.norm();

        Eigen::Matrix<T, 3, 1> homogeneous_u1, homogeneous_u2;
        T denominator1 = undistortionDenominator(r1d, hyp_lambdas);
        T denominator2 = undistortionDenominator(r2d, hyp_lambdas);

        is_correct = (denominator1 > 0.0 && denominator2 > 0.0);
        u1 = u1d / denominator1;
        u2 = u2d / denominator2;
        homogeneous_u1 = u1.homogeneous();
        homogeneous_u2 = u2.homogeneous();
        Eigen::Matrix<T, 3, 1> l1 = hyp_F * homogeneous_u1;
        Eigen::Matrix<T, 3, 1> l2 = hyp_F.transpose() * homogeneous_u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(homogeneous_u2);

        std::pair<T, T> residuals;
        residuals.first = alpha * r * err / n1;
        residuals.second = alpha * r * err / n2;
        return residuals;
    }

    template<typename T>
    T findRoot(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 2, 1> &u1) {
        size_t n_lambda = hyp_lambdas.size();
        size_t deg = 2 * n_lambda;
        T r_u = u1.norm();
        Eigen::Matrix<T, Eigen::Dynamic, 1> coeff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
                deg + 1);

        for (size_t k = n_lambda; k > 0; --k) {
            if (std::abs(hyp_lambdas[k - 1]) < 1e-8)
                deg -= 2;
            coeff(2 * k, 0) = r_u * hyp_lambdas[k - 1];
        }
        coeff(1, 0) = T(-1);
        coeff(0, 0) = r_u;
        T r_d = std::numeric_limits<T>::max();


        coeff /= coeff[deg];
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> companion(deg, deg);
        companion.setZero();
        companion.col(deg - 1) = -1 * coeff.topRows(deg).template cast<double>();
        companion.block(1, 0, deg - 1, deg - 1).setIdentity();


        auto r = companion.eigenvalues();

        for (int kk = 0; kk < r.rows(); ++kk) {
            T real = T(r[kk].real());
            T imag = T(r[kk].imag());
            if (ceres::abs(imag) < 1e-9 && real > T(0) && r_d > real) {
                r_d = real;
            }

        }
        return r_d;
    }

    template<typename T>
    std::pair<T, T>
    computeEpipolarCurveDistanceError(const T &hyp_lambda,
                                      const Eigen::Matrix<T, 3, 3> &hyp_F,
                                      const Eigen::Matrix<T, 2, 1> &u1d,
                                      const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h,
                                      Eigen::Matrix<T, 2, 1> &u1,
                                      Eigen::Matrix<T, 2, 1> &u2, Eigen::Matrix<T, 2, 1> &curve_point1,
                                      Eigen::Matrix<T, 2, 1> &curve_point2, bool &is_correct) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> hyp_lambdas(1);

        hyp_lambdas(0) = hyp_lambda;
        T wT(w), hT(h);
        T d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T alpha = undistortionDenominator(dr, hyp_lambdas);

        T r1d, r2d;
        r1d = u1d.norm();
        r2d = u2d.norm();

        Eigen::Matrix<T, 3, 1> homogeneous_u1, homogeneous_u2;
        Eigen::Matrix<T, 2, 1> line_point1, line_point2;
        T denominator1 = undistortionDenominator(r1d, hyp_lambdas);
        T denominator2 = undistortionDenominator(r2d, hyp_lambdas);

        is_correct = (denominator1 > 0.0 && denominator2 > 0.0);
        u1 = u1d / denominator1;
        u2 = u2d / denominator2;
        homogeneous_u1 = u1.homogeneous();
        homogeneous_u2 = u2.homogeneous();
        Eigen::Matrix<T, 3, 1> l1 = hyp_F * homogeneous_u1;
        Eigen::Matrix<T, 3, 1> l2 = hyp_F.transpose() * homogeneous_u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(homogeneous_u2);


        std::pair<T, T> residuals;
        residuals.first = err / n1;
        residuals.second = err / n2;
        line_point1 = u1 + residuals.first * l1.template block<2, 1>(0, 0) / n1;
        line_point2 = u2 + residuals.second * l1.template block<2, 1>(0, 0) / n1;

        T r1u = line_point1.norm();
        T r2u = line_point2.norm();
        T root_r1d = (T(1) - ceres::sqrt(T(1.0) - T(4.0) * hyp_lambda * r1u * r1u)) / (T(2.0) * hyp_lambda * r1u);
        T root_r2d = (T(1) - ceres::sqrt(T(1.0) - T(4.0) * hyp_lambda * r2u * r2u)) / (T(2.0) * hyp_lambda * r2u);

        T dd1 = undistortion_utils::undistortionDenominator<T>(root_r1d, hyp_lambdas);
        T dd2 = undistortion_utils::undistortionDenominator<T>(root_r2d, hyp_lambdas);

        curve_point1 = dd1 * line_point1;
        curve_point2 = dd2 * line_point2;


        residuals.first = r * (u1d - curve_point1).norm();
        residuals.second = r * (u2d - curve_point2).norm();
        return residuals;
    }

    template<typename T>
    std::pair<T, T>
    computeEpipolarCurveDistanceError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas,
                                      const Eigen::Matrix<T, 3, 3> &hyp_F,
                                      const Eigen::Matrix<T, 2, 1> &u1d,
                                      const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h,
                                      Eigen::Matrix<T, 2, 1> &u1,
                                      Eigen::Matrix<T, 2, 1> &u2, Eigen::Matrix<T, 2, 1> &curve_point1,
                                      Eigen::Matrix<T, 2, 1> &curve_point2, bool &is_correct) {

        T wT(w), hT(h);
        T d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T alpha = undistortionDenominator(dr, hyp_lambdas);

        T r1d, r2d;
        r1d = u1d.norm();
        r2d = u2d.norm();

        Eigen::Matrix<T, 3, 1> homogeneous_u1, homogeneous_u2;
        Eigen::Matrix<T, 2, 1> line_point1, line_point2;
        T denominator1 = undistortionDenominator(r1d, hyp_lambdas);
        T denominator2 = undistortionDenominator(r2d, hyp_lambdas);

        is_correct = (denominator1 > 0.0 && denominator2 > 0.0);
        u1 = u1d / denominator1;
        u2 = u2d / denominator2;
        homogeneous_u1 = u1.homogeneous();
        homogeneous_u2 = u2.homogeneous();
        Eigen::Matrix<T, 3, 1> l1 = hyp_F * homogeneous_u1;
        Eigen::Matrix<T, 3, 1> l2 = hyp_F.transpose() * homogeneous_u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(homogeneous_u2);


        std::pair<T, T> residuals;
        residuals.first = err / n1;
        residuals.second = err / n2;
        line_point1 = u1 + residuals.first * l1.template block<2, 1>(0, 0) / n1;
        line_point2 = u2 + residuals.second * l1.template block<2, 1>(0, 0) / n1;


        T root_r1d = (findRoot<T>(hyp_lambdas, line_point1));
        T root_r2d = (findRoot<T>(hyp_lambdas, line_point2));

        T dd1 = undistortion_utils::undistortionDenominator<T>(root_r1d, hyp_lambdas);
        T dd2 = undistortion_utils::undistortionDenominator<T>(root_r2d, hyp_lambdas);

        curve_point1 = dd1 * line_point1;
        curve_point2 = dd2 * line_point2;


        residuals.first = r * (u1d - curve_point1).norm();
        residuals.second = r * (u2d - curve_point2).norm();
        return residuals;
    }

    template<typename T>
    std::pair<T, T>
    computeEpipolarCurveDistanceError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas,
                                      const Eigen::Matrix<T, 3, 3> &hyp_F,
                                      const Eigen::Matrix<T, 2, 1> &u1d,
                                      const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h,
                                      Eigen::Matrix<T, 2, 1> &u1,
                                      Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {

        Eigen::Matrix<T, 2, 1>
                curve_point1;
        Eigen::Matrix<T, 2, 1>
                curve_point2;
        return computeEpipolarCurveDistanceError<T>(hyp_lambdas,
                                                    hyp_F,
                                                    u1d,
                                                    u2d, r,
                                                    w,
                                                    h,
                                                    u1,
                                                    u2, curve_point1,
                                                    curve_point2, is_correct);
    }

    template<typename T>
    std::pair<T, T>
    computeEpipolarCurveDistanceError(const T &hyp_lambda,
                                      const Eigen::Matrix<T, 3, 3> &hyp_F,
                                      const Eigen::Matrix<T, 2, 1> &u1d,
                                      const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h,
                                      Eigen::Matrix<T, 2, 1> &u1,
                                      Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {

        Eigen::Matrix<T, 2, 1>
                curve_point1;
        Eigen::Matrix<T, 2, 1>
                curve_point2;
        return computeEpipolarCurveDistanceError<T>(hyp_lambda,
                                                    hyp_F,
                                                    u1d,
                                                    u2d, r,
                                                    w,
                                                    h,
                                                    u1,
                                                    u2, curve_point1,
                                                    curve_point2, is_correct);
    }
    class UndistortionProblemHelper {
    public:
        static constexpr double EPS = 1e-10;
    private:
        double expected_percent_of_inliers_ = 0.2;
        double confidence_interval_;
        double w_, h_, r_;
        double alpha_;
        Points u1d_, u2d_, u1_, u2_, nearest_u1d_, nearest_u2d_;
        long number_of_points_;

    private:
        Eigen::Matrix<double, Eigen::Dynamic, 1> hyp_lambdas_;
        Eigen::Matrix<double, 3, 3> hyp_F_;
        double quantile_;
        std::vector<double> left_residuals_, right_residuals_, errors_;
        std::vector<int> inliers_ind_;

    public:
        UndistortionProblemHelper(double w, double h, double r, const Points &u1d, const Points &u2d,
                                  double percent = 0.1);


        void computeErrors();

        void estimateQuantile();


        size_t findInliers(const std::string &out_name);


        void setHypLambdas(const Eigen::Matrix<double, -1, 1> &hyp_lambdas_);

        void setHypF(const Eigen::Matrix<double, 3, 3> &hyp_F_);


        double getExpectedPercentOfInliers() const;

        double getConfidenceInterval() const;

        double getW() const;

        double getH() const;

        double getR() const;

        double getAlpha() const;

        const Points &getU1d() const;

        const Points &getU2d() const;

        const Points &getU1() const;

        const Points &getU2() const;

        const Points &getNearestU1d() const;

        const Points &getNearestU2d() const;

        long getNumberOfPoints() const;

        const Eigen::Matrix<double, -1, 1> &getLambdas() const;

        const Eigen::Matrix<double, 3, 3> &getF() const;

        double getQuantile() const;

        const std::vector<double> &getLeftResiduals() const;

        const std::vector<double> &getRightResiduals() const;

        const std::vector<double> &getErrors() const;

        const std::vector<int> &getInliersIndices() const;

    };
}
#endif //AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H
