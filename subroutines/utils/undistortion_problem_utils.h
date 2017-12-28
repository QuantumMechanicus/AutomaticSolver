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
    bool checkUndistortionInvertibility(const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) {
        T k1 = lambdas(0);
        T k2 = T(0);
        if (lambdas.rows() > 1) {
            k2 = lambdas(1);
        }
        return (k1 > T(-2) and ((k1 > T(2) and (T(-1) - k1 < k2 and k2 < -(k1 * k1 / T(12)))) or
                                (k1 <= T(2) and (T(-1) - k1 < k2 and k2 < (T(1) - k1) / T(3.0)))));

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
    T undistortionDenominatorDerivative(const T &r_distorted, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) {
        T denominator(0.0);
        T c(2.0);
        T r_distorted2 = r_distorted * r_distorted;
        T r_distorted_pow = r_distorted;
        for (int i = 0; i < lambdas.rows(); ++i) {
            denominator += c * lambdas[i] * r_distorted_pow;
            c = c + T(2.0);
            r_distorted_pow *= r_distorted2;
        }

        return denominator;
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

    //FindRoot --- depricated
    template<typename T>
    T findRoot(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 2, 1> &u1) {
        size_t n_lambda = hyp_lambdas.size();
        size_t deg = 2 * n_lambda;
        T r_u = u1.norm();
        Eigen::Matrix<T, Eigen::Dynamic, 1> coeff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
                deg + 1);

        for (size_t k = n_lambda; k > 0; --k) {
            if (ceres::abs(hyp_lambdas[k - 1]) < 1e-8)
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


        T root_r1d_estimation, root_r2d_estimation;

        T epsilon1, epsilon2;
        epsilon1 = std::min(ceres::abs(residuals.first), T(0.5));
        epsilon2 = std::min(ceres::abs(residuals.second), T(0.5));

        root_r1d_estimation = r1d - epsilon1;
        root_r2d_estimation = r2d - epsilon2;
        if (root_r1d_estimation < T(0))
            root_r1d_estimation += T(2) * epsilon1;
        if (root_r2d_estimation < T(0))
            root_r2d_estimation += T(2) * epsilon2;


        T r1u = line_point1.norm();
        T r2u = line_point2.norm();

        if (hyp_lambdas.rows() > 1) {

            int kkk = 0;
            while (kkk < 200) {
                kkk++;

                root_r1d_estimation = root_r1d_estimation -
                                      (r1u *
                                       undistortion_utils::undistortionDenominator(root_r1d_estimation, hyp_lambdas) -
                                       root_r1d_estimation) / (

                                              r1u *
                                              undistortion_utils::undistortionDenominatorDerivative(root_r1d_estimation,
                                                                                                    hyp_lambdas) - T(1)
                                      );

                root_r2d_estimation = root_r2d_estimation -
                                      (r2u *
                                       undistortion_utils::undistortionDenominator(root_r2d_estimation, hyp_lambdas) -
                                       root_r2d_estimation) / (

                                              r2u *
                                              undistortion_utils::undistortionDenominatorDerivative(root_r2d_estimation,
                                                                                                    hyp_lambdas) - T(1)
                                      );

            }
        } else {
            auto tmp1 = T(4.0) * hyp_lambdas(0) * r1u * r1u;
            auto tmp2 = T(4.0) * hyp_lambdas(0) * r2u * r2u;
            if (T(1) - tmp1 < T(0) or T(1) - tmp2 < T(0)) {
                root_r1d_estimation = T(std::numeric_limits<double>::max());
                root_r2d_estimation = T(std::numeric_limits<double>::max());
                residuals.first = T(std::numeric_limits<double>::max());
                residuals.second = T(std::numeric_limits<double>::max());
                is_correct = false;
                return residuals;
            } else {
                root_r1d_estimation = (T(1) - ceres::sqrt(T(1.0) - tmp1)) /
                                      (T(2.0) * hyp_lambdas(0) * r1u);
                root_r2d_estimation = (T(1) - ceres::sqrt(T(1.0) - tmp2)) /
                                      (T(2.0) * hyp_lambdas(0) * r2u);
            }

        }
        if (r1u < T(1e-9))
            root_r1d_estimation = T(0);
        if (r2u < T(1e-9))
            root_r2d_estimation = T(0);

        T dd1 = undistortion_utils::undistortionDenominator<T>(root_r1d_estimation, hyp_lambdas);
        T dd2 = undistortion_utils::undistortionDenominator<T>(root_r2d_estimation, hyp_lambdas);

        curve_point1 = dd1 * line_point1;
        curve_point2 = dd2 * line_point2;


        residuals.first = r * (u1d - curve_point1).norm();
        residuals.second = r * (u2d - curve_point2).norm();
        if (residuals.first != residuals.first or residuals.second != residuals.second) {
            residuals.first = T(std::numeric_limits<double>::max());
            residuals.second = T(std::numeric_limits<double>::max());
            is_correct = false;
            std::fstream f("wtf", std::fstream::out | std::fstream::app);
            f << root_r1d_estimation << " " << root_r2d_estimation << " " << curve_point1.transpose() << " "
              << curve_point2.transpose() << " ! " << u1d.transpose() << " ! " << " r: " << r << " dd1: " << dd1
              << " dd2: " << dd2 << " " << r1u << " " << r2u << residuals.first << " "
              << residuals.second << "err and n1" << err << " " << n1 << " " << n2 << std::endl;
        }
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

    class UndistortionProblemHelper {
    public:
        static constexpr double EPS = 1e-10;
    private:
        double expected_percent_of_inliers_;
        double confidence_interval_;
        double w_, h_, r_;
        double alpha_;
        Points u1d_, u2d_, u1_, u2_, nearest_u1d_, nearest_u2d_;
        long number_of_points_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> hyp_lambdas_;
        Eigen::Matrix<double, 3, 3> hyp_F_;
        double quantile_;
        std::vector<double> left_residuals_, right_residuals_, errors_;
        std::vector<int> inliers_ind_;

    public:
        UndistortionProblemHelper(double w, double h, double r, const Points &u1d, const Points &u2d,
                                  double percent = 0.2);

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
