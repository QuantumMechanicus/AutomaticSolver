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

namespace undistortion_utils {

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
    T findRoot(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 2, 1> &u1d) {
        size_t n_lambda = hyp_lambdas.size();
        size_t deg = 2 * n_lambda;
        T r_u = u1d.norm();
        Eigen::Matrix<T, Eigen::Dynamic, 1> coeff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
                deg + 1);

        for (size_t k = n_lambda; k > 0; --k)
            coeff(2 * k, 0) = r_u * hyp_lambdas[k - 1];

        coeff(1, 0) = -1;
        coeff(0, 0) = r_u;
        T r_d = std::numeric_limits<T>::max();


        coeff /= coeff[deg];
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> companion(deg, deg);
        companion.setZero();
        companion.col(deg - 1) = -1 * coeff.topRows(deg);
        companion.block(1, 0, deg - 1, deg - 1).setIdentity();


        auto r = companion.eigenvalues();
        for (int kk = 0; kk < r.rows(); ++kk) {
            T real = r[kk].real();
            T imag = r[kk].imag();
            if (std::abs(imag) < 1e-9 && real > 0 && r_d > real) {
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


        T root_r1d = findRoot<T>(hyp_lambdas, line_point2);
        T root_r2d = findRoot<T>(hyp_lambdas, line_point2);
        T dd1 = undistortion_utils::undistortionDenominator<T>(r1d, hyp_lambdas);
        T dd2 = undistortion_utils::undistortionDenominator<T>(r2d, hyp_lambdas);

        curve_point1 = dd1 * line_point1;
        curve_point2 = dd2 * line_point2;


        residuals.first = r * (u1d - curve_point1).norm();
        residuals.second = r * (u2d - curve_point2).norm();
        return residuals;
    }


    template<typename T>
    class UndistortionProblemHelper {
    public:
        typedef Eigen::Matrix<T, 2, Eigen::Dynamic> Points;
        using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;
        static constexpr double EPS = 1e-10;
    private:
        double expected_percent_of_inliers_ = 0.2;
        double confidence_interval_;
        T w_, h_, r_;
        T alpha_;
        Points u1d_, u2d_, u1_, u2_, nearest_u1d_, nearest_u2d_;
        long number_of_points_;

    private:
        Eigen::Matrix<T, Eigen::Dynamic, 1> hyp_lambdas_;
        Eigen::Matrix<T, 3, 3> hyp_F_;
        T quantile_;
        std::vector<T> left_residuals_, right_residuals_, errors_;
        std::vector<int> inliers_ind_;

    public:
        UndistortionProblemHelper(T w, T h, T r, const Points &u1d, const Points &u2d, double percent = 0.1) : w_(w),
                                                                                                               h_(h),
                                                                                                               r_(r),
                                                                                                               u1d_(u1d),
                                                                                                               u2d_(u2d),
                                                                                                               expected_percent_of_inliers_(
                                                                                                                       percent) {
            assert(u1d.cols() == u2d.cols() && "Points correspondences must be the same size");
            hyp_lambdas_.resize(1, Eigen::NoChange);
            hyp_lambdas_.setZero();
            hyp_F_.setZero();
            number_of_points_ = u1d.cols();
            quantile_ = T(0);
            normalizePoints<T>(u1d_, w, h, r);
            normalizePoints<T>(u2d_, w, h, r);
            confidence_interval_ = 0;


            assert(number_of_points_ >= 0 && "Number of points less zero?");
            left_residuals_.resize(number_of_points_);
            right_residuals_.resize(number_of_points_);
            errors_.resize(number_of_points_);
            u1_.resize(Eigen::NoChange, number_of_points_);
            u2_.resize(Eigen::NoChange, number_of_points_);
            nearest_u1d_.resize(Eigen::NoChange, number_of_points_);
            nearest_u2d_.resize(Eigen::NoChange, number_of_points_);
            u1_.setZero();
            u2_.setZero();
            nearest_u1d_.setZero();
            nearest_u2d_.setZero();
            inliers_ind_.resize(0);
            T d = std::max(h_, w_) / T(2.0);
            T dr = d / r_;
            T alpha_ = undistortionDenominator(dr, hyp_lambdas_);

        }


        void computeErrors() {


            for (size_t k = 0; k < number_of_points_; ++k) {
                Eigen::Matrix<T, 2, 1> u1k, cu1k;
                Eigen::Matrix<T, 2, 1> u2k, cu2k;
                u1k.setZero();
                u2k.setZero();
                bool is_correct;

                auto errors = undistortion_utils::computeEpipolarCurveDistanceError<T>(hyp_lambdas_, hyp_F_,
                                                                                       u1d_.col(k), u2d_.col(k),
                                                                                       r_, w_,
                                                                                       h_, u1k, u2k, cu1k, cu2k,
                                                                                       is_correct);
                u1_.col(k) = u1k;
                u2_.col(k) = u2k;
                nearest_u1d_.col(k) = cu1k;
                nearest_u2d_.col(k) = cu2k;
                left_residuals_[k] = errors.first;
                right_residuals_[k] = errors.second;
                errors_[k] = std::abs(errors.first) + std::abs(errors.second);

            }
        }


        void estimateQuantile() {
            computeErrors();
            std::nth_element(errors_.begin(), errors_.begin() + int(errors_.size() * expected_percent_of_inliers_),
                             errors_.end());
            quantile_ = errors_[int(errors_.size() * expected_percent_of_inliers_)];
            confidence_interval_ =
                    quantile_ * boost::math::erfc_inv((0.95 + 1.0)) /
                    boost::math::erfc_inv((expected_percent_of_inliers_ + 1.0));
        }


        size_t findInliers(const std::string &out_name) {

            if (quantile_ == 0) {
                estimateQuantile();
            }


            size_t goods = 0;

            std::fstream minimal_summary(out_name + "_summary", std::ios_base::out);
            std::fstream f_lambdas(out_name + "_lambda", std::ios_base::out);
            std::fstream f_left_inliers_original_picture(out_name + "_left", std::ios_base::out);
            std::fstream f_right_inliers_original_picture(out_name + "_right", std::ios_base::out);
            std::fstream f_left_inliers_undistorted_picture(out_name + "_undistorted_left", std::ios_base::out);
            std::fstream f_right_inliers_undistorted_picture(out_name + "_undistorted_right", std::ios_base::out);
            std::fstream f_recomputed_fundamental_matrix(out_name + "_recomputed_F", std::ios_base::out);
            std::fstream f_errors(out_name + "_all_errors", std::ios_base::out);
            std::fstream f_left_inliers_unscaled_undistorted_picture(out_name + "_unscaled_undistorted_left",
                                                                     std::ios_base::out);
            std::fstream f_right_inliers_unscaled_undistorted_picture(out_name + "_unscaled_undistorted_right",
                                                                      std::ios_base::out);
            std::fstream f_left_nearest_backprojected_point(out_name + "_curve_left", std::ios_base::out);
            std::fstream f_right_nearest_backprojected_point(out_name + "_curve_right", std::ios_base::out);
            Eigen::Matrix<T, 1, 2> center;
            center(0, 0) = w_ / T(2.0);
            center(0, 1) = h_ / T(2.0);


            double sum_squared_err = 0;

            minimal_summary << "Distortion coefficients: \n";
            for (size_t k = 0; k < hyp_lambdas_.size(); ++k) {
                f_errors << hyp_lambdas_[k] << " ";
                f_lambdas << hyp_lambdas_[k] << " ";
            }

            minimal_summary << "\nFundamental matrix: \n" << hyp_F_ << "\n";
            minimal_summary << "Quantile and confidence interval: " << quantile_ << " " << confidence_interval_ << "\n";
            inliers_ind_.resize(0);

            for (size_t k = 0; k < u1d_.cols(); ++k) {

                f_errors << errors_[k] << "\n";
                if (std::abs(errors_[k]) - confidence_interval_ < EPS) {
                    sum_squared_err += (left_residuals_[k] * left_residuals_[k] +
                                        right_residuals_[k] * right_residuals_[k]);
                    inliers_ind_.push_back(k);
                    ++goods;
                    f_left_inliers_original_picture << r_ * u1d_.col(k).transpose().leftCols(2) + center << "\n";
                    f_right_inliers_original_picture << r_ * u2d_.col(k).transpose().leftCols(2) + center << "\n";
                    f_left_inliers_undistorted_picture << r_ * alpha_ * u1_.col(k).transpose() + center << "\n";
                    f_right_inliers_undistorted_picture << r_ * alpha_ * u2_.col(k).transpose() + center << "\n";
                    f_left_inliers_unscaled_undistorted_picture << u1_.col(k).transpose() << "\n";
                    f_right_inliers_unscaled_undistorted_picture << u2_.col(k).transpose() << "\n";
                    f_left_nearest_backprojected_point << r_ * nearest_u1d_.col(k).transpose() + center << "\n";
                    f_right_nearest_backprojected_point << r_ * nearest_u2d_.col(k).transpose() + center << "\n";

                }

            }

            Eigen::Matrix3d shift;
            Eigen::Matrix3d scale;
            shift << 1, 0, -w_ / 2.0,
                    0, 1, -h_ / 2.0,
                    0, 0, 1;
            scale << 1.0 / (alpha_ * r_), 0, 0,
                    0, 1.0 / (alpha_ * r_), 0,
                    0, 0, 1;

            Eigen::Matrix3d recompute_F = shift.transpose() * scale.transpose() * hyp_F_ * scale * shift;
            f_recomputed_fundamental_matrix << recompute_F;
            minimal_summary << "Number of inliers: " << goods << "\nSquared error: " << sum_squared_err << std::endl;
            return goods;
        }


        void setHypLambdas(const Eigen::Matrix<T, -1, 1> &hyp_lambdas_) {
            UndistortionProblemHelper::hyp_lambdas_ = hyp_lambdas_;
        }

        void setHypF(const Eigen::Matrix<T, 3, 3> &hyp_F_) {
            UndistortionProblemHelper::hyp_F_ = hyp_F_;
        }


        double getExpectedPercentOfInliers() const {
            return expected_percent_of_inliers_;
        }

        double getConfidenceInterval() const {
            return confidence_interval_;
        }

        T getW() const {
            return w_;
        }

        T getH() const {
            return h_;
        }

        T getR() const {
            return r_;
        }

        T getAlpha() const {
            return alpha_;
        }

        const Points &getU1d() const {
            return u1d_;
        }

        const Points &getU2d() const {
            return u2d_;
        }

        const Points &getU1() const {
            return u1_;
        }

        const Points &getU2() const {
            return u2_;
        }

        const Points &getNearestU1d() const {
            return nearest_u1d_;
        }

        const Points &getNearestU2d() const {
            return nearest_u2d_;
        }

        long getNumberOfPoints() const {
            return number_of_points_;
        }

        const Eigen::Matrix<T, -1, 1> &getLambdas() const {
            return hyp_lambdas_;
        }

        const Eigen::Matrix<T, 3, 3> &getF() const {
            return hyp_F_;
        }

        T getQuantile() const {
            return quantile_;
        }

        const std::vector<T> &getLeftResiduals() const {
            return left_residuals_;
        }

        const std::vector<T> &getRightResiduals() const {
            return right_residuals_;
        }

        const std::vector<T> &getErrors() const {
            return errors_;
        }

        const std::vector<int> &getInliersIndices() const {
            return inliers_ind_;
        }
    };
}
#endif //AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H
