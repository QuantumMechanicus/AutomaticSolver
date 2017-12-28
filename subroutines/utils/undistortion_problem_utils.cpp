//
// Created by danielbord on 12/12/17.
//

#include "undistortion_problem_utils.h"

namespace undistortion_utils {
    UndistortionProblemHelper::UndistortionProblemHelper(double w, double h, double r, const Points &u1d,
                                                         const Points &u2d,
                                                         double percent) : w_(w),
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
        quantile_ = double(0);
        normalizePoints<double>(u1d_, w, h, r);
        normalizePoints<double>(u2d_, w, h, r);
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

        alpha_ = 0;

    }


    void UndistortionProblemHelper::computeErrors() {
        errors_.resize(number_of_points_, std::numeric_limits<double>::max());
        for (size_t k = 0; k < number_of_points_; ++k) {
            Eigen::Matrix<double, 2, 1> u1k, cu1k;
            Eigen::Matrix<double, 2, 1> u2k, cu2k;
            u1k.setZero();
            u2k.setZero();
            bool is_correct;

            auto errors = undistortion_utils::computeEpipolarCurveDistanceError<double>(hyp_lambdas_, hyp_F_,
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


    void UndistortionProblemHelper::estimateQuantile() {
        computeErrors();

        std::vector<double> errors(errors_);
        std::nth_element(errors.begin(), errors.begin() + int(errors.size() * expected_percent_of_inliers_),
                         errors.end());
        quantile_ = errors[int(errors.size() * expected_percent_of_inliers_) + 1];
        confidence_interval_ =
                quantile_ * boost::math::erfc_inv((0.95 + 1.0)) /
                boost::math::erfc_inv((expected_percent_of_inliers_ + 1.0));
    }


    size_t UndistortionProblemHelper::findInliers(const std::string &out_name) {

        alpha_ = undistortionDenominator<double>(( std::max(h_, w_) / (r_*double(2.0))), hyp_lambdas_);
        estimateQuantile();
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
        Eigen::Matrix<double, 1, 2> center;
        center(0, 0) = w_ / double(2.0);
        center(0, 1) = h_ / double(2.0);


        double sum_squared_err = 0;

        minimal_summary << "Distortion coefficients: \n";
        for (size_t k = 0; k < hyp_lambdas_.size(); ++k) {
            minimal_summary << hyp_lambdas_[k] << " ";
            f_lambdas << hyp_lambdas_[k] << " ";
        }
        f_lambdas << "\n" << alpha_ << std::endl;
        minimal_summary << "\nFundamental matrix: \n" << hyp_F_ << "\n";
        minimal_summary << "Quantile and confidence interval: " << quantile_ << " " << confidence_interval_ << " " << errors_.size() << " " << expected_percent_of_inliers_ << "\n";
        inliers_ind_.resize(0);

        for (size_t k = 0; k < u1d_.cols(); ++k) {

            f_errors << errors_[k] << "\n";
            if (errors_[k] - confidence_interval_ < EPS) {
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


    void UndistortionProblemHelper::setHypLambdas(const Eigen::Matrix<double, -1, 1> &hyp_lambdas_) {
        UndistortionProblemHelper::hyp_lambdas_ = hyp_lambdas_;
    }

    void UndistortionProblemHelper::setHypF(const Eigen::Matrix<double, 3, 3> &hyp_F_) {
        UndistortionProblemHelper::hyp_F_ = hyp_F_;
    }


    double UndistortionProblemHelper::getExpectedPercentOfInliers() const {
        return expected_percent_of_inliers_;
    }

    double UndistortionProblemHelper::getConfidenceInterval() const {
        return confidence_interval_;
    }

    double UndistortionProblemHelper::getW() const {
        return w_;
    }

    double UndistortionProblemHelper::getH() const {
        return h_;
    }

    double UndistortionProblemHelper::getR() const {
        return r_;
    }

    double UndistortionProblemHelper::getAlpha() const {
        return alpha_;
    }

    const Points &UndistortionProblemHelper::getU1d() const {
        return u1d_;
    }

    const Points &UndistortionProblemHelper::getU2d() const {
        return u2d_;
    }

    const Points &UndistortionProblemHelper::getU1() const {
        return u1_;
    }

    const Points &UndistortionProblemHelper::getU2() const {
        return u2_;
    }

    const Points &UndistortionProblemHelper::getNearestU1d() const {
        return nearest_u1d_;
    }

    const Points &UndistortionProblemHelper::getNearestU2d() const {
        return nearest_u2d_;
    }

    long UndistortionProblemHelper::getNumberOfPoints() const {
        return number_of_points_;
    }

    const Eigen::Matrix<double, -1, 1> &UndistortionProblemHelper::getLambdas() const {
        return hyp_lambdas_;
    }

    const Eigen::Matrix<double, 3, 3> &UndistortionProblemHelper::getF() const {
        return hyp_F_;
    }

    double UndistortionProblemHelper::getQuantile() const {
        return quantile_;
    }

    const std::vector<double> &UndistortionProblemHelper::getLeftResiduals() const {
        return left_residuals_;
    }

    const std::vector<double> &UndistortionProblemHelper::getRightResiduals() const {
        return right_residuals_;
    }

    const std::vector<double> &UndistortionProblemHelper::getErrors() const {
        return errors_;
    }

    const std::vector<int> &UndistortionProblemHelper::getInliersIndices() const {
        return inliers_ind_;
    }

}