//
// Created by danielbord on 12/12/17.
//

#ifndef AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H
#define AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <fstream>

namespace undistortion_utils {

    template<typename T>
    T undistortionDenominator(const T &r_distorted2, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) {
        T denominator(1.0);
        T r_distorted2_pow = r_distorted2;
        for (int i = 0; i < lambdas.rows(); ++i) {
            denominator += lambdas[i] * r_distorted2_pow;
            r_distorted2_pow *= r_distorted2;
        }
        return denominator;
    }

    template<typename T>
    std::pair<T, T>
    computeError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 3, 3> &hyp_F,
                 const Eigen::Matrix<T, 2, 1> &u1d,
                 const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h, Eigen::Matrix<T, 2, 1> &u1,
                 Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {
        T wT(w), hT(h);
        T d = std::max(hT, wT)/2.0;
        T dr = d / r;
        T dr2 = dr * dr;
        T alpha = undistortionDenominator(dr2, hyp_lambdas);
        //T alpha = 4 / (4 + hyp_lambda * (d / r_) * (d / r_));

        T r1d, r2d;
        r1d = u1d.squaredNorm();
        r2d = u2d.squaredNorm();

        Eigen::Matrix<T, 3, 1> homogeneous_u1, homogeneous_u2;
        T denominator1 = undistortionDenominator(r1d, hyp_lambdas);
        T denominator2 = undistortionDenominator(r2d, hyp_lambdas);

        is_correct = (denominator1 > T(0) && denominator2 > T(0));
        u1 = u1d / denominator1;
        u2 = u2d / denominator2;
        homogeneous_u1 = u1.homogeneous();
        homogeneous_u2 = u2.homogeneous();
        Eigen::Matrix<T, 3, 1> l1 = hyp_F * homogeneous_u1;
        Eigen::Matrix<T, 3, 1> l2 = hyp_F.transpose() * homogeneous_u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(homogeneous_u1);
        T err2 = err * err;
        std::pair<T, T> residuals;
        residuals.first = alpha * r * err / n1;
        residuals.second = alpha * r * err / n2;
        return residuals;
    }


    template<typename T>
    class UndistortionProblemHelper {
    public:
        typedef Eigen::Matrix<T, 2, Eigen::Dynamic> Points;
        using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;
        static constexpr double EPS = 1e-8;
        static constexpr double confidence_interval = 3.36752;
    private:
        T w_, h_, r_;
        Points u1d_, u2d_, u1_, u2_;
        int number_of_points_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> hyp_lambdas_;
        Eigen::Matrix<T, 3, 3> hyp_F_;
        T quantile_;

    public:
        UndistortionProblemHelper(T w, T h, T r, const Points &u1d, const Points &u2d) : w_(w), h_(h), r_(r),
                                                                                         u1d_(u1d), u2d_(u2d) {
            assert(u1d.cols() == u2d.cols() && "Points correspondences must be the same size");
            hyp_lambdas_.setZero();
            hyp_F_.setZero();
            number_of_points_ = u1d.cols();
            quantile_ = T(0);
            u1d_.resize(Eigen::NoChange, u1d.cols());
            u2d_.resize(Eigen::NoChange, u1d.cols());
            auto ones = Eigen::VectorXd::Ones(u1d_.cols()).transpose();

            u1d_.row(0) = u1d_.row(0) - ones * w_ / 2.0;

            u1d_.row(1) = u1d_.row(1) - ones * h_ / 2.0;

            u2d_.row(0) = u2d_.row(0) - ones * w_ / 2.0;

            u2d_.row(1) = u2d_.row(1) - ones * h_ / 2.0;
            u1d_ = u1d_ / r_;
            u2d_ = u2d_ / r_;
        }


        void computeErrors(std::vector<T> &left_errors,
                           std::vector<T> &right_errors) {


            left_errors.resize(number_of_points_);
            right_errors.resize(number_of_points_);
            u1_.resize(Eigen::NoChange, number_of_points_);
            u2_.resize(Eigen::NoChange, number_of_points_);


            for (size_t k = 0; k < number_of_points_; ++k) {
                Eigen::Matrix<double, 2, 1> u1k;
                Eigen::Matrix<double, 2, 1> u2k;
                bool is_correct;
                auto errors = undistortion_utils::computeError<double>(hyp_lambdas_, hyp_F_, u1d_.col(k), u2d_.col(k),
                                                                       r_, w_,
                                                                       h_, u1k, u2k, is_correct);
                u1_.col(k) = u1k;
                u2_.col(k) = u2k;
                left_errors[k] = std::abs(errors.first);
                right_errors[k] = std::abs(errors.second);

            }
        }


        void estimateQuantile() {
            std::vector<T> err1, err2;
            computeErrors(err1, err2);
            std::vector<T> err_sample;
            err_sample.resize(number_of_points_);
            for (size_t k = 0; k < u1d_.cols(); ++k) {
                err_sample[k] = (err1[k] + err2[k]) / 2.0;
            }
            std::nth_element(err_sample.begin(), err_sample.begin() + err_sample.size() / 10, err_sample.end());
            quantile_ = err_sample[err_sample.size() / 10];

        }


        size_t findInliers(const std::string &out_name) {
            std::vector<T> err1, err2;
            computeErrors(err1, err2);
            size_t goods = 0;

            std::fstream errf(out_name, std::ios_base::out);
            std::fstream errf1(out_name + "_left", std::ios_base::out);
            std::fstream errf2(out_name + "_right", std::ios_base::out);
            std::fstream errf3(out_name + "_undistorted_left", std::ios_base::out);
            std::fstream errf4(out_name + "_undistorted_right", std::ios_base::out);
            std::fstream errf5(out_name + "_recomp_F", std::ios_base::out);
            std::fstream errf6(out_name + "_all_errs", std::ios_base::out);

            Eigen::Matrix<double, 1, 2> center;
            center(0, 0) = w_ / 2.0;
            center(0, 1) = h_ / 2.0;
            T d = std::max(h_, w_);
            T dr = d / r_;
            T dr2 = dr * dr;
            T alpha = undistortionDenominator(dr2, hyp_lambdas_);

            double full_err = 0;
            errf6 << "Left distorted points: " << u1d_ << "\n\n";
            errf6 << "Right distorted points: " << u2d_ << "\n\n";
            errf6 << "Distortion coefficients: ";
            for (size_t k = 0; k < hyp_lambdas_.size(); ++k)
                errf6 << hyp_lambdas_[k] << " ";
            errf6 << "\n Fundamental matrix: " << hyp_F_ << "\n";
            errf6 << "Quantile: " << quantile_ << "\n";

            for (size_t k = 0; k < u1d_.cols(); ++k) {
                double err = err1[k] + err2[k];
                errf6 << err << "\n";
                if (std::abs(err) < quantile_ * confidence_interval) {
                    full_err += (err1[k] * err1[k] + err2[k] * err2[k]);
                    ++goods;
                    errf1 << r_ * u1d_.col(k).transpose().leftCols(2) + center << "\n";
                    errf2 << r_ * u2d_.col(k).transpose().leftCols(2) + center << "\n";
                    errf3 << r_ * alpha * u1_.col(k).transpose() + center << "\n";
                    errf4 << r_ * alpha * u2_.col(k).transpose() + center << "\n";

                }

            }

            Eigen::Matrix3d shift;
            Eigen::Matrix3d scale;
            shift << 1, 0, -w_ / 2.0,
                    0, 1, -h_ / 2.0,
                    0, 0, 1;
            scale << 1.0 / (alpha * r_), 0, 0,
                    0, 1.0 / (alpha * r_), 0,
                    0, 0, 1;

            Eigen::Matrix3d recompute_F = shift.transpose() * scale.transpose() * hyp_F_ * scale * shift;
            errf5 << quantile_ << "\n" << recompute_F;
            errf << "Number of inliers: " << goods << "\nSquared error: " << full_err << std::endl;

            return goods;
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

        int getNumberOfPoints() const {
            return number_of_points_;
        }

        const Eigen::Matrix<T, -1, 1> &getHypLambdas() const {
            return hyp_lambdas_;
        }

        const Eigen::Matrix<T, 3, 3> &getHypF() const {
            return hyp_F_;
        }

        T getQuantile() const {
            return quantile_;
        }

        void setHypLambdas(const Eigen::Matrix<T, -1, 1> &hyp_lambdas_) {
            UndistortionProblemHelper::hyp_lambdas_ = hyp_lambdas_;
        }

        void setHypF(const Eigen::Matrix<T, 3, 3> &hyp_F_) {
            UndistortionProblemHelper::hyp_F_ = hyp_F_;
        }
    };
}
#endif //AUTOMATICSOLVER_UNDISTORTION_PROBLEM_UTILS_H
