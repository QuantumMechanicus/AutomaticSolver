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
    computeError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 3, 3> &hyp_F,
                 const Eigen::Matrix<T, 2, 1> &u1d,
                 const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h, Eigen::Matrix<T, 2, 1> &u1,
                 Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {
        T wT(w), hT(h);
        T d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T alpha = undistortionDenominator(dr, hyp_lambdas);
        //T alpha = 4 / (4 + hyp_lambda * (d / r_) * (d / r_));

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
        T err2 = err * err;
        //std::cout << "err " << err2 << " " << denominator1 << " !! " << denominator2  << "\n" << hyp_F << "\n\n" << homogeneous_u1 << "\n" << l1 << "\n" << l2 << std::endl;
        //exit(0);
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
                break;
            }

        }
        return r_d;
    }

    template<typename T>
    std::pair<T, T>
    computeCurveError(const Eigen::Matrix<T, Eigen::Dynamic, 1> &hyp_lambdas, const Eigen::Matrix<T, 3, 3> &hyp_F,
                      const Eigen::Matrix<T, 2, 1> &u1d,
                      const Eigen::Matrix<T, 2, 1> &u2d, const T &r, double w, double h, Eigen::Matrix<T, 2, 1> &u1,
                      Eigen::Matrix<T, 2, 1> &u2, bool &is_correct) {


        T r1d, r2d;
        r1d = u1d.norm();
        r2d = u2d.norm();

        Eigen::Matrix<T, 3, 1> homogeneous_u1, homogeneous_u2;
        Eigen::Matrix<T, 2, 1> line_point1, line_point2, curve_point1, curve_point2;
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


        //std::cout << "err " << err2 << " " << denominator1 << " !! " << denominator2  << "\n" << hyp_F << "\n\n" << homogeneous_u1 << "\n" << l1 << "\n" << l2 << std::endl;
        //exit(0);
        std::pair<T, T> residuals;
        residuals.first = err / n1;
        residuals.second = err / n2;
        line_point1 = u1 + residuals.first * l1.template block<2, 1>(0, 0) / n1;
        line_point2 = u2 + residuals.second * l1.template block<2, 1>(0, 0) / n1;

        T r1u = u1.squaredNorm();
        T r2u = u2.squaredNorm();

        T root_r1d = findRoot<T>(hyp_lambdas, u1);
        T root_r2d = findRoot<T>(hyp_lambdas, u1);
        T dd = undistortion_utils::undistortionDenominator<T>(r1d, hyp_lambdas);

        curve_point1 = r * dd * u1;
        curve_point2 = r * dd * u2;
        residuals.first = (u1d - curve_point1).norm();
        residuals.second = (u2d - curve_point2).norm();
        return residuals;
    }


    template<typename T>
    class UndistortionProblemHelper {
    public:
        typedef Eigen::Matrix<T, 2, Eigen::Dynamic> Points;
        using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;
        static constexpr double EPS = 1e-10;
    private:
        double prcnt_ = 0.2;
        double confidence_interval_;
        T w_, h_, r_;
        Points u1d_, u2d_, u1_, u2_;
        int number_of_points_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> hyp_lambdas_;
        Eigen::Matrix<T, 3, 3> hyp_F_;
        T quantile_;

    public:
        UndistortionProblemHelper(T w, T h, T r, const Points &u1d, const Points &u2d, double prcnt = 0.1) : w_(w),
                                                                                                             h_(h),
                                                                                                             r_(r),
                                                                                                             u1d_(u1d),
                                                                                                             u2d_(u2d),
                                                                                                             prcnt_(prcnt) {
            assert(u1d.cols() == u2d.cols() && "Points correspondences must be the same size");
            hyp_lambdas_.resize(1, Eigen::NoChange);
            hyp_lambdas_.setZero();
            hyp_F_.setZero();
            number_of_points_ = u1d.cols();
            quantile_ = T(0);
            normalizePoints<T>(u1d_, w, h, r);
            normalizePoints<T>(u2d_, w, h, r);
            confidence_interval_ = 3.36752;

        }


        void computeErrors(std::vector<T> &left_errors,
                           std::vector<T> &right_errors) {


            left_errors.resize(number_of_points_);
            right_errors.resize(number_of_points_);
            u1_.resize(Eigen::NoChange, number_of_points_);
            u2_.resize(Eigen::NoChange, number_of_points_);
            u1_.setZero();
            u2_.setZero();

            for (size_t k = 0; k < number_of_points_; ++k) {
                Eigen::Matrix<T, 2, 1> u1k;
                Eigen::Matrix<T, 2, 1> u2k;
                u1k.setZero();
                u2k.setZero();
                bool is_correct;
                auto errors = undistortion_utils::computeError<T>(hyp_lambdas_, hyp_F_, u1d_.col(k), u2d_.col(k),
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
                err_sample[k] = (err1[k] + err2[k]);
            }

            std::nth_element(err_sample.begin(), err_sample.begin() + int(err_sample.size() * prcnt_),
                             err_sample.end());
            quantile_ = err_sample[int(err_sample.size() * prcnt_)];

        }


        size_t findInliers(const std::string &out_name, std::vector<int> &inliers_ind) {
            if (quantile_ == 0) {
                estimateQuantile();
            }
            confidence_interval_ =
                    quantile_ * boost::math::erfc_inv((0.95 + 1.0)) / boost::math::erfc_inv((prcnt_ + 1.0));

            std::cout << "conf:" << confidence_interval_ << " " << quantile_ << std::endl;

            std::vector<T> err1, err2;
            computeErrors(err1, err2);
            size_t goods = 0;

            std::fstream errf(out_name, std::ios_base::out);
            std::fstream errf1(out_name + "_left", std::ios_base::out);
            std::fstream errf2(out_name + "_right", std::ios_base::out);
            std::fstream errf3(out_name + "_undistorted_left", std::ios_base::out);
            std::fstream errf4(out_name + "_undistorted_right", std::ios_base::out);
            std::fstream errf5(out_name + "_recomp_F", std::ios_base::out);
            std::fstream errf6(out_name + "_all_errors", std::ios_base::out);
            std::fstream errf7(out_name + "_unscaled_undistorted_left", std::ios_base::out);
            std::fstream errf8(out_name + "_unscaled_undistorted_right", std::ios_base::out);

            Eigen::Matrix<T, 1, 2> center;
            center(0, 0) = w_ / T(2.0);
            center(0, 1) = h_ / T(2.0);
            T d = std::max(h_, w_)/T(2.0);
            T dr = d / r_;

            T alpha = undistortionDenominator(dr, hyp_lambdas_);

            double full_err = 0;
            errf6 << "Left distorted points: \n" << u1d_ << "\n\n";
            errf6 << "Right distorted points: \n" << u2d_ << "\n\n";
            errf6 << "Distortion coefficients: ";
            for (size_t k = 0; k < hyp_lambdas_.size(); ++k)
                errf6 << hyp_lambdas_[k] << " ";
            errf6 << "\nFundamental matrix: \n" << hyp_F_ << "\n";
            errf6 << "Quantile: " << quantile_ << " " << confidence_interval_ << "\n";
            inliers_ind.resize(0);
            //errf3 << "TEST:\n" << center << "\n" << u1_ << "\n";
            for (size_t k = 0; k < u1d_.cols(); ++k) {
                double err = (err1[k] + err2[k]);
                errf6 << err << "\n";
                if (std::abs(err) - quantile_ * confidence_interval_ < EPS) {
                    full_err += (err1[k] * err1[k] + err2[k] * err2[k]);
                    inliers_ind.push_back(k);
                    ++goods;
                    errf1 << r_ * u1d_.col(k).transpose().leftCols(2) + center << "\n";
                    errf2 << r_ * u2d_.col(k).transpose().leftCols(2) + center << "\n";
                    errf3 << r_ * alpha * u1_.col(k).transpose() + center << "\n";
                    errf4 << r_ * alpha * u2_.col(k).transpose() + center << "\n";
                    errf7 << u1_.col(k).transpose() << "\n";
                    errf8 << u2_.col(k).transpose() << "\n";

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
            errf5 << recompute_F;
            errf << "Number of inliers: " << goods << "\nSquared error: " << full_err << std::endl;
            errf7 << std::endl << alpha << std::endl;
            errf8 << std::endl << alpha << std::endl;

            return goods;
        }

        size_t findInliers(const std::string &out_name) {
            std::vector<int> dumb_vec;
            return findInliers(out_name, dumb_vec);
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
