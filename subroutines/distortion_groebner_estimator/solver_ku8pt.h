//
// Created by danielbord on 10/30/17.
//

#ifndef AUTOMATICSOLVER_SOLVER_KU8PT_H
#define AUTOMATICSOLVER_SOLVER_KU8PT_H


#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <limits>
#include "undistortion_problem_utils.h"

namespace eight_points_problem {
    class AutomaticEstimator {
    public:
        typedef Eigen::Matrix<double, 7, 1> GPolynomial;
        typedef Eigen::Matrix<double, 2, 8> EightPoints;
        typedef Eigen::Matrix<double, 2, Eigen::Dynamic> Points;
        typedef std::pair<std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>, std::vector<double> > FudnamentalMatricesAndDistrotionCoefficients;

    private:
        undistortion_utils::UndistortionProblemHelper helper_;

        FudnamentalMatricesAndDistrotionCoefficients
        solver_ku8pt(const GPolynomial &g1, const GPolynomial &g2, const GPolynomial &g3,
                     const GPolynomial &g4, const GPolynomial &g5, const GPolynomial &g6,
                     const GPolynomial &g7, const GPolynomial &g8);

        FudnamentalMatricesAndDistrotionCoefficients run_solver8pt(EightPoints u1d, EightPoints u2d);


        double estimateQuantile(double hyp_lambda, const Eigen::Matrix3d &hyp_F);

        size_t
        findInliers(double hyp_lambda, const Eigen::Matrix3d &hyp_F, const std::string &out_name);

    public:

        AutomaticEstimator(double w, double h, const Points &u1d, const Points &u2d, double prcnt = 0.1);

        double estimate(Eigen::Matrix3d &F, double &Lambda, const std::string &lambdas_distribution_file,
                        const std::string &inliers_output_file,
                        int number_of_RANSAC_iterations, double lower_threshold = 0.25, double upper_threshold = -1);

    };

}
#endif //AUTOMATICSOLVER_SOLVER_KU8PT_H
