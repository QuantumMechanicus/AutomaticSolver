//
// Created by danielbord on 10/30/17.
//

#ifndef AUTOMATICSOLVER_SOLVER_KU8PT_H
#define AUTOMATICSOLVER_SOLVER_KU8PT_H

#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <limits>


class AutomaticEstimator {
public:
    typedef Eigen::Matrix<double, 7, 1> G_polynomial;
    typedef Eigen::Matrix<double, 2, 8> EightPoints;
    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> Points;
    typedef std::pair<std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>, std::vector<double> > PairOfMatricesFandLambdas;

private:
    static constexpr double EPS = 1e-8;
    double w_, h_;
    Points u1d_, u2d_;

    PairOfMatricesFandLambdas solver_ku8pt(const G_polynomial &g1, const G_polynomial &g2, const G_polynomial &g3,
                                           const G_polynomial &g4, const G_polynomial &g5, const G_polynomial &g6,
                                           const G_polynomial &g7, const G_polynomial &g8);

    PairOfMatricesFandLambdas run_solver8pt(EightPoints u1d, EightPoints u2d);

    void computeErrors(double hyp_lambda, const Eigen::Matrix3d &hyp_F, std::vector<double> &left_errors,
                       std::vector<double> &right_errors);

    double estimateQuantile(double hyp_lambda, const Eigen::Matrix3d &hyp_F);

    size_t findInliers(double hyp_lambda, const Eigen::Matrix3d &hyp_F, double quantile, const std::string &out_name);

public:

    AutomaticEstimator(double w, double h, const Points &u1d, const Points &u2d);

    double estimate(Eigen::Matrix3d &F, double &Lambda, const std::string &name_f,
                    const std::string &inliers_f,
                    int numberOfIterations, double threshold = 0.25, double threshold2 = -1);

};


#endif //AUTOMATICSOLVER_SOLVER_KU8PT_H
