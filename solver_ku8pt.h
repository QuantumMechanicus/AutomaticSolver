//
// Created by danielbord on 10/30/17.
//

#ifndef AUTOMATICSOLVER_SOLVER_KU8PT_H
#define AUTOMATICSOLVER_SOLVER_KU8PT_H

#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <limits>

const double EPS = 1e-8;
typedef Eigen::Matrix<double, 7, 1> G_polynomial;
typedef Eigen::Matrix<double, 2, 8> EightPoints;

typedef std::pair<std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>, std::vector<double> > PairOfMatricesFandLambdas;

PairOfMatricesFandLambdas solver_ku8pt(G_polynomial &g1, G_polynomial &g2, G_polynomial &g3,
                                       G_polynomial &g4, G_polynomial &g5, G_polynomial &g6,
                                       G_polynomial &g7, G_polynomial &g8);

PairOfMatricesFandLambdas run_solver8pt(EightPoints &u1d, EightPoints &u2d);

double getFundamentalMatrixAndLambda(Eigen::Matrix<double, 2, Eigen::Dynamic>
                                     &u1d,
                                     Eigen::Matrix<double, 2, Eigen::Dynamic> &u2d, double w, double h,
                                     Eigen::Matrix3d &F, double &Lambda, const std::string &name_f,
                                     int numberOfIterations, double threshold = 0.25, double threshold2 = -1);

std::size_t updateNumberofIters(double confidence, double error_prob, std::size_t n_points, std::size_t n_iters);

double estimateQuantile(Eigen::Matrix<double, 2, Eigen::Dynamic> &u1d, Eigen::Matrix<double, 2, Eigen::Dynamic> &u2d,
                        double w, double h,
                        double hyp_lambda, Eigen::Matrix3d &hyp_F);

#endif //AUTOMATICSOLVER_SOLVER_KU8PT_H
