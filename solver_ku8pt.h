//
// Created by danielbord on 10/30/17.
//

#ifndef AUTOMATICSOLVER_SOLVER_KU8PT_H
#define AUTOMATICSOLVER_SOLVER_KU8PT_H

#include "Eigen/Dense"
#include <vector>
#include <iostream>

const double EPS = 1e-8;
typedef Eigen::Matrix<double, 7, 1> G_polynomial;
typedef Eigen::Matrix<double, 2, 8> EightPoints;

typedef std::pair<std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>, std::vector<double> > PairOfMatricesFandLambdas;

PairOfMatricesFandLambdas solver_ku8pt(G_polynomial &g1, G_polynomial &g2, G_polynomial &g3,
                                       G_polynomial &g4, G_polynomial &g5, G_polynomial &g6,
                                       G_polynomial &g7, G_polynomial &g8);

PairOfMatricesFandLambdas run_solver8pt(EightPoints &u1d, EightPoints &u2d);

size_t getFundamentalMatrixAndLambda(Eigen::Matrix<double, 2, Eigen::Dynamic>
                                     &u1d,
                                     Eigen::Matrix<double, 2, Eigen::Dynamic> &u2d,
                                     Eigen::Matrix3d &F, double &Lambda,
                                     int numberOfIterations, double threshold = 1);

#endif //AUTOMATICSOLVER_SOLVER_KU8PT_H
