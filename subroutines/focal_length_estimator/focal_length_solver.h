//
// Created by danielbord on 12/18/17.
//

#ifndef AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H
#define AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H

#include <Eigen/Dense>

class FocalLengthEstimator
{
    Eigen::Matrix3d F;
    void findEpipoles(Eigen::Vector3d &left_epipole, Eigen::Vector3d &right_epipole)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        left_epipole = fmatrix_svd.matrixU().col(2);
        right_epipole = fmatrix_svd.matrixV().col(2);

    }

    void rotateEpipoles(Eigen::Vector3d &left_epipole, Eigen::Vector3d &right_epipole)
    {

    }
public:

    double estimate()
    {
        Eigen::Matrix3d recompute_F, left_e, right_e;
        left_e.setZero();
        right_e.setZero();
        left_e(0, 0) = right_epipole(2);
        left_e(0, 0) = 1;
        left_e(2, 2) = -right_epipole(0);

        right_e(0, 0) = left_epipole(2);
        right_e(0, 0) = 1;
        right_e(2, 2) = -left_epipole(0);
        //TODO ROTATE EPIPOLES AND F
        recompute_F = left_e.inverse()*F*right_e.inverse();
        double a, b, c, d;
        a = recompute_F(0,0);
        b = recompute_F(0, 1);
        c = recompute_F(1, 0);
        d = recompute_F(1, 1);
        double k2 = -a*c*left_epipole(0)*left_epipole(0)/(a*c*left_epipole(2)*left_epipole(2) + b*d);

        return sqrt(k2);
    }
};

#endif //AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H
