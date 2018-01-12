//
// Created by danielbord on 12/18/17.
//

#ifndef AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H
#define AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H

#include <Eigen/Dense>
#include <iostream>

class FocalLengthEstimator {
    static double constexpr scale_focal0 = 1e3;
    Eigen::Matrix3d F;

private:
    Eigen::Vector3d makeQuadricEquation() {
        Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        //std::cout << "SVD:\n" << fmatrix_svd.singularValues().transpose() << std::endl;

        double a, b, u13, u23, v13, v23;
        a = fmatrix_svd.singularValues()[0];
        b = fmatrix_svd.singularValues()[1];

        u13 = fmatrix_svd.matrixU()(2, 0);
        u23 = fmatrix_svd.matrixU()(2, 1);
        v13 = fmatrix_svd.matrixV()(2, 0);
        v23 = fmatrix_svd.matrixV()(2, 1);

        Eigen::Vector3d res;
        res[2] =   a * a * (1 - u13 * u13) * (1 - v13 * v13)
                 - b * b * (1 - u23 * u23) * (1 - v23 * v23);
        res[1] =   a * a * (u13 * u13 + v13 * v13 - 2 * u13 * u13 * v13 * v13)
                 - b * b * (u23 * u23 + v23 * v23 - 2 * u23 * u23 * v23 * v23);
        res[0] = a * a * u13 * u13 * v13 * v13 - b * b * u23 * u23 * v23 * v23;
        res = res / res[2];


        //double ff2 = (-u23*v13*(a*u13*v13+b*u23*v23))/(a*u13*u23*(1-v13*v13)+b*v13*v23*(1-u23*u23));
        //std::cout << "FF2: "<< ff2 << " " << 1.0/ff2 << std::endl;
        return res;


    }

public:
    FocalLengthEstimator() {F.setZero();}
    explicit FocalLengthEstimator(const Eigen::Matrix3d &F) : F(F) {}

    void setF(const Eigen::Matrix3d &F) {
        FocalLengthEstimator::F = F;
    }


    double estimate() {
        Eigen::Vector3d coeff = makeQuadricEquation();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> companion(2, 2);
        companion.setZero();
        companion.col(1) = -1 * coeff.topRows(2);
        companion(1, 0) = 1;
        //std::cout << coeff.transpose() << std::endl;
        //std::cout << companion << std::endl;
        double f2 = 1;
        auto r = companion.eigenvalues();
        for (int kk = 0; kk < r.rows(); ++kk) {
            double real = r[kk].real();
            double imag = r[kk].imag();

            std::cout << r[kk] << std::endl;
            if (std::abs(imag) < 1e-9 && real > 0) {
                f2 = real;
                std::cout << "real root: " << real << std:: endl;
            }

        }
        return std::sqrt(f2);
    }
};

#endif //AUTOMATICSOLVER_FOCAL_LENGTH_SOLVER_H
