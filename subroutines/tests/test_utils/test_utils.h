//
// Created by danielbord on 9/11/17.
//

#ifndef ROTATIONAVERAGING_TEST_UTILS_H
#define ROTATIONAVERAGING_TEST_UTILS_H

#include <random>
#include <fstream>
#include "Eigen/Dense"
#include <sophus/se3.hpp>
#include <ceres/ceres.h>

namespace test_utils {
    void generateSmallRotationAndBias(Eigen::Matrix3d &matrixR,
                                      Eigen::Vector3d &vectorT);

    void generateRotation(Sophus::SO3d &omega);

    void generateCalibrationMatrix(Eigen::Matrix3d &matrixK);

    void addRotationNoise(double angle, double length, Sophus::SO3d &omega);

    Sophus::SO3d generateRotationPerturbation();

    Eigen::Vector3d generateTranslationPerturbation();

    Eigen::Vector3d generatePointOnSphere(const Eigen::Vector3d &shift = Eigen::Vector3d::Zero(), double raidus = 1);

    void generateWorldSphere(Eigen::Matrix3d &matrixR, Eigen::Vector3d &vectorT,
                             Eigen::Matrix<double, 3, Eigen::Dynamic> &worldPoints1,
                             Eigen::Matrix<double, 3, Eigen::Dynamic> &worldPoints2,
                             std::size_t size);

    void generateTest(Eigen::Matrix3d &matrixR, Eigen::Vector3d &vectorT,
                      Eigen::Matrix3d &matrixK,
                      Eigen::Matrix<double, 3, Eigen::Dynamic> &imagePoints1,
                      Eigen::Matrix<double, 3, Eigen::Dynamic> &imagePoints2,
                      std::size_t numberOfPoints, bool noise = false);
}


#endif //ROTATIONAVERAGING_TEST_UTILS_H
