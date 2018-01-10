//
// Created by danielbord on 9/11/17.
//

#include "test_utils.h"

namespace test_utils {
    void generateSmallRotationAndBias(Eigen::Matrix3d &matrixR,
                                      Eigen::Vector3d &vectorT) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(-0.5, 1);
        matrixR = (Sophus::SO3d::rotX(distribution(gen)) * Sophus::SO3d::rotY(distribution(gen)) *
                   Sophus::SO3d::rotZ(distribution(gen))).matrix();
        vectorT = Eigen::Vector3d::Random();
        vectorT.normalize();
    }

    void generateRotation(Sophus::SO3d &omega)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        srand((unsigned int) gen());
        omega = Sophus::SO3d(Eigen::Quaterniond::UnitRandom());
    }

    void addRotationNoise(double angle, double length, Sophus::SO3d &omega) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(angle - length / 2.0, angle + length / 2.0);
        std::uniform_int_distribution<int> coin(0, 1);
        std::srand((unsigned int) gen());

        if (coin(gen) == 1)
        {
            auto tmp = omega.unit_quaternion();
            tmp.coeffs() *= -1;
            omega.setQuaternion(tmp);
        }

        Eigen::Quaterniond shift;
        Eigen::Vector3d ax;

        double random_angle = distribution(gen);
        ax = Eigen::Vector3d::Random();
        shift.w() = std::cos(random_angle/2.0);
        shift.vec() = std::sin(random_angle/2.0)*ax;

        omega *= Sophus::SO3d(shift);

    }

    Eigen::Vector3d generateTranslationPerturbation()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(-0.003, 0.003);
        Eigen::Vector3d res;
        res << distribution(gen), distribution(gen), distribution(gen);
        return res;
    }

    Sophus::SO3d generateRotationPerturbation()
    {
        return Sophus::SO3d::exp(generateTranslationPerturbation());
    }


    void generateCalibrationMatrix(Eigen::Matrix3d &matrixK) {
        matrixK.setZero();
        matrixK(0, 0) = 1;
        matrixK(1, 1) = 1;
        matrixK(2, 2) = 1;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> urd(10, 50);


        double focal = urd(gen);
        //double p_x = urd(gen) / 5 - 5;
        //double p_y = urd(gen) / 5 - 5;



        matrixK(0, 0) = focal;
        matrixK(1, 1) = focal;
        matrixK(2, 2) = 1;


    }

    Eigen::Vector3d generatePointOnSphere(const Eigen::Vector3d &shift, double raidus)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0, 1);
        double x, y, z, norm;
        x = distribution(gen);
        y = distribution(gen);
        z = distribution(gen);
        norm = std::sqrt(x * x + y * y + z * z);
        return (Eigen::Vector3d(raidus*x / norm, raidus*y / norm, raidus*z / norm) + shift);
    }

    void generateWorldSphere(Eigen::Matrix3d &matrixR, Eigen::Vector3d &vectorT,
                             Eigen::Matrix<double, 3, Eigen::Dynamic> &worldPoints1,
                             Eigen::Matrix<double, 3, Eigen::Dynamic> &worldPoints2,
                             std::size_t size) {

        worldPoints1.resize(Eigen::NoChange, size);
        worldPoints2.resize(Eigen::NoChange, size);



        for (size_t k = 0; k < size; ++k) {
            Eigen::Vector3d shift;
            shift << 10, 10, 10;
            worldPoints1.col(k) = generatePointOnSphere(shift);
        }
        Sophus::SO3d d(matrixR);
        worldPoints2 = matrixR * worldPoints1;
        for (size_t k = 0; k < size; ++k)
            worldPoints2.col(k) += vectorT;

    }

    void addNoise( Eigen::Matrix<double, 3, Eigen::Dynamic> &worldPoints)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0, 0.0005);
        for (size_t j = 0; j < worldPoints.cols(); ++j)
        {
            Eigen::Vector3d noise;
            noise << distribution(gen), distribution(gen), distribution(gen);
            worldPoints.col(j) += noise;
        }
    }

    void generateTest(Eigen::Matrix3d &matrixR, Eigen::Vector3d &vectorT,
                      Eigen::Matrix3d &matrixK,
                      Eigen::Matrix<double, 3, Eigen::Dynamic> &imagePoints1,
                      Eigen::Matrix<double, 3, Eigen::Dynamic> &imagePoints2,
                      std::size_t numberOfPoints, bool noise) {

        Eigen::Matrix<double, 3, Eigen::Dynamic> worldPoints1, worldPoints2;
        generateCalibrationMatrix(matrixK);
        generateSmallRotationAndBias(matrixR, vectorT);
        generateWorldSphere(matrixR, vectorT, worldPoints1, worldPoints2,
                            numberOfPoints);
        if (noise) {
            addNoise(worldPoints1);
            addNoise(worldPoints2);
        }

        imagePoints1 = matrixK * worldPoints1;
        imagePoints2 = matrixK * worldPoints2;


    }

}
