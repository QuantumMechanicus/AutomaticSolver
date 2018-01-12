//
// Created by danielbord on 12/18/17.
//
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include "focal_length_solver.h"
#include "Sophus/sophus/so3.hpp"
#include "undistortion_problem_utils.h"
#include <boost/math/special_functions/erf.hpp>
int main(int argc, char *argv[]) {

    /*std::string f_f;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("fund_f", po::value<std::string>(&f_f)->required(), "File with %n_pic estimated fundamental matrices");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);


        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return -1;
        }
        boost::program_options::notify(vm);

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -2;
    }*/
    std::cout << boost::math::erfc_inv((0.95 + 1.0)) / boost::math::erfc_inv((0.3 + 1.0));
    Eigen::Matrix<long double, 2, 1> lmbd;
    lmbd(0) = -0.989201;
    lmbd(1) = 0.00703376;
    double rr =1.663544010175267/2;

    std::cout << "al:" <<   undistortion_utils::undistortionDenominator<long double>(rr, lmbd) << std::endl;
    Eigen::Matrix3d f, k, k2, r, t, f2, kkk;

    r = Sophus::SO3d::rotX(0.4).matrix() * Sophus::SO3d::rotY(0.9).matrix();
    t.setZero();
    Eigen::Vector3d tt;
    kkk.setZero();
    kkk(0, 0) = kkk(1,1) = 3115.0/1000.0;
    kkk(2,2) = 1;
    tt(0) = 1;
    tt(1) = 1;
    tt(2) = -1;
    tt.normalize();
    t(0, 1) = -tt(2);
    t(0, 2) = tt(1);
    t(1, 0) = tt(2);
    t(1, 2) = -tt(0);
    t(2, 0) = -tt(1);
    t(2, 1) = tt(0);

    k.setZero();
    k(0, 0) = 1500;
    k(1, 1) = 1500;
    k(2, 2) = 1;

    k2.setZero();
    k2(0, 0) = 1500;
    k2(1, 1) = 1500;
    k2(2, 2) = 1;
    f = k.inverse() * t * r * k2.inverse();
    double fnorm = f.norm();
    f = f / fnorm;
    std::cout << std::endl << "F: \n" << f << std::endl << std::endl;
    FocalLengthEstimator est(f);
    auto ans = est.estimate();
    std::cout << "Estimated: " << ans << std::endl;

    std::fstream ff("/home/danielbord/CLionProjects/AutomaticSolver/subroutines/focal_length_estimator/testF.txt", std::fstream::in);
    double allf = 0;
    int c = 0;
    while (!ff.eof()) {
        std::cout <<"GO\n";
        double f11, f12, f13, f21, f22, f23, f31, f32, f33;

        ff >> f11 >> f12 >> f13 >> f21 >> f22 >> f23 >> f31 >> f32 >> f33;

        Eigen::Matrix3d mf, me;

        mf << f11, f12, f13, f21, f22, f23, f31, f32, f33;



        std::cout << "Fund:\n" << mf << std::endl;
        mf = mf/mf.norm();
        est.setF(mf);
        double est_f = est.estimate();

        std::cout << "est f: " << est_f << std::endl;

        allf += est_f;
        ++c;
        k.setZero();
        k(0, 0) = est_f;
        k(1, 1) =est_f;
        k(2, 2) = 1;
        me = mf;
        std::cout << "Ess:\n" << mf << std::endl;

    }
    std::cout << "mean: " << allf/c << " " << 180*2*std::atan(c/(allf))/M_PI<< std::endl;
    /*f << 27.948, -2381.82 , 26.8546,
    2378.93 , 33.8998 , 1015.57,
                      -12.1748 ,-988.609 ,0.994689;
    f = f/f.norm();
    est.setF(f);
    std::cout << est.estimate() << std::endl;

    f2 <<-6.69767 , -562.042 ,  6.48043,
    590.018  , -5.211 , 230.666,
                      -8.12151, -257.849 , 0.99994;
    f2 = f2/f2.norm();
    est.setF(f2);
    std::cout << est.estimate() << std::endl;

    f2 << -43.1166 , 3894.88, -58.0881,
                    -3970.98, -14.2969, -1651.18,
    44.7582 , 1704.44 , 1.03902;
    f2 = f2/f2.norm();
    est.setF(f2);
    std::cout << est.estimate() << std::endl;

    f2 << -5.92595 ,-331.757,  5.19802,
    348.025 ,-5.81519 ,  142.641,
                      -6.37953, -157.177,  1.00003;

    f2 = f2/f2.norm();
    est.setF(f2);
    std::cout << est.estimate() << std::endl;

    f2*/
    return 0;
}