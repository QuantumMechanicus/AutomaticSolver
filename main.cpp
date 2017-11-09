#include "solver_ku8pt.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;

void readPointsFromFile(std::string &name, Eigen::Matrix<double, 2, Eigen::Dynamic> &m) {
    std::fstream f(name, std::fstream::in);
    std::vector<std::pair<double, double>> v;
    while (!f.eof()) {
        double x, y;
        f >> x;
        f >> y;
        if (f.eof())
            break;
        v.emplace_back(std::make_pair(x, y));
    }
    m.resize(Eigen::NoChange, v.size());
    for (size_t k = 0; k < v.size(); ++k) {
        m(0, k) = v[k].first;
        m(1, k) = v[k].second;
    }
    f.close();
}

size_t goodPoints(Eigen::Matrix<double, 2, Eigen::Dynamic> &u1d, Eigen::Matrix<double, 2, Eigen::Dynamic> &u2d,
                         double hyp_lambda, Eigen::Matrix3d &hyp_F, double threshold,
                         double &modelErr) {
    //hyp_F.transposeInPlace();
    Eigen::Matrix<double, 2, Eigen::Dynamic> r1d;
    r1d.resize(Eigen::NoChange, u1d.cols());
    r1d.row(0) = (u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));
    r1d.row(1) = (u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));

    Eigen::Matrix<double, 2, Eigen::Dynamic> r2d;
    r2d.resize(Eigen::NoChange, u2d.cols());

    r2d.row(0) = (u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1)));
    r2d.row(1) = (u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1)));


    Eigen::Matrix<double, 2, Eigen::Dynamic> ones;
    ones.resize(Eigen::NoChange, r1d.cols());
    ones.setOnes();
    modelErr = 0;
    size_t goods = 0;
    auto u1 = u1d.cwiseProduct((ones + hyp_lambda * r1d).cwiseInverse());
    auto u2 = u2d.cwiseProduct((ones + hyp_lambda * r2d).cwiseInverse());

    Eigen::Matrix<double, 3, Eigen::Dynamic> uu1, uu2;
    uu1.resize(Eigen::NoChange, u1.cols());
    uu2.resize(Eigen::NoChange, u2.cols());
    uu1.row(0) = u1.row(0);
    uu2.row(0) = u2.row(0);

    uu1.row(1) = u1.row(1);
    uu2.row(1) = u2.row(1);


    uu1.row(2).setOnes();
    uu2.row(2).setOnes();

    Eigen::Matrix<double, 3, Eigen::Dynamic> l1 = (hyp_F * uu1);
    Eigen::Matrix<double, 3, Eigen::Dynamic> l2 = (hyp_F.transpose() * uu2);

    std::fstream f_errs("inliers.txt", std::fstream::out | std::ofstream::trunc);
    f_errs << "F:\n" << hyp_F << "\nLambda: " << hyp_lambda << "\nTresh: " << threshold << "\nInliers:\n"
           << std::endl;

    for (size_t k = 0; k < u1d.cols(); ++k) {
        double c1 = l1.col(k).template topRows<2>().norm();
        double c2 = l2.col(k).template topRows<2>().norm();
        double err = std::abs(uu2.col(k).dot(l1.col(k)) / c1) + std::abs(uu1.col(k).dot(l2.col(k)) / c2);

        if (std::abs(err) < threshold) {
            Eigen::Matrix<double, 1, 2> shift;
            shift << 7360.0 / 2.0, 4912.0 / 2.0;
            f_errs << goods << " " << err << "\n";
            f_errs << u1d.col(k).transpose() + shift << "\n";
            f_errs << u2d.col(k).transpose() +shift << "\n";

            ++goods;
        }
        modelErr += err;
    }
    f_errs.close();
    return goods;
}

int main(int argc, char *argv[]) {

    bool doublePrecision;
    int iter;
    std::string input1, input2;

    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1)->required(),
                 "Input filename1 (matches and camera parameters)")
                ("input2", po::value<std::string>(&input2)->required(),
                 "Input filename2 (matches and camera parameters)");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        // clang-format on
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return -1;
        }
        boost::program_options::notify(vm);

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -2;
    }
    std::cout << input1 << std::endl;
    std::cout << input2 << std::endl;
    Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
    readPointsFromFile(input1, u1d);
    readPointsFromFile(input2, u2d);
    std::cout << u1d << std::endl;
    std::cout << u2d << std::endl;
    u1d.row(0) = u1d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * 7360.0 / 2.0;

    u1d.row(1) = u1d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * 4912.0 / 2.0;

    u2d.row(0) = u2d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * 7360.0 / 2.0;

    u2d.row(1) = u2d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * 4912.0 / 2.0;
    u1d = u1d/1000;
    u2d = u2d/1000;
    Eigen::Matrix3d F;
    double Lambda;

    size_t inl = getFundamentalMatrixAndLambda(u1d, u2d, F, Lambda, 10000, 1);

    std::fstream output_f(input1 + "_results", std::fstream::app);
    output_f << "\n" << inl << " - number of inliesrs \n" << Lambda << " - lambda\n\n F:\n\n" << F << std::endl;
    output_f.close();
    double f_err = 0;
    goodPoints(u1d, u2d, Lambda, F, 1, f_err);
    return 0;
}