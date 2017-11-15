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

int main(int argc, char *argv[]) {

    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2, distr_f;
    int iter = 10000;
    double tr = 0.25;
    double tr2 = -1;
    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1)->required(),
                 "Input filename1 (matches and camera parameters)")
                ("input2", po::value<std::string>(&input2)->required(),
                 "Input filename2 (matches and camera parameters)")
                ("f", po::value<std::string>(&distr_f)->required(),
                        "Output file for lambdas")
                ("threshold1", po::value<double>(&tr), "Lambda threshold")
                ("threshold2", po::value<double>(&tr2), "Lambda threshold")
                ("iters", po::value<int>(&iter), "Number of iterations")
                ("w", po::value<double>(&w), "Width")
                ("h", po::value<double>(&h), "Height");

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

    Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
    readPointsFromFile(input1, u1d);
    readPointsFromFile(input2, u2d);

    Eigen::Matrix3d F;
    double Lambda;

    size_t inl = getFundamentalMatrixAndLambda(u1d, u2d, w, h, F, Lambda, distr_f, iter, tr, tr2);
    std::cout << std::endl;
    /*std::fstream output_f(input1 + "_results", std::fstream::app);
    output_f << "\n" << inl << " - number of inliesrs \n" << Lambda << " - lambda\n\nF:\n\n" << F << std::endl;

    output_f.close();
    double f_err = 0;
    goodPoints(u1d, u2d, Lambda, F, 5, f_err, w, h);*/
    return 0;
}