#include "solver_ku8pt.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;

bool readPointsFromFile(std::string &name, Eigen::Matrix<double, 2, Eigen::Dynamic> &m) {
    std::fstream f(name, std::fstream::in);
    if (f.good()) {
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
        return true;
    } else
        return false;
}


int main(int argc, char *argv[]) {
    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2, distr_f = "./automatic_solver_results/lambdas_distribution", f_l = "./automatic_solver_results/estimated_lambda", f_f = "./automatic_solver_results/estimated_f", inliers_f = "./automatic_solver_results/inliers";
    int iter = 10000;
    double lambda_upper_bound = 0.25;
    double lambda_lower_bound = -1;

    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1)->required(),
                 "Input filename1  --- first camera points")
                ("input2", po::value<std::string>(&input2)->required(),
                 "Input filename2 --- second camera points")
                ("distr_f", po::value<std::string>(&distr_f), "Output file for lambdas distribution")
                ("low_threshold", po::value<double>(&lambda_upper_bound), "Lambda upper threshold")
                ("up_threshold", po::value<double>(&lambda_lower_bound), "Lambda lower threshold")
                ("fund_f", po::value<std::string>(&f_f), "Output file for fundamental matrix estimation")
                ("lambd_f", po::value<std::string>(&f_l), "Output file for lambda estimation")
                ("iters", po::value<int>(&iter), "Number of iterations")
                ("w", po::value<double>(&w), "Width")
                ("h", po::value<double>(&h), "Height")
                ("inliers_f", po::value<std::string>(&inliers_f), "Output file for inliers");

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
    }
    Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
    readPointsFromFile(input1, u1d);
    readPointsFromFile(input2, u2d);
    Eigen::Matrix3d F;
    double Lambda;
    AutomaticEstimator estimator(w, h, u1d, u2d);
    double q = estimator.estimate(F, Lambda, distr_f, inliers_f, iter, lambda_upper_bound, lambda_lower_bound);


    std::fstream output_f1(f_f, std::fstream::out);
    std::fstream output_f2(f_l, std::fstream::out);


    F /= F(2, 2);

    output_f2 << Lambda << "\n";
    output_f1 << F << "\n";

    output_f1.close();
    output_f2.close();
    return 0;
}