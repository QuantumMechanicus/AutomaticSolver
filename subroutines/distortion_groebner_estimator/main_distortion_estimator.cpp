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
    double w, h, r;
    std::string input1, input2, f_lambda_distribution, f_estimated_lambda, f_estimated_fundamental_matrix, f_inliers_list;
    int iter;
    double lambda_upper_bound;
    double lambda_lower_bound;
    double prcnt_inl;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1)->required(),
                 "Input filename1  --- first camera correspondent points")
                ("input2", po::value<std::string>(&input2)->required(),
                 "Input filename2 --- second camera correspondent points")
                ("distr_f", po::value<std::string>(&f_lambda_distribution)->default_value(
                        "./automatic_solver_results/lambdas_distribution"),
                 "Output file for estimated distortion parameter distribution")
                ("up_threshold", po::value<double>(&lambda_upper_bound)->default_value(0.25), "Lambda upper threshold")
                ("low_threshold", po::value<double>(&lambda_lower_bound)->default_value(-2), "Lambda lower threshold")
                ("fund_f", po::value<std::string>(&f_estimated_fundamental_matrix)->default_value(
                        "./automatic_solver_results/estimated_f"), "Output file for fundamental matrix estimation")
                ("lambd_f", po::value<std::string>(&f_estimated_lambda)->default_value(
                        "./automatic_solver_results/estimated_lambda"), "Output file for lambda estimation")
                ("iters", po::value<int>(&iter)->default_value(10000), "Number of iterations")
                ("w", po::value<double>(&w)->default_value(7360), "Width")
                ("h", po::value<double>(&h)->default_value(4912), "Height")
                ("q", po::value<double>(&prcnt_inl)->default_value(0.1), "quantile to minimize")
                ("inliers_f",
                 po::value<std::string>(&f_inliers_list)->default_value("./automatic_solver_results/inliers"),
                 "Output file for inliers");

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
    Eigen::Matrix3d fundamental_matrix;
    double lambda;
    eight_points_problem::AutomaticEstimator estimator(w, h, u1d, u2d, prcnt_inl);
    double q = estimator.estimate(fundamental_matrix, lambda, f_lambda_distribution, f_inliers_list, iter, lambda_lower_bound,
                                  lambda_upper_bound);


    std::fstream output_f1(f_estimated_fundamental_matrix, std::fstream::out);
    std::fstream output_f2(f_estimated_lambda, std::fstream::out);


    fundamental_matrix /= fundamental_matrix(2, 2);

    output_f2 << lambda << "\n";
    output_f1 << fundamental_matrix << "\n";

    output_f1.close();
    output_f2.close();
    return 0;
}