#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
#include <iomanip>
#include <tbb/tbb.h>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <boost/program_options.hpp>
#include "undistortion_problem_utils.h"


int main(int argc, char *argv[]) {
    std::string input_image_name;
    std::vector<double> v_lambdas;
    std::string out_path = "./";
    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("img_in", po::value<std::string>(&input_image_name)->required(),
                 "Input image name")
                ("out_path", po::value<std::string>(&out_path),
                 "Output directory")
                ("lambdas", po::value<std::vector<double> >(&v_lambdas)->multitoken()->required(), "description");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);

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

    cv::Mat in, out, mx, my;
    in = cv::imread(input_image_name);
    int rows = in.rows;
    int cols = in.cols;


    double r_img = std::sqrt(cols * cols / 4.0 + rows * rows / 4.0);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> map_x(rows, cols), map_y(rows, cols);

    double d = std::max(rows, cols) / 2.0;
    double dr = d / r_img;

    size_t first_non_zero = v_lambdas.size() - 1;
    while (v_lambdas[first_non_zero] == 0)
        --first_non_zero;
    v_lambdas.resize(first_non_zero + 1);
    Eigen::VectorXd lambdas = Eigen::Map<Eigen::VectorXd>(v_lambdas.data(), v_lambdas.size());
    std::cout << lambdas.transpose() << " " << cols << " " << rows << std::endl;

    double alpha = undistortion_utils::undistortionDenominator<double>(dr, lambdas.cast<double>());

    size_t n_lambda = v_lambdas.size();
    size_t deg = 2 * n_lambda;

    tbb::parallel_for(tbb::blocked_range<int>(0, cols), [&](auto range) {

        for (int i = range.begin(); i != range.end(); ++i) {
            for (int j = 0; j < rows; ++j) {
                double ii = ((i - cols / 2.0) / r_img) / alpha;
                double jj = ((j - rows / 2.0) / r_img) / alpha;
                double r_u = std::sqrt(ii * ii + jj * jj + 1e-7);
                Eigen::Matrix<double, Eigen::Dynamic, 1> coeff = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(
                        deg + 1);

                for (size_t k = n_lambda; k > 0; --k)
                    coeff(2 * k, 0) = r_u * lambdas[k - 1];

                coeff(1, 0) = -1;
                coeff(0, 0) = r_u;
                double r_d = std::numeric_limits<double>::max();


                coeff /= coeff[deg];
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> companion(deg, deg);
                companion.setZero();
                companion.col(deg - 1) = -1 * coeff.topRows(deg);
                companion.block(1, 0, deg - 1, deg - 1).setIdentity();


                auto r = companion.eigenvalues();
                for (int kk = 0; kk < r.rows(); ++kk) {
                    double real = r[kk].real();
                    double imag = r[kk].imag();
                    if (std::abs(imag) < 1e-9 && real > 0 && r_d > real) {
                        r_d = real;
                        break;
                    }

                }
                double dd = undistortion_utils::undistortionDenominator<double>(r_d, lambdas.cast<double>());
                map_x(j, i) = r_img * ii * dd + cols / 2.0;
                map_y(j, i) = r_img * jj * dd + rows / 2.0;
            }
        }

    });
    cv::eigen2cv(map_x.cast<float>().eval(), mx);
    cv::eigen2cv(map_y.cast<float>().eval(), my);
    cv::remap(in, out, mx, my, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0, 0, 0));
    cv::imwrite(out_path, out);
    std::cout << "Undistortion've done\n";

    return 0;
}
