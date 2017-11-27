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

double undistortionDenominator(const double &r_distorted2, const Eigen::Matrix<double, Eigen::Dynamic, 1> &lambdas) {
    double denominator(1.0);
    double r_distorted2_pow = r_distorted2;
    for (int i = 0; i < lambdas.cols(); ++i) {
        denominator += lambdas[i] * r_distorted2_pow;
        r_distorted2_pow *= r_distorted2;
    }
    return denominator;
}

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

    Eigen::MatrixXd map_x(rows, cols), map_y(rows, cols);

    double d = std::max(rows, cols) / 2.0;
    double dr = d / r_img;
    double dr2 = dr * dr;
    Eigen::VectorXd lambdas = Eigen::Map<Eigen::VectorXd>(v_lambdas.data(), v_lambdas.size());
    std::cout << lambdas.transpose() << " " << cols << " " << rows << std::endl;
    double alpha = undistortionDenominator(dr2, lambdas);
    size_t n_lambda = v_lambdas.size();

#ifdef PARALLEL
    tbb::parallel_for(tbb::blocked_range<int>(0, cols), [&](auto range) {
        for (int i = range.begin(); i != range.end(); ++i) {
#else
    for (int i = 0; i < cols; ++i) {
#endif
        for (int j = 0; j < rows; ++j) {
            double ii = (i - cols / 2.0) / r_img / alpha;
            double jj = (j - rows / 2.0) / r_img / alpha;
            double r_u = std::sqrt(ii * ii + jj * jj + 1e-7);

            Eigen::VectorXd coeff = Eigen::VectorXd::Zero(2 * n_lambda + 1);

            for (size_t k = 2 * n_lambda; k >= 2; k -= 2)
                coeff(k, 0) = r_u * lambdas[k / 2 - 1];

            coeff(1, 0) = -1;
            coeff(0, 0) = r_u;
           // std::cout << coeff.transpose() << i << " " << j << std::endl;

            double r_d = std::numeric_limits<double>::max();
            size_t deg = 0;

            Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
            for (size_t kk = 2 * n_lambda; kk >= 0; --kk) {
                if (coeff[kk] != 0) {
                    deg = kk;
                    break;
                }
            }

            if (!Eigen::isfinite(coeff.array()).all() or std::abs(coeff[deg]) < 1e-9)
                std::cout << "A " << coeff;
            solver.compute(coeff.block(0, 0, deg + 1, 1).eval());
            Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType r = solver.roots();
            for (int iii = 0; iii < r.rows(); ++iii) {
                double real = r[iii].real();
                double imag = r[iii].imag();
                if (std::abs(imag) < 1e-9 && real > 0 && r_d > real) {
                    r_d = real;
                }

            }
            double dd = undistortionDenominator(r_d * r_d, lambdas);
            map_x(j, i) = r_img * ii * dd + cols / 2.0;
            map_y(j, i) = r_img * jj * dd + rows / 2.0;
        }
    }
#ifdef PARALLEL
    });
#endif

    cv::eigen2cv(map_x.cast<float>().eval(), mx);
    cv::eigen2cv(map_y.cast<float>().eval(), my);
    cv::remap(in, out, mx, my, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0, 0, 0));
    cv::imwrite(out_path, out);

#if 0
    f1 << std::setprecision(16) << std::showpos << map_x;
    f2 << std::setprecision(16) << std::showpos << map_y;
#endif
    return 0;
}
