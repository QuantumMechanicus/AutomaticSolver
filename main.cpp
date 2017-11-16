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
                         double w, double h,
                         double hyp_lambda, Eigen::Matrix3d &hyp_F, double threshold) {

    double r = std::sqrt((w / 2.0) * (w / 2.0) + (h / 2.0) * (h / 2.0));
    if (r == 0) {
        r = 1;
    }
    u1d.row(0) = u1d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * w / 2.0;

    u1d.row(1) = u1d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * h / 2.0;

    u2d.row(0) = u2d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * w / 2.0;

    u2d.row(1) = u2d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * h / 2.0;
    u1d = u1d / r;
    u2d = u2d / r;
    double d = std::max(h, w);
    double alpha = 4/(4+hyp_lambda*(d/r)*(d/r));
    alpha = 1.0/alpha;
    Eigen::Matrix3d shift;
    Eigen::Matrix3d scale;
    shift << 1, 0, -w / 2.0,
            0, 1, -h / 2.0,
            0, 0, 1;
    scale << 1.0/(alpha*r), 0, 0,
            0, 1.0/(alpha*r), 0,
            0, 0, 1;
    Eigen::Matrix3d recompute_F = shift.transpose()*scale.transpose()*hyp_F*scale*shift;
    hyp_F = recompute_F;
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

    size_t goods = 0;
    auto u1 = r*alpha*u1d.cwiseProduct((ones + hyp_lambda * r1d).cwiseInverse());
    auto u2 = r*alpha*u2d.cwiseProduct((ones + hyp_lambda * r2d).cwiseInverse());

    Eigen::Matrix<double, 3, Eigen::Dynamic> uu1, uu2;
    uu1.resize(Eigen::NoChange, u1.cols());
    uu2.resize(Eigen::NoChange, u2.cols());
    uu1.setOnes();
    uu2.setOnes();
    uu1.row(0) = u1.row(0) + Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * w / 2.0;
    uu2.row(0) = u2.row(0) + Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * w / 2.0;

    uu1.row(1) = u1.row(1) + Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * h / 2.0;
    uu2.row(1) = u2.row(1) + Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1, u1d.cols()) * h / 2.0;


    Eigen::Matrix<double, 3, Eigen::Dynamic> l1 = (recompute_F * uu1);
    Eigen::Matrix<double, 3, Eigen::Dynamic> l2 = (recompute_F.transpose() * uu2);

    std::fstream errf1("inlp1.txt", std::fstream::app);
    std::fstream errf2("inlp2.txt", std::fstream::app);
    for (size_t k = 0; k < u1d.cols(); ++k) {
        double c1 = l1.col(k).template topRows<2>().norm();
        double c2 = l2.col(k).template topRows<2>().norm();
        double err = std::abs(uu2.col(k).dot(l1.col(k)) / c1) + std::abs(uu1.col(k).dot(l2.col(k)) / c2);

        if (std::abs(err) < threshold) {
            ++goods;
            errf1 << uu1.col(k).transpose() << "\n";
            errf2 << uu2.col(k).transpose() << "\n";
        }

    }
    errf1.close();
    errf2.close();
    std::cout << goods << std::endl;
    return goods;
}
int main(int argc, char *argv[]) {

    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2, distr_f;
    int iter = 10000;
    double tr = 0.25;
    double tr2 = -1;
    double threshold = 1;
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
                ("threshold", po::value<double>(&threshold),
                 "RANSAC threshold")
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
    std::cout << u1d << std::endl;
    size_t inl = getFundamentalMatrixAndLambda(u1d, u2d, w, h, F, Lambda, distr_f, iter, tr, tr2, threshold);
    std::cout << std::endl;

    std::cout << u1d << std::endl;
    std::fstream output_f("inl.txt", std::fstream::app);

    goodPoints(u1d, u2d, w, h, Lambda, F, threshold);
    F /= F(2,2);
    output_f << "\n" << inl << " - number of inliesrs \n" << Lambda << " - lambda\n\nF:\n\n" << F << std::endl;
    output_f.close();
    //double f_err = 0;
    //goodPoints(u1d, u2d, Lambda, F, 5, f_err, w, h);*/
    return 0;
}