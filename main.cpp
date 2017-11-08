#include "solver_ku8pt.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;
void readPointsFromFile(std::string &name, Eigen::Matrix<double, 2, Eigen::Dynamic> &m)
{
    std::fstream f(name, std::fstream::in);
    std::vector<std::pair<double, double>> v;
    while (!f.eof())
    {
        double x, y;
        f >> x;
        f >> y;
        if (f.eof())
            break;
        v.emplace_back(std::make_pair(x, y));
    }
    m.resize(Eigen::NoChange, v.size());
    for (size_t k = 0; k < v.size(); ++k)
    {
        m(0, k) = v[k].first;
        m(1, k) = v[k].second;
    }
    f.close();
}
int main(int argc, char *argv[]) {

    bool doublePrecision;
    int iter;
    std::string output, plyOutput, input1, input2, poses;

    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1)->required(), "Input filename1 (matches and camera parameters)")
                ("input2", po::value<std::string>(&input2)->required(), "Input filename2 (matches and camera parameters)")
                ;

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
    u1d.row(0) = u1d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1,u1d.cols())*7360.0/2.0;

    u1d.row(1) = u1d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1,u1d.cols())*4912.0/2.0;

    u2d.row(0) = u2d.row(0) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1,u1d.cols())*7360.0/2.0;

    u2d.row(1) = u2d.row(1) - Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(1,u1d.cols())*4912.0/2.0;
    //u1d = u1d/1000;
    //u2d = u2d/1000;
    Eigen::Matrix3d F;
    double Lambda;

    size_t inl = getFundamentalMatrixAndLambda(u1d, u2d, F, Lambda, 900, 4);



    return 0;
}