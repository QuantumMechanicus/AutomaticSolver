#include "solver_ku8pt.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include "ceres/ceres.h"
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
                  double hyp_lambda, Eigen::Matrix3d &hyp_F, double quantile) {

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
    double alpha = 4 / (4 + hyp_lambda * (d / r) * (d / r));
    alpha = 1.0 / alpha;
    Eigen::Matrix3d shift;
    Eigen::Matrix3d scale;
    shift << 1, 0, -w / 2.0,
            0, 1, -h / 2.0,
            0, 0, 1;
    scale << 1.0 / (alpha * r), 0, 0,
            0, 1.0 / (alpha * r), 0,
            0, 0, 1;
    Eigen::Matrix3d recompute_F = shift.transpose() * scale.transpose() * hyp_F * scale * shift;
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
    auto u1 = r * alpha * u1d.cwiseProduct((ones + hyp_lambda * r1d).cwiseInverse());
    auto u2 = r * alpha * u2d.cwiseProduct((ones + hyp_lambda * r2d).cwiseInverse());

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

    std::fstream errf1("inlp1.txt", std::ios_base::out);
    std::fstream errf2("inlp2.txt", std::ios_base::out);
    for (size_t k = 0; k < u1d.cols(); ++k) {
        double c1 = l1.col(k).template topRows<2>().norm();
        double c2 = l2.col(k).template topRows<2>().norm();
        double err = std::abs(uu2.col(k).dot(l1.col(k)) / c1) + std::abs(uu1.col(k).dot(l2.col(k)) / c2);

        if (std::abs(err) < quantile * 5.36752) {
            ++goods;
            errf1 << uu1.col(k).transpose().leftCols(2) << "\n";
            errf2 << uu2.col(k).transpose().leftCols(2) << "\n";
        }

    }
    errf1.close();
    errf2.close();
    std::cout << goods << std::endl;
    return goods;
}

class ErrorFunctor {
private:
    Eigen::Vector3d left_point, right_point;
    double w, h;
public:
    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h) : left_point(left_point),
                                                                                          right_point(right_point), w(w), h(h) {}

    template<typename T>
    bool operator()(const T* lambda_ptr, const T* f_ptr, T* residuals) const
    {
        using Vector2T = Eigen::Matrix<T, 2, 1>;
        using Vector3T = Eigen::Matrix<T, 3, 1>;
        using Matrix3T = Eigen::Matrix<T, 3, 3>;
        Vector3T left_point_T = left_point.cast<T>();
        Vector3T right_point_T = right_point.cast<T>();
        Vector2T center;
        center[0] = T(w / 2.0);
        center[1] = T(h / 2.0);
        left_point_T.template block<2, 1>(0, 0) -= center;
        right_point_T.template block<2, 1>(0, 0) -= center;


        T lambda = *lambda_ptr; //(atan(*lambda_ptr) / T(M_PI) + T(0.5)) * T(5.0 / 4.0) - T(1.0);
        if (lambda < T(-1.0) || lambda > T(0.25))
            return false;
        Matrix3T F;
        for (size_t k = 0; k < 3; ++k)
            for (size_t j = 0; j < 3; ++j)
            {
                if (k < 2 || j < 2)
                {
                    F(k,j) = f_ptr[3*j + k];
                }
            }
        F(2,2) = T(1.0);
        T wT(w), hT(h);


        T r = ceres::sqrt((wT / 2.0) * (wT / 2.0) + (hT / 2.0) * (hT / 2.0));
        T d = std::max(hT, wT);
        T alpha = (T(4.0) + lambda * (d / r) * (d / r))/ T(4.0);


        Matrix3T scale;
        scale.setZero();
        scale(0, 0) = scale(1, 1) = T(1.0) / alpha / r;
        scale(2, 2) = T(1.0);

        Matrix3T recompute_F = scale.transpose() * F * scale;

        T r1d, r2d;
        r1d = (left_point_T.template block<2,1>(0,0)).squaredNorm();
        r2d = (right_point_T.template block<2,1>(0,0)).squaredNorm();



        Vector3T u1, u2;
        u1 = r * alpha * left_point_T/(T(1.0) + lambda * r1d);
        u2 = r * alpha * right_point_T/(T(1.0) + lambda * r2d);
        u1[2] = u2[2] = T(1.0);

        Vector3T l1 = recompute_F * u1;
        Vector3T l2 = recompute_F.transpose() * u2;

        T n1 = l1.template block<2, 1>(0, 0).squaredNorm();
        T n2 = l2.template block<2, 1>(0, 0).squaredNorm();
        T err = l1.dot(u2);
        T err2 = err * err;

        residuals[0] = err / sqrt(n1);
        residuals[1] = err / sqrt(n2);
        return true;
    }
};

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
#if 1
    Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
    readPointsFromFile(input1, u1d);
    readPointsFromFile(input2, u2d);

    Eigen::Matrix3d F;
    double Lambda;
    std::cout << u1d << std::endl;
    double q = getFundamentalMatrixAndLambda(u1d, u2d, w, h, F, Lambda, distr_f, iter, tr, tr2);
    std::cout << std::endl;
    std::cout << u1d << std::endl;
    std::fstream output_f("inl.txt", std::fstream::app);

    //goodPoints(u1d, u2d, w, h, Lambda, F, threshold);
    goodPoints(u1d, u2d, w, h, Lambda, F, q);
    F /= F(2, 2);
    output_f << "\n" << q << " - quantile \n" << Lambda << " - lambda\n\nF:\n\n" << F << std::endl;
    output_f.close();
    //double f_err = 0;

#else
    Eigen::Matrix3d F;
    double Lambda=-0.898384;
    F << -5.54865e-10,  5.38068e-08, -0.000100788,
    2.13784e-08 ,  1.4244e-07,   0.00231868,
   -3.09032e-05,  -0.00314068   ,         1;
#endif
    ceres::Problem problem;
    double *lambda_ptr = &Lambda;
    double *f_ptr = F.data();

    problem.AddParameterBlock(lambda_ptr, 1);
    problem.AddParameterBlock(f_ptr, 8);


    Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
    std::string name1 = "inlp1.txt";
    std::string name2 = "inlp2.txt";
    readPointsFromFile(name1, i1d);
    readPointsFromFile(name2, i2d);

    for (size_t k = 0; k < i1d.cols(); ++k) {
        Eigen::Vector3d left, right;
        left.template block<2,1>(0,0) = i1d.col(k);
        right.template block<2,1>(0,0) = i2d.col(k);

        problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ErrorFunctor, 2, 1, 8>(new ErrorFunctor(left, right, w, h)), /*new ceres::HuberLoss(15)*/ nullptr,lambda_ptr, f_ptr);
    }
    ceres::Solver::Options options;
    options.max_num_iterations = 500;
    options.linear_solver_type = ceres::DENSE_QR;
    options.function_tolerance = 1e-15;
    options.parameter_tolerance = 1e-15;
    options.minimizer_progress_to_stdout = true;

    // Solve
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;
    std::cout << F << std::endl;
    std::cout << Lambda << std::endl;
    return 0;
}