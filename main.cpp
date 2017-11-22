#include "solver_ku8pt.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include "ceres/ceres.h"
#include <iostream>

namespace po = boost::program_options;

template<typename T>
using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;

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
                  double hyp_lambda, Eigen::Matrix3d &hyp_F, double quantile, std::string out_name) {

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

    std::fstream errf1(out_name + "1", std::ios_base::out);
    std::fstream errf2(out_name + "2", std::ios_base::out);
    double full_err = 0;
    for (size_t k = 0; k < u1d.cols(); ++k) {
        double c1 = l1.col(k).template topRows<2>().norm();
        double c2 = l2.col(k).template topRows<2>().norm();
        double err1 = std::abs(uu2.col(k).dot(l1.col(k)) / c1);
        double err2 = std::abs(uu1.col(k).dot(l2.col(k)) / c2);
        double err = err1 + err2;

        if (std::abs(err) < quantile * 5.36752) {
            full_err += (err1*err1 + err2*err2);
            ++goods;
            errf1 << u1d.col(k).transpose().leftCols(2) << "\n";
            errf2 << u2d.col(k).transpose().leftCols(2) << "\n";
        }

    }
    std::cout << recompute_F << std::endl;
    std::cout << "Sq err: " << full_err << std::endl;
    errf1.close();
    errf2.close();
    std::cout << goods << std::endl;
    return goods;
}


static int invocations = 0;
bool firstInvocation = true;

template<typename T>
void foo(const Eigen::Matrix<T, 3, 3> &m) {
    if (++invocations <= 3) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << m(i, j).a << " ";
            }
            std::cout << std::endl;
        }
    }
}

template<>
void foo<double>(const Eigen::Matrix3d &n) {
    if (++invocations <= 3) {
        std::cout << std::endl << n << std::endl;
    }
}

template<typename T>
double double_cast(const T& t) { return t.a; }

template<>
double double_cast<double>(const double &d) { return d; }


class ErrorFunctor {
private:
    Eigen::Vector3d left_point, right_point;
    double w, h;
    int nLambda;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h, const int nLambda)
            : left_point(left_point),
              right_point(right_point), w(w), h(h), nLambda(nLambda) {
        std::cout << "Residual@ L:" << left_point.transpose() << " " << right_point.transpose() << std::endl;
    }

    template<typename T>
    bool operator()(T const * const * parameters, T *residuals) const {
        const T *lambda_ptr = parameters[0];
        const T *f_ptr = parameters[1];

        using Vector2T = Eigen::Matrix<T, 2, 1>;
        using Vector3T = Eigen::Matrix<T, 3, 1>;
        using Matrix3T = Eigen::Matrix<T, 3, 3>;
        using VectorNL = Eigen::Matrix<T, Eigen::Dynamic, 1>;
        Vector3T left_point_T = left_point.cast<T>();
        Vector3T right_point_T = right_point.cast<T>();
        Vector2T center;
        center[0] = T(w / 2.0);
        center[1] = T(h / 2.0);
        T wT(w), hT(h);
        T r = ceres::sqrt((wT / T(2.0)) * (wT / T(2.0)) + (hT / T(2.0)) * (hT / T(2.0)));

        VectorNL lambdas = Eigen::Map<const VectorNL>(lambda_ptr, nLambda);
        Eigen::VectorXd lambdas_d;
        double r_d = double_cast(r);
        lambdas_d.resize(nLambda, 1);
        /*for (int i = 0; i < nLambda; ++i) {
            lambdas_d[i] = double_cast(lambdas[i]);
            if (i && std::abs(lambdas_d[i]) > 0.1)
                return false;
        }*/




        T lambda = lambdas[0]; //(atan(*lambda_ptr) / T(M_PI) + T(0.5)) * T(5.0 / 4.0) - T(1.0);
        //if (lambda < T(-1.0) || lambda > T(0.25))
       //     return false;
        Matrix3T F;
        for (size_t k = 0; k < 3; ++k)
            for (size_t j = 0; j < 3; ++j) {
                if (k < 2 || j < 2) {
                    F(k, j) = f_ptr[3 * j + k];
                }
            }
        F(2, 2) = T(1.0);

        const T &d = std::max(hT, wT);
        T dr2 = (d / (T(2.0)*r));
        dr2 *= dr2;
        T dr_pow = dr2;
        T nominator (1.0);
        for (int i = 0; i < nLambda; ++i) {
            nominator += dr_pow * lambdas[i];
            dr_pow *= dr2;
        }

        T alpha = nominator / T(1.0);
        if (nominator <= .0)
            return false;


        Matrix3T scale;
        Matrix3T shift;
        scale.setZero();
        shift.setZero();
        shift(0, 0) = shift(1, 1) = shift(2, 2) = T(1.0);
        shift(0, 2) = -wT / T(2.0);
        shift(1, 2) = -hT / T(2.0);
        scale(0, 0) = scale(1, 1) = T(1.0) / (alpha * r);
        scale(2, 2) = T(1.0);

        Matrix3T recompute_F =
                        scale.transpose()*F*scale
        ;
        foo<T>(recompute_F);

        T r1d, r2d;
        r1d = (left_point_T.template block<2, 1>(0, 0)).squaredNorm();
        r2d = (right_point_T.template block<2, 1>(0, 0)).squaredNorm();


        Vector3T u1, u2;
        T denominator1(1.0);
        T denominator2(1.0);
        T r1d_pow = r1d;
        T r2d_pow = r2d;
        for (int i = 0; i < nLambda; ++i) {
            denominator1 += lambdas[i] * r1d_pow;
            denominator2 += lambdas[i] * r2d_pow;
            r1d_pow *= r1d;
            r2d_pow *= r2d;
        }
        u1 = r * alpha * left_point_T / denominator1;
        u2 = r * alpha * right_point_T / denominator2;
        u1[2] = u2[2] = T(1.0);

        Vector3T l1 = recompute_F * u1;
        Vector3T l2 = recompute_F.transpose() * u2;

        T n1 = l1.template block<2, 1>(0, 0).squaredNorm();
        T n2 = l2.template block<2, 1>(0, 0).squaredNorm();
        T err = l1.dot(u2);
        T err2 = err * err;

        residuals[0] = err / sqrt(n1);
        residuals[1] = err / sqrt(n2);
        return denominator1 > 0.0 && denominator2 > 0.0;
    }
};

int main(int argc, char *argv[]) {
    std::cout << argv[0] << std::endl;
    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2, distr_f;
    int iter = 10000;
    double tr = 0.25;
    double tr2 = -1;
    double threshold = 1;
    bool non_linear_opt = 0;
    int nLambda;
    std::string f_l;
    std::string f_f;
    std::vector<std::string> vec_names;
    std::vector<std::string> vec_names2;
    std::string inliers_f;
    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("input1", po::value<std::string>(&input1),
                 "Input filename1 (matches and camera parameters)")
                ("input2", po::value<std::string>(&input2),
                 "Input filename2 (matches and camera parameters)")
                ("f", po::value<std::string>(&distr_f),
                 "Output file for lambdas")
                ("threshold", po::value<double>(&threshold),
                 "RANSAC threshold")
                ("threshold1", po::value<double>(&tr), "Lambda threshold")
                ("threshold2", po::value<double>(&tr2), "Lambda threshold")
                ("iters", po::value<int>(&iter), "Number of iterations")
                ("w", po::value<double>(&w), "Width")
                ("h", po::value<double>(&h), "Height")
                ("optim", po::value<bool>(&non_linear_opt))
                ("inliers_f", po::value<std::string>(&inliers_f))
                ("f_f", po::value<std::string>(&f_f)->required())
                ("f_l", po::value<std::string>(&f_l)->required())
                ("f_names", po::value<std::vector<std::string> >(&vec_names)->multitoken(), "description")
                ("f_names2", po::value<std::vector<std::string> >(&vec_names2)->multitoken(), "description")
                ("nlambda", po::value<int>(&nLambda)->default_value(1));

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
    std::cout << vec_names.size() << " --- names:\n";
    for (size_t kk =0; kk < vec_names.size(); ++kk)
        std::cout << vec_names[kk] << std::endl;

    if (!non_linear_opt) {
        Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
        readPointsFromFile(input1, u1d);
        readPointsFromFile(input2, u2d);

        Eigen::Matrix3d F;
        double Lambda;
        std::cout << u1d << std::endl;
        double q = getFundamentalMatrixAndLambda(u1d, u2d, w, h, F, Lambda, distr_f, iter, tr, tr2);
        std::cout << std::endl;
        std::cout << u1d << std::endl;
        std::fstream output_f1(f_f, std::fstream::out);
        std::fstream output_f2(f_l, std::fstream::out);

        //goodPoints(u1d, u2d, w, h, Lambda, F, threshold);
        goodPoints(u1d, u2d, w, h, Lambda, F, q, inliers_f);
        F /= F(2, 2);
        output_f2 << Lambda << "\n";
        output_f1 << F << "\n";

        output_f1.close();
        output_f2.close();
        //double f_err = 0;
    } else {

        int n_f = 9;
        std::cout << f_f << " --- f_f\n";
        std::cout << f_l << " --- f_l\n";

        std::fstream input_f1(f_f, std::fstream::in);
        std::fstream input_f2(f_l, std::fstream::in);
        Eigen::Matrix<double, Eigen::Dynamic, 1> Lambda(nLambda, 1);
        Lambda.setZero();
        StdVector<Eigen::Matrix3d> Fvec(n_f);
        for (size_t kk = 0; kk < n_f; ++kk)
        {
            for (int ii = 0; ii < 3; ++ii)
                for (int jj = 0; jj < 3; ++jj)
                    input_f1 >> Fvec[kk](ii, jj);
            double foo;
            input_f2  >> foo;
            Lambda[0] += foo;
        }
        Lambda[0] /= n_f;


        std::cout << Lambda.transpose() << " --- lambda" << std::endl;
        //F.transposeInPlace();


        ceres::Problem problem;
        double *lambda_ptr = Lambda.data();


        problem.AddParameterBlock(lambda_ptr, nLambda);


        int residuals = 0;
        for (size_t kk = 0; kk < n_f; ++kk) {
            double *f_ptr = Fvec[kk].data();
            Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
            std::string name1 = vec_names[kk];
            std::string name2 = vec_names2[kk];
            readPointsFromFile(name1, i1d);
            readPointsFromFile(name2, i2d);

            problem.AddParameterBlock(f_ptr, 8);

            for (size_t k = 0; k < i1d.cols(); ++k) {
                Eigen::Vector3d left, right;
                left.template block<2, 1>(0, 0) = i1d.col(k);
                right.template block<2, 1>(0, 0) = i2d.col(k);
                left[2] = right[2] = 1.0;

                auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                        new ErrorFunctor(left, right, w, h, nLambda));
                fun->AddParameterBlock(nLambda);
                fun->AddParameterBlock(8);
                fun->SetNumResiduals(2);
                problem.AddResidualBlock(fun, /*new ceres::HuberLoss(15)*/ nullptr, lambda_ptr,
                                         f_ptr);
                ++residuals;
            }

        }
        std::cout << "Problem size: " << residuals << "points" << std::endl;


        ceres::Solver::Options options;
        //options.max_trust_region_radius = 0.01;
        options.max_num_iterations = 500;
        options.linear_solver_type = ceres::DENSE_QR;
        options.function_tolerance = 1e-15;
        options.parameter_tolerance = 1e-15;
        options.minimizer_progress_to_stdout = true;
        //options.preconditioner_type = ceres::IDENTITY;
        //options.jacobi_scaling = false;

#if 0
        const double * params[] = {lambda_ptr, f_ptr};
        double residual = 0.0;
        for (auto& adcf: adcfs) {
            double currRes = 0.0;
            adcf->Evaluate(params, &currRes, nullptr);
            residual += currRes * currRes;
        }
        std::cout << "Residual @ initial point: " << residual << std::endl;
#endif
        // Solve
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << std::endl;
        //std::cout << F << std::endl;
        std::cout << Lambda.transpose() << std::endl;
        std::cout << "final residual " << std::sqrt(summary.final_cost / residuals) << " (per point)" << std::endl;
    }
    return 0;
}