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
        v.emplace_back(x, y);
    }
    m.resize(Eigen::NoChange, v.size());
    for (size_t k = 0; k < v.size(); ++k) {
        m(0, k) = v[k].first;
        m(1, k) = v[k].second;
    }
    f.close();
}


template<typename T>
double double_cast(const T &t) { return t.a; }

template<>
double double_cast<double>(const double &t) { return t; }


class ErrorFunctor {
private:
    Eigen::Vector3d left_point, right_point;
    double w, h;
    int nLambda;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h,
                 const int nLambda)
            : left_point(left_point),
              right_point(right_point), w(w), h(h), nLambda(nLambda) {
        std::cout << "Residual@ L:" << left_point.transpose() << " " << right_point.transpose() << std::endl;
    }

    template<typename T>
    T undistortionDenominator(const T &r_distorted2, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) const {
        T denominator(1.0);
        T r_distorted2_pow = r_distorted2;
        for (int i = 0; i < nLambda; ++i) {
            denominator += lambdas[i] * r_distorted2_pow;
            r_distorted2_pow *= r_distorted2;
        }
        return denominator;
    }

    template<typename T>
    bool operator()(T const *const *parameters, T *residuals) const {
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

        Matrix3T F;
        for (size_t k = 0; k < 3; ++k)
            for (size_t j = 0; j < 3; ++j) {
                    F(k, j) = f_ptr[3 * j + k];
                }
            }

        /*Eigen::JacobiSVD<Matrix3T> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vector3T singularValues = svd.singularValues();
        Matrix3T S = Matrix3T::Zero();
        singularValues /= singularValues(0);
        singularValues(2) = T(0);
        S.diagonal() = singularValues;
        F = svd.matrixU()*S*svd.matrixV().conjugate();*/

        const T &d = std::min(hT, wT) / T(2.0);
        T dr = d / r;
        T dr2 = dr * dr;
        T alpha = undistortionDenominator(dr2, lambdas);


        if (alpha <= .0)
            return false;


        /*Matrix3T scale;
        Matrix3T shift;
        scale.setZero();
        shift.setZero();
        shift(0, 0) = shift(1, 1) = shift(2, 2) = T(1.0);
        shift(0, 2) = -wT / T(2.0);
        shift(1, 2) = -hT / T(2.0);
        scale(0, 0) = scale(1, 1) = T(1.0) / (alpha * r);
        scale(2, 2) = T(1.0);*/

        Matrix3T recompute_F = F;


        T r1d, r2d;
        r1d = (left_point_T.template block<2, 1>(0, 0)).squaredNorm();
        r2d = (right_point_T.template block<2, 1>(0, 0)).squaredNorm();


        Vector3T u1, u2;
        T denominator1 = undistortionDenominator(r1d, lambdas);
        T denominator2 = undistortionDenominator(r2d, lambdas);
        u1 =left_point_T / denominator1;
        u2 = right_point_T / denominator2;
        u1[2] = u2[2] = T(1.0);

        Vector3T l1 = recompute_F * u1;
        Vector3T l2 = recompute_F.transpose() * u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(u2);
        T err2 = err * err;

        residuals[0] = alpha*r*err / n1;
        residuals[1] = alpha*r*err / n2;
        return denominator1 > 0.0 && denominator2 > 0.0;
    }
};

void normalizePoints(auto &points, double w, double h, double r) {
    points.row(0).array() -= w / 2.0;
    points.row(1).array() -= h / 2.0;
    points /= r;
}

int main(int argc, char *argv[]) {
    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2;
    int nLambda;
    int nPictures;
    int n_f;
    std::string f_l;
    std::string f_f;
    std::vector<std::string> vec_names;
    std::vector<std::string> vec_names2;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("w", po::value<double>(&w)->required(), "Width")
                ("h", po::value<double>(&h)->required(), "Height")
                ("n_pic", po::value<int>(&n_f)->required(), "Number of pirctures")
                ("lambda_f", po::value<std::string>(&f_l)->required(), "File with %n_pic estimated lambdas")
                ("fund_f", po::value<std::string>(&f_f)->required(), "File with %n_pic estimated fundamental matrices")
                ("left_inl_f", po::value<std::vector<std::string> >(&vec_names)->multitoken()->required(),
                 "Path to left inliers")
                ("right_inl_f", po::value<std::vector<std::string> >(&vec_names2)->multitoken()->required(),
                 "Path to right inliers")
                ("nlambda", po::value<int>(&nLambda)->default_value(1), "Number of parameters in model");

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


    std::fstream input_f1(f_f, std::fstream::in);
    std::fstream input_f2(f_l, std::fstream::in);
    Eigen::Matrix<double, Eigen::Dynamic, 1> Lambda(nLambda, 1);
    Lambda.setZero();
    StdVector<Eigen::Matrix3d> Fvec(n_f);
    for (size_t kk = 0; kk < n_f; ++kk) {
        for (int ii = 0; ii < 3; ++ii)
            for (int jj = 0; jj < 3; ++jj)
                input_f1 >> Fvec[kk](ii, jj);
        double cur_l;
        input_f2 >> cur_l;
        Lambda[0] += cur_l;
    }
    Lambda[0] /= n_f;
    r = std::sqrt(w * w + h * h) / 2.0;


    ceres::Problem problem;
    double *lambda_ptr = Lambda.data();


    problem.AddParameterBlock(lambda_ptr, nLambda);
    for (size_t k = 0; k < Fvec.size(); ++k) {
        //std::cout << Fvec[k] << std::endl;
        Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
        singular_values[2] = 0.0;

        Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                              fmatrix_svd.matrixV().transpose();
        //std::cout << Fvec[k] << std::endl;
        //std::cout << Fvec[k].determinant() << "\n"<< Fvec[k] << "\n\n";
    }

    int residuals = 0;
    for (size_t kk = 0; kk < n_f; ++kk) {
        double *f_ptr = Fvec[kk].data();
        Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
        std::string name1 = vec_names[kk];
        std::string name2 = vec_names2[kk];
        readPointsFromFile(name1, i1d);
        readPointsFromFile(name2, i2d);

        normalizePoints(i1d, w, h, r);
        normalizePoints(i2d, w, h, r);
        std::cout << i1d << std::endl;
        std::cout << Lambda[0] << std::endl;
        problem.AddParameterBlock(f_ptr, 9);

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
            const double* ptrs[] = {lambda_ptr, f_ptr};
            double ptrs_res[2];
            fun->Evaluate(ptrs, ptrs_res, nullptr);
            std::cout << ptrs_res[0] << " " << ptrs_res[1] << " " << left.transpose() << " " << right.transpose() << std::endl;
            ++residuals;
        }

    }
    std::cout << "Problem size: " << residuals << "points" << std::endl;


    ceres::Solver::Options options;
    //options.max_trust_region_radius = 0.01;
    options.max_num_iterations = 500;
    options.linear_solver_type = ceres::DENSE_QR;
    options.function_tolerance = 1e-16;
    options.parameter_tolerance = 1e-16;
    options.minimizer_progress_to_stdout = true;
    options.preconditioner_type = ceres::IDENTITY;
    options.jacobi_scaling = false;

    // Solve
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;
    //std::cout << F << std::endl;
    std::cout << Lambda.transpose() << std::endl;
    std::cout << "final residual " << std::sqrt(summary.final_cost / residuals) << " (per point)" << std::endl;
    for (size_t k = 0; k < Fvec.size(); ++k) {
        std::cout << Fvec[k].determinant() << "\n"<< Fvec[k] << "\n\n";
    }

    return 0;
}