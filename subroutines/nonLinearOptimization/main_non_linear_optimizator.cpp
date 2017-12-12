#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include "ceres/ceres.h"
#include <iostream>
#include "undistortion_problem_utils.h"

namespace po = boost::program_options;

template<typename T>
using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;

void readPointsFromFile(const std::string &name, Eigen::Matrix<double, 2, Eigen::Dynamic> &m) {
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
    int nLambda, nLambda2;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h,
                 const int nLambda, const int nLambda2)
            : left_point(left_point),
              right_point(right_point), w(w), h(h), nLambda(nLambda), nLambda2(nLambda2) {
        //std::cout << "Residual@ L:" << left_point.transpose() << " " << right_point.transpose() << std::endl;
    }


    template<typename T>
    T undistortionFraction(const T &r_distorted2, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) const {

        return undistortion_utils::undistortionDenominator<T>(r_distorted2, lambdas.topRows(nLambda));
    }

    template<typename T>
    bool operator()(T const *const *parameters, T *residuals) const {
        const T *lambda_ptr = parameters[0];
        const T *f_ptr = parameters[1];

        using Vector2T = Eigen::Matrix<T, 2, 1>;
        using Vector3T = Eigen::Matrix<T, 3, 1>;
        using Matrix3T = Eigen::Matrix<T, 3, 3>;
        using VectorNL = Eigen::Matrix<T, Eigen::Dynamic, 1>;
        Vector2T left_point_T = left_point.topRows(2).cast<T>();
        Vector2T right_point_T = right_point.topRows(2).cast<T>();
        Vector2T center;
        center[0] = T(w / 2.0);
        center[1] = T(h / 2.0);
        T wT(w), hT(h);
        T r = ceres::sqrt((wT / T(2.0)) * (wT / T(2.0)) + (hT / T(2.0)) * (hT / T(2.0)));

        VectorNL lambdas = Eigen::Map<const VectorNL>(lambda_ptr, nLambda);
        Eigen::VectorXd lambdas_d;
        double r_d = double_cast(r);
        lambdas_d.resize(nLambda, 1);

        Matrix3T F; //=  Eigen::Map<const Matrix3T >(f_ptr);
        for (size_t k = 0; k < 3; ++k)
            for (size_t j = 0; j < 3; ++j) {
                if (k < 2 || j < 2)
                    F(k, j) = f_ptr[3 * j + k];
            }

        F(2, 2) = T(1);
        Eigen::JacobiSVD<Matrix3T> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vector3T singular_values = fmatrix_svd.singularValues();
        //std::cout << double_cast(singular_values[0]) << " " <<  double_cast(singular_values[1]) <<  " "<< double_cast(singular_values[2]) << "\n";
        singular_values[2] = T(0.0);

        F = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
            fmatrix_svd.matrixV().transpose();

        F = F / F(2, 2);
        //std::cout <<  double_cast(F.determinant()) << " det\n";

        const T &d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T dr2 = dr * dr;
        T alpha = undistortionFraction(dr2, lambdas);


        if (alpha <= .0)
            return false;

        Vector2T u1, u2;
        bool is_correct;
        auto res = undistortion_utils::computeError<T>(lambdas, F, left_point_T, right_point_T, r, w, h, u1, u2,
                                                       is_correct);

        residuals[0] = res.first;
        residuals[1] = res.second;
        return is_correct;
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

    int non_linear_iter = 2;
    std::vector<std::string> vec_names, left_cor_vec_names;
    std::vector<std::string> vec_names2, right_cor_vec_names;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("w", po::value<double>(&w)->required(), "Width")
                ("h", po::value<double>(&h)->required(), "Height")
                ("n_pic", po::value<int>(&n_f)->required(), "Number of pirctures")
                ("n_iters", po::value<int>(&non_linear_iter), "Number of iterations")
                ("lambda_f", po::value<std::string>(&f_l)->required(), "File with %n_pic estimated lambdas")
                ("fund_f", po::value<std::string>(&f_f)->required(), "File with %n_pic estimated fundamental matrices")
                ("left_inl_f", po::value<std::vector<std::string> >(&vec_names)->multitoken()->required(),
                 "Path to left inliers")
                ("right_inl_f", po::value<std::vector<std::string> >(&vec_names2)->multitoken()->required(),
                 "Path to right inliers")
                ("left_corr_f", po::value<std::vector<std::string> >(&left_cor_vec_names)->multitoken()->required(),
                 "Path to left correspondeces")
                ("right_corr_f", po::value<std::vector<std::string> >(&right_cor_vec_names)->multitoken()->required(),
                 "Path to right correspondeces")
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

    int nLambda2 = 0;
    Eigen::Matrix<double, Eigen::Dynamic, 1> Lambda(nLambda + nLambda2, 1);
    StdVector<Eigen::Matrix3d> Fvec(n_f);
    std::vector<double> lmbds(n_f);
    r = std::sqrt(w * w + h * h) / 2.0;
    for (size_t iters = 0; iters < non_linear_iter; ++iters) {

        if (iters == 0) {

            std::fstream input_f1(f_f, std::fstream::in);
            std::fstream input_f2(f_l, std::fstream::in);
            Lambda.setZero();
            for (size_t kk = 0; kk < n_f; ++kk) {
                for (int ii = 0; ii < 3; ++ii)
                    for (int jj = 0; jj < 3; ++jj)
                        input_f1 >> Fvec[kk](ii, jj);
                double cur_l;
                input_f2 >> cur_l;
                Lambda[0] += cur_l;
                lmbds[kk] = cur_l;
            }
            Lambda[0] /= n_f;
        }


        std::cout << "Start lambdas is: " << Lambda.transpose() << std::endl;
        ceres::Problem problem;
        double *lambda_ptr = Lambda.data();
        std::cout << "Start Fs:\n";

        problem.AddParameterBlock(lambda_ptr, nLambda + nLambda2);
        for (size_t k = 0; k < Fvec.size(); ++k) {
            if (abs(Fvec[k].determinant()) > 1e-5) {
                Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
                singular_values[2] = 0.0;

                Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                          fmatrix_svd.matrixV().transpose();
                Fvec[k] = Fvec[k] / Fvec[k](2, 2);
                std::cout << Fvec[k] << " \n";
            }
        }

        std::cout << "Find inliers...\n";


        int residuals = 0;
        for (size_t kk = 0; kk < n_f; ++kk) {
            double *f_ptr = Fvec[kk].data();
            Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
            std::string name1;
            std::string name2;
            if (iters < 1) {
                name1 = vec_names[kk];
                name2 = vec_names2[kk];
            } else {
                name1 = std::to_string(iters) + std::to_string(kk + 1) + "_found_inliers_left";
                name2 = std::to_string(iters) + std::to_string(kk + 1) + "_found_inliers_right";
            }
            readPointsFromFile(name1, i1d);
            readPointsFromFile(name2, i2d);

            normalizePoints(i1d, w, h, r);
            normalizePoints(i2d, w, h, r);
            std::string out_name = std::to_string(iters + 1) + std::to_string(kk + 1) + "_found_inliers";
            int inl;
            if (iters == 0) {
                Lambda[0] = lmbds[kk];
            }

            //inl = findInliers(w, h, i1d, i2d, Lambda, Fvec[kk], out_name);
            //readPointsFromFile(out_name+"_left", i1d);
            //readPointsFromFile(out_name+"_right", i2d);
            //normalizePoints(i1d, w, h, r);
            //normalizePoints(i2d, w, h, r);

            problem.AddParameterBlock(f_ptr, 8);

            for (size_t k = 0; k < i1d.cols(); ++k) {
                Eigen::Vector3d left, right;
                left.template block<2, 1>(0, 0) = i1d.col(k);
                right.template block<2, 1>(0, 0) = i2d.col(k);
                left[2] = right[2] = 1.0;

                auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                        new ErrorFunctor(left, right, w, h, nLambda, nLambda2));
                fun->AddParameterBlock(nLambda + nLambda2);
                fun->AddParameterBlock(8);
                fun->SetNumResiduals(2);
                problem.AddResidualBlock(fun, nullptr, lambda_ptr,
                                         f_ptr);

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
        std::cout << Lambda.transpose() << std::endl;
        std::cout << "final residual " << std::sqrt(summary.final_cost / residuals) << " (per point)" << std::endl;
        for (size_t k = 0; k < Fvec.size(); ++k) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;

            Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                      fmatrix_svd.matrixV().transpose();
            Fvec[k] = Fvec[k] / Fvec[k](2, 2);
        }
        std::vector<double> residualsf;
        problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &residualsf, NULL, NULL);
        std::fstream fm("res", std::fstream::out);
        std::cout << "SZ: " << residualsf.size() << std::endl;
        for (size_t mm = 0; mm < residualsf.size(); ++mm)
            fm << residualsf[mm] << "\n";


    }
    std::fstream fm("fundamental_matrices", std::fstream::out);
    for (size_t k = 0; k < Fvec.size(); ++k) {
        Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
        singular_values[2] = 0.0;

        Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                  fmatrix_svd.matrixV().transpose();
        Fvec[k] = Fvec[k] / Fvec[k](2, 2);
        Eigen::Matrix3d scale;
        Eigen::Matrix3d shift;
        scale.setZero();
        shift.setZero();
        shift(0, 0) = shift(1, 1) = shift(2, 2) = (1.0);
        shift(0, 2) = -w / (2.0);
        shift(1, 2) = -h / (2.0);

        double d = std::max(h, w) / (2.0);
        double dr = d / r;
        double dr2 = dr * dr;
        double alpha = undistortion_utils::undistortionDenominator(dr2, Lambda);
        scale(0, 0) = scale(1, 1) = (1.0) / (alpha * r);
        scale(2, 2) = (1.0);
        fm << shift.transpose() * scale.transpose() * Fvec[k] * scale * shift << "\n";
    }
    /*for (size_t kk = 0; kk < n_f; ++kk) {
        double *f_ptr = Fvec[kk].data();
        Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
        std::string name1 = vec_names[kk];
        std::string name2 = vec_names2[kk];
        readPointsFromFile(name1, i1d);
        readPointsFromFile(name2, i2d);

        normalizePoints(i1d, w, h, r);
        normalizePoints(i2d, w, h, r);
        std::cout << i1d << std::endl;
        for (size_t k = 0; k < i1d.cols(); ++k) {
            Eigen::Vector3d left, right;
            left.template block<2, 1>(0, 0) = i1d.col(k);
            right.template block<2, 1>(0, 0) = i2d.col(k);
            left[2] = right[2] = 1.0;

            auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                    new ErrorFunctor(left, right, w, h, nLambda));
            fun->AddParameterBlock(nLambda);
            fun->AddParameterBlock(9);
            fun->SetNumResiduals(2);

            const double* ptrs[] = {lambda_ptr, f_ptr};
            double ptrs_res[2];
            fun->Evaluate(ptrs, ptrs_res, nullptr);
            errr = errr + ptrs_res[0]*ptrs_res[0] +  ptrs_res[1]*ptrs_res[1];
            //std::cout << ptrs_res[0] << " " << ptrs_res[1]  << std::endl;

        }
    }
    std::cout << "Check res: " << errr << std::endl;*/
    return 0;
}