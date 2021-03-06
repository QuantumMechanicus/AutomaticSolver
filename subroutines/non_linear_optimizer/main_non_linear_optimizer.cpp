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


class ErrorFunctor {
private:
    Eigen::Vector3d left_point, right_point;
    double w, h;
    int nLambda;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h,
                 const int nLambda, const int nLambda2)
            : left_point(left_point),
              right_point(right_point), w(w), h(h), nLambda(nLambda) {
    }


    template<typename T>
    T undistortionDenominator(const T &r_distorted, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) const {

        return undistortion_utils::undistortionDenominator<T>(r_distorted, lambdas);
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

        bool is_invertible = undistortion_utils::checkUndistortionInvertibility<T>(lambdas);

        if (!is_invertible)
            return false;

        Matrix3T F;
        for (size_t k = 0; k < 3; ++k)
            for (size_t j = 0; j < 3; ++j) {
                if (k < 2 || j < 2)
                    F(k, j) = f_ptr[3 * j + k];
            }

        F(2, 2) = T(1);
        Eigen::JacobiSVD<Matrix3T> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vector3T singular_values = fmatrix_svd.singularValues();
        singular_values[2] = T(0.0);
        F = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
            fmatrix_svd.matrixV().transpose();
        F = F / F(2, 2);

        const T &d = std::max(hT, wT) / T(2.0);
        T dr = d / r;
        T alpha = undistortionDenominator<T>(dr, lambdas);

        if (alpha <= .0)
            return false;

        Vector2T u1, u2;
        bool is_correct = true;
        std::pair<T, T> res = undistortion_utils::computeEpipolarCurveDistanceError<T>(lambdas, F, left_point_T,
                                                                                       right_point_T,
                                                                                       r,
                                                                                       w, h,
                                                                                       u1, u2,
                                                                                       is_correct);


        residuals[0] = res.first;
        residuals[1] = res.second;
        return is_correct;
    }
};


int main(int argc, char *argv[]) {
    double w, h, prcnt_inl, r;
    std::string input1, input2;
    int number_of_distortion_coefficients;
    int number_of_pictures;
    std::string f_lambdas;
    std::string f_fundamental_matrices;
    std::string f_optimzer_results;
    int non_linear_iter;
    std::vector<std::string> left_cor_vec_names;
    std::vector<std::string> right_cor_vec_names;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("w", po::value<double>(&w)->required(), "Width")
                ("h", po::value<double>(&h)->required(), "Height")
                ("n_pic", po::value<int>(&number_of_pictures)->required(), "Number of pirctures")
                ("n_iters", po::value<int>(&non_linear_iter)->default_value(1), "Number of iterations")
                ("lambda_f", po::value<std::string>(&f_lambdas)->required(), "File with %n_pic estimated lambdas")
                ("fund_f", po::value<std::string>(&f_fundamental_matrices)->required(),
                 "File with %n_pic estimated fundamental matrices")
                ("left_corr_f", po::value<std::vector<std::string> >(&left_cor_vec_names)->multitoken()->required(),
                 "Path to left correspondeces")
                ("right_corr_f", po::value<std::vector<std::string> >(&right_cor_vec_names)->multitoken()->required(),
                 "Path to right correspondeces")
                ("q", po::value<double>(&prcnt_inl)->default_value(0.1), "quantile to minimize")
                ("results_f", po::value<std::string>(&f_optimzer_results)->default_value("./OptimizerResults/"),
                 "Path to results directory")
                ("nlambda", po::value<int>(&number_of_distortion_coefficients)->default_value(1),
                 "Number of parameters in denominator of model");

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

    std::fstream input_f1(f_fundamental_matrices, std::fstream::in);
    std::fstream input_f2(f_lambdas, std::fstream::in);
    Eigen::Matrix<double, Eigen::Dynamic, 1> Lambda(number_of_distortion_coefficients, 1);

    Lambda.setZero();
    StdVector<Eigen::Matrix3d> Fvec(number_of_pictures);


    r = std::sqrt(w * w + h * h) / 2.0;
    double minimal_residuals_per_points = std::numeric_limits<double>::max();
    for (size_t iters = 0; iters < non_linear_iter; ++iters) {

        if (iters == 0) {
            for (size_t kk = 0; kk < number_of_pictures; ++kk) {
                for (int ii = 0; ii < 3; ++ii)
                    for (int jj = 0; jj < 3; ++jj)
                        input_f1 >> Fvec[kk](ii, jj);
                double cur_l;
                input_f2 >> cur_l;
                Lambda[0] += cur_l;
            }
            Lambda[0] /= number_of_pictures;
        }

        ceres::Problem problem;
        double *lambda_ptr = Lambda.data();

        for (size_t ll = 0; ll < Lambda.size(); ++ll)
            std::cout << Lambda[ll] << " ";
        std::cout << std::endl;
        problem.AddParameterBlock(lambda_ptr, number_of_distortion_coefficients);


        for (size_t k = 0; k < Fvec.size(); ++k) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;
            Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                      fmatrix_svd.matrixV().transpose();
        }

        int residuals = 0;
        for (size_t kk = 0; kk < number_of_pictures; ++kk) {
            double *f_ptr = Fvec[kk].data();
            Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d, i1d, i2d;
            std::string name1 = left_cor_vec_names[kk];
            std::string name2 = right_cor_vec_names[kk];
            readPointsFromFile(name1, u1d);
            readPointsFromFile(name2, u2d);


            undistortion_utils::UndistortionProblemHelper helper(w, h, r, u1d, u2d, prcnt_inl);
            helper.setHypLambdas(Lambda);
            helper.setHypF(Fvec[kk]);
            helper.findInliers(f_optimzer_results + "non_linear_optimizer_info_for" + std::to_string(kk) + "_image_on" +
                               std::to_string(iters) + "_iteration");
            auto &inliers_ind = helper.getInliersIndices();

            i1d.resize(Eigen::NoChange, inliers_ind.size());
            i2d.resize(Eigen::NoChange, inliers_ind.size());
            for (size_t ll = 0; ll < inliers_ind.size(); ++ll) {
                i1d.col(ll) = u1d.col(inliers_ind[ll]);
                i2d.col(ll) = u2d.col(inliers_ind[ll]);
            }


            undistortion_utils::normalizePoints(i1d, w, h, r);
            undistortion_utils::normalizePoints(i2d, w, h, r);


            problem.AddParameterBlock(f_ptr, 8);

            for (size_t k = 0; k < i1d.cols(); ++k) {
                Eigen::Vector3d left, right;
                left.template block<2, 1>(0, 0) = i1d.col(k);
                right.template block<2, 1>(0, 0) = i2d.col(k);
                left[2] = right[2] = 1.0;

                auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                        new ErrorFunctor(left, right, w, h, number_of_distortion_coefficients, 0));
                fun->AddParameterBlock(number_of_distortion_coefficients);
                fun->AddParameterBlock(8);
                fun->SetNumResiduals(2);
                problem.AddResidualBlock(fun, /*new ceres::HuberLoss(15)*/ nullptr, lambda_ptr,
                                         f_ptr);
                /*const double *ptrs[] = {lambda_ptr, f_ptr};
                 double ptrs_res[2];
                 fun->Evaluate(ptrs, ptrs_res, nullptr);
                 std::cout << ptrs_res[0] << " " << ptrs_res[1] << " " << left.transpose() << " " << right.transpose()
                           << std::endl;*/
                ++residuals;
            }

        }
        std::cout << "Problem size: " << residuals << "points" << std::endl;


        ceres::Solver::Options options;
        //options.max_trust_region_radius = 0.01;
        options.max_num_iterations = 500;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.num_threads = 8;
        options.function_tolerance = 1e-16;
        options.parameter_tolerance = 1e-16;
        options.minimizer_progress_to_stdout = true;
        options.preconditioner_type = ceres::IDENTITY;
        options.jacobi_scaling = false;

        // Solve
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << std::endl;
        std::cout << Lambda.transpose() << " --- estimated lambdas" << std::endl;
        double curr_residual_per_point = std::sqrt(summary.final_cost / residuals);
        std::cout << "final residual " << curr_residual_per_point << " (per point)" << std::endl;
        if (curr_residual_per_point < minimal_residuals_per_points) {
            minimal_residuals_per_points = curr_residual_per_point;
        }
        std::fstream all_fund_matrices(f_optimzer_results + "on_" + std::to_string(iters) + "_iter_all. ff",
                                       std::fstream::out);
        for (size_t k = 0; k < Fvec.size(); ++k) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;

            Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                      fmatrix_svd.matrixV().transpose();
            std::cout << Fvec[k].determinant() << "\n" << Fvec[k] / Fvec[k](2, 2) << "\n\n";
            all_fund_matrices << Fvec[k] / Fvec[k](2, 2) << std::endl;
        }
    }
    std::cout << "minimal residual " << minimal_residuals_per_points << " (per point)" << std::endl;
    return 0;
}