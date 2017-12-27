#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include "ceres/ceres.h"
#include <iostream>
#include "undistortion_problem_utils.h"
#include <boost/math/special_functions/erf.hpp>
#include <unsupported/Eigen/Polynomials>

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
    int nLambda;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point, const double w, const double h,
                 const int nLambda, const int nLambda2)
            : left_point(left_point),
              right_point(right_point), w(w), h(h), nLambda(nLambda) {
        //std::cout << "Residual@ L:" << left_point.transpose() << " " << right_point.transpose() << std::endl;
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
        //Eigen::VectorXd lambdas_d;




        double r_d = double_cast(r);
        //lambdas_d.resize(nLambda, 1);
        //lambdas2_d.resize(nLambda2, 1);

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
        T alpha = undistortionDenominator<T>(dr, lambdas);
        //for (size_t kk = 0; kk < lambdas.size(); ++kk)
        //     std::cout << lambdas[kk] << " !!! ";
        //std::cout << "\n" << alpha << "A" << std::endl;
        // std::cout <<  dr2 << "dr2" << std::endl;
        if (alpha <= .0)
            return false;

        Vector2T u1, u2;
        bool is_correct;
        std::pair<T, T> res;
        if (nLambda == 1) {
            res = undistortion_utils::computeEpipolarCurveDistanceError<T>(lambdas(0), F, left_point_T, right_point_T,
                                                                           r,
                                                                           w, h,
                                                                           u1, u2,
                                                                           is_correct);
        } else {
            res = undistortion_utils::computeEpipolarLineDistanceError<T>(lambdas, F, left_point_T, right_point_T, r,
                                                                          w, h,
                                                                          u1, u2,
                                                                          is_correct);
        }
        residuals[0] = res.first;
        residuals[1] = res.second;
        return is_correct;
    }
};


int main(int argc, char *argv[]) {
    double w = 7360.0, h = 4912.0, prcnt_inl = 0.2, r;
    std::string input1, input2;
    int nLambda, nLambda2;
    int nPictures;
    int n_f;
    std::string f_l;
    std::string f_f;

    int non_linear_iter = 1;
    std::vector<std::string> left_cor_vec_names;
    std::vector<std::string> right_cor_vec_names;
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
                ("left_corr_f", po::value<std::vector<std::string> >(&left_cor_vec_names)->multitoken()->required(),
                 "Path to left correspondeces")
                ("right_corr_f", po::value<std::vector<std::string> >(&right_cor_vec_names)->multitoken()->required(),
                 "Path to right correspondeces")
                ("q", po::value<double>(&prcnt_inl), "quantile to minimize")
                ("nlambda", po::value<int>(&nLambda)->default_value(1), "Number of parameters in denominator of model");

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


    r = std::sqrt(w * w + h * h) / 2.0;
    double minimal_residuals_per_points = std::numeric_limits<double>::max();
    for (size_t iters = 0; iters < non_linear_iter; ++iters) {

        if (iters == 0) {
            for (size_t kk = 0; kk < n_f; ++kk) {
                for (int ii = 0; ii < 3; ++ii)
                    for (int jj = 0; jj < 3; ++jj)
                        input_f1 >> Fvec[kk](ii, jj);
                double cur_l;
                input_f2 >> cur_l;
                Lambda[0] += cur_l;
            }
            Lambda[0] /= n_f;
        }

        ceres::Problem problem;
        double *lambda_ptr = Lambda.data();

        for (size_t ll = 0; ll < Lambda.size(); ++ll)
            std::cout << Lambda[ll] << " ";
        std::cout << std::endl;
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
            Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d, i1d, i2d;
            std::string name1 = left_cor_vec_names[kk];
            std::string name2 = right_cor_vec_names[kk];
            readPointsFromFile(name1, u1d);
            readPointsFromFile(name2, u2d);


            std::cout << "Prcnt:inl " << prcnt_inl << std::endl;
            std::cout << "Prcnt:inl " << boost::math::erfc_inv((0.95 + 1.0)) / boost::math::erfc_inv((0.3 + 1.0))
                      << std::endl;
            std::cout << "Prcnt:inl " << boost::math::erfc_inv((0.95 + 1.0)) / boost::math::erfc_inv((prcnt_inl + 1.0))
                      << std::endl;
            undistortion_utils::UndistortionProblemHelper helper(w, h, r, u1d, u2d, prcnt_inl);
            helper.setHypLambdas(Lambda);
            helper.setHypF(Fvec[kk]);
            helper.findInliers("OptimizerResults/non_linear_optimizer_inl_comp" + std::to_string(kk) + "_img_" +
                               std::to_string(iters) + "_iter");
            auto &inliers_ind = helper.getInliersIndices();

            i1d.resize(Eigen::NoChange, inliers_ind.size());
            i2d.resize(Eigen::NoChange, inliers_ind.size());
            for (size_t ll = 0; ll < inliers_ind.size(); ++ll) {
                i1d.col(ll) = u1d.col(inliers_ind[ll]);
                i2d.col(ll) = u2d.col(inliers_ind[ll]);
            }

            // readPointsFromFile(vec_names[kk], i1d);
            //readPointsFromFile(vec_names2[kk], i2d);


            undistortion_utils::normalizePoints(i1d, w, h, r);
            undistortion_utils::normalizePoints(i2d, w, h, r);


            problem.AddParameterBlock(f_ptr, 8);

            for (size_t k = 0; k < i1d.cols(); ++k) {
                Eigen::Vector3d left, right;
                left.template block<2, 1>(0, 0) = i1d.col(k);
                right.template block<2, 1>(0, 0) = i2d.col(k);
                left[2] = right[2] = 1.0;

                auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                        new ErrorFunctor(left, right, w, h, nLambda, 0));
                fun->AddParameterBlock(nLambda);
                fun->AddParameterBlock(8);
                fun->SetNumResiduals(2);
                problem.AddResidualBlock(fun, /*new ceres::HuberLoss(15)*/ nullptr, lambda_ptr,
                                         f_ptr);
                const double *ptrs[] = {lambda_ptr, f_ptr};
                double ptrs_res[2];
                fun->Evaluate(ptrs, ptrs_res, nullptr);
                /*std::cout << ptrs_res[0] << " " << ptrs_res[1] << " " << left.transpose() << " " << right.transpose()
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
        //std::cout << F << std::endl;
        std::cout << Lambda.transpose() << ": lambdas" << std::endl;
        double curr_residual_per_point = std::sqrt(summary.final_cost / residuals);
        std::cout << "final residual " << curr_residual_per_point << " (per point)" << std::endl;
        if (curr_residual_per_point < minimal_residuals_per_points) {
            minimal_residuals_per_points = curr_residual_per_point;
        }
        std::fstream all_fund_matrices("./OptimizerResults/all_f" + std::to_string(iters) + "_iters",
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