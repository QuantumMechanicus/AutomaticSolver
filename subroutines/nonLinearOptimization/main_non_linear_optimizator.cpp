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


double undistortionDenominator(const double &r_distorted2, const Eigen::Matrix<double, Eigen::Dynamic, 1> &lambdas) {
    double denominator(1.0);
    double r_distorted2_pow = r_distorted2;
    for (int i = 0; i < lambdas.rows(); ++i) {
        denominator += lambdas[i] * r_distorted2_pow;
        r_distorted2_pow *= r_distorted2;
    }
    return denominator;
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
    T undistortionDenominator(const T &r_distorted2, const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) const {
        T denominator(1.0);
        T nominator(1.0);
        T r_distorted2_pow = r_distorted2;
        for (int i = 0; i < nLambda; ++i) {
            denominator += lambdas[i] * r_distorted2_pow;
            r_distorted2_pow *= r_distorted2;
        }
        r_distorted2_pow = r_distorted2;
        for (int i = 0; i < nLambda2; ++i) {
            nominator += lambdas[nLambda + i] * r_distorted2_pow;
            r_distorted2_pow *= r_distorted2;
        }
        return denominator / nominator;
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
        u1 = left_point_T / denominator1;
        u2 = right_point_T / denominator2;
        u1[2] = u2[2] = T(1.0);

        Vector3T l1 = recompute_F * u1;
        Vector3T l2 = recompute_F.transpose() * u2;

        T n1 = l1.template block<2, 1>(0, 0).norm();
        T n2 = l2.template block<2, 1>(0, 0).norm();
        T err = l1.dot(u2);
        T err2 = err * err;

        residuals[0] = alpha * r * err / n1;
        residuals[1] = alpha * r * err / n2;
        return denominator1 > 0.0 && denominator2 > 0.0;
    }
};

void normalizePoints(auto &points, double w, double h, double r) {
    points.row(0).array() -= w / 2.0;
    points.row(1).array() -= h / 2.0;
    points /= r;
}


size_t findInliers(double w, double h, Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, Eigen::Matrix<double, 2, Eigen::Dynamic> u2d, const Eigen::Matrix<double, Eigen::Dynamic, 1> &lambdas, const Eigen::Matrix3d &hyp_F,
            const std::string &out_name) {

    Eigen::Matrix<double, 2, Eigen::Dynamic> u1, u2;
    u1.resize(Eigen::NoChange, u1d.cols());
    u2.resize(Eigen::NoChange, u2d.cols());
    double r  = std::sqrt((h/2.0)*(h/2.0) + (w/2.0)*(w/2.0));
    double d = std::max(h, w) / (2.0);
    double dr = d / r;
    double dr2 = dr * dr;
    double alpha = undistortionDenominator(dr2, lambdas);



    Eigen::VectorXd r1d(u1d.cols());

    r1d.col(0) = (u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1))).transpose();


    Eigen::VectorXd r2d(u1d.cols());


    r2d.col(0) = (u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1))).transpose();


    for (size_t k = 0; k < u1d.cols();++k)
    {
        u1.col(k) = u1d.col(k)/undistortionDenominator(r1d[k], lambdas);
        u2.col(k) = u2d.col(k)/undistortionDenominator(r2d[k], lambdas);

    }


    Eigen::Matrix<double, 3, Eigen::Dynamic> uu1, uu2;
    uu1.resize(Eigen::NoChange, u1.cols());
    uu2.resize(Eigen::NoChange, u2.cols());
    uu1.setOnes();
    uu2.setOnes();
    uu1.row(0) = u1.row(0);
    uu2.row(0) = u2.row(0);

    uu1.row(1) = u1.row(1);
    uu2.row(1) = u2.row(1);


    Eigen::Matrix<double, 3, Eigen::Dynamic> l1 = (hyp_F * uu1);
    Eigen::Matrix<double, 3, Eigen::Dynamic> l2 = (hyp_F.transpose() * uu2);

    size_t goods = 0;

    std::fstream errf1(out_name + "_left", std::ios_base::out);
    std::fstream errf2(out_name + "_right", std::ios_base::out);
    Eigen::Matrix<double, 1, 2> center;
    center(0, 0) = w / 2.0;
    center(0, 1) = h / 2.0;
    std::vector<double> err_sample(u1d.cols());


    for (size_t k = 0; k < u1d.cols(); ++k) {
        double c1 = l1.col(k).template topRows<2>().norm();
        double c2 = l2.col(k).template topRows<2>().norm();
        double err1 = alpha * r* std::abs(uu2.col(k).dot(l1.col(k)) / c1);
        double err2 = alpha * r*std::abs(uu1.col(k).dot(l2.col(k)) / c2);
        err_sample[k] = err1 + err2;
    }
    std::cout << "W\n";
    std::nth_element(err_sample.begin(), err_sample.begin() + err_sample.size() / 10, err_sample.end());
    double quantile = err_sample[err_sample.size() / 10];
    const double confidence_interval = 3.36752;
    for (size_t k = 0; k < u1d.cols(); ++k) {
        double err = err_sample[k];

        if (std::abs(err) < quantile * confidence_interval) {
            ++goods;
            errf1 << r * u1d.col(k).transpose().leftCols(2) + center << "\n";
            errf2 << r * u2d.col(k).transpose().leftCols(2) + center << "\n";
        }

    }

    errf1.close();
    errf2.close();

    return goods;
}

int main(int argc, char *argv[]) {
    double w = 7360.0, h = 4912.0, r;
    std::string input1, input2;
    int nLambda;
    int nPictures;
    int n_f;
    std::string f_l;
    std::string f_f;


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
    const int non_linear_iter = 1;
    int nLambda2 = 0;
    Eigen::Matrix<double, Eigen::Dynamic, 1> Lambda(nLambda + nLambda2, 1);
    StdVector<Eigen::Matrix3d> Fvec(n_f);
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
            }
            Lambda[0] /= n_f;
        }



        ceres::Problem problem;
        double *lambda_ptr = Lambda.data();


        problem.AddParameterBlock(lambda_ptr, nLambda + nLambda2);
        for (size_t k = 0; k < Fvec.size(); ++k) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;

            Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                      fmatrix_svd.matrixV().transpose();
            Fvec[k] = Fvec[k] / Fvec[k](2, 2);
        }

        int residuals = 0;
        for (size_t kk = 0; kk < n_f; ++kk) {
            double *f_ptr = Fvec[kk].data();
            Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
            std::string name1;
            std::string name2;
            if (iters == 0) {
                name1 = vec_names[kk];
                name2 = vec_names2[kk];
            } else{
                name1 = std::to_string(kk +1)+"_found_inliers_left";
                name2 = std::to_string(kk +1)+"_found_inliers_right";
            }
            readPointsFromFile(name1, i1d);
            readPointsFromFile(name2, i2d);

            normalizePoints(i1d, w, h, r);
            normalizePoints(i2d, w, h, r);

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
                problem.AddResidualBlock(fun, new ceres::HuberLoss(15), lambda_ptr,
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
        std::cout << "Find inliers...\n";
        for (size_t image_n = 0; image_n < n_f; ++image_n) {
            Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
            readPointsFromFile(left_cor_vec_names[image_n],i1d);
            readPointsFromFile(right_cor_vec_names[image_n],i2d);
            normalizePoints(i1d, w, h, r);
            normalizePoints(i2d, w, h, r);

            std::string out_name = std::to_string(image_n+1)+"_found_inliers";
            findInliers(w, h, i1d, i2d, Lambda, Fvec[image_n],out_name);
        }

        for (size_t k = 0; k < Fvec.size(); ++k) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(Fvec[k], Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;

            Fvec[k] = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                      fmatrix_svd.matrixV().transpose();
            Fvec[k] = Fvec[k] / Fvec[k](2, 2);
        }

    }
    std::fstream fm("fundamental_matrices", std::fstream::out);
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