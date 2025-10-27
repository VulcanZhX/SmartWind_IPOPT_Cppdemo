#include <cmath>
#include <algorithm>
#include "Turbine.hpp"


Eigen::VectorXd velocity_function(
    const Eigen::VectorXd& delta_u_coff,
    double coeff,
    const Eigen::VectorXd& sigma_square,
    const Eigen::VectorXd& exp_term,
    const Eigen::VectorXd & x_d,
    const Eigen::VectorXd & y_d,
    const Eigen::VectorXd & z_d,
    const Turbine& turbine,
    const Eigen::VectorXd & deflection,
    double u_initial
);

Eigen::VectorXd deflection_function(
    const Eigen::VectorXd& x_d,
    const Turbine& turbine,
    const std::vector<int>& idx,
    const std::vector<Eigen::VectorXd>& y_model
);

Eigen::VectorXd turbulence_function(
    const Eigen::VectorXd& sigma_tm,
    const Eigen::VectorXd& R_half,
    const Eigen::VectorXd& x_d,
    const Eigen::VectorXd& y_d,
    const Eigen::VectorXd& z_d,
    const Turbine& turbine,
    const Eigen::VectorXd& deflection,
    const Eigen::VectorXd& turbulence_init);

// combination function
Eigen::VectorXd combination_function(
    const Eigen::VectorXd& u_wake,
    const Eigen::VectorXd& u_wake_induced
);

Eigen::VectorXd turbulence_combination_function(
    const Eigen::VectorXd& turbulence_wake,
    const Eigen::VectorXd& turbulence_wake_induced
);

Eigen::MatrixXd rot_func(Eigen::MatrixXd& pos, Eigen::VectorXd& center, double angle);

//// 2d-vector to MatrixXd
//Eigen::MatrixXd vector2matrix(const std::vector<std::vector<double>>& vec) {
//    int rows = vec.size();
//    int cols = vec[0].size();
//    Eigen::MatrixXd mat(rows, cols);
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            mat(i, j) = vec[i][j];
//        }
//    }
//    return mat;
//}
//
//// Alternative way of 2d-vector to MatrixXd
//Eigen::MatrixXd vector2matrix_alt(const std::vector<std::vector<double>>& vec) {
//    int rows = vec.size();
//    int cols = vec[0].size();
//    std::vector<double> flat_vec;
//	flat_vec.reserve(rows * cols);
//    for (const auto& row : vec) {
//        flat_vec.insert(flat_vec.end(), row.begin(), row.end());
//    }
//    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mat_map(flat_vec.data(), rows, cols);
//	return Eigen::MatrixXd(mat_map);
//}
//
//// MatrixXd to 2d-vector
//std::vector<std::vector<double>> matrix2vector(const Eigen::MatrixXd& mat) {
//    std::vector<std::vector<double>> vec(mat.rows(), std::vector<double>(mat.cols()));
//    for (int i = 0; i < mat.rows(); ++i) {
//        for (int j = 0; j < mat.cols(); ++j) {
//            vec[i][j] = mat(i, j);
//        }
//    }
//    return vec;
//}