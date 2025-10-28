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
