#include "Toolset.hpp"
using namespace Eigen;

VectorXd turbulence_function(
    const VectorXd& sigma_tm,
    const VectorXd& R_half,
    const VectorXd& x_d,
    const VectorXd& y_d,
    const VectorXd& z_d,
    const Turbine& turbine,
    const VectorXd& deflection,
    const VectorXd& turbulence_init)
{
    int n = x_d.size();
    VectorXd I_d = VectorXd::Zero(n);

    // mask: 只计算下游点
    VectorXi mask(n);
    for (int i = 0; i < n; ++i) mask(i) = (x_d(i) >= 0) ? 1 : 0;

    // y_d(abs(y_d) < 0.0001) = 0.0001;
    VectorXd y_d_mod = y_d;
    for (int i = 0; i < n; ++i) {
        if (std::abs(y_d_mod(i)) < 0.0001) y_d_mod(i) = 0.0001;
    }

    // dr = sqrt((y_d - deflection).^2 + z_d.^2);
    VectorXd dr = ((y_d_mod - deflection).array().square() + z_d.array().square()).sqrt();

    // 参数
    double D = turbine.rotor_diameter;
    double gamma_deg = turbine.yaw_angle;
    double gamma_rad = gamma_deg * M_PI / 180.0;
    double Ct = turbine.getCt();
    double Ia = turbine.turbulence;

    // 推力系数接近0，附加湍流强度为0
    double d = 2.3 * std::pow(Ct * std::cos(gamma_rad), -1.2);
    double e = std::pow(Ia, 0.1);
    VectorXd q = 0.7 * std::pow(Ct * std::cos(gamma_rad), -3.2) * std::pow(Ia, -0.45) * (1.0 + x_d.array() / D).array().pow(-2);

    VectorXd dImax = (d + e * x_d.array() / D + q.array()).array().inverse();

    // 横向分布
    VectorXd shape_tm = VectorXd::Zero(n);
    VectorXd k1 = VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        if (dr(i) < R_half(i)) {
            shape_tm(i) = 1 - 0.15 * (1 + std::cos(M_PI * dr(i) / R_half(i)));
            k1(i) = std::sin(M_PI * 0.5 * dr(i) / R_half(i)); 
        }
        else {
            shape_tm(i) = std::exp(-std::pow(dr(i) - R_half(i), 2) / (2 * std::pow(sigma_tm(i), 2)));
            k1(i) = 1.0;
        }
    }

    // alpha_tm = atan2d(abs(z_d), abs(y_d))
    VectorXd alpha_tm(n);
    for (int i = 0; i < n; ++i) {
        alpha_tm(i) = std::atan2(std::abs(z_d(i)), std::abs(y_d_mod(i))) * 180.0 / M_PI;
    }

    // delta_tm
    VectorXd delta_tm = VectorXd::Zero(n);
    VectorXd exp_part = ((dr - R_half).array().square() / (2 * sigma_tm.array().square())).array().unaryExpr([](double v) { return std::exp(-v); });

    for (int i = 0; i < n; ++i) {
        double sind_alpha = std::sin(alpha_tm(i) * M_PI / 180.0);
        if (z_d(i) >= 0) {
            delta_tm(i) = 0.23 * Ia * sind_alpha * k1(i) * exp_part(i);
        }
        else {
            delta_tm(i) = -1.23 * Ia * sind_alpha * k1(i) * exp_part(i);
        }
    }

    // 结果
    for (int i = 0; i < n; ++i) {
        if (mask(i)) {
            I_d(i) = dImax(i) * shape_tm(i) + delta_tm(i);
        }
    }

    VectorXd wake_turbulence = I_d.array() * turbulence_init.array();
    return wake_turbulence;
}