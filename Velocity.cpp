//#include <vector>
//#include <cmath>
//#include <algorithm>
//#include "Turbine.hpp"
//
//struct Turbine {
//    double Ct;
//    double yaw_angle;
//    double rotor_diameter;
//};
//
//void velocity_function(
//    const std::vector<double>& delta_u_coff,
//    const std::vector<double>& coeff,
//    const std::vector<double>& sigma_square,
//    const std::vector<double>& exp_term,
//    const std::vector<double>& x_d,
//    const std::vector<double>& y_d,
//    const std::vector<double>& z_d,
//    const Turbine& turbine,
//    double deflection,
//    const std::vector<double>& u_initial,
//    std::vector<double>& turb_u_wake,
//    std::vector<double>& turb_v_wake,
//    std::vector<double>& turb_w_wake
//) {
//    size_t n = x_d.size();
//    turb_u_wake.resize(n);
//    turb_v_wake.resize(n, 0.0);
//    turb_w_wake.resize(n, 0.0);
//
//    double cos_yaw = std::cos(turbine.yaw_angle * M_PI / 180.0);
//
//    for (size_t i = 0; i < n; ++i) {
//        double delta_U = (1.0 - std::sqrt(1.0 - turbine.Ct * cos_yaw * cos_yaw)) * delta_u_coff[i];
//        double percent_deficit = delta_U * coeff[i] * std::exp(-std::pow(y_d[i] - deflection, 2) / sigma_square[i]) * exp_term[i];
//        double deficit = percent_deficit * u_initial[i];
//
//        // MATLAB: deficit(x_d < 0) = 0;
//        //         deficit(x_d > 12000) = 0;
//        if (x_d[i] < 0.0 || x_d[i] > 12000.0) {
//            deficit = 0.0;
//        }
//
//        turb_u_wake[i] = deficit;
//        // v, w ·ÖÁ¿ÎªÁã
//        // turb_v_wake[i] = 0.0;
//        // turb_w_wake[i] = 0.0;
//    }
//}
