#include "Toolset.hpp"
using namespace Eigen;

inline static double cosd(double deg) { return std::cos(deg * M_PI / 180.0); }
inline static double sind(double deg) { return std::sin(deg * M_PI / 180.0); }

//����˵����
// x_d: std::vector<double>�����λ�ã�rotate_x-x_sorted��
// turbine: Turbine���󣬰����������
// idx: std::vector<int>���߼�������0/1��(1x159, ��ʾ��Щ������compact)
// y_model: std::vector<std::vector<VectorXd>>��ƫתģ�����������ʽ��Ͻ��
VectorXd deflection_function(
    const VectorXd& x_d,
    const Turbine& turbine,
    const std::vector<int>& idx,
    const std::vector<VectorXd>& y_model
) {
    // ��ϰ汾
    double rotor_D = turbine.rotor_diameter;
	int n = x_d.size();

    // ��Eigen VectorXd��Ϊ�м���
    VectorXd y_c = VectorXd::Zero(n);

    int index;
    if (rotor_D - 228 == 0) 
        index = 1;
    else if (rotor_D - 158 == 0) 
        index = 2;
    else
        index = 3;


    // Ԥ��������
    double factor = turbine.getCt() * std::pow(cosd(turbine.yaw_angle), 2) * sind(turbine.yaw_angle);

    // idxתEigen����
    std::vector<int> idx_vec = idx;
    bool any_idx = false;
	for (int i : idx_vec) if (i) { any_idx = true; break; } // ����Ƿ����κ����е㣬����Ϊtrue
    if (any_idx) {
        // y_model[0][index] �� std::vector<double>
		y_c = factor * y_model[index];
    }

	return y_c;
}
