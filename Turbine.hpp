#ifndef TURBINE_HPP
#define TURBINE_HPP

#include <Dense>
#include <vector>
#include <array>
const double M_PI = 3.1416;
//using namespace Eigen;

class Turbine {
public:
    // ���캯��
    Turbine(
        const std::vector<std::vector<double>>& PT,
        const std::array<std::vector<int>, 3>& t_qz,
        const std::vector<double>& turbulence_sheer_veer,
        const std::vector<double>& turbine_diameter_vector,
        const std::vector<double>& turbine_hub_height_vector,
        const std::vector<double>& rated_power_vector,
        const std::vector<double>& life_total_vector,
        const std::vector<double>& repair_c_vector,
        const std::vector<double>& fatigue,
        const std::vector<double>& fatigue_p,
		int count_tn // �����ţ���1��ʼ
    );

    // ����
    double getAverageVelocity() const;
    double getCp() const;
    double getCt() const;
    double getPower() const;
    double getComprehensiveFatigueCoefficient() const;
    double getSingleTurbineObjective() const;
    double getSingleTurbineGeneration() const;

    // �����������
    double fCp(double at_wind_speed) const;
    double fCt(double at_wind_speed) const;

    // ���·��١�������
    //u_initial: 1x1 scalar indicating the global natural windspeed
	void updateVelocities(const std::vector<double>& u_wake, const double u_initial, int nt);
    void updateTurbulenceIntensity(const std::vector<double>& u_turbulence_wake, int nt);

    // ���¶���ת�Ӱ뾶��ֱ��
    void setRotorRadius(double r);
    void setRotorDiameter(double d);

    // ��������
    void updateGrid();

    // ������Ա�������ɸ�����Ҫ��Ϊprivate/protected����getter/setter��
    double rotor_diameter;
    double hub_height;
    double rated_power;
    double life_total;
    double repair_c;
    double optimization_period = 1;
    double annual_average_power;
    double past_comprehensive_fatigue_coefficient;
    double ref_turbulence = 0.1;
    double dis_coefficient = 0.7;
    double consumed_comprehensive_fatigue_coefficient = 0;
    double life_work_coeff;
    double life_turbulence_coeff;
    double generator_efficiency = 1;

	// optim variable
    double yaw_angle = 0;

    double air_density = 1.225;
    double turbulence_ambient;
    double rotor_radius;
    std::vector<double> velocities_u; // wind speed at grid points, to be updated
    double turbulence;
    std::vector<double> efficiency_value;
    std::vector<std::vector<double>> power_thrust_table_11;
    std::vector<std::vector<double>> power_thrust_table_68;
    std::vector<std::vector<double>> power_thrust_table_83;
    std::vector<double> grid;
    int turbine_idx_org;
    static constexpr int points_turbine_grid = 1;

private:
    // �ڲ���������
    std::vector<double> calculateTurbineVelocities(const std::vector<double>& local_wind_speed, int nt) const;
    void updateRadius();
};

#endif // TURBINE_HPP
