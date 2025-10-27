#include "Turbine.hpp"
#include "Toolset.hpp"

/*ע�⣺������/�����������������£�����ֵ�����飩ʱʹ��vector,����һ��ʹ��Eigen��*/

typedef std::vector<std::tuple<Eigen::MatrixXd, Turbine, int>> TURB_CELLARR;
static constexpr std::size_t MaxTurbines = 159; //����������������һ��

struct Turbine_cell {
    Eigen::MatrixXd layout;
    std::vector<Turbine> turbine_chart;
    std::array<int, MaxTurbines> idx;
};

//// ����ṹ�壺����������Ҫ���ص� rotated ��Ϣ
// �ṹ�嶨��
struct RotatedResult {
    // �����ķ������ (N x 3)
    Eigen::MatrixXd sorted_coords;      // (N, 3)
    std::array<int, MaxTurbines> sorted_indexes;     // (N,)

    // ���з��������������ת����
    Eigen::VectorXd rot_x_grid;
    Eigen::VectorXd rot_y_grid;
    Eigen::VectorXd rot_z_grid;

    // x_d, y_d, z_d
    std::vector<Eigen::VectorXd> x_d;
    std::vector<Eigen::VectorXd> y_d;
    std::vector<Eigen::VectorXd> z_d;

    // β������ģ��
    std::vector<Eigen::VectorXd> delta_u;    
    std::vector<Eigen::VectorXd> sigma_square;
    std::vector<double> coeff; 
    std::vector<Eigen::VectorXd> exp_term; 

    // ����ģ��
    std::vector<Eigen::VectorXd> sigma_tm; 
    std::vector<Eigen::VectorXd> R_half;

    // ƫתģ�� y_model: ÿ̨���3�����ʽ���
    std::vector<std::vector<Eigen::VectorXd>> y_model;

    // mask: (P, N) �߼����루0/1��
    std::vector<std::vector<int>> idx_mask;
};

class WindFarmOptimization {
public:
    // ����Դ����
    double wind_speed;
    double wind_direction;
    double turbulence_intensity;
    double added_turbulence_intensity = 0.0;
    double wind_shear;
    double wind_veer;
    int n_turbines;
    std::array<int, MaxTurbines> idx_org;

    // �糡����(originally included in turbinechart class)
    Eigen::MatrixXd layout; // ÿ��һ�������[x, y, z]
    std::vector<Turbine> turbine_chart;
    TURB_CELLARR turbine_tuple;

    // �Ż���ز���
    double yaw_lower = -30.0;
    double yaw_upper = 30.0;
    // �������(92+67=159)
    std::vector<int> qz_12;
    std::vector<int> qz_3;

    // �糡�������
    Eigen::VectorXd u; // ���ٳ�
    Eigen::VectorXd turbulence; // ������
    
    // ���캯��
    WindFarmOptimization(
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
        const std::vector<std::vector<double>>& layout_farm,
        double wind_speed,
        double wind_direction
        );

    // �糡�߽�
    Eigen::VectorXd getBounds() const;

    // �������
    std::vector<int> getIndexes();

    //��ȡ����xyz����
    Eigen::VectorXd getX() const;
    Eigen::VectorXd getY() const;
    Eigen::VectorXd getZ() const;

    // rotated���Ԥ����
    RotatedResult compute_rotated();

    // ��ȡ���з��ƫ����
    Eigen::VectorXd getYawAngles() const;

    // �������з��ƫ����
    void setYawAngles(const Eigen::VectorXd& yaw_angles);

    // ��ȡ���з������
    Eigen::VectorXd getTurbinesPower() const;

    // ��ȡ���з������
    Eigen::VectorXd getTurbinesTurbulence() const;

    // ��ȡ�糡�ܹ���
    double getFarmPower() const;

    // ��ȡ����12�糡����
    double getFarmQingzhou12Power() const;

    // ��ȡ����3�糡����
    double getFarmQingzhou3Power() const;

    // ��ȡ��Ⱥ���������

	// ��ȡ��Ⱥ���з��ƽ������

    // ��ȡ��Ⱥ������Ż�ָ��

    // ��ȡ��Ⱥ��ָ��

    // �Ż�����ڣ�ʾ���������Ż���
    Eigen::VectorXd optimizeYawAngles();

    // ����β��Ӱ��
    void calculateWake();

    // �����ӿںͳ�Ա����������Ҫ����...

    Turbine_cell rotated_turbine(Eigen::VectorXd& center, double angle);

private:
    // Ԥ���㻺���
    // ...
};

Turbine_cell sortInX(const Turbine_cell& turbine_arr);
Eigen::VectorXd polyval(const Eigen::VectorXd& x, const Eigen::VectorXd& coeffs);