#include "mech.h"
#include <Eigen\Dense>
using namespace std;
using namespace Eigen;

namespace INS
{
    typedef struct ST_Kalman_Mat
    {
        Eigen::Matrix<double, 21, 21> F_21_21;
        Eigen::Matrix<double, 21, 18> G_21_18;
        Eigen::Matrix<double, 18, 18> q_18_18;
        Eigen::Matrix<double, 21, 1> x_21_1;
        Eigen::Matrix<double, 21, 21> P_21_21;
        ST_Kalman_Mat() : F_21_21(Eigen::Matrix<double, 21, 21>::Zero()),
                          G_21_18(Eigen::Matrix<double, 21, 18>::Zero()),
                          q_18_18(Eigen::Matrix<double, 18, 18>::Zero()),
                          x_21_1(Eigen::Matrix<double, 21, 1>::Zero()),
                          P_21_21(Eigen::Matrix<double, 21, 21>::Zero())
        {
            ;
        }

    } Kalman_Mat_T;

    typedef struct ST_IMU_parameter
    {
        Eigen::Vector3d g_b;
        Eigen::Vector3d g_s;
        Eigen::Vector3d a_b;
        Eigen::Vector3d a_s;
        ST_IMU_parameter() : g_b(Eigen::Vector3d::Zero()),
                             g_s(Eigen::Vector3d::Zero()),
                             a_b(Eigen::Vector3d::Zero()),
                             a_s(Eigen::Vector3d::Zero())
        {
            ;
        }
    } IMU_parameter_T;

    // GNSS ���
    typedef struct ST_GNSS_data
    {
        double timestamp;               //ʱ���
        Eigen::Vector3d pos_blh;        //λ��-γ����
        Eigen::Vector3d pos_blh_std;    //λ��-γ���߱�׼��
    } GNSS_data_T;

    // �ο��������
    typedef struct ST_NAV_data
    {
        double week;
        double second;
        Eigen::Vector3d pos_blh;
        Eigen::Vector3d vel_ned;
        Eigen::Vector3d att_ypr;
    } NAV_data_T;

    class intergrated : public mechanization
    {
    public:
        intergrated(m_State state, IMU_data data, Vector3d lb, Vector3d pos_std, Vector3d vel_std, Vector3d att_std);
        ~intergrated();

    private:
        Kalman_Mat_T kal_mat_t1;            // k-1ʱ�̵���չ�������˲��ľ���
        Kalman_Mat_T kal_mat_t2;            // kʱ�̵���չ�������˲��ľ���
        Eigen::Vector3d l_b;                // �˱���Ϣ �ɳ�ʼ������
        GNSS_data_T GNSS_data;              // GNSS rtk ����
        IMU_parameter_T imu_parameter;      // IMU�ı궨����, ��������¶�����
        // Eigen::Vector3d z_k;
        // Eigen::Matrix3d R_3_3;

    public:
        void EKFpredictUpdate();            // Ԥ�����
        void EKFmeasurementUpdate();        // �������

        void ErrorFeedback();               // �����ܺ���
        void AtitudeErrorFeedback();        // ��̬����
        void PosErrorFeedback();            // λ������
        void VelErrorFeedback();            // �ٶ�����
        void IMUErrorFeedback();            // IMU�궨��������

        int IntergtatedNavUpdate(const IMU_data &imudata, const GNSS_data_T &gnssdata); // ��ϵ��������ܺ���
        void stateExtrapolation();          // λ������

        IMU_data IMUdataFix(const IMU_data &data);  // IMUԤ����
    };
}