#include "intergratedNav.h"
#include "Eigen\Dense"

namespace INS
{
    intergrated::intergrated(m_State state, IMU_data data, Eigen::Vector3d lb, Vector3d pos_std, Vector3d vel_std, Vector3d att_std) : mechanization(state, data), l_b(lb)
    {
        Eigen::Matrix<double, 21, 1> P_21_1;
        Eigen::Matrix<double, 3, 3> C_np;
        Eigen::Matrix<double, 3, 3> phi_cross;
        Eigen::Vector3d phi_std;

        Eigen::Matrix<double, 3, 3> I_3_3;
        Eigen::Matrix<double, 3, 1> I_3_1;

        I_3_3 = Eigen::Matrix<double, 3, 3>::Identity();
        I_3_1.setOnes();

        C_np = Euler2Quart(att_std).matrix();
        phi_cross = C_np - I_3_3;
        phi_std << phi_cross(1, 2), phi_cross(2, 0), phi_cross(0, 1);

        P_21_1.block<3, 1>(0, 0) = (pos_std.array().square()).matrix();
        P_21_1.block<3, 1>(3, 0) = (vel_std.array().square()).matrix();
        P_21_1.block<3, 1>(6, 0) = (phi_std.array().square()).matrix();

        P_21_1.block<3, 1>(9, 0) = pow2Func(gb_std) * I_3_1;
        P_21_1.block<3, 1>(12, 0) = pow2Func(ab_std) * I_3_1;
        P_21_1.block<3, 1>(15, 0) = pow2Func(gs_std) * I_3_1;
        P_21_1.block<3, 1>(18, 0) = pow2Func(as_std) * I_3_1;

        kal_mat_t2.P_21_21 = P_21_1.asDiagonal();
        kal_mat_t1 = kal_mat_t2;
    }
    intergrated::~intergrated()
    {
        ;
    }
    void intergrated::EKFpredictUpdate()
    {
        // 计算F阵
        double R_M;
        double R_N;
        Eigen::Vector3d f_b;
        Eigen::Vector3d omega_ib_b;
        Eigen::Vector3d omega_ie_n;
        Eigen::Vector3d omega_en_n;
        Eigen::Vector3d omega_in_n;
        Eigen::Vector3d gp;
        Eigen::Matrix3d F_rr;
        Eigen::Matrix3d F_vr;
        Eigen::Matrix3d F_vv;
        Eigen::Matrix3d F_phir;
        Eigen::Matrix3d F_phiv;
        Eigen::Matrix3d I_3_3;
        Eigen::Matrix3d C_bn;

        Eigen::Matrix<double, 21, 21> F_21_21;
        Eigen::Matrix<double, 21, 21> Phi_21_21;

        Eigen::Matrix<double, 21, 21> I_21_21;
        Eigen::Matrix<double, 21, 18> G_21_18;
        Eigen::Matrix<double, 18, 18> q_18_18;
        Eigen::Matrix<double, 21, 21> Q_21_21;

        R_M = getRM(state_t2.GEO_eb(Latitude));
        R_N = getRN(state_t2.GEO_eb(Latitude));
        f_b = imudata_t2.accel / INS_T;
        omega_ib_b = imudata_t2.gyro / INS_T;
        omega_ie_n = Pos2AngleVelocity_ie(state_t2.GEO_eb);
        omega_en_n = PosVec2AngleVelocity_en(state_t2.GEO_eb, state_t2.NED_vec);
        omega_in_n = omega_ie_n + omega_en_n;
        gp = getG(state_t2.GEO_eb);

        getFrr(state_t2.GEO_eb, state_t2.NED_vec, R_N, R_M, F_rr);
        getFvr(state_t2.GEO_eb, state_t2.NED_vec, R_N, R_M, gp(D), F_vr);
        getFvv(state_t2.GEO_eb, state_t2.NED_vec, R_N, R_M, F_vv);
        getFphir(state_t2.GEO_eb, state_t2.NED_vec, R_N, R_M, F_phir);
        getFphiv(state_t2.GEO_eb, R_N, R_M, F_phiv);
        I_3_3 = Eigen::Matrix3d::Identity();
        C_bn = state_t2.q_bn.matrix();

        F_21_21 = F_21_21.setZero();
        F_21_21.block<3, 3>(0, 0) = F_rr;
        F_21_21.block<3, 3>(0, 3) = I_3_3;

        F_21_21.block<3, 3>(3, 0) = F_vr;
        F_21_21.block<3, 3>(3, 3) = F_vv;
        F_21_21.block<3, 3>(3, 6) = Vector2CrossMatrix(C_bn * f_b);
        F_21_21.block<3, 3>(3, 12) = C_bn;
        F_21_21.block<3, 3>(3, 18) = C_bn * (f_b.asDiagonal());

        F_21_21.block<3, 3>(6, 0) = F_phir;
        F_21_21.block<3, 3>(6, 3) = F_phiv;
        F_21_21.block<3, 3>(6, 6) = -(Vector2CrossMatrix(omega_in_n));
        F_21_21.block<3, 3>(6, 9) = -C_bn;
        F_21_21.block<3, 3>(6, 15) = -C_bn * (omega_ib_b.asDiagonal());
        
        
        F_21_21.block<3, 3>(9, 9) = (-1.0 / Tgb) * I_3_3;
        F_21_21.block<3, 3>(12, 12) = (-1.0 / Tab) * I_3_3;
        F_21_21.block<3, 3>(15, 15) = (-1.0 / Tgs) * I_3_3;
        F_21_21.block<3, 3>(18, 18) = (-1.0 / Tas) * I_3_3;
        //
        kal_mat_t2.F_21_21 = F_21_21;
        // 陈. lecture 8 式(20)
        I_21_21 = Eigen::Matrix<double, 21, 21>::Identity();
        Phi_21_21 = I_21_21 + kal_mat_t1.F_21_21 * INS_T;

        G_21_18.setZero();
        G_21_18.block<3, 3>(3, 0) = C_bn;
        G_21_18.block<3, 3>(6, 3) = C_bn;
        G_21_18.block<3, 3>(9, 6) = I_3_3;
        G_21_18.block<3, 3>(12, 9) = I_3_3;
        G_21_18.block<3, 3>(15, 12) = I_3_3;
        G_21_18.block<3, 3>(18, 15) = I_3_3;
        kal_mat_t2.G_21_18 = G_21_18;

        q_18_18.setZero();
        q_18_18.block<3, 3>(0, 0) = pow2Func(VRW_std) * I_3_3;
        q_18_18.block<3, 3>(3, 3) = pow2Func(ARW_std) * I_3_3;
        q_18_18.block<3, 3>(6, 6) = 2.0 * pow2Func(gb_std) * I_3_3 / Tgb;
        q_18_18.block<3, 3>(9, 9) = 2.0 * pow2Func(ab_std) * I_3_3 / Tab;
        q_18_18.block<3, 3>(12, 12) = 2.0 * pow2Func(gs_std) * I_3_3 / Tgs;
        q_18_18.block<3, 3>(15, 15) = 2.0 * pow2Func(as_std) * I_3_3 / Tas;
        kal_mat_t2.q_18_18 = q_18_18;
        // 陈. lecture 8 式(23)
        Q_21_21 = 0.5 * (Phi_21_21 * kal_mat_t1.G_21_18 * kal_mat_t1.q_18_18 * (kal_mat_t1.G_21_18.transpose()) * (Phi_21_21.transpose()) + kal_mat_t2.G_21_18 * kal_mat_t2.q_18_18 * (kal_mat_t2.G_21_18.transpose())) * INS_T;
        // 陈. lecture 8 式(4)
        kal_mat_t2.x_21_1 = Phi_21_21 * kal_mat_t1.x_21_1;
        // 陈. lecture 8 式(5)
        kal_mat_t2.P_21_21 = Phi_21_21 * kal_mat_t1.P_21_21 * (Phi_21_21.transpose()) + Q_21_21;
    }
    void intergrated::EKFmeasurementUpdate()
    {
        Eigen::Matrix<double, 3, 21> H_3_21;
        Eigen::Matrix3d C_bn;
        Eigen::Matrix<double, 21, 3> K_21_3;
        Eigen::Matrix3d I_3_3;
        Eigen::Matrix<double, 21, 21> I_21_21;

        double R_M;
        double R_N;
        Eigen::Vector3d DR_inv_3_1;
        Eigen::Matrix3d DR_inv_3_3;                 // DR^-1用于

        Eigen::Vector3d r_G_n;                      // imu推算的GNSS天线相位中心的位置

        Eigen::Vector3d deltaz_r;                   // 观测向量
        Eigen::Matrix3d Rk_3_3;                     // 测量噪声方差阵

        I_3_3 = Eigen::Matrix3d::Identity();
        I_21_21 = Eigen::Matrix<double, 21, 21>::Identity();
        C_bn = state_t2.q_bn.toRotationMatrix();

        Rk_3_3 = (((GNSS_data.pos_blh_std.array()).square()).matrix()).asDiagonal();// 根据gnss位置标准差求测量方差阵
        // 陈. lecture 8 式(33) H阵
        H_3_21.setZero();
        H_3_21.block<3, 3>(0, 0) = I_3_3;
        H_3_21.block<3, 3>(0, 6) = Vector2CrossMatrix(C_bn * l_b);
        // 陈. lecture 8 式(25)
        R_M = getRM(state_t2.GEO_eb(Latitude));
        R_N = getRN(state_t2.GEO_eb(Latitude));
        DR_inv_3_1 << 1.0 / (R_M + state_t2.GEO_eb(Height)), 1.0 / ((R_N + state_t2.GEO_eb(Height)) * cos(state_t2.GEO_eb(Latitude))), -1.0;
        DR_inv_3_3 = DR_inv_3_1.asDiagonal();
        // 陈. lecture 8 式(26)
        r_G_n = state_t2.GEO_eb + DR_inv_3_3 * C_bn * l_b;
        // 陈. lecture 8 式(31)
        deltaz_r = DR_inv_3_3.inverse() * (r_G_n - GNSS_data.pos_blh);
        // 陈. lecture 8 式(6)
        K_21_3 = kal_mat_t2.P_21_21 * (H_3_21.transpose()) * ((H_3_21 * kal_mat_t2.P_21_21 * (H_3_21.transpose()) + Rk_3_3).inverse());
        // 陈. lecture 8 式(7)
        kal_mat_t2.x_21_1 = kal_mat_t2.x_21_1 + K_21_3 * (deltaz_r - H_3_21 * kal_mat_t2.x_21_1);
        // 陈. lecture 8 式(8)
        kal_mat_t2.P_21_21 = (I_21_21 - K_21_3 * H_3_21) * kal_mat_t2.P_21_21 * (I_21_21 - K_21_3 * H_3_21).transpose() + K_21_3 * Rk_3_3 * K_21_3.transpose();
    }

    void intergrated::ErrorFeedback()
    {
        AtitudeErrorFeedback();
        PosErrorFeedback();
        VelErrorFeedback();
        IMUErrorFeedback();
    }
    // void intergrated::AtitudeErrorFeedback()
    // {
    //     // 赖.
    //     Eigen::Vector3d phi;
    //     Eigen::Matrix3d I_3_3;
    //     Eigen::Matrix3d C_bp;
    //     Eigen::Matrix3d C_tp;
    //     Eigen::Matrix3d C_bt;
    //     Eigen::Quaterniond q_bt;

    //     I_3_3 = Eigen::Matrix3d::Identity();
    //     phi = kal_mat_t2.x_21_1.block<3, 1>(6, 0);
    //     C_bp = state_t2.q_bn.matrix();
    //     C_tp = I_3_3 - Vector2CrossMatrix(phi);
    //     C_bt = C_tp.inverse() * C_bp;
    //     q_bt = C_bt;

    //     state_t2.q_bn = q_bt;
    //     state_t2.q_bn.normalize(); //归一化，否则旋转矩阵可能不正交
    //     state_t2.e_bn = DCM2Euler(((state_t2.q_bn).matrix()));
    //     // 反馈后清零
    //     kal_mat_t2.x_21_1.block<3, 1>(6, 0).setZero();
    // }

    void intergrated::PosErrorFeedback()
    {
        double R_M;
        double R_N;

        Eigen::Vector3d delta_r;
        Eigen::Vector3d DR_inv_3_1;
        Eigen::Matrix3d DR_inv_3_3;                 // DR^-1用于

        delta_r = kal_mat_t2.x_21_1.block<3, 1>(0, 0);

        R_M = getRM(state_t2.GEO_eb(Latitude));
        R_N = getRN(state_t2.GEO_eb(Latitude));
        DR_inv_3_1 << 1.0 / (R_M + state_t2.GEO_eb(Height)), 1.0 / ((R_N + state_t2.GEO_eb(Height)) * cos(state_t2.GEO_eb(Latitude))), -1.0;
        DR_inv_3_3 = DR_inv_3_1.asDiagonal();

        state_t2.GEO_eb = state_t2.GEO_eb - DR_inv_3_3 * delta_r;
        state_t2.q_ne = NED2Quart(state_t2.GEO_eb);

        kal_mat_t2.x_21_1.block<3, 1>(0, 0).setZero();
    }
    // void intergrated::VelErrorFeedback()
    // {
    //     Eigen::Vector3d delta_v;
    //     Eigen::Matrix3d C_tp;
    //     Eigen::Vector3d phi;
    //     Eigen::Matrix3d I_3_3;

    //     I_3_3 = Eigen::Matrix3d::Identity();
    //     phi = kal_mat_t2.x_21_1.block<3, 1>(6, 0);
    //     C_tp = I_3_3 - Vector2CrossMatrix(phi);

    //     delta_v = kal_mat_t2.x_21_1.block<3, 1>(3, 0);
    //     state_t2.NED_vec = C_tp.inverse() * state_t2.NED_vec - delta_v;
    //     kal_mat_t2.x_21_1.block<3, 1>(3, 0).setZero();
    // }

    void intergrated::AtitudeErrorFeedback()
    {
        Eigen::Vector3d phi;
        Eigen::Matrix3d I_3_3;
        Eigen::Quaterniond q_pt;
        Eigen::Matrix3d C_bt;

        I_3_3 = Eigen::Matrix3d::Identity();
        phi = kal_mat_t2.x_21_1.block<3, 1>(6, 0);
        q_pt = SpinVector2Quart(phi);
        state_t2.q_bn = q_pt * state_t2.q_bn;
        state_t2.q_bn.normalize(); //归一化，否则旋转矩阵可能不正交
        state_t2.e_bn = DCM2Euler(((state_t2.q_bn).matrix()));
        // 反馈后清零
        kal_mat_t2.x_21_1.block<3, 1>(6, 0).setZero();
    }
    // void intergrated::PosErrorFeedback()
    // {
    //     Eigen::Vector3d delta_theta;
    //     Eigen::Quaterniond q_nc;
    //     double R_M;
    //     double R_N;
    //     double latitude_fix;
    //     double longitude_fix;
    //     double height_fix;
    //     R_M = getRM(state_t2.GEO_eb(Latitude));
    //     R_N = getRN(state_t2.GEO_eb(Latitude));
    //     // Eun-Hwan. 式(2.35)
    //     delta_theta << kal_mat_t2.x_21_1(1, 0) / (R_N + state_t2.GEO_eb(Height)),
    //         -kal_mat_t2.x_21_1(0, 0) / (R_M + state_t2.GEO_eb(Height)),
    //         -kal_mat_t2.x_21_1(1, 0) * tan(state_t2.GEO_eb(Latitude)) / (R_N + state_t2.GEO_eb(Height));
    //     q_nc = SpinVector2Quart(delta_theta).conjugate();// Eun-Hwan. 式(3.80b)
    //     state_t2.q_ne = state_t2.q_ne * q_nc;// Eun-Hwan. 式(3.80a)
    //     Quart2GEO(state_t2.q_ne, &latitude_fix, &longitude_fix);
    //     height_fix = state_t2.GEO_eb(Height) + kal_mat_t2.x_21_1(2, 0);
    //     state_t2.GEO_eb << latitude_fix, longitude_fix, height_fix;

    //     kal_mat_t2.x_21_1.block<3, 1>(0, 0).setZero();
    // }
    void intergrated::VelErrorFeedback()
    {
        Eigen::Vector3d v_temp;
        v_temp = kal_mat_t2.x_21_1.block<3, 1>(3, 0);
        state_t2.NED_vec = state_t2.NED_vec - v_temp;
        //delta_v = delta_v - v_temp;
        kal_mat_t2.x_21_1.block<3, 1>(3, 0).setZero();
    }
    void intergrated::IMUErrorFeedback()
    {
        Eigen::Vector3d I_3_1;
        Eigen::Vector3d gb_3_1;
        Eigen::Vector3d gs_3_1;
        Eigen::Vector3d ab_3_1;
        Eigen::Vector3d as_3_1;

        gb_3_1 = kal_mat_t2.x_21_1.block<3, 1>(9, 0);
        ab_3_1 = kal_mat_t2.x_21_1.block<3, 1>(12, 0);
        gs_3_1 = kal_mat_t2.x_21_1.block<3, 1>(15, 0);
        as_3_1 = kal_mat_t2.x_21_1.block<3, 1>(18, 0);

        I_3_1.setOnes();
        // imu_parameter.g_b = imu_parameter.g_b + ((I_3_1 + imu_parameter.g_b).array() * gb_3_1.array()).matrix();
        // imu_parameter.a_b = imu_parameter.a_b + ((I_3_1 + imu_parameter.a_b).array() * ab_3_1.array()).matrix();
        // imu_parameter.g_s = ((I_3_1 + imu_parameter.g_s).array() * (I_3_1 + gs_3_1).array()).matrix() - I_3_1;
        // imu_parameter.a_s = ((I_3_1 + imu_parameter.a_s).array() * (I_3_1 + as_3_1).array()).matrix() - I_3_1;

        imu_parameter.g_b += gb_3_1;
        imu_parameter.a_b += ab_3_1;
        imu_parameter.g_s += gs_3_1;
        imu_parameter.a_s += as_3_1;


        kal_mat_t2.x_21_1.block<3, 1>(9, 0).setZero();
        kal_mat_t2.x_21_1.block<3, 1>(12, 0).setZero();
        kal_mat_t2.x_21_1.block<3, 1>(15, 0).setZero();
        kal_mat_t2.x_21_1.block<3, 1>(18, 0).setZero();
    }

    int intergrated::IntergtatedNavUpdate(const IMU_data &imudata, const GNSS_data_T &gnssdata)
    {
        int GNSS_flag = 0;
        if (imudata.timestamp == gnssdata.timestamp)
        {
            GNSS_flag = 1;
            GNSS_data = gnssdata;
            mechanizationinterUpdate(IMUdataFix(imudata));
            EKFpredictUpdate();
            EKFmeasurementUpdate();
            ErrorFeedback();

            imudata_t1 = imudata_t2;
            state_t1 = state_t2;
            kal_mat_t1 = kal_mat_t2;
        }
        else if ((imudata_t2.timestamp < gnssdata.timestamp) && (imudata.timestamp > gnssdata.timestamp))
        {
            GNSS_flag = 1;
            GNSS_data = gnssdata;
            double ratio;
            IMU_data imudatamid;
            ratio = (gnssdata.timestamp - imudata_t1.timestamp) / (imudata.timestamp - imudata_t1.timestamp);
            imudatamid.timestamp = gnssdata.timestamp;
            imudatamid.accel = ratio * imudata.accel;
            imudatamid.gyro = ratio * imudata.gyro;
            INS_T = imudatamid.timestamp - imudata_t2.timestamp;
            mechanizationinterUpdate(IMUdataFix(imudatamid));
            EKFpredictUpdate();
            EKFmeasurementUpdate();
            ErrorFeedback();
            imudata_t1 = imudata_t2;
            state_t1 = state_t2;
            kal_mat_t1 = kal_mat_t2;
            imudatamid.timestamp = imudata.timestamp;
            imudatamid.accel = (1.0 - ratio) * imudata.accel;
            imudatamid.gyro = (1.0 - ratio) * imudata.gyro;
            INS_T = imudatamid.timestamp - imudata_t2.timestamp;
            mechanizationinterUpdate(IMUdataFix(imudatamid));
            EKFpredictUpdate();
            imudata_t1 = imudata_t2;
            state_t1 = state_t2;
            kal_mat_t1 = kal_mat_t2;

        }
        else
        {
            mechanizationinterUpdate(IMUdataFix(imudata));
            EKFpredictUpdate();
            imudata_t1 = imudata_t2;
            state_t1 = state_t2;
            kal_mat_t1 = kal_mat_t2;
        }
        // else if ((imudata.timestamp < gnssdata.timestamp) && (imudata.timestamp + 0.005 > gnssdata.timestamp))
        // {
        //     GNSS_flag = 1;
        //     GNSS_data = gnssdata;
        //     mechanizationinterUpdate(IMUdataFix(imudata));
        //     EKFpredictUpdate();
        //     kal_mat_t1 = kal_mat_t2;
        //     stateExtrapolation(GNSS_data.timestamp);
        //     EKFpredictUpdate();
        //     EKFmeasurementUpdate();
        //     ErrorFeedback();

        //     state_t1 = state_t2;
        //     kal_mat_t1 = kal_mat_t2;
        // }
        // else if (imudata.timestamp + 0.005 <= gnssdata.timestamp)
        // {
        //     if((imudata_t2.timestamp == gnssdata.timestamp - 1)&&(imudata.timestamp - imudata_t2.timestamp < 0.0049))
        //     {
        //         stateExtrapolation1(imudata);
        //         imudata_t1 = imudata_t2;
        //         state_t1 = state_t2;
        //     }
        //     else
        //     {
        //         mechanizationinterUpdate(IMUdataFix(imudata));
        //         EKFpredictUpdate();

        //         imudata_t1 = imudata_t2;
        //         state_t1 = state_t2;
        //         kal_mat_t1 = kal_mat_t2;
        //     }

        // }

        return GNSS_flag;
    }
    void intergrated::stateExtrapolation(const double &timestamp)
    {
        m_State state_ext;
        IMU_data imudata_ext;
        IMU_data imudata_ext1;
        double GNSS_imu_time;
        double ratio;
        double height_ext;
        double latitude_ext;
        double longitude_ext;
        Eigen::Quaterniond q_ne_delta_theta;
        Eigen::Quaterniond q_ne_delta_theta_ratio;
        Vector3d ne_delta_theta;
        Eigen::Quaterniond q_bn_delta_theta;
        Eigen::Quaterniond q_bn_delta_theta_ratio;
        Vector3d bn_delta_theta;

        GNSS_imu_time = timestamp - imudata_t2.timestamp;
        ratio = GNSS_imu_time / INS_T;
        imudata_ext.timestamp = GNSS_data.timestamp;
        imudata_ext.accel = ratio * imudata_t2.accel;
        imudata_ext.gyro = ratio * imudata_t2.gyro;

        state_ext.timestamp = GNSS_data.timestamp;
        q_ne_delta_theta = state_t1.q_ne.inverse() * state_t2.q_ne;
        ne_delta_theta = Quart2SpinVector(q_ne_delta_theta);
        q_ne_delta_theta_ratio = SpinVector2Quart(ratio * ne_delta_theta);
        state_ext.q_ne = state_t2.q_ne * q_ne_delta_theta_ratio;
        Quart2GEO(state_ext.q_ne, &latitude_ext, &longitude_ext);
        height_ext = state_t2.GEO_eb(Height) - state_t2.NED_vec(D) * GNSS_imu_time;
        state_ext.GEO_eb << latitude_ext, longitude_ext, height_ext;
        state_ext.NED_vec = state_t2.NED_vec + GNSS_imu_time * a_t2;
        q_bn_delta_theta = state_t1.q_bn.inverse() * state_t2.q_bn;
        bn_delta_theta = Quart2SpinVector(q_bn_delta_theta);
        q_bn_delta_theta_ratio = SpinVector2Quart(ratio * bn_delta_theta);
        state_ext.q_bn = state_t2.q_bn * q_bn_delta_theta_ratio;
        state_ext.e_bn = DCM2Euler(state_ext.q_bn.toRotationMatrix());

        state_t1 = state_t2;
        imudata_t1 = imudata_t2;
        a_t2 = (state_ext.NED_vec - state_t2.NED_vec) / (GNSS_imu_time);

        state_t2 = state_ext;
        imudata_t2 = imudata_ext;
        INS_T = GNSS_imu_time;
        // imudata_ext1 = imudata_t2;
        // mechanizationinterUpdate(imudata_ext);
        // imudata_t1 = imudata_ext1;
        // state_t1 = state_t2;

        
    }
    void intergrated::stateExtrapolation1(const IMU_data &imudata)
    {
        m_State state_ext;
        IMU_data imudata_ext;
        double GNSS_imu_time;
        double ratio;
        double height_ext;
        double latitude_ext;
        double longitude_ext;
        Eigen::Quaterniond q_ne_delta_theta;
        Eigen::Quaterniond q_ne_delta_theta_ratio;
        Vector3d ne_delta_theta;
        Eigen::Quaterniond q_bn_delta_theta;
        Eigen::Quaterniond q_bn_delta_theta_ratio;
        Vector3d bn_delta_theta;

        GNSS_imu_time = imudata.timestamp - imudata_t2.timestamp;
        ratio = GNSS_imu_time / (imudata.timestamp - imudata_t1.timestamp);
        imudata_ext.timestamp = GNSS_data.timestamp;
        imudata_ext.accel = ratio * imudata.accel;
        imudata_ext.gyro = ratio * imudata.gyro;

        state_ext.timestamp = GNSS_data.timestamp;
        q_ne_delta_theta = state_t1.q_ne.inverse() * state_t2.q_ne;
        ne_delta_theta = Quart2SpinVector(q_ne_delta_theta);
        q_ne_delta_theta_ratio = SpinVector2Quart(ratio * ne_delta_theta);
        state_ext.q_ne = state_t2.q_ne * q_ne_delta_theta_ratio;
        Quart2GEO(state_ext.q_ne, &latitude_ext, &longitude_ext);
        height_ext = state_t2.GEO_eb(Height) - state_t2.NED_vec(D) * GNSS_imu_time;
        state_ext.GEO_eb << latitude_ext, longitude_ext, height_ext;
        state_ext.NED_vec = state_t2.NED_vec + GNSS_imu_time * a_t2;
        q_bn_delta_theta = state_t1.q_bn.inverse() * state_t2.q_bn;
        bn_delta_theta = Quart2SpinVector(q_bn_delta_theta);
        q_bn_delta_theta_ratio = SpinVector2Quart(ratio * bn_delta_theta);
        state_ext.q_bn = state_t2.q_bn * q_bn_delta_theta_ratio;
        state_ext.e_bn = DCM2Euler(state_ext.q_bn.toRotationMatrix());

        // state_t1 = state_t2;
        // imudata_t1 = imudata_t2;
        // //delta_v = state_ext.NED_vec - state_t2.NED_vec;

        // state_t2 = state_ext;
        // imudata_t2 = imudata;
        // INS_T = GNSS_imu_time;

        mechanizationinterUpdate(IMUdataFix(imudata_ext));
        imudata_t1 = imudata_t2;
        state_t1 = state_t2;

        
    }
    IMU_data intergrated::IMUdataFix(const IMU_data &data)
    {
        IMU_data data_fixed;
        Eigen::Vector3d I_3_1;
        INS_T = data.timestamp - imudata_t2.timestamp;
        I_3_1.setOnes();
        data_fixed.timestamp = data.timestamp;
        data_fixed.accel = (((data.accel - imu_parameter.a_b * INS_T).array()) / ((I_3_1 + imu_parameter.a_s).array())).matrix();
        data_fixed.gyro = (((data.gyro - imu_parameter.g_b * INS_T).array()) / ((I_3_1 + imu_parameter.g_s).array())).matrix();
        return data_fixed;
    }
}
