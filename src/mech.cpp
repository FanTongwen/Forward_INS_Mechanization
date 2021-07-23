#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include "mech.h"
#include "mech_utils.hpp"
namespace INS {
	dataReader::dataReader()
	{
		imu_file_name = "../data/IMU.bin";
		ref_file_name = "../data/Reference.bin";

		imu_file_handle.open(imu_file_name, ios::binary);
		ref_file_handle.open(ref_file_name, ios::binary);
	}

	dataReader::~dataReader()
	{
		imu_file_handle.close();
		ref_file_handle.close();
	}

	//void dataReader::imuRead1ms()
	//{
	//	double tempdata[35000];
	//	file_handle.read((char *)&tempdata, 35000 * sizeof(double));
	//	for (int i = 0; i < 5000; i++)
	//	{
	//		imu_time[i] = tempdata[7 * i];
	//		gyroscope_x[i] = tempdata[7 * i + 1];
	//		gyroscope_y[i] = tempdata[7 * i + 2];
	//		gyroscope_z[i] = tempdata[7 * i + 3];
	//		accelerometer_x[i] = tempdata[7 * i + 4];
	//		accelerometer_y[i] = tempdata[7 * i + 5];
	//		accelerometer_z[i] = tempdata[7 * i + 6];
	//	}
	//}

	void dataReader::imuRead()
	{
		imu_file_handle.read((char *)&imudata, sizeof(imudata));
	}

	void dataReader::refRead()
	{
		ref_file_handle.read((char *)&refdata, sizeof(refdata));
	}

	mechanization::mechanization(m_State state, IMU_data data) : state_t2(state), imudata_t1(data), imudata_t2(data)
	{
		//��ʼʱ�䣬��ʼλ�ã���ʼ�ٶȣ���ʼ��̬
		//��ʼʱ���imu���
		state_t2.q_ne = NED2Quart(state_t2.GEO_eb);
		state_t2.q_bn = Euler2Quart(state_t2.e_bn);
		state_t1 = state_t2;
		a_t2 << 0.0, 0.0, 0.0;
	}

	mechanization::~mechanization()
	{
		;
	}

	void mechanization::attitudeUpdate()
	{
		//double R_M;		//����Ȧ�뾶
		//double R_N;		//î��Ȧ�뾶
		//Eigen::Vector3d Omega_ie_n;
		//Eigen::Vector3d Omega_en_n;
		//Eigen::Vector3d Zeta_tk;
		//Eigen::Vector3d Phi_tk;
		//Eigen::Quaterniond Quat_b_t1t2;
		//R_M = EarthLongAxis * (1.0 - pow2Func(EarthECC)) / sqrt(pow3Func(1.0 - pow2Func(EarthECC*sin(GEO_eb(0)*D2R))));
		//R_N = EarthLongAxis / sqrt(1.0 - pow2Func(EarthECC * sin(GEO_eb(0)*D2R)));
		//Omega_ie_n << EarthAngleVelocity * cos(GEO_eb(0)*D2R), 0.0, -EarthAngleVelocity * sin(GEO_eb(0)*D2R);
		//Omega_en_n << NED_vec(1) / (R_N + GEO_eb(2)),
		//	-NED_vec(0) / (R_M + GEO_eb(2)),
		//	-NED_vec(1)*tan(GEO_eb(0)*D2R) / (R_N + GEO_eb(2));
		//Zeta_tk = (Omega_ie_n + Omega_en_n)*INS_T;
		//Quat_n_t2t1.w() = cos(0.5*Zeta_tk.norm());
		//Quat_n_t2t1.vec() = (-sin(0.5*Zeta_tk.norm()) / (0.5*Zeta_tk.norm()))*0.5*Zeta_tk;

		//Phi_tk = gyro_t2 + (gyro_t1.cross(gyro_t2)) / 12.0;
		//Quat_b_t1t2.w() = cos(0.5*Phi_tk.norm());
		//Quat_b_t1t2.vec() = (sin(0.5*Phi_tk.norm()) / (0.5*Phi_tk.norm()))*0.5*Phi_tk;
		//Quat_nb = Quat_n_t2t1 * Quat_nb * Quat_b_t1t2;

		//Eigen::Vector3d xi_t2;
		//Eigen::Vector3d zeta1_t2;
		//Eigen::Quaterniond q_nn_t2_t1;
		//Eigen::Quaterniond q_ee_t1_t2;
		Eigen::Quaterniond q_delta_theta;		//��tk-1��tk��neϵλ�ñ仯
		Eigen::Vector3d delta_theta;			//q_delta_theta��Ӧ�ĵ�Ч��תʸ��
		Eigen::Quaterniond q_delta_halftheta;	//0.5theta����תʸ����Ӧ����Ԫ��
		Eigen::Quaterniond q_ne_mid_mid;		//tk-0.5��Ӧ��ne��Ԫ��
		double height_mid;						//����λ���ڲ���м�ʱ��λ��
		double longitude_mid;					//ͨ����Ԫ���ڲ���м�ʱ�̾���
		double latitude_mid;					//ͨ����Ԫ���ڲ���м�ʱ��γ��
		Eigen::Vector3d GEO_eb_mid;				//nϵ���м�ʱ�̵�λ����������
		Eigen::Vector3d NED_vec_mid;			//nϵ���м�ʱ�̵��ٶ���������
		Eigen::Vector3d omega_ie_n_mid;			//nϵ���м�ʱ�̼����ie���ٶ�
		Eigen::Vector3d omega_en_n_mid;			//nϵ���м�ʱ�̼����en���ٶ�
		Eigen::Vector3d zeta_t2;				//q_nn_t1_t2�Ĺ�������Ӧ�ĵ�Ч��תʸ��
		Eigen::Quaterniond q_nn_t1_t2;			//tk-1��tk��nϵ�仯��Ԫ��
		Eigen::Vector3d phi_t2;					//q_bb_t2_t1��Ӧ�ĵ�Ч��תʸ��
		Eigen::Quaterniond q_bb_t2_t1;			//tk��tk-1��bϵ�仯��Ԫ��


		// ���� state_t2.q_ne
		//NED_vec_mid = (state_t2.NED_vec + state_t1.NED_vec) / 2.0;
		//omega_en_n_mid = PosVec2AngleVelocity_en(GEO_eb_mid, NED_vec_mid);
		//xi_t2 = omega_ie_e * INS_T;
		//zeta1_t2 = (omega_ie_n_mid + omega_en_n_mid)*INS_T;
		//q_nn_t2_t1 = SpinVector2Quart(zeta1_t2);
		//q_ee_t1_t2 = (SpinVector2Quart(xi_t2)).conjugate();
		//state_t2.q_ne = q_ee_t1_t2 * state_t1.q_ne * q_nn_t2_t1;
		// ����λ�ø��¼����м�ʱ��λ��
		//q_delta_theta = q_ne_t1_t1.inverse() * q_ne_t2_t2;//Eun-Hwan. ʽ(2.62)
		q_delta_theta = state_t1.q_ne.inverse() * state_t2.q_ne;//Eun-Hwan. ʽ(2.62)
		delta_theta = Quart2SpinVector(q_delta_theta);
		q_delta_halftheta = SpinVector2Quart(0.5 * delta_theta);
		//q_ne_mid_mid = q_ne_t1_t1 * q_delta_halftheta;//Eun-Hwan. ʽ(2.63)
		q_ne_mid_mid = state_t1.q_ne * q_delta_halftheta;//Eun-Hwan. ʽ(2.63)
		// ��t=k-0.5ʱ�̵�λ�á��ٶ�
		Quart2GEO(q_ne_mid_mid, &latitude_mid, &longitude_mid);
		height_mid = (state_t2.GEO_eb(Height) + state_t1.GEO_eb(Height)) / 2.0;
		GEO_eb_mid << latitude_mid, longitude_mid, height_mid;
		NED_vec_mid = (state_t2.NED_vec + state_t1.NED_vec) / 2.0;
		omega_ie_n_mid = Pos2AngleVelocity_ie(GEO_eb_mid);
		omega_en_n_mid = PosVec2AngleVelocity_en(GEO_eb_mid, NED_vec_mid);

		zeta_t2 = (omega_ie_n_mid + omega_en_n_mid)*INS_T;
		q_nn_t1_t2 = (SpinVector2Quart(zeta_t2)).conjugate();//Eun-Hwan. ʽ(2.61)

		phi_t2 = imudata_t2.gyro + (imudata_t1.gyro.cross(imudata_t2.gyro)) / 12.0;//Eun-Hwan. ʽ(2.60)
		q_bb_t2_t1 = SpinVector2Quart(phi_t2);//Eun-Hwan. ʽ(2.57)
		// ��̬����
		//q_bn_t1_t1 = q_bn_t2_t2;
		state_t2.q_bn = q_nn_t1_t2 * state_t1.q_bn * q_bb_t2_t1;//Eun-Hwan. ʽ(2.56a),(2.56b)
		state_t2.q_bn.normalize();//��һ����������ת������ܲ�����
		state_t2.e_bn = DCM2Euler(((state_t2.q_bn).matrix()));
	}

	void mechanization::velocityUpdate()
	{
		Eigen::Vector3d Delta_v_fk_bk1;

		Eigen::Vector3d omega_ie_n;				//eϵ�����iϵ�Ľ��ٶ���nϵ��ͶӰ
		Eigen::Vector3d omega_en_n;				//eϵ�����iϵ�Ľ��ٶ���nϵ��ͶӰ

		double h_mid;							//�м�ʱ�̸߶�
		double longitude_mid;					//�м�ʱ�̾���
		double latitude_mid;					//�м�ʱ��γ��
		Eigen::Vector3d NED_vec_mid;			//nϵ���м�ʱ�̵��ٶ���������
		Eigen::Vector3d GEO_eb_mid;				//nϵ���м�ʱ�̵�λ����������
		Eigen::Vector3d omega_ie_n_mid;			//nϵ���м�ʱ�̼����ie���ٶ�
		Eigen::Vector3d omega_en_n_mid;			//nϵ���м�ʱ�̼����en���ٶ�

		Eigen::Vector3d zeta_mid;
		Eigen::Vector3d xi_mid;
		Eigen::Quaterniond q_nn_mid_t1;
		Eigen::Quaterniond q_ee_t1_mid;
		Eigen::Quaterniond q_ne_mid_mid;		//�м�ʱ�̵�neϵ��Ԫ�������ڼ����м�ʱ�̾�γ��
		Eigen::Vector3d zeta_t2;
		Eigen::Matrix3d C_bn_t1;				//��һʱ�̵�bn��̬�������Ҿ���
		Eigen::Vector3d delta_v_f;				//����������ٶ�����
		// ���������Ĳ��� �Ϲ��� �����ߵ��㷨����ϵ���ԭ����
		//ʽ(3-2-19) 
		//ʽ(3-2-23) 
		//ʽ(3-2-24) 
		//ʽ(3-2-25)
		double g_e = 9.780325;					//�������
		double beta = 5.30240E-3;				//������������
		double beta_1 = 5.82E-6;
		double beta_2 = 3.08279729293E-6;
		double beta_3 = 8.08E-9;
		double g_Lh;							//������ģֵ
		double true_geo_angle;					//�洹�ߺ͵����ߵļн�
		Eigen::Vector3d g_n_mid;				//����ʸ����nϵͶӰ

		Eigen::Vector3d delta_v_gcor;			//�����Ϳ�����������ٶ�����

		// Eun-Hwan. ʽ(2.48b) �����������ٶȸı���м����
		//Delta_v_fk_bk1 = accelerometer_t2 + 0.5*gyro_t2.cross(accelerometer_t2) + 
		//	(gyro_t1.cross(accelerometer_t2) + accelerometer_t1.cross(gyro_t2))/12.0;
		Delta_v_fk_bk1 = imudata_t2.accel + 0.5*imudata_t2.gyro.cross(imudata_t2.accel) +
			(imudata_t1.gyro.cross(imudata_t2.accel) + imudata_t1.accel.cross(imudata_t2.gyro)) / 12.0;
		// Eun-Hwan. ʽ(2.28) ʽ(2.30) ��t=k-1ʱ�̵�ie��en���ٶ� omega_ie_n omega_en_n
		omega_ie_n = Pos2AngleVelocity_ie(state_t1.GEO_eb);
		omega_en_n = PosVec2AngleVelocity_en(state_t1.GEO_eb, state_t1.NED_vec);
		// ��t=k-0.5ʱ�̵�λ�á��ٶ�
		zeta_mid = (omega_ie_n + omega_en_n) * 0.5 * INS_T;
		q_nn_mid_t1 = SpinVector2Quart(zeta_mid);//Eun-Hwan. ʽ(2.51c)
		xi_mid = omega_ie_e * 0.5 * INS_T;
		q_ee_t1_mid = (SpinVector2Quart(xi_mid)).conjugate();//Eun-Hwan. ʽ(2.51d)
		// q_ne_t1_t1 = NED2Quart(GEO_eb);
		q_ne_mid_mid = q_ee_t1_mid * state_t1.q_ne * q_nn_mid_t1;//Eun-Hwan. ʽ(2.51a),ʽ(2.51b)
		Quart2GEO(q_ne_mid_mid, &latitude_mid, &longitude_mid);//�����м�ʱ�̵ľ��Ⱥ�γ��
		h_mid = state_t1.GEO_eb(Height) - state_t1.NED_vec(D)*INS_T*0.5;//�����м�ʱ�̵ĸ߶�
		NED_vec_mid = state_t1.NED_vec + 0.5 * INS_T * a_t2;//���Ƶ��м�ʱ���ٶ� Eun-Hwan. ʽ(2.52a)
		GEO_eb_mid << latitude_mid, longitude_mid, h_mid;//���м�ʱ��λ�ñ���Ϊvector3d��ʽ
		// ��t=k-0.5ʱ�̶�Ӧλ���ٶȵ�ie���ٶȡ�en���ٶ�
		omega_ie_n_mid = Pos2AngleVelocity_ie(GEO_eb_mid);
		omega_en_n_mid = PosVec2AngleVelocity_en(GEO_eb_mid, NED_vec_mid);
		// ��ʽ2.48a
		zeta_t2  = (omega_ie_n_mid + omega_en_n_mid)*INS_T;
		//C_bn_t1 = q_bn_t2_t2.toRotationMatrix();
		C_bn_t1 = state_t1.q_bn.toRotationMatrix();
		delta_v_f = (Eigen::Matrix3d::Identity() - 0.5*(Vector2CrossMatrix(zeta_t2)))*
			C_bn_t1*Delta_v_fk_bk1;//Eun-Hwan. ʽ(2.49a)
		// ��t=k-0.5������ g_n_mid
		g_Lh = g_e * (1.0 + beta * pow2Func(sin(latitude_mid)) - beta_1 * pow2Func(sin(2.0*latitude_mid))) - beta_2 * h_mid;
		true_geo_angle = (beta_3 * h_mid * sin(2.0 * latitude_mid)) / g_Lh;
		//g_n_mid << -g_Lh * sin(true_geo_angle), 0, g_Lh * cos(true_geo_angle);//�Ϲ���. ʽ(3-2-25)
		//g_n_mid << -g_Lh * sin(true_geo_angle), 0, g_Lh;//�Ϲ���. ʽ(3-2-25)
		//g_n_mid = 0.99999999973*getG(Eigen::Vector3d(latitude_mid, longitude_mid, h_mid));
		g_n_mid = getG(Eigen::Vector3d(latitude_mid, longitude_mid, h_mid));
		// �����м�ʱ��������ʽ(2.53)
		delta_v_gcor = (g_n_mid - (2.0 * omega_ie_n_mid + omega_en_n_mid).cross(NED_vec_mid)) * INS_T;//Eun-Hwan. ʽ(2.53)
		// �ٶȸ��� ʽ2.47
		//NED_vec_t1 = NED_vec;
		a_t2 = (delta_v_f + delta_v_gcor) / INS_T;//Eun-Hwan. ʽ(2.52b)
		state_t2.NED_vec = state_t1.NED_vec + delta_v_f + delta_v_gcor;//Eun-Hwan. ʽ(2.47)

	}

	void mechanization::positionUpdate()
	{
		/*
		Eigen::Vector3d xi_t2;
		Eigen::Vector3d zeta_t2;
		Eigen::Quaterniond q_nn_t2_t1;
		Eigen::Quaterniond q_ee_t1_t2;
		double longitude_t2, latitude_t2, height_t2;
		NED_vec_mid = (state_t2.NED_vec + state_t1.NED_vec) / 2.0;
		omega_en_n_mid = PosVec2AngleVelocity_en(GEO_eb_mid, NED_vec_mid);
		xi_t2 = omega_ie_e * INS_T;
		zeta_t2 = (omega_ie_n_mid + omega_en_n_mid)*INS_T;
		q_nn_t2_t1 = SpinVector2Quart(zeta_t2);
		q_ee_t1_t2 = (SpinVector2Quart(xi_t2)).conjugate();
		//q_ne_t1_t1 = q_ne_t2_t2;
		//q_ne_t2_t2 = q_ee_t1_t2 * q_ne_t1_t1 * q_nn_t2_t1;
		state_t2.q_ne = q_ee_t1_t2 * state_t1.q_ne * q_nn_t2_t1;
		//Quart2GEO(q_ne_t2_t2, &latitude_t2, &longitude_t2);
		Quart2GEO(state_t1.q_ne, &latitude_t2, &longitude_t2);
		height_t2 = state_t1.GEO_eb(Height) - NED_vec_mid(D)*INS_T;
		// λ�ø���
		//GEO_eb_t1 = GEO_eb;
		state_t2.GEO_eb << latitude_t2, longitude_t2, height_t2;
		*/
		double R_M_t1;
		double R_N_mid;
		double height_mid;
		double latitude_mid;

		state_t2.GEO_eb(Height) = state_t1.GEO_eb(Height) - 0.5 * (state_t2.NED_vec(D) + state_t1.NED_vec(D)) * INS_T;// �����. lecture6 (54)
		height_mid = 0.5*(state_t2.GEO_eb(Height) + state_t1.GEO_eb(Height));
		R_M_t1 = EarthLongAxis * (1.0 - EarthECC2) / sqrt(pow3Func(1.0 - EarthECC2 * pow2Func(sin(state_t1.GEO_eb(Latitude)))));
		state_t2.GEO_eb(Latitude) = state_t1.GEO_eb(Latitude) + (state_t2.NED_vec(N) + state_t1.NED_vec(N)) / (2.0*(R_M_t1 + height_mid)) * INS_T;// �����. lecture6 (55)
		latitude_mid = 0.5*(state_t2.GEO_eb(Latitude) + state_t1.GEO_eb(Latitude));
		R_N_mid = EarthLongAxis / sqrt(1.0 - EarthECC2 * pow2Func(sin(latitude_mid)));
		state_t2.GEO_eb(Longitude) = state_t1.GEO_eb(Longitude) + (state_t2.NED_vec(E) + state_t1.NED_vec(E)) / (2.0*(R_N_mid + height_mid) * cos(latitude_mid)) * INS_T;// �����. lecture6 (56)
		state_t2.q_ne = NED2Quart(state_t2.GEO_eb);
	}

	void mechanization::mechanizationUpdate(const IMU_data& data)
	{
		imudata_t2 = data;
		INS_T = imudata_t2.timestamp - imudata_t1.timestamp;
		velocityUpdate();
		positionUpdate();
		attitudeUpdate();
		state_t2.timestamp = imudata_t2.timestamp;

		imudata_t1 = imudata_t2;
		state_t1 = state_t2;
	}

	void mechanization::mechanizationinterUpdate(const IMU_data& data)
	{
		imudata_t2 = data;
		INS_T = imudata_t2.timestamp - imudata_t1.timestamp;
		velocityUpdate();
		positionUpdate();
		attitudeUpdate();
		state_t2.timestamp = imudata_t2.timestamp;
	}
	m_State mechanization::getstate()
	{
		return state_t2;
	}
}
//void main(int a)
//{
//	
//	Eigen::MatrixX2d m(2, 2);
//	Eigen::MatrixX2d n(2, 2);
//	Eigen::MatrixX2d mn(2, 2);
//	Eigen::Vector3d v3d(1, 2, 3);
//	Eigen::Vector3d v3da(1, 2, 3);
//	Eigen::Vector3d v3db(1, 3, 3);
//	Eigen::Quaterniond qd(1, 2, 3, 4);
//	Eigen::Quaterniond qda(1, 1, 1, 2);
//	Eigen::Quaterniond qdb(1, 1, 1, 1);
//	m << 1, 2, 3, 4;
//	n << 5, 6, 7, 8;
//	mn = m * n;
//	std::cout << v3d << std::endl << m << std::endl << n << std::endl;
//	std::cout << mn(0, 0) << std::endl;
//	std::cout << v3d(0) << v3d(1) << v3d(2) << v3d(2) << std::endl;
//	v3d << 3, 7, 56;
//	std::cout << v3d << std::endl;
//	v3d << 1, 5, 7; //ֻ��ͬʱ��ֵ3��
//	std::cout << v3d << std::endl;
//	// ��Ԫ������
//	std::cout << "��Ԫ������" << std::endl;
//	std::cout << qd.coeffs() << std::endl;
//	std::cout << (qd.conjugate()).coeffs() << std::endl;//����
//	qd.coeffs() << 3, 2, 1, 4;
//	std::cout << qd.coeffs() << std::endl;
//	qd.vec() << v3d;
//	v3d << 1, 2, 3;
//	std::cout << qd.coeffs() << std::endl;
//	qd.vec() = v3d;
//	std::cout << qd.coeffs() << std::endl;
//
//	std::cout << std::endl;
//	std::cout << qda.coeffs() << std::endl;
//	std::cout << qdb.coeffs() << std::endl;
//	qd = qda * qdb;
//	std::cout << qd.coeffs() << std::endl;
//	// ��˲���
//	std::cout << "��˲���" << std::endl;
//	std::cout << v3da << std::endl;
//	std::cout << v3db << std::endl;
//	v3d = v3da.cross(v3db);
//	std::cout << v3d << std::endl;
//	std::cout << v3da << std::endl;
//	std::cout << v3db << std::endl;
//	std::cout << std::endl;
//	std::cout << Eigen::Matrix3d::Identity() << std::endl;
//	// ���Ǻ���
//	std::cout << "��˲���" << std::endl;
//	std::cout << atan(-1.0) << std::endl;
//	std::cout << atan(1.0) << std::endl;
//	system("pause");
//	
//	//imudataReader imureader;
//	//
//	////for (int j = 0; j < 20; j++)
//	////{
//	////	imureader.imuRead1ms();
//	////	for (int i = 0; i < 5000; i++)
//	////	{
//	////		//cout << setprecision(8) << imureader.imu_time[i] << endl;
//	////	}
//	////}
//
//	//
//	//system("pause");
//
//}

