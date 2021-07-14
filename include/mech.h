#pragma once
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include "mech_utils.hpp"

using namespace std;


namespace INS {

	class dataReader
	{
	public:
		dataReader();
		~dataReader();
	private:
		string imu_file_name;
		string ref_file_name;
	public:
		//double imu_time[5000];
		//double gyroscope_x[5000];
		//double gyroscope_y[5000];
		//double gyroscope_z[5000];
		//double accelerometer_x[5000];
		//double accelerometer_y[5000];
		//double accelerometer_z[5000];
		IMU_data imudata;
		Ref_data refdata;
		ifstream ref_file_handle;
		ifstream imu_file_handle;
		//void imuRead1ms();	
		void imuRead();
		void refRead();
	};



	class mechanization
	{
	public:
		mechanization(m_State state, IMU_data data);
		~mechanization();
	private:
		void velocityUpdate();
		void positionUpdate();
		void attitudeUpdate();
	public:
		void mechanizationUpdate(IMU_data& data);
		m_State getstate();
	private:
		//Eigen::Vector3d gyro_t2;				//������tk���
		//Eigen::Vector3d gyro_t1;				//������tk-1���
		//Eigen::Vector3d accelerometer_t2;		//���ٶȼ�tk���
		//Eigen::Vector3d accelerometer_t1;		//���ٶȼ�tk-1���
		//	��̬
		//Eigen::Quaterniond q_bn_t1_t1;			//tk��tk-1ʱnϵ�任����Ԫ��
		//Eigen::Quaterniond q_bn_t2_t2;			//tk��tk-1ʱnϵ�任����Ԫ��
		//	λ��
		//Eigen::Vector3d GEO_eb;					//λ�õĵ�������bϵ�����eϵ������ γ�� ���� �߶�
		//Eigen::Vector3d GEO_eb_t1;				//λ�õĵ�������bϵ�����eϵ������ γ�� ���� �߶�
		Eigen::Vector3d GEO_eb_mid;				//nϵ���м�ʱ�̵�λ����������
		//Eigen::Quaterniond q_ne_t2_t2;			//tkʱ�̵�ne��Ԫ��
		//Eigen::Quaterniond q_ne_t1_t1;			//tk-1ʱ�̵�ne��Ԫ��
		// �ٶ�
		//Eigen::Vector3d NED_vec;				//��������ڵ�����ٶ���NED����ϵ��ͶӰ
		//Eigen::Vector3d NED_vec_t1;				//nϵ��tk-1ʱ�̵��ٶ���������
		Eigen::Vector3d NED_vec_mid;			//nϵ���м�ʱ�̵��ٶ���������
		Eigen::Vector3d delta_v;
		Eigen::Vector3d omega_ie_n_mid;			//nϵ���м�ʱ�̼����ie���ٶ�
		Eigen::Vector3d omega_en_n_mid;			//nϵ���м�ʱ�̼����en���ٶ�
		// 
		m_State state_t2;
		m_State state_t1;
		//
		IMU_data imudata_t1;
		IMU_data imudata_t2;
	};
};