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
		//Eigen::Vector3d gyro_t2;				//陀螺仪tk输出
		//Eigen::Vector3d gyro_t1;				//陀螺仪tk-1输出
		//Eigen::Vector3d accelerometer_t2;		//加速度计tk输出
		//Eigen::Vector3d accelerometer_t1;		//加速度计tk-1输出
		//	姿态
		//Eigen::Quaterniond q_bn_t1_t1;			//tk到tk-1时n系变换的四元数
		//Eigen::Quaterniond q_bn_t2_t2;			//tk到tk-1时n系变换的四元数
		//	位置
		//Eigen::Vector3d GEO_eb;					//位置的地理坐标b系相对于e系的向量 纬度 经度 高度
		//Eigen::Vector3d GEO_eb_t1;				//位置的地理坐标b系相对于e系的向量 纬度 经度 高度
		Eigen::Vector3d GEO_eb_mid;				//n系下中间时刻的位置向量坐标
		//Eigen::Quaterniond q_ne_t2_t2;			//tk时刻的ne四元数
		//Eigen::Quaterniond q_ne_t1_t1;			//tk-1时刻的ne四元数
		// 速度
		//Eigen::Vector3d NED_vec;				//载体相对于地球的速度在NED坐标系的投影
		//Eigen::Vector3d NED_vec_t1;				//n系下tk-1时刻的速度向量坐标
		Eigen::Vector3d NED_vec_mid;			//n系下中间时刻的速度向量坐标
		Eigen::Vector3d delta_v;
		Eigen::Vector3d omega_ie_n_mid;			//n系下中间时刻计算的ie角速度
		Eigen::Vector3d omega_en_n_mid;			//n系下中间时刻计算的en角速度
		// 
		m_State state_t2;
		m_State state_t1;
		//
		IMU_data imudata_t1;
		IMU_data imudata_t2;
	};
};