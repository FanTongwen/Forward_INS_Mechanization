#pragma once
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

namespace INS{
	// WGS-84 spheroid parameters:
	#define		EarthLongAxis				6378137.0E0				//地球长半轴（即赤道横截面半径a）/m
	#define		EarthECC					0.0818191908426E0		//地球偏心率
	#define		EarthECC2					0.00669437999014E0		//地球椭球偏心率平方
	#define		EarthAngleVelocity			7.2921151467E-5			//地球自转角速度 /rad/s

	//常用的数学定义
	#define		PI							3.1415926535897932		/* pi */
	#define		pow2Func(a)					((a)*(a))
	#define		pow3Func(a)					((a)*(a)*(a))
	#define		D2R							(PI/180.0)				/* deg to rad */
	#define		R2D							(180.0/PI)				/* rad to deg */
	//INS参数
	#define		INS_SampleRate				200.0
	#define		INS_T						(1.0/INS_SampleRate)
	//INS初始参数
	#define		TIME_BEGIN					91620.0

	const double Ra = 6378137.0;            // 长轴
	const double Rb = 6356752.3142;         // 短轴
	const double E2 = 0.00669437999013;     // 第一偏心率e的平方 e=sqrt(a^2-b^2)/a
	const double GA = 9.7803267715;
	const double GB = 9.8321863685;
	const double F = 1.0 / 298.257223563;
	const double GM = 3.986005e14; // 万有引力常量G * 地球质量M
	enum
	{
		Latitude,
		Longitude,
		Height
	};

	enum
	{
		N,
		E,
		D
	};

	enum
	{
		Yaw,
		Pitch,
		Roll
	};
	typedef struct
	{
		double timestamp = 0;
		Eigen::Vector3d gyro;
		Eigen::Vector3d accel;
	}IMU_data;

	typedef struct
	{
		double timestamp;
		Eigen::Vector3d pos;
		Eigen::Vector3d vel;
		Eigen::Vector3d att;
	}Ref_data;

	typedef struct
	{
		double timestamp;
		Eigen::Quaterniond q_bn;//姿态
		Eigen::Vector3d GEO_eb;//位置
		Eigen::Vector3d NED_vec;//速度
		Eigen::Quaterniond q_ne;
		Eigen::Vector3d e_bn;//姿态 euler
	}m_State;

	const Eigen::Vector3d omega_ie_e(0.0, 0.0, EarthAngleVelocity);

	static Eigen::Matrix3d Vector2CrossMatrix(Eigen::Vector3d a)
	{
		Eigen::Matrix3d a_crossMatrix;
		a_crossMatrix << 0, -a(2), a(1),
						 a(2), 0, -a(0),
						 -a(1), a(0), 0;
		return a_crossMatrix;
	}

	static Eigen::Vector3d Pos2AngleVelocity_ie(Eigen::Vector3d Pos)
	{
		Eigen::Vector3d omega_ie;
		omega_ie << EarthAngleVelocity * cos(Pos(Latitude)),
			0.0,
			-EarthAngleVelocity * sin(Pos(Latitude));
		return omega_ie;
	}

	static Eigen::Vector3d PosVec2AngleVelocity_en(Eigen::Vector3d Pos, Eigen::Vector3d Vec)
	{
		Eigen::Vector3d omega_en;
		double R_M;
		double R_N;
		R_M = EarthLongAxis * (1.0 - EarthECC2) / sqrt(pow3Func(1.0 - EarthECC2 * pow2Func(sin(Pos(Latitude)))));
		R_N = EarthLongAxis / sqrt(1.0 - EarthECC2 * pow2Func(sin(Pos(Latitude))));
		omega_en << Vec(E) / (R_N + Pos(Height)),
			-Vec(N) / (R_M + Pos(Height)),
			-Vec(E)*tan(Pos(Latitude)) / (R_N + Pos(Height));
		return omega_en;
	}

	static Eigen::Quaterniond SpinVector2Quart(Eigen::Vector3d vec)
	{
		Eigen::Quaterniond quat;
		quat.w() = cos((0.5*vec).norm());
		quat.vec() = (sin((0.5*vec).norm()) / ((0.5*vec).norm()))*(0.5*vec);
		return quat;
	}

	static Eigen::Vector3d Quart2SpinVector(Eigen::Quaterniond quat)
	{
		Eigen::Vector3d vec;
		if (quat.w() == 0)
		{
			vec = PI * quat.vec();
		}
		else
		{
			double phi_half;
			double f;
			phi_half = atan2((quat.vec()).norm(), quat.w());
			f = sin(phi_half) / (2.0 * phi_half);
			vec = quat.vec() / f;
		}
		return vec;
	}

	static Eigen::Quaterniond NED2Quart(Eigen::Vector3d vec)
	{
		Eigen::Quaterniond quat;
		double temp_angle1, temp_angle2;

		temp_angle1 = -PI / 4.0 - vec(Latitude) / 2.0;
		temp_angle2 = vec(Longitude) / 2.0;
		quat.coeffs() <<	-sin(temp_angle1)*sin(temp_angle2),
							sin(temp_angle1)*cos(temp_angle2),
							cos(temp_angle1)*sin(temp_angle2),
							cos(temp_angle1)*cos(temp_angle2);
		return quat;
	}

	static void Quart2GEO(const Eigen::Quaterniond quat, double * lat, double * longi)
	{
		*longi = atan2(quat.z(), quat.w()) * 2.0;
		*lat = (atan2(quat.y(), quat.w()) * -2.0) - (PI * 0.5);
	}

	static Eigen::Quaterniond Euler2Quart(Eigen::Vector3d vec)
	{
		Eigen::Quaterniond quat;
		Eigen::AngleAxisd rollAngle(vec(Roll), Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd pitchAngle(vec(Pitch), Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd yawAngle(vec(Yaw), Eigen::Vector3d::UnitZ());
		quat = yawAngle * pitchAngle * rollAngle;
		return quat;
	}

	static Eigen::Vector3d EulerNorm(Eigen::Vector3d vec)
	{
		Eigen::Vector3d vec_norm;
		if (vec(Yaw) > PI / 2.0)
		{
			vec_norm(Yaw) = vec(Yaw) - PI;
			vec_norm(Pitch) = -vec(Pitch) - PI;
			vec_norm(Roll) = vec(Roll) + PI;
		}
		else
		{
			vec_norm = vec;
		}
		return vec_norm;
	}

	static Eigen::Vector3d DCM2Euler(Eigen::MatrixX3d dcm)
	{
		Eigen::Vector3d vec;
		vec <<	atan2(dcm(1, 0), dcm(0, 0)), 
				atan2(-dcm(2, 0), sqrt(dcm(2, 1) * dcm(2, 1) + dcm(2, 2) * dcm(2, 2))), 
				atan2(dcm(2, 1), dcm(2, 2));
		return vec;
	}

	// 计算当前位置重力加速度
	static Eigen::Vector3d getG(const Eigen::Vector3d &pos)
	{
		// 计算当前纬度的重力加速度
		double g_lat = (EarthLongAxis * GA * pow(cos(pos[Latitude]), 2) + Rb * GB * pow(sin(pos[Latitude]), 2)) / sqrt(pow(EarthLongAxis, 2) * pow(cos(pos[Latitude]), 2) + pow(Rb, 2) * pow(sin(pos[Latitude]), 2));
		double m = pow(EarthAngleVelocity, 2) * pow(Ra, 2) * Rb / GM;
		// 计算当前位置的重力加速度
		Eigen::Vector3d g_n;
		g_n << 0, 0, g_lat * (1 - 2.0 / Ra * (1 + F + m - 2 * F * pow(sin(pos[Latitude]), 2)) * pos[Height] + 3 * pow(pos[Height], 2) / pow(Ra, 2));
		// g_n <<  0, 0,  9.80665;
		return g_n;
	}
};