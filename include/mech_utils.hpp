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
	// #define		INS_T						(1.0/INS_SampleRate)
	//INS初始参数
	#define		TIME_BEGIN					91620.0

	const double Ra = 6378137.0;            // 长轴
	const double Rb = 6356752.3142;         // 短轴
	const double E2 = 0.00669437999013;     // 第一偏心率e的平方 e=sqrt(a^2-b^2)/a
	const double GA = 9.7803267715;
	const double GB = 9.8321863685;
	const double F = 1.0 / 298.257223563;
	const double GM = 3.986005e14; // 万有引力常量G * 地球质量M
	// 组合导航单位换算
	const double hour = 3600.0;
	const double sqrt_hour = 60.0;
	const double mGal = 1.0E-5;
	const double ppm = 1.0E-6;
	// 组合导航初始化参数

	// 组合导航参数
	const double Tgs = 4.0 * hour;
	const double Tgb = 4.0 * hour;
	const double Tas = 4.0 * hour;
	const double Tab = 4.0 * hour;
	const double VRW_std = 0.03 / sqrt_hour;
	const double ARW_std = 0.003 * D2R /sqrt_hour;
	const double gb_std = 0.027 * D2R / hour;
	const double ab_std = 15.0 * mGal;
	const double gs_std = 300.0 * ppm;
	const double as_std = 300.0 * ppm;

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

	static double getRM(const double &latitude_in)
	{
		return EarthLongAxis * (1.0 - EarthECC2) / sqrt(pow3Func(1.0 - EarthECC2 * pow2Func(sin(latitude_in))));
	}
	// 计算RN
	static double getRN(const double &latitude_in)
	{
		return EarthLongAxis / sqrt(1.0 - EarthECC2 * pow2Func(sin(latitude_in)));
	}

	// 计算当前位置重力加速度
	static Eigen::Vector3d getG(const Eigen::Vector3d &pos)
	{
		// 计算当前纬度的重力加速度
		double g_lat = (EarthLongAxis * GA * pow(cos(pos[Latitude]), 2) + Rb * GB * pow(sin(pos[Latitude]), 2)) / sqrt(pow(EarthLongAxis, 2) * pow(cos(pos[Latitude]), 2) + pow(Rb, 2) * pow(sin(pos[Latitude]), 2));
		double m = pow(EarthAngleVelocity, 2) * pow(Ra, 2) * Rb / GM;
		//double R_M = getRM(pos(Latitude));
		//double R_N = getRN(pos(Latitude));

		// 计算当前位置的重力加速度
		Eigen::Vector3d g_n;
		//g_n << 0, 0, g_lat * R_M * R_N / (pow2Func(sqrt(R_M * R_N) + pos(Height)));
		g_n << 0, 0, g_lat * (1 - 2.0 / Ra * (1 + F + m - 2 * F * pow(sin(pos[Latitude]), 2)) * pos[Height] + 3 * pow(pos[Height], 2) / pow(Ra, 2));
		// g_n <<  0, 0,  9.80665;
		return g_n;
	}
	// 计算h=0处重力
	static double getGp(const Eigen::Vector3d &pos)
	{
		double gp;
		gp = (EarthLongAxis * GA * pow(cos(pos[Latitude]), 2) + Rb * GB * pow(sin(pos[Latitude]), 2)) / sqrt(pow(EarthLongAxis, 2) * pow(cos(pos[Latitude]), 2) + pow(Rb, 2) * pow(sin(pos[Latitude]), 2));
		return gp;
	}

	// 计算RM
	
	//
	static void getFrr(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, const double &R_N, const double &R_M, Eigen::Matrix3d &F_rr)
	{
		F_rr << -vel(D) / (R_M + pos(Height)), 							0, 																vel(N) / (R_M + pos(Height)),
				(vel(E) * tan(pos(Latitude))) / (R_N + pos(Height)), 	-(vel(D) + vel(N) * tan(pos(Latitude))) / (R_N + pos(Height)), 	vel(E) / (R_N + pos(Height)),
				0, 														0, 																0;
	}
	static void getFvr(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, const double &R_N, const double &R_M, const double &gp, Eigen::Matrix3d &F_vr)
	{
		double vN, vE, vD;
		double phi, h;
		double omega_e;
		vN = vel(N);
		vE = vel(E);
		vD = vel(D);
		phi = pos(Latitude);
		h = pos(Height);
		omega_e = EarthAngleVelocity;
		F_vr << -(2.0 * vE * omega_e * cos(phi)) / (R_M + h) - (pow2Func(vE) * pow2Func(1.0 / cos(phi))) / ((R_M + h) * (R_N + h)),
			0.0,
			(vN * vD) / (pow2Func(R_M + h)) - (pow2Func(vE) * tan(phi)) / (pow2Func(R_N + h)),
			2.0 * omega_e * (vN * cos(phi) - vD * sin(phi)) / (R_M + h) + (vN * vE * pow2Func(1.0 / cos(phi))) / ((R_M + h) * (R_N + h)),
			0.0,
			(vE * vD + vN * vE * tan(phi)) / (pow2Func(R_N + h)),
			2.0 * omega_e * vE * sin(phi) / (R_M + h),
			0.0,
			-pow2Func(vE) / pow2Func(R_N + h) - pow2Func(vN) / pow2Func(R_M + h) + 2.0 * gp / (sqrt(R_M * R_N) + h);
	}

	static void getFvv(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, const double &R_N, const double &R_M, Eigen::Matrix3d &F_vv)
	{
		double vN, vE, vD;
		double phi, h;
		double omega_e;
		vN = vel(N);
		vE = vel(E);
		vD = vel(D);
		phi = pos(Latitude);
		h = pos(Height);
		omega_e = EarthAngleVelocity;
		F_vv << 	vD / (R_M + h), 											-2.0 * (omega_e * sin(phi) + (vE * tan(phi) / (R_N + h))), 	vN / (R_M + h),
					2.0 * omega_e * sin(phi) + (vE * tan(phi) / (R_N + h)), 	(vD + vN * tan(phi)) / (R_N + h), 							2.0 * omega_e * cos(phi) + vE / (R_N + h),
					-(2.0 * vN) / (R_M + h), 									-2.0 * (omega_e * cos(phi) + vE / (R_N + h)), 				0.0;
	}

	static void getFphir(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, const double &R_N, const double &R_M, Eigen::Matrix3d &F_phir)
	{
		double vN, vE, vD;
		double phi, h;
		double omega_e;
		vN = vel(N);
		vE = vel(E);
		vD = vel(D);
		phi = pos(Latitude);
		h = pos(Height);
		omega_e = EarthAngleVelocity;

		F_phir << 	-omega_e * sin(phi) / (R_M + h), 																0.0, 	vE / pow2Func(R_N + h),
					0.0, 																							0.0, 	-vN / pow2Func(R_M + h),
					-omega_e * cos(phi) / (R_M + h) - (vE * pow2Func(1.0 / cos(phi)) / ((R_M + h) * (R_N + h))), 	0.0, 	-vE * tan(phi) / pow2Func(R_N + h);
	}
	static void getFphiv(const Eigen::Vector3d &pos, const double &R_N, const double &R_M, Eigen::Matrix3d &F_phiv)
	{
		double phi, h;
		phi = pos(Latitude);
		h = pos(Height);

		F_phiv << 	0.0, 			1.0/(R_N + h), 			0.0,
					-1.0/(R_M + h), 0.0, 					0.0,
					0.0, 			-tan(phi)/(R_N + h), 	0.0;

	}
};