#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include "mech.h"
#include "mech_utils.hpp"

int main(int argc, char** argv)
{
	
	// 变量声明
	INS::dataReader datareader;
	INS::m_State state_in;
	FILE *fout_result;

	fopen_s(&fout_result, "../data/mech_result.txt", "w+");
	int i = 0;

	do {
		datareader.imuRead();
	} while (datareader.imudata.timestamp < TIME_BEGIN);

	// 初始化
	state_in.timestamp = datareader.imudata.timestamp;
	state_in.NED_vec	<<	0.0,					0.0,					0.0;
	state_in.GEO_eb		<<	23.1373950708*D2R,		113.3713651222*D2R,		2.175;
	state_in.e_bn		<<	-75.7498049314083*D2R,	-2.14251290749072*D2R,	0.0107951084511778*D2R;
	INS::mechanization mechan(state_in, datareader.imudata);
	while (!datareader.imu_file_handle.eof())
	{
		datareader.imuRead();
		mechan.mechanizationUpdate(datareader.imudata);
		datareader.refRead();
		fprintf(fout_result, "%10.7lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %10.7lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf, %7.10lf\n",	
			datareader.refdata.timestamp, 
			datareader.refdata.pos(INS::Latitude),
			datareader.refdata.pos(INS::Longitude),
			datareader.refdata.pos(INS::Height),
			datareader.refdata.vel(INS::N),
			datareader.refdata.vel(INS::E),
			datareader.refdata.vel(INS::D),
			datareader.refdata.att(2),
			datareader.refdata.att(1),
			datareader.refdata.att(0),
			datareader.imudata.timestamp,
			(mechan.getstate()).GEO_eb(INS::Latitude) * R2D,
			(mechan.getstate()).GEO_eb(INS::Longitude) * R2D,
			(mechan.getstate()).GEO_eb(INS::Height),
			(mechan.getstate()).NED_vec(INS::N),
			(mechan.getstate()).NED_vec(INS::E),
			(mechan.getstate()).NED_vec(INS::D),
			(mechan.getstate()).e_bn(INS::Yaw) * R2D,
			(mechan.getstate()).e_bn(INS::Pitch) * R2D,
			(mechan.getstate()).e_bn(INS::Roll) * R2D
			);
		i++;
		/*
		if (i == 10000)
		{
			i = 0;
			std::cout.precision(12);
			std::cout << "ref: ";
			std::cout << datareader.refdata.timestamp << " ";
			std::cout << datareader.refdata.pos(INS::Latitude) << " "
				<< datareader.refdata.pos(INS::Longitude) << " "
				<< datareader.refdata.pos(INS::Height) << " "
				<< std::endl;
			std::cout << datareader.refdata.vel(INS::N) << " "
				<< datareader.refdata.vel(INS::E) << " "
				<< datareader.refdata.vel(INS::D) << " "
				<< std::endl;
			std::cout << datareader.refdata.att(2) << " "
				<< datareader.refdata.att(1) << " "
				<< datareader.refdata.att(0) << " "
				<< std::endl;
			std::cout << "calc: ";
			std::cout << datareader.imudata.timestamp << " ";
			std::cout << (mechan.getstate()).GEO_eb(INS::Latitude) * R2D << " "
				<< (mechan.getstate()).GEO_eb(INS::Longitude) * R2D << " "
				<< (mechan.getstate()).GEO_eb(INS::Height) << " "
				<< std::endl;
			std::cout << (mechan.getstate()).NED_vec(INS::N) << " "
				<< (mechan.getstate()).NED_vec(INS::E) << " "
				<< (mechan.getstate()).NED_vec(INS::D) << " "
				<< std::endl;
			std::cout << (mechan.getstate()).e_bn(INS::Yaw) * R2D << " "
				<< (mechan.getstate()).e_bn(INS::Pitch) * R2D << " "
				<< (mechan.getstate()).e_bn(INS::Roll) * R2D << " "
				<< std::endl;
			std::cout << std::endl;
		}
		*/
	}
	
	system("pause");
	return 0;
	//INS::mechanization()

}