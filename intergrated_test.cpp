#include <Eigen\Dense>
#include <iostream>
#include "intergratedNav.h"
using namespace std;

const double inter_time_begin = 456300.0;

void matrix_test()
{
    Eigen::Matrix3d mat3d;
    Eigen::Vector3d vec3d;
    Eigen::Matrix<double, 5, 2> mat_5_2;
    Eigen::Matrix<double, 2, 5> mat_2_5;
    mat3d.setZero();
    cout << mat3d << endl;
    mat3d.setOnes();
    cout << mat3d << endl;
    cout << Eigen::Matrix3d::Identity() << endl;
    vec3d.setOnes();
    mat3d = vec3d.asDiagonal();
    cout << vec3d << endl;
    cout << mat3d << endl;
    mat_5_2.setOnes();
    mat_2_5 = mat_5_2.transpose();
    cout << mat_5_2 << endl;
    cout << mat_2_5 << endl;
    mat3d << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    cout << mat3d << endl;
    cout << (mat3d.array().square()).matrix() << endl;
    cout << mat3d << endl; 
}

int main(int argc, char **argv)
{

    std::ifstream imudata_fs("../internavdata/A15_imu.bin", std::ios::binary);
    std::ifstream gnssdata_fs("../internavdata/GNSS_RTk.txt");
    std::ifstream result_fs("../internavdata/truth.nav");

    //matrix_test();
    INS::m_State state;
    INS::IMU_data imudata;
    INS::GNSS_data_T gnssdata;
    Eigen::Vector3d lb;
    Eigen::Vector3d pos_std;
    Eigen::Vector3d vel_std;
    Eigen::Vector3d att_std;

    INS::NAV_data_T navdata;
    FILE *fout_result;

	fopen_s(&fout_result, "../internavdata/combresult.nav", "w+");

    double data_done = 1;

    do{
        imudata_fs.read((char *)(&imudata), sizeof(imudata));
    } while (imudata.timestamp < inter_time_begin);

    do{
        gnssdata_fs >> gnssdata.timestamp >> gnssdata.pos_blh(INS::Latitude) >> gnssdata.pos_blh(INS::Longitude) >> gnssdata.pos_blh(INS::Height) >>
        gnssdata.pos_blh_std(INS::Latitude) >> gnssdata.pos_blh_std(INS::Longitude) >> gnssdata.pos_blh_std(INS::Height);
    } while (gnssdata.timestamp < inter_time_begin);
    gnssdata.pos_blh(INS::Latitude) = gnssdata.pos_blh(INS::Latitude) * D2R;
    gnssdata.pos_blh(INS::Longitude) = gnssdata.pos_blh(INS::Longitude) * D2R;
    // INS状态初始化
    state.timestamp = inter_time_begin;
	state.NED_vec	    <<	0.0,					0.0,					0.0;
	state.GEO_eb		<<	30.444787369*D2R,		114.471863247*D2R,		20.910;
	state.e_bn		    <<	185.696*D2R,	        -2.0345*D2R,	        0.854*D2R;
    // 杆臂信息
    lb                  <<  0.136,                  -0.301,                 -0.184;
    // 初始导航状态std
    pos_std             <<  0.005,                  0.004,                  0.008;
    vel_std             <<  0.003,                  0.004,                  0.004;
    att_std             <<  0.023*D2R,              0.003*D2R,              0.003*D2R;

    INS::intergrated interNav(state, imudata, lb, pos_std, vel_std, att_std);

    while(!imudata_fs.eof())
    {
        imudata_fs.read((char *)(&imudata), sizeof(imudata));
        if (data_done == 1)
        {
            gnssdata_fs >> gnssdata.timestamp >> gnssdata.pos_blh(INS::Latitude) >> gnssdata.pos_blh(INS::Longitude) >> gnssdata.pos_blh(INS::Height) >>
                gnssdata.pos_blh_std(INS::Latitude) >> gnssdata.pos_blh_std(INS::Longitude) >> gnssdata.pos_blh_std(INS::Height);
            gnssdata.pos_blh(INS::Latitude) = gnssdata.pos_blh(INS::Latitude) * D2R;
            gnssdata.pos_blh(INS::Longitude) = gnssdata.pos_blh(INS::Longitude) * D2R;
        }
        data_done = interNav.IntergtatedNavUpdate(imudata, gnssdata);
        state = interNav.getstate();
        fprintf(fout_result, "%5d %10.7lf %7.10lf %7.10lf %7.10lf %7.10lf %7.10lf %7.10lf %7.10lf %7.10lf, %7.10lf\n",	
            2017,
			state.timestamp,
			state.GEO_eb(INS::Latitude) * R2D,
			state.GEO_eb(INS::Longitude) * R2D,
			state.GEO_eb(INS::Height),
			state.NED_vec(INS::N),
			state.NED_vec(INS::E),
			state.NED_vec(INS::D),
			state.e_bn(INS::Roll) * R2D,
			state.e_bn(INS::Pitch) * R2D,
			state.e_bn(INS::Yaw) * R2D
			);
        if (0)
        {
            do
            {
                result_fs >> navdata.week >> navdata.second >>
                    navdata.pos_blh(INS::Latitude) >> navdata.pos_blh(INS::Longitude) >> navdata.pos_blh(INS::Height) >>
                    navdata.vel_ned(INS::N) >> navdata.vel_ned(INS::E) >> navdata.vel_ned(INS::D) >>
                    navdata.att_ypr(INS::Roll) >> navdata.att_ypr(INS::Pitch) >> navdata.att_ypr(INS::Yaw);
            } while ((state.timestamp - navdata.second) > 0.005);
            std::cout.precision(12);
			std::cout << "ref: ";
            cout << navdata.pos_blh(INS::Latitude) << " " << navdata.pos_blh(INS::Longitude) << " " << navdata.pos_blh(INS::Height) << endl;
            cout << navdata.vel_ned(INS::N) << " " << navdata.vel_ned(INS::E) << " " << navdata.vel_ned(INS::D) << endl;
            cout << navdata.att_ypr(INS::Yaw) << " " << navdata.att_ypr(INS::Pitch) << " " << navdata.att_ypr(INS::Roll) << endl;
            std::cout << "compute :";
            cout << state.GEO_eb(INS::Latitude) * R2D << " " << state.GEO_eb(INS::Longitude) * R2D << " " << state.GEO_eb(INS::Height) << endl;
            cout << state.NED_vec(INS::N) << " " << state.NED_vec(INS::E) << " " << state.NED_vec(INS::D) << endl;
            cout << state.e_bn(INS::Yaw) * R2D << " " << state.e_bn(INS::Pitch) * R2D << " " << state.e_bn(INS::Roll) * R2D << endl;
        }
    }

    imudata_fs.close();
    gnssdata_fs.close();
    result_fs.close();
    fclose(fout_result);

    system("pause");
    return 0;
}