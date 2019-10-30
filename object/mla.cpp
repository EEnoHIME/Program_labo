#include "mla.h"
#include "metadata.h"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

using namespace std;

struct mla_center_array mla[700][700];
struct mla_center_array mla_rot[700][700];

double m2pixel(double scale, double pixelPitch) {
	double pixel_num = scale / pixelPitch;
	return pixel_num;
}

int roundUp(double n) {
	if (n >= 0) {
		return (int)n + 1;
	}
	else {
		return (int)n - 1;
	}
}

int roundDown(double n) {
	return (int)n;
}

int roundOff(double n) {
	double decimal = 0;

	decimal = n - (int)n;

	if (decimal >= 0.5 || decimal <= -0.5) {
		return roundUp(n);
	}
	else {
		return roundDown(n);
	}
}



/*回転考慮なし*/
int mla_center_func(double default_x, double default_y,double lensDiameter,int m_max,int n_max) {
	//declaration
	double X = default_x;
	double Y = default_y;
	double lDiameter = lensDiameter;
	int M_max = m_max;
	int N_max = n_max;
	//counter
	int i = 0;
	int j = 0;

	for (j = 0; j < N_max; ++j) {
		for (i = 0; i < M_max; ++i) {
			X += lensDiameter;
			mla[j][i].x = X;
			mla[j][i].y = Y;
		}
		X = default_x;
		Y += lensDiameter*sin(M_PI/ 3);
	}
	for (j = 0; j < N_max; ++j) {
		for (i = 0; i < M_max; ++i) {
			if(j%2==0){ mla[j][i].x += 7+lensDiameter/2; 
			}
			else mla[j][i].x += 7;
			mla[j][i].y -= 1;
		}
	}

	return 0;
}

/*回転考慮あり*/
int  mla_center_func_rot(double d_X, double d_Y, metadata d, int m_max, int n_max) {
	double X;
	double Y;
	double img_c_x = d_X + d.w/2;
	double img_c_y = d_Y + d.h/2;
	int M_max = m_max;
	int N_max = n_max;
	int i = 0;
	int j = 0;

	for (j = 0; j < N_max; ++j) {
		for (i = 0; i < M_max; ++i) {
			X = (mla[j][i].x-img_c_x)*cos(d.rot) - (mla[j][i].y-img_c_y)*sin(d.rot) + img_c_x;
			Y = (mla[j][i].x - img_c_x)*sin(d.rot) + (mla[j][i].y - img_c_y)*cos(d.rot) + img_c_y;
			mla_rot[j][i].x = X;
			mla_rot[j][i].y = Y;
		}
	}
	return 0;
}
