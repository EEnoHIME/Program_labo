#include "metadata.h"
#include "mla.h"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

extern metadata md;
extern mla_center_array mla[700][700];
extern mla_center_array mla_rot[700][700];


using namespace std;

int get_mla_ary(char *argv[])
{
	json_parse(argv[0]);
	md.lensDiameter = m2pixel(md.lP, md.pP);
	md.default_x = m2pixel(md.x_off, md.pP);
	md.default_y = m2pixel(md.y_off, md.pP);
	md.R = (lensDiameter + 1) / 2;
	md.m_max = roundOff((md.w + R )/lensDiameter);
	md.n_max = roundOff((md.h+R)/(lensDiameter*sin(M_PI / 3)));

	mla_center_func(default_x, default_y, lensDiameter, m_max, n_max);
	mla_center_func_rot(default_x, default_y, md, m_max, n_max);

	//cout << md.n_max << endl;
	//cout << md.m_max << endl;
	//cout << mla_rot[0][0].x << "," << mla_rot[0][0].y << endl;

	return 0;
}
