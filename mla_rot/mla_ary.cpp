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

int get_mla_ary(const char *inputfile)
{
	json_parse(inputfile);
	md.lensDiameter = m2pixel(md.lP, md.pP);
	md.default_x = m2pixel(md.x_off, md.pP);
	md.default_y = m2pixel(md.y_off, md.pP);
	md.R = (md.lensDiameter + 1) / 2;
	md.m_max = roundOff((md.w + md.R )/md.lensDiameter);
	md.n_max = roundOff((md.h+md.R)/(md.lensDiameter*sin(M_PI / 3)));

	mla_center_func(md.default_x, md.default_y, md.lensDiameter, md.m_max, md.n_max);
	mla_center_func_rot(md.default_x, md.default_y, md, md.m_max, md.n_max);

	cout << md.n_max << endl;
	cout << md.m_max << endl;
	cout << mla_rot[0][0].x << "," << mla_rot[0][0].y << endl;

	return 0;
}
