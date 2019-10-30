#include <sstream>
#include <iostream>
#include <math.h>

#ifndef METADATA_H
#define METADATA_H
typedef struct metadata {
	int h;
	int w;
	double pP;
	double x_off;
	double y_off;
	double z_off;
	double lP;
	double rot;
	double lensDiameter;
	double default_x;
	double default_y;
	double R ;
	int m_max;
	int n_max;
}; 
extern int json_parse(const char *inputfile);
#endif
