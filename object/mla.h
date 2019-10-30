#include "metadata.h"
#ifndef MLA_H
#define MLA_H
typedef struct mla_center_array {
	double x;
	double y;
};
typedef struct image_data {
	int x;
	int y;
	int value;
	int m;
	int n;
	double d;
	int e;
};
#endif

extern double m2pixel(double scale, double pixelPitch);
extern int roundUp(double n);
extern int roundDown(double n);
extern int roundOff(double n);
extern int json_parse(metadata md);
extern int mla_center_func(double default_x, double default_y,double lensDiameter,int m_max,int n_max);
extern int mla_center_func_rot(double d_X, double d_Y, metadata d, int m_max, int n_max);
extern int get_mla_ary(const char *inputfile);
