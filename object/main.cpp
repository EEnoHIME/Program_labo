#include "image_data.h"
#include "metadata.h"
#include "equation.h"
#include "mla.h"
#include "voronoi.h"
#include "weight.h"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

extern img_data data[8000][8000];
extern mla_center_array mla[700][700];
extern mla_center_array mla_rot[700][700];
extern value_for_eq coeff;
extern double a,b;
//extern double compensation_func(int value,int d,double a,double b);


using namespace std;

int main(int argc,char *argv[])
{
	readPGM(argv[1]);
	get_mla_ary(argv[2]);
	cout << mla_rot[100][100].x << endl;
	cout << mla_rot[100][100].y<< endl;
	voronoi();

	calc_coeff(md.h, md.h);
	a = solve_equation();
	cout << "a is " << a << endl;
	cout << "b is " << b << endl;


	//output_csv();
	return 0;
}
