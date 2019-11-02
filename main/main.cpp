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

extern metadata md;
extern img_data data[8000][8000];
extern mla_center_array mla[700][700];
extern mla_center_array mla_rot[700][700];
extern value_for_eq coeff;
extern double a,b;
extern double A, B, C;

using namespace std;

int main(int argc,char *argv[])
{
	readPGM(argv[1]);
	get_mla_ary(argv[2]);
	voronoi();
	calc_coeff(md.h, md.h);

	cout << "A is " << A << endl;
	cout << "B is " << B << endl;
	cout << "C is " << C << endl;

	a = solve_equation();
	cout << "a is " << a << endl;
	cout << "b is " << b << endl;

	output_csv(md.w, md.h, argv[1]);
	return 0;
}
