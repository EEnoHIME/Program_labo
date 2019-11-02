#include "equation.h"
#include <iostream>  // for cout
#include <math.h>    // for fabs()
#include <stdio.h>   // for printf()

using namespace std;
extern double A,B,C;

//double A = 1;
//double B = 6;
//double C = 9;


double solve_equation(){
	double ans_a;
    
	double D = B*B - 4 * A*C;

	if (D < 0) {
		cout << "虚数解" << endl;
		return 0;
	}

	else{
		double bunsi = sqrt(D) - B;
		double bunbo = 2 *A;
		ans_a = bunsi / bunbo;
		return ans_a;
	}
}

 /*
 int main()
{
    double a = solve_equation();
	cout << a << endl;
	return 0;
}
*/

 
