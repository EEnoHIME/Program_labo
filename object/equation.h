#include "image_data.h"
#include <iostream>  // for cout
#include <math.h>    // for fabs()
#include <stdio.h>   // for printf()


#ifndef EQUATION
#define EQUATION
typedef struct value_for_eq{
    double sqq[4];
    double sqd[4];
    double c[4];
};
#endif			

extern double F(double x);
extern double  G(double x);
extern int calc_coeff(int m,int n);
extern double solve_equation();