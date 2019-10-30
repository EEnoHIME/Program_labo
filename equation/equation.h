#include "image_data.h"
#include <iostream>  // for cout
#include <math.h>    // for fabs()
#include <stdio.h>   // for printf()


#ifndef EQUATION
#define EQUATION
typedef struct value_for_eq{
    double 4_sq[4];
    double 2_sq[4];
    double const[4];
};
#endif			

extern double F(double x);
extern double  G(double x);
extern int calc_coeff(img_data **img_d,value_for_eq coeff,int m,int n);
extern double solve_equation();