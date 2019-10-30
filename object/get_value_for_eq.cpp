#include "equation.h"
#include "image_data.h"
#include <iostream>  // for cout
#include <math.h>    // for fabs()
#include <stdio.h>   // for printf()

//extern struct img_data data[8000][8000];
//struct value_for_eq coeff;

value_for_eq coeff;
extern img_data data[8000][8000];
extern double b;

int calc_coeff(int m,int n){
    int i,j;
    for(j=1; j<n+1; j++){
        for(i=1;i<m+1;i++){
            coeff.sqq[0] += (data[j][i].value-b)*data[j][i].d*data[j][i].d*data[j][i].d*data[j][i].d;
            coeff.sqq[1] += (data[j][i-1].value-b)*data[j][i-1].d*data[j][i-1].d*data[j][i-1].d*data[j][i-1].d;
            coeff.sqq[2] += (data[j-1][i].value-b)*data[j-1][i].d*data[j-1][i].d*data[j-1][i].d*data[j-1][i].d;
            coeff.sqq[3] += (data[j-1][i-1].value-b)*data[j-1][i-1].d*data[j-1][i-1].d*data[j-1][i-1].d*data[j-1][i-1].d;

            coeff.sqd[0] += 2*(data[j][i].value-b)*data[j][i].d*data[j][i].d;
            coeff.sqd[1] += 2*(data[j][i-1].value-b)*data[j][i-1].d*data[j][i-1].d;
            coeff.sqd[2] += 2*(data[j-1][i].value-b)*data[j-1][i].d*data[j-1][i].d;
            coeff.sqd[3] += 2*(data[j-1][i-1].value-b)*data[j-1][i-1].d*data[j-1][i-1].d;

			coeff.c[0] += data[j][i].value - b;
			coeff.c[1] += data[j][i - 1].value - b;
			coeff.c[2] += data[j - 1][i].value - b;
			coeff.c[3] += data[j - 1][i - 1].value - b;
        }
    }
    return 0;
}