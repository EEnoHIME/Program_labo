#include "image_data.h"
#include "metadata.h"
#include "weight.h"
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>

extern double a;
extern double b;
extern img_data data[8000][8000];
extern metadata md;

using namespace std;
//int data[5][5];

double compensation_func(int value,int d,double a,double b){
    double div = (d*d*a*a+1)*(d*d*a*a+1);
    double noise = (1-div)*b/div;
    double compensation = value/div - noise;
    return compensation/value;
}

int output_csv(int w,int h,char *output_filename){
	ofstream ofs(output_filename);
    for (int j=0;j<h;j++){
		for(int i=1;i<w;i++){
			int value = data[j][i].value;
            double d = data[j][i].d;
            double ans = compensation_func(value,d,a,b);
			ofs << ans << ",";
		}
		ofs << endl;
	}
    return 0;
}
/*
int main(int argc,char *argv[]) {
	for (int j = 0; j < 5; j++) {
		for (int i = 0; i < 5; i++) {
			data[j][i] = j*i;
		}
	}
	output_csv(5,5,argv[1]);
	return 0;
}
*/
