#include "image_data.h"
#include "metadata.h"
#include "weight.h"
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>

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

int output_csv(int w,int h,char *input_pgmname){
	string RAW = input_pgmname;
	string G1 = input_pgmname;
	string G2 = input_pgmname;
	string B = input_pgmname;
	string R = input_pgmname;

	//RAW = G1 = G2 = B = R = input_pgmname;
	int l = RAW.size();

	RAW.erase(l-5);
	G1.erase(l-5);
	G2.erase(l-5);
	B.erase(l-5);
	R.erase(l-5);
	RAW += ".csv";
	G1 += "_g1.csv";
	G2 += "_g1.csv";
	B += "_b.csv";
	G1 += "_r.csv";

	//char *output_csv_g1, *output_csv_r, *output_csv_b, *output_csv_g2;
	/*
	while (1) {
		if (input_pgmname[num] == '.'&&input_pgmname[num+1] == 'p' && input_pgmname[num + 2] == 'g' && ouput_filename[num + 3] == 'm') {
			input_pgmname[num + 1] = 'c';
			input_pgmname[num + 2] = 's';
			input_pgmname[num + 3] = 'v'; 
			break;
		}
		num++;
	}*/

	ofstream ofs;
	ofstream ofs_G1;
	ofstream ofs_R;
	ofstream ofs_B;
	ofstream ofs_G2;

	ofs.open(RAW.c_str());
	ofs_G1.open(G1.c_str());
	ofs_G2.open(G2.c_str());
	ofs_R.open(R.c_str());
	ofs_B.open(B.c_str());


    for (int j=0;j<h;j++){
		for(int i=1;i<w;i++){
			int value = data[j][i].value;
            double d = data[j][i].d;
            double ans = compensation_func(value,d,a,b);

			ofs << ans << ",";

			if (i % 2 == 1 && j % 2 == 1) { ofs_G1 << ans << ","; } //G1
			else if (i % 2 == 0 && j % 2 == 1) { ofs_R << ans << ","; } //R 
			else if (i % 2 == 1 && j % 2 == 0) { ofs_B << ans << ","; } //B
			else if (i % 2 == 0 && j % 2 == 0) { ofs_G2 << ans << ","; } //G2
		}
		ofs_G1 << endl;
		ofs_R << endl;
		ofs_B << endl;
		ofs_G2 << endl;
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
