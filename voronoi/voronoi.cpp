#include "voronoi.h"
#include "mla.h"
#include "image_data.h"
#include "metadata.h"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

extern metadata md;
extern mla_center_array mla[700][700];
extern mla_center_array mla_rot[700][700];
extern img_data data[8000][8000];

struct ID id;

using namespace std;

void voronoi(){
	double min_d = 15;
    cout << "( i , j ) belongs to ( m , n )" << endl;
    for(int j=1;j<md.h+1;j++){
        for(int i=1;i<md.w+1;i++){
            id.d = 14;
            id.m = -1;
            id.n = -1;
            for(int n=0;n<md.n_max;n++){
                for(int m=0;m<md.m_max;m++){
                    double dx = mla_rot[n][m].x-i;
                    double dy = mla_rot[n][m].y-j;
                    double distance = sqrt(dx*dx+dy*dy);
                    if(distance < min_d){
                        id.m = m;
                        id.n = n;
                        id.d = distance;
                    }
                }
            } 
            data[j][i].m = id.m;
            data[j][i].n = id.n;
			data[j][i].d = id.d;
            //cout << "(" << i << "," << j <<")" << " belongs to " << "(" << data[j][i].m << "," << data[j][i].n <<")" << endl;
        }
		cout << "Line_" << j  << " is done " <<  endl;
    }
	cout << "Voronoi is done " << endl;
}