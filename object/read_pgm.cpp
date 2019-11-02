#include "image_data.h"
#include "metadata.h"
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <bits/stdc++.h>

using namespace cv;
using namespace std;

extern metadata md;
struct img_data data[8000][8000];

void readPGM(const char* filename) {

	int height, width, maxval;
	int i, j;

	//ファイルポインタ用変数finを宣言
	FILE *fin;
	char buf[256];

	//ファイルオープン
	//ファイルが見つからなかったらプログラム終了
	if ((fin = fopen(filename, "rb")) == NULL) {
		printf("%sが見つかりません\n", filename);
		exit(1);
	}

	//ファイルの形式がP5じゃなかったらプログラム終了
	fgets(buf, 256, fin);
	if (buf[0] != 'P' || buf[1] != '5') {
		exit(1);
	}

	//printf("%s", buf);

	//ファイルの幅と高さを読み込む
	//int width,height;
	while (*(fgets(buf, 256, fin)) == '#');
	sscanf(buf, "%d %d\n", &width, &height);
	cout << width << endl;
	cout << height << endl;;

	//ファイルの最大輝度値を読み込む
	//int maxvalue;
	while (*(fgets(buf, 256, fin)) == '#');
	sscanf(buf, "%d\n", &maxval);
	
	for (int y=0; y < 8000;y++) {
		for (int x = 0; x< 8000; x++) {
			data[y][x].value = 0;
		}
	}

	for (int j = 1; j< height+1; j++) {
		for (int i = 1; i< width+1; i++) {
			data[j][i].y = j;
			data[j][i].x = i;
			data[j][i].value = fgetc(fin) * 256 + fgetc(fin);
			data[j][i].e = data[j][i - 1].value + data[j - 1][i].value - data[j - 1][i - 1].value - data[j][i].value;
		}
	}
		
		fclose(fin);
	}


   

/*int main(int argc,char *argv[]) {
	readPGM(argv[1]);
	cout << data[100][100].value << endl;
}*/
