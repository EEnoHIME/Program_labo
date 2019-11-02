#ifndef IMGD_H
#define IMGD_H
typedef struct img_data {
	int x;
	int y;
	int value;
	int m;
	int n;
	double d;
	int e;
};
extern void readPGM(const char* filename);
#endif