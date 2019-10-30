#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mrp.h"

//extern POINT_REF dyx[NUM_OF_VGEN][126];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];
int ***roff_up;
int ***roff_left;
int ***roff_right;
int mindy, mindx, maxdx;
int mindy2, mindx2, maxdx2;
//int roff_const[NUM_OF_VGEN][MAX_PRD_ORDER];
int ***roff_up2;
int ***roff_left2;
int ***roff_right2;
int roff_const2[MAX_PRD_ORDER];
extern POINT_REF dyx2[];

void init_ref_offset5_2(int prd_order, int height, int width, DECODER *dec)
{	
	printf("init_ref_offset5_2->");
	int k, dy, dx, y, x, flag, temp_k, temp_dy, temp_dx, yy, xx;
	mindy2 = 0;
	mindx2 = 0;
	maxdx2 = 0;
	for (k = 0; k < prd_order; k++) {
		dy = dyx2[k].y;
		dx = dyx2[k].x;
		if(dy > 0) printf("you must change dyx2[] in common.c\n");
		if(dy == 0 && dx >= 0) printf("you must change dyx2[] in common.c\n");
		if (dy < mindy2) mindy2 = dy;
		if (dx < mindx2) mindx2 = dx;
		if (dx > maxdx2) maxdx2 = dx;
		roff_const2[k] = dy * width + dx;
	}
	roff_up2 = (int ***)alloc_3d_array(abs(mindy2), width, prd_order, sizeof(int));
	roff_left2 = (int ***)alloc_3d_array(height-abs(mindy2), abs(mindx2), prd_order, sizeof(int));
	roff_right2 = (int ***)alloc_3d_array(height-abs(mindy2), abs(maxdx2), prd_order, sizeof(int));
	//printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
	printf("roff_up->");
	for(y = 0; y < abs(mindy2); y++){
		for(x = 0; x < width; x++){
			for(k = 0; k < prd_order; k++){
				if(y == 0 && x == 0){
					dx = 0;
					dy = height;
					roff_up2[y][x][k] = dy * width + dx;
				}else{
					dy = dyx2[k].y;
					dx = dyx2[k].x;
					if(y + dy < 0 || x + dx < 0 || x + dx >= width){
						flag = 0;
						temp_k = k;
						temp_k--;
						if(color(y, x) == 1 || color(y, x) == 2){
							while(temp_k >= 0){
								temp_dy = dyx2[temp_k].y;
								temp_dx = dyx2[temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(dyx2[k].p == dyx2[temp_k].p){
										roff_up2[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
							}
						}else{
							while(temp_k >= 0){
								temp_dy = dyx2[temp_k].y;
								temp_dx = dyx2[temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(dyx2[k].p == 0 || dyx2[k].p == 2){
										if(dyx2[temp_k].p == 0 || dyx2[temp_k].p == 2){
											roff_up2[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}else{
										if(dyx2[k].p == dyx2[temp_k].p){
											roff_up2[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}
									
								}
								temp_k--;
							}
						}
						
						if(flag == 0){
							if(x != 0){
								
								temp_dy = 0;
								temp_dx = -1;
							}else{
								temp_dy = -1;
								temp_dx = 0;
							}
						}
						roff_up2[y][x][k] = temp_dy * width + temp_dx;
					}else{
						roff_up2[y][x][k] = dy * width + dx;
					}
				}
				//printf("roff_up[%d][%d][%d] = %d\n", y, x, k, roff_up[y][x][k]);
			}
		}
	}
	printf("roff_left2->");
	for(y = 0; y < height-abs(mindy2); y++){
		yy = y + abs(mindy2);
		for(x = 0; x < abs(mindx2); x++){
			for(k = 0; k < prd_order; k++){
				dy = dyx2[k].y;
				dx = dyx2[k].x;
				if(yy + dy < 0 || x + dx < 0 || x + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, x) == 1 || color(yy, x) == 2){
						while(temp_k >= 0){
							temp_dy = dyx2[temp_k].y;
							temp_dx = dyx2[temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(dyx2[k].p == dyx2[temp_k].p){
									roff_left2[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = dyx2[temp_k].y;
							temp_dx = dyx2[temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(dyx2[k].p == 0 || dyx2[k].p == 2){
									if(dyx2[temp_k].p == 0 || dyx2[temp_k].p == 2){
										roff_left2[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(dyx2[k].p == dyx2[temp_k].p){
										roff_left2[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								
							}
							temp_k--;
						}
					}
					
					if(flag == 0){
						if(x != 0){
							temp_dy = 0;
							temp_dx = -1;
						}else{
							temp_dy = -1;
							temp_dx = 0;
						}
					}
					roff_left2[y][x][k] = temp_dy * width + temp_dx;
				}else{
					roff_left2[y][x][k] = dy * width + dx;
				}
			}
		}
	}
	printf("roff_right2->");
	for(y = 0; y < height-abs(mindy2); y++){
		yy = y + abs(mindy2);
		for(x = 0; x < abs(maxdx2); x++){
			xx = x + (width-abs(maxdx2));
			for(k = 0; k < prd_order; k++){
				dy = dyx2[k].y;
				dx = dyx2[k].x;
				if(yy + dy < 0 || xx + dx < 0 || xx + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, xx) == 1 || color(yy, xx) == 2){
						while(temp_k >= 0){
							temp_dy = dyx2[temp_k].y;
							temp_dx = dyx2[temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(dyx2[k].p == dyx2[temp_k].p){
									roff_right2[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = dyx2[temp_k].y;
							temp_dx = dyx2[temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(dyx2[k].p == 0 || dyx2[k].p == 2){
									if(dyx2[temp_k].p == 0 || dyx2[temp_k].p == 2){
										roff_right2[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(dyx2[k].p == dyx2[temp_k].p){
										roff_right2[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								
							}
							temp_k--;
						}
					}
					
					if(flag == 0){
						if(xx != 0){
							temp_dy = 0;
							temp_dx = -1;
						}else{
							temp_dy = -1;
							temp_dx = 0;
						}
					}
					roff_right2[y][x][k] = temp_dy * width + temp_dx;
				}else{
					roff_right2[y][x][k] = dy * width + dx;
				}
			}
		}
	}
	printf("ok\n");
}

/*void init_ref_offset5(int *prd_order, int height, int width, DECODER *dec)
{	
	int k, dy, dx, y, x, flag, temp_k, temp_dy, temp_dx, yy, xx, ar;
	mindy = 0;
	mindx = 0;
	maxdx = 0;
	//roff_const = (int *)alloc_mem(prd_order * sizeof(int));
	for(ar = 0; ar < NUM_OF_VGEN; ar++){
		for (k = 0; k < prd_order[ar]; k++) {
			dy = dyx[ar][k].y;
			dx = dyx[ar][k].x;
			if (dy < mindy) mindy = dy;
			if (dx < mindx) mindx = dx;
			if (dx > maxdx) maxdx = dx;
			roff_const[ar][k] = dy * width + dx;
		}
	}
	roff_up = (int ***)alloc_3d_array(abs(mindy), width, dec->max_prd_order_all, sizeof(int));
	roff_left = (int ***)alloc_3d_array(height-abs(mindy), abs(mindx), dec->max_prd_order_all, sizeof(int));
	roff_right = (int ***)alloc_3d_array(height-abs(mindy), abs(maxdx), dec->max_prd_order_all, sizeof(int));
	for(k = 0; k < dec->max_prd_order_all; k++){
		for(y = 0; y < abs(mindy); y++){
			for(x = 0; x < width; x++){
				roff_up[y][x][k] = -1;
			}
		}
		for(y = 0; y < height-abs(mindy); y++){
			for(x = 0; x < abs(mindx); x++){
				roff_left[y][x][k] = -1;
			}
		}
		for(y = 0; y < height-abs(mindy); y++){
			for(x = 0; x < abs(maxdx); x++){
				roff_right[y][x][k] = -1;
			}
		}
	}
	printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
	printf("roff_up\n");
	for(y = 0; y < abs(mindy); y++){
		for(x = 0; x < width; x++){
			ar = dec->area[y][x];
			for(k = 0; k < prd_order[ar]; k++){
				if(y == 0 && x == 0){
					dx = 0;
					dy = height;
					roff_up[y][x][k] = dy * width + dx;
				}else{
					dy = dyx[ar][k].y;
					dx = dyx[ar][k].x;
					if(y + dy < 0 || x + dx < 0 || x + dx >= width){
						flag = 0;
						temp_k = k;
						temp_k--;
						if(color(y, x) == 1 || color(y, x) == 2){
							while(temp_k >= 0){
								temp_dy = dyx[ar][temp_k].y;
								temp_dx = dyx[ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(dyx[ar][k].p == dyx[ar][temp_k].p){
										roff_up[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
							}
						}else{
							while(temp_k >= 0){
								temp_dy = dyx[ar][temp_k].y;
								temp_dx = dyx[ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(dyx[ar][k].p == 0 || dyx[ar][k].p == 2){
										if(dyx[ar][temp_k].p == 0 || dyx[ar][temp_k].p == 2){
											roff_up[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}else{
										if(dyx[ar][k].p == dyx[ar][temp_k].p){
											roff_up[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}
									
								}
								temp_k--;
							}
						}
						if(flag == 0){
							if(x != 0){
								
								temp_dy = 0;
								temp_dx = -1;
							}else{
								temp_dy = -1;
								temp_dx = 0;
							}
						}
						roff_up[y][x][k] = temp_dy * width + temp_dx;
					}else{
						roff_up[y][x][k] = dy * width + dx;
					}
				}
			}
		}
	}
	printf("roff_left\n");
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(mindx); x++){
			ar = dec->area[yy][x];
			for(k = 0; k < prd_order[ar]; k++){
				dy = dyx[ar][k].y;
				dx = dyx[ar][k].x;
				if(yy + dy < 0 || x + dx < 0 || x + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, x) == 1 || color(yy, x) == 2){
						while(temp_k >= 0){
							temp_dy = dyx[ar][temp_k].y;
							temp_dx = dyx[ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(dyx[ar][k].p == dyx[ar][temp_k].p){
									roff_left[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = dyx[ar][temp_k].y;
							temp_dx = dyx[ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(dyx[ar][k].p == 0 || dyx[ar][k].p == 2){
									if(dyx[ar][temp_k].p == 0 || dyx[ar][temp_k].p == 2){
										roff_left[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(dyx[ar][k].p == dyx[ar][temp_k].p){
										roff_left[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								
							}
							temp_k--;
						}
					}
					if(flag == 0){
						if(x != 0){
							temp_dy = 0;
							temp_dx = -1;
						}else{
							temp_dy = -1;
							temp_dx = 0;
						}
					}
					roff_left[y][x][k] = temp_dy * width + temp_dx;
				}else{
					roff_left[y][x][k] = dy * width + dx;
				}
			}
		}
	}
	printf("roff_right\n");
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(maxdx); x++){
			xx = x + (width-abs(maxdx));
			ar = dec->area[yy][xx];
			for(k = 0; k < prd_order[ar]; k++){
				dy = dyx[ar][k].y;
				dx = dyx[ar][k].x;
				if(yy + dy < 0 || xx + dx < 0 || xx + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, xx) == 1 || color(yy, xx) == 2){
						while(temp_k >= 0){
							temp_dy = dyx[ar][temp_k].y;
							temp_dx = dyx[ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(dyx[ar][k].p == dyx[ar][temp_k].p){
									roff_right[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = dyx[ar][temp_k].y;
							temp_dx = dyx[ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(dyx[ar][k].p == 0 || dyx[ar][k].p == 2){
									if(dyx[ar][temp_k].p == 0 || dyx[ar][temp_k].p == 2){
										roff_right[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(dyx[ar][k].p == dyx[ar][temp_k].p){
										roff_right[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								
							}
							temp_k--;
						}
					}
					if(flag == 0){
						if(xx != 0){
							temp_dy = 0;
							temp_dx = -1;
						}else{
							temp_dy = -1;
							temp_dx = 0;
						}
					}
					roff_right[y][x][k] = temp_dy * width + temp_dx;
				}else{
					roff_right[y][x][k] = dy * width + dx;
				}
			}
		}
	}
}
*/

/*inline int ref_offset5(int height, int width, int y, int x, int ar, int k)
{
	if(y + mindy >= 0 && x + mindx >= 0 && x + maxdx < width){
		return(roff_const[ar][k]);
	}else if(y + mindy < 0){
		return(roff_up[y][x][k]);
	}else if(y + mindy >= 0 && x + mindx < 0){
		return(roff_left[y+mindy][x][k]);
	}else{
		return(roff_right[y+mindy][x-width+maxdx][k]);
	}
	
		
}*/

inline int ref_offset5_2(int height, int width, int y, int x, int k)
{
	//int roff = 0;
	//if(k >=  MAX_PRD_ORDER) {printf("over\n"); exit(1);}
	if(y + mindy2 >= 0){
		if(x + mindx2 >= 0){
			if(x + maxdx2 < width){
				return(roff_const2[k]);
				//roff = roff_const[k];
			}else{
				return(roff_right2[y+mindy2][x-width+maxdx2][k]);
				//roff = roff_right[y+mindy][x-width+maxdx][k];
			}
		}else{
			return(roff_left2[y+mindy2][x][k]);
			//roff = roff_left[y+mindy][x][k];
		}
	}else{
		return(roff_up2[y][x][k]);
		//roff = roff_up[y][x][k];
	}
}

uint getbits(FILE *fp, int n)
{
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);
	while (n > bitpos) {
		n -= bitpos;
		x = (x << bitpos) | bitbuf;
		bitbuf = getc(fp) & 0xff;
		bitpos = 8;
	}
	bitpos -= n;
	x = (x << n) | (bitbuf >> bitpos);
	bitbuf &= ((1 << bitpos) - 1);
	return (x);
}

void dec_make_label(DECODER *dec)//enc_make_labelと同じ
{
		FILE *fp1, *fp2;
		int x, y;
		int xg, yg, area;
		int height, width;
		height = dec->height;
		width = dec->width;

		//POINT tl, br;
		fp1 = fileopen("outfile_label.txt", "rb");//label取り込み
		fp2 = fileopen("outfile_range.txt", "rb");//range取り込み
		fscanf(fp1, "%d\n", &dec->num_of_vgen);
		fscanf(fp1, "%d\n", &dec->num_of_vgenc);
		printf("vgen,vgenc = %d,%d\n", dec->num_of_vgen, dec->num_of_vgenc);
		dec->num_of_geny = (int *)alloc_mem(dec->num_of_vgen * sizeof(int));
		dec->num_of_genx = (int *)alloc_mem(dec->num_of_vgen * sizeof(int));
		for(area = 0; area < dec->num_of_vgen; area++){
			fscanf(fp1, "%d\n",&dec->num_of_geny[area]);
			fscanf(fp1, "%d\n",&dec->num_of_genx[area]);
			printf("geny,genx = %d,%d\n",dec->num_of_geny[area], dec->num_of_genx[area]);
			if(dec->max_geny < dec->num_of_geny[area]) dec->max_geny = dec->num_of_geny[area];
			if(dec->max_genx < dec->num_of_genx[area]) dec->max_genx = dec->num_of_genx[area];
		}
		
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				fscanf(fp1, "%d\n", &dec->label[y][x].x);
				fscanf(fp1, "%d\n", &dec->label[y][x].y);
				fscanf(fp1, "%d\n", &dec->label[y][x].area);
			}
		}
		dec->gen_tl = (POINT ***)alloc_3d_array(dec->num_of_vgen, dec->max_geny, dec->max_genx, sizeof(POINT));//vseg
		dec->gen_br = (POINT ***)alloc_3d_array(dec->num_of_vgen, dec->max_geny, dec->max_genx, sizeof(POINT));//vseg
		for(area = 0; area < dec->num_of_vgen; area++){
			for(yg = 0; yg < dec->num_of_geny[area]; yg++){
				for(xg = 0; xg < dec->num_of_genx[area]; xg++){
					dec->gen_br[area][yg][xg].x = 0;
					dec->gen_br[area][yg][xg].y = 0;
					dec->gen_tl[area][yg][xg].x = width-1;
					dec->gen_tl[area][yg][xg].y = height-1;
					fscanf(fp2, "%d\n", &dec->gen_tl[area][yg][xg].x);
					fscanf(fp2, "%d\n", &dec->gen_tl[area][yg][xg].y);
					fscanf(fp2, "%d\n", &dec->gen_br[area][yg][xg].x);
					fscanf(fp2, "%d\n", &dec->gen_br[area][yg][xg].y);
					if(dec->gen_tl[area][yg][xg].x < 0){dec->gen_tl[area][yg][xg].x = 0;
					}else if(dec->gen_tl[area][yg][xg].x > dec->width - 1){dec->gen_tl[area][yg][xg].x = dec->width - 1;}
					if(dec->gen_tl[area][yg][xg].y < 0) {dec->gen_tl[area][yg][xg].y = 0;
					}else if(dec->gen_tl[area][yg][xg].y > dec->height - 1){dec->gen_tl[area][yg][xg].y = dec->height - 1;}
					if(dec->gen_br[area][yg][xg].x < 0) {dec->gen_br[area][yg][xg].x = 0;
					}else if(dec->gen_br[area][yg][xg].x > dec->width - 1){dec->gen_br[area][yg][xg].x = dec->width - 1;}
					if(dec->gen_br[area][yg][xg].y < 0) {dec->gen_br[area][yg][xg].y = 0;
					}else if(dec->gen_br[area][yg][xg].y > dec->height - 1){dec->gen_br[area][yg][xg].y = dec->height - 1;}
				}
			}
		}
		
		fclose(fp1);
		fclose(fp2);
		for(y = 0; y < dec->height; y++){
			for(x = 0; x < dec->width; x++){
				dec->area[y][x] = (char)(dec->label[y][x].area);
			}
		}
		
		return;
}

DECODER *init_decoder(FILE *fp)
{
	DECODER *dec;
	int i;
	int max;
	int cl, co, ar, gr, k, x, y;
	int height, width, num_kind_prd, max_class, maxval, max_prd_order_all, coef_precision;
	dec = (DECODER *)alloc_mem(sizeof(DECODER));
	if (getbits(fp, 16) != MAGIC_NUMBER) {
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}
	dec->version = getbits(fp, 8);
	width = dec->width = getbits(fp, 16);
	height = dec->height = getbits(fp, 16);
	maxval = dec->maxval = getbits(fp, 16);
	dec->num_comp = getbits(fp, 4);
	dec->label = (POINT **)alloc_2d_array(height, width, sizeof(POINT));//vseg
	dec->area = (char **)alloc_2d_array(height, width, sizeof(char));//追加
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			dec->label[y][x].area = 0;
			dec->area[y][x] = 0;

		}
	}
	dec_make_label(dec);//vseg//dec->label構造体各メンバーにデータ読み込み
	dec->num_class = (int *)alloc_mem(dec->num_of_vgen * sizeof(int));
	dec->max_prd_order = (int *)alloc_mem(dec->num_of_vgen * sizeof(int));
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		dec->max_prd_order[ar] = 0;
	}
	
	
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		dec->num_class[ar] = getbits(fp, 8);
	}
  max = 0;
  for(ar = 0; ar < dec->num_of_vgen; ar++){
    if(max < dec->num_class[ar]){
      max = dec->num_class[ar];
    }
  }
	max_class = dec->max_class = max;
	//dec->num_center_class = getbits(fp, 7);//vseg2class
	dec->num_group = getbits(fp, 6);
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		dec->max_prd_order[ar] = getbits(fp, 8);
		printf("max_prd_order[%d] = %d\n", ar, dec->max_prd_order[ar]);
	}
	max = 0;
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		if(max < dec->max_prd_order[ar]){
			max = dec->max_prd_order[ar];
		}
	}
	max_prd_order_all = dec->max_prd_order_all = max;
	
	//max_prd_order = dec->max_prd_order = getbits(fp, 8);
	num_kind_prd = dec->num_kind_prd = getbits(fp, 3);
	dec->num_pmodel = getbits(fp, 6) + 1;
	coef_precision = dec->coef_precision = getbits(fp, 4) + 1;
	dec->pm_accuracy = getbits(fp, 3) - 1;
	dec->f_huffman = getbits(fp, 1);
	dec->quadtree_depth = (getbits(fp, 1)) ? QUADTREE_DEPTH : -1;
	dec->maxprd = maxval << coef_precision;
	dec->max_coef = (2 << coef_precision);
	dec->predictor = (int ****)alloc_4d_array(max_class, num_kind_prd,
		dec->num_of_vgen, max_prd_order_all, sizeof(int));
	dec->nzconv = (int ****)alloc_4d_array(max_class, num_kind_prd,
		dec->num_of_vgen, max_prd_order_all, sizeof(int));
	dec->th = (int ****)alloc_4d_array(max_class, num_kind_prd,
			dec->num_of_vgen, dec->num_group - 1, sizeof(int));
	for(cl = 0; cl < max_class; cl++){
		for(co = 0; co < num_kind_prd; co++){
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for(k = 0; k < max_prd_order_all; k++){
					dec->predictor[cl][co][ar][k] = 0;
					dec->nzconv[cl][co][ar][k] = k;
				}
				for(gr = 0; gr < (dec->num_group - 1); gr++){
						dec->th[cl][co][ar][gr] = 0;
				}
			}
		}
	}
	dec->num_nzcoef = (int ***)alloc_3d_array(max_class, num_kind_prd,dec->num_of_vgen, sizeof(int));
	init_3d_array(dec->num_nzcoef, max_class, num_kind_prd, dec->num_of_vgen, 0);
	//dec->prd = (int **)alloc_2d_array(dec->height, dec->width, sizeof(int));
	dec->err = (int **)alloc_2d_array(height + 1, width, sizeof(int));
	dec->org = (int **)alloc_2d_array(height + 1, width, sizeof(int));//lfc
	init_2d_array(dec->org, height+1,width, 0);
	init_2d_array(dec->err, height+1,width, 0);
	dec->org[height][0] = (maxval + 1) >> 1; //lfc
	dec->err[height][0] = (maxval + 1) >> 2; //lfc
	dec->ctx_weight = init_ctx_weight();//lfc
	
	
	dec->num_pixel = (int ***)alloc_3d_array(dec->num_of_vgen, dec->max_geny, dec->max_genx, sizeof(int));//vseg
	dec->gen_class = (int ***)alloc_3d_array(dec->num_of_vgen, dec->max_geny, dec->max_genx, sizeof(int));//vseg

	
	if (dec->quadtree_depth > 0) {
		int xx, yy;
		yy = (height + MAX_BSIZE - 1) / MAX_BSIZE;
		xx = (width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = dec->quadtree_depth - 1; i >= 0; i--) {
			dec->qtmap[i] = (char **)alloc_2d_array(yy, xx, sizeof(char));
			for (y = 0; y < yy; y++) {
				for (x = 0; x < xx; x++) {
					dec->qtmap[i][y][x] = 0;
				}
			}
			yy <<= 1;
			xx <<= 1;
		}
	}
	dec->class = (int **)alloc_2d_array(height, width, sizeof(int));
	init_2d_array(dec->class, height, width, 0);
	if (dec->num_pmodel > 1) {
		dec->pm_idx = (int *)alloc_mem(dec->num_group * sizeof(int));
	}
	else {
		dec->pm_idx = NULL;
	}
	dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);
	if (dec->f_huffman == 1) {
		dec->sigma = sigma_h;
	}
	else {
		dec->sigma = sigma_a;
	}
	dec->mtfbuf = (int *)alloc_mem(max_class * sizeof(int));
	dec->ord2mhd = (int *)alloc_mem(max_prd_order_all * sizeof(int));
	for (i = 0; i < max_prd_order_all; i++) {
		dec->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
	}
	dec->prd_mhd = dec->ord2mhd[max_prd_order_all - 1] + 1;
	dec->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for (i = 0; i < NUM_ZMODEL; i++) {
		dec->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
	return (dec);
}

int decode_vlc(FILE *fp, VLC *vlc)
{
    int i, k, min, off;
    uint code;

    code = min = off = k = 0;
    for (i = 0; i < vlc->max_len; i++) {
	code = (code << 1) | getbits(fp, 1);
	k = vlc->off[i];
	if (k < 0) {
	    min <<= 1;
	} else {
	    if (code <= vlc->code[vlc->index[k]]) break;
	    min = (vlc->code[vlc->index[k]] + 1) << 1;
	    off = k + 1;
	}
    }
    i = off + code - min;
    return (vlc->index[i]);
}

int decode_golomb(FILE *fp, int m)
{
    int v = 0;
    while (getbits(fp, 1) == 0) {
	v++;
    }
    v = (v << m) | getbits(fp, m);
    return (v);
}



void decode_predictor2(FILE *fp, DECODER *dec)
{
    int k, m, cl, co, coef, sgn, d, zero_m, coef_m, ar;
printf("***decode_predictor***\n");
	if (dec->f_huffman == 1) {
		for (k = 0; k < dec->max_prd_order_all; k++) {
			m = getbits(fp, 4);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
					for (co = 0; co <dec->num_kind_prd; co++) {
						if(k < dec->max_prd_order[ar]){
							coef = decode_golomb(fp, m);
							if (coef > 0) {
								sgn = getbits(fp, 1);
								if (sgn) {
									coef = -coef;
								}
							}
							dec->predictor[cl][co][ar][k] = coef;
						}
					}
				}
			}
		}
	} else {
		//printf("check\n");
		PMODEL *pm;
		int b;

		pm = &dec->spm;
		pm->size = dec->max_coef + NUM_ZMODEL + 5;
		pm->cumfreq[dec->max_coef + 5] = 0;
		for(k = dec->max_coef + 5; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}
		b = dec->max_coef + 2;
		for (k = 0; k < dec->max_prd_order_all; k++) {
			//printf("check0\n");
			zero_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5,
			dec->max_coef + NUM_ZMODEL + 5) - (dec->max_coef + 5);
			//printf("zero_m[%d] = %d\n",k ,zero_m);
			coef_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + 13)
				- (dec->max_coef + 5);
			//printf("coef_m_m[%d] = %d\n",k ,coef_m);
			//printf("check2\n");
			pm->cumfreq[b] = 0;
			pm->freq[b] = TOT_ZEROFR - dec->zero_fr[zero_m];
			pm->freq[b + 1] = dec->zero_fr[zero_m];
			pm->cumfreq[b + 1] = pm->freq[b];
			pm->cumfreq[b + 2] = TOT_ZEROFR;
			set_spmodel(pm, dec->max_coef + 1, coef_m);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
					for (co = 0; co < dec->num_kind_prd; co++) {
						//printf("check3\n");
						if(k < dec->max_prd_order[ar]){
							coef = rc_decode(fp, dec->rc, pm, dec->max_coef + 2,
								dec->max_coef + 4) - (dec->max_coef + 2);
							//printf("check4\n");
							if (coef == 1) {
								//printf("check5\n");
								coef = rc_decode(fp, dec->rc, pm, 1, dec->max_coef + 1);
								//printf("check6\n");
								sgn = rc_decode(fp, dec->rc, pm, dec->max_coef + 5,
									dec->max_coef + 7) - (dec->max_coef + 5);
								//printf("check7\n");
								if (sgn) {
									coef = -coef;
								}
								dec->predictor[cl][co][ar][k] = coef;
							} else {
								dec->predictor[cl][co][ar][k] = 0;
							}
						}
					}
				}
			}
		}
	}
	
	
	//非ゼロ係数のフラグ格納
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		for (cl = 0; cl < dec->num_class[ar]; cl++) {
			for (co = 0; co < dec->num_kind_prd; co++) {
				d = 0;
				for (k = 0; k < dec->max_prd_order[ar]; k++) {
					if (dec->predictor[cl][co][ar][k] != 0) {
						dec->nzconv[cl][co][ar][d++] = k;
					}
				}
				dec->num_nzcoef[cl][co][ar] = d;
			}
		}
	}
    return;
}



void decode_predictor(FILE *fp, DECODER *dec)
{
    int k, m, cl, co, ar, coef, sgn, d;
	printf("***decode_predictor***\n");
	if (dec->f_huffman == 1) {//dec->f_huhhman = 0なので下のelse処理しか行わない
		for (k = 0; k < dec->max_prd_order_all; k++) {
			m = getbits(fp, 4);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
					for (co = 0; co < dec->num_kind_prd; co++) {
						if(k < dec->max_prd_order[ar]){
							coef = decode_golomb(fp, m);
							if (coef > 0) {
								sgn = getbits(fp, 1);
								if (sgn) {
									coef = -coef;
								}
							}
							dec->predictor[cl][co][ar][k] = coef;
						}
					}
				}
			}
		}
	} else {//算術符号化はここから
		PMODEL *pm;

		pm = &dec->spm;
		pm->size = dec->max_coef + 18;
		pm->cumfreq[dec->max_coef + 2] = 0;
		for(k = dec->max_coef + 2; k < pm->size; k++) {
	    pm->freq[k] = 1;
        pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}
		for (k = 0; k < dec->max_prd_order_all; k++) {
			m = rc_decode(fp, dec->rc, pm, dec->max_coef + 2,
                          dec->max_coef + 18) - (dec->max_coef + 2);
			set_spmodel(pm, dec->max_coef + 1, m);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
					for (co = 0; co < dec->num_kind_prd; co++) {
						if(k < dec->max_prd_order[ar]){
							coef = rc_decode(fp, dec->rc, pm, 0, dec->max_coef + 1);
							if (coef > 0) {
								sgn = rc_decode(fp, dec->rc, pm, dec->max_coef+2,
								  dec->max_coef + 4) - (dec->max_coef + 2);
								if (sgn) {
									coef = -coef;
								}
							}
							dec->predictor[cl][co][ar][k] = coef;
						}
					}
				}
			}
		}
	}
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		for (cl = 0; cl < dec->num_class[ar]; cl++) {
			for (co = 0; co < dec->num_kind_prd; co++) {
				d = 0;
				for (k = 0; k < dec->max_prd_order[ar]; k++) {
					if (dec->predictor[cl][co][ar][k] != 0) {
						dec->nzconv[cl][co][ar][d++] = k;
					}
				}
				dec->num_nzcoef[cl][co][ar] = d;
			}
		}
	}

  return;
}



void decode_threshold(FILE *fp, DECODER *dec)
{
    int cl, co, ar, gr, m, k;
		printf("***decode_threshold***\n");
    if (dec->f_huffman == 1) {
			m = getbits(fp, 4);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
	    		for (co = 0; co < dec->num_kind_prd; co++) {
          	k = 0;
	        	for (gr = 1; gr < dec->num_group; gr++) {
		    			if (k <= MAX_UPARA) {
		        		if (getbits(fp, 1)) k += decode_golomb(fp, m) + 1;
		    			}
		    			dec->th[cl][co][ar][gr - 1] = k;
						}
	    		}
				}
			}
			if (dec->num_pmodel > 1) {
	    	for (k = 1; (1 << k) < dec->num_pmodel; k++);
	    	for (gr = 0; gr < dec->num_group; gr++) {
					dec->pm_idx[gr] = getbits(fp, k);
	    	}
			}
    } else {// Arithmetic
			PMODEL *pm;
			pm = &dec->spm;
			pm->size = 16;
			pm->cumfreq[0] = 0;
			for (k = 0; k < pm->size; k++) {
	    	pm->freq[k] = 1;
	    	pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
			}
			m = rc_decode(fp, dec->rc, pm, 0, pm->size);
			set_spmodel(pm, MAX_UPARA + 2, m);
			for(ar = 0; ar < dec->num_of_vgen; ar++){
				for (cl = 0; cl < dec->num_class[ar]; cl++) {
	    		for (co = 0; co < dec->num_kind_prd; co++) {
          	k = 0;
						for (gr = 1; gr < dec->num_group; gr++) {
							if (k <= MAX_UPARA) {
								k += rc_decode(fp, dec->rc, pm, 0, pm->size - k);
							}
							dec->th[cl][co][ar][gr - 1] = k;
						}
	    		}
				}
			}
			if (dec->num_pmodel > 1) {
	    	pm->size = dec->num_pmodel;
	    	pm->freq[0] = 0;
	    	for (k = 0; k < pm->size; k++) {
					pm->freq[k] = 1;
					pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	    	}
	    	for (gr = 0; gr < dec->num_group; gr++) {
					dec->pm_idx[gr] = rc_decode(fp, dec->rc, pm, 0, pm->size);
	    	}
			}
    }
    return;
}

//改良
void decode_class_vseg(FILE *fp, DECODER *dec)
{
    int i, j, x, y, cl;
    VLC *vlc;
    PMODEL *pm, cpm[1];
		int xg, yg, area;
		printf("***decode_class_vseg***\n");
		for (yg = 0; yg < dec->max_geny; yg++) { //vsegﾂ渉可甘ｺﾂ値ﾂ静敖津ｨ
			for (xg = 0; xg < dec->max_genx; xg++) {
				for (area = 0; area < dec->num_of_vgen; area++) {
					dec->num_pixel[area][yg][xg] = 0;
				}
			}
		}

		for(y = 0; y < dec->height; y++){
			for(x = 0; x < dec->width; x++){
				dec->num_pixel[dec->label[y][x].area][dec->label[y][x].y][dec->label[y][x].x]++;
			}
		}
		cpm->cumfreq = NULL;
	  // mtf_code
		double p;
		int mtf_code[dec->max_class];

		pm = &dec->spm;
		printf("dec->num_class[0] = %d\n", dec->num_class[0]);
		
		for(area = 0; area < dec->num_of_vgen; area++){
			cpm->size = dec->num_class[area];
			cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
			for(i = 0; i < (cpm->size * 2 + 1); i++){
				cpm->freq[i] = 0;
			}
			set_spmodel(pm, PMCLASS_LEVEL, -1);
			cpm->cumfreq = &cpm->freq[cpm->size];
			cpm->cumfreq[0] = 0;
			for (i = 0; i < dec->num_class[area]; i++) {
				mtf_code[i] = rc_decode(fp, dec->rc, pm, 0, pm->size);
				if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
					for (j = 0; j < pm->size; j++) {
						if (j < mtf_code[i]) {
							pm->freq[j] /= 2;
						} else {
							pm->freq[j] *= 2;
						}
						if (pm->freq[j] <= 0) pm->freq[j] = 1;
						pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
					}
				}
			}

			for (i = 0; i < dec->num_class[area]; i++) {
				p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
													* PMCLASS_MAX/PMCLASS_LEVEL);
				cpm->freq[i] = p * (1 << 10);
				if (cpm->freq[i] <= 0) cpm->freq[i] = 1;
				cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
			}
			vlc = NULL;
		//
			for (i = 0; i < dec->max_class; i++) {
				dec->mtfbuf[i] = i;
			}

			for (yg = 0; yg < dec->num_of_geny[area]; yg++) { //vsegﾂ渉可甘ｺﾂ値ﾂ静敖津ｨ
				for (xg = 0; xg < dec->num_of_genx[area]; xg++) {
					if(dec->num_pixel[area][yg][xg] == 0) {
						dec->gen_class[area][yg][xg] = -1;
					}else{
						i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
						mtf_classlabel_vseg(dec->max_genx, dec->gen_class, dec->mtfbuf, xg, yg, area, dec->num_class);

						for (cl = 0; cl < dec->num_class[area]; cl++) {
							if (dec->mtfbuf[cl] == i) break;
						}
						dec->gen_class[area][yg][xg] = cl;
						for (y = dec->gen_tl[area][yg][xg].y; y < dec->gen_br[area][yg][xg].y + 1; y++) {//ﾂクﾂδ可スﾂ渉堕つｫﾂ債楪づ?
							for (x = dec->gen_tl[area][yg][xg].x; x < dec->gen_br[area][yg][xg].x + 1; x++) {
								if(dec->label[y][x].x == xg && dec->label[y][x].y == yg && dec->label[y][x].area == area){
									dec->class[y][x] = cl;
								}
							}
						}
					}
				}
			}
			free(cpm->freq);
		}

    if (dec->f_huffman == 1) {
				free_vlc(vlc);
    }
    return;
}

inline int calc_udec(DECODER *dec, int y, int x)
{
	int k, u;//rx, ry;// lfc
	int *wt_p, *err_p, height, width/*, ar*/;// **err;//lfc
	height = dec->height;
	width = dec->width;
	//ar = dec->area[y][x];
    u = 0;
    wt_p = dec->ctx_weight;
	err_p = &dec->err[y][x];//lfc
	for (k = 0; k < NUM_UPELS; k++){
		u += err_p[ref_offset5_2(height, width, y, x, k)] * (*wt_p++);
	}
    u >>= 6;
    if (u > MAX_UPARA) u = MAX_UPARA;
    return (u);
}

inline int calc_prd(IMAGE *img, DECODER *dec, int cl, int co, int ar, int y, int x)//lfc
{//引数で与えられた(x,y)座標の予測値を計算してprdに保存。そのprdを返す。
	int k, prd, prd_order, *coef_p, *nzc_p, l, *org_p, height, width; //rx, ry, i;//lfc
    prd_order = dec->num_nzcoef[cl][co][ar];
    prd = 0;
    coef_p = dec->predictor[cl][co][ar];
    nzc_p = dec->nzconv[cl][co][ar];
		org_p = &dec->org[y][x];//lfc
	height = dec->height;
	width = dec->width;
		for (k = 0; k < prd_order; k++){//追加
		l = nzc_p[k];
		prd += org_p[ref_offset5_2(height, width, y, x, l)] * (coef_p[l]);
		}//lfc
	
    if (prd < 0) prd = 0;
    else if (prd > dec->maxprd) prd = dec->maxprd;
		/*if(y == 0){
			printf("prd[%d][%d] = %d\n", y, x, prd);
		}*/
    return (prd);
}

IMAGE *decode_image(FILE *fp, DECODER *dec)
{
  int x, y, cl, co, ar, gr, prd, u, e, E, p;
  int *th_p, height, width, *class_p, *org_p, *err_p;
  char *area_p;
  img_t *val_p;
  int maxprd, maxval, num_group, pm_accuracy, coef_precision;
  maxprd = dec->maxprd;
  maxval = dec->maxval;
  num_group = dec->num_group;
  pm_accuracy = dec->pm_accuracy;
  coef_precision = dec->coef_precision;
  height = dec->height;
  width = dec->width;
  IMAGE *img;
	printf("***decode_image***\n");
  img = alloc_image(width, height, dec->maxval);
  if (dec->f_huffman == 1) {
		//printf("x");
		VLC *vlc;
		dec->vlcs = init_vlcs(dec->pmodels, dec->num_group, 1);
		for (y = 0; y < height; y++) {
	  	for (x = 0; x < width; x++) {
				cl = dec->class[y][x];
				co = color(y, x);
				ar = (int)dec->area[y][x];
    		u = calc_udec(dec, y, x);
				th_p = dec->th[cl][co][ar];
				for (gr = 0; gr < dec->num_group - 1; gr++) {
		    	if (u < *th_p++) break;
				}
				prd = calc_prd(img, dec, cl, co, ar, y, x);
				prd >>= (dec->coef_precision - 1);
				p = (prd + 1) >> 1;
				vlc = &dec->vlcs[gr][0];
				dec->err[y][x] = E = decode_vlc(fp, vlc);
				e = E2e(E, p, prd & 1, dec->maxval);
				img->val[y][x] = p + e;
				dec->org[y][x] = img->val[y][x]; //lfc
	  	}
		}
  } else {
		PMODEL *pm;
	if (pm_accuracy < 0) {
			//printf("y");
	  	for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
		    	cl = dec->class[y][x];
				co = color(y, x);
				ar = (int)dec->area[y][x];
		    	u = calc_udec(dec, y, x);
		    	th_p = dec->th[cl][co][ar];
		    	for (gr = 0; gr < dec->num_group - 1; gr++) {
						if (u < *th_p++) break;
		    	}
		    	prd = calc_prd(img, dec, cl, co, ar, y, x);
		    	prd >>= (dec->coef_precision - 1);
		    	p = (prd + 1) >> 1;
		    	pm = dec->pmodels[gr][0];
		    	dec->err[y][x] = E = rc_decode(fp, dec->rc, pm, 0, pm->size);
		    	e = E2e(E, p, prd & 1, dec->maxval);
		    	img->val[y][x] = p + e;
					dec->org[y][x] = img->val[y][x]; //lfc
				}
	  	}
	} else {
		//printf("z");
	    int mask, shift, base;
	    mask = (1 << pm_accuracy) - 1;
	    shift = coef_precision - pm_accuracy;
		class_p = dec->class[0];
		area_p = dec->area[0];
		org_p = dec->org[0];
		val_p = img->val[0];
		err_p = dec->err[0];
	    for (y = 0; y < height; y++) {
			//printf("y = %d\n", y);
				for (x = 0; x < width; x++) {
					//printf("(x,y) = (%d, %d)\n", x, y);
					//printf("a");
					cl = *class_p++;
					co = color(y, x);
					ar = (int)*area_p++;
					u = calc_udec(dec, y, x);
					//printf("dec->upara[%d][%d] = %d\n", y, x, u);
					th_p = dec->th[cl][co][ar];
					for (gr = 0; gr < num_group - 1; gr++) {
						//printf("dec->th[%d][%d][%d][%d] = %d\n", cl, co, ar, gr, dec->th[cl][co][ar][gr]);
						if (u < *th_p++) break;
					}
					prd = calc_prd(img, dec, cl, co, ar, y, x);
						
					base = (maxprd - prd + (1 << shift) / 2) >> shift;
					pm = dec->pmodels[gr][0] + (base & mask);
					base >>= pm_accuracy;
					//予測誤差を符号化する際、予測誤差が-255～255まであたいをとるのでレベルシフトさせて符号化している
					p = rc_decode(fp, dec->rc, pm, base, base+maxval+1) - base;//baseっていうのはレベルシフト分と予測値をあわせたもの
					//printf("p[%d][%d] = %d\n", y, x, p);
					*val_p++ = p;
					prd >>= (coef_precision - 1);
					e = (p << 1) - prd - 1;
					if (e < 0) e = -(e + 1);
					*err_p++ = e;
					*org_p++ = p; //lfc
				}
			}
		}
	}
	//exit(1);
	return (img);
}

void write_pgm(IMAGE *img, char *filename)
{
		printf("***write_pgm***\n");
    int i, j;
    FILE *fp;
    fp = fileopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n%d\n", img->width, img->height, img->maxval);
    for (i = 0; i < img->height; i++) {
	for (j = 0; j < img->width; j++) {
	    putc(img->val[i][j], fp);
	}
    }
    fclose(fp);
		printf("write pgm success.\n");
    return;
}


void print_nzcoef(DECODER *dec)
{
	int cl, co, ar, count;
	double sum, ave;
	sum = count = 0;
	ave = 0;
	for(ar = 0; ar < dec->num_of_vgen; ar++){
		for(cl = 0; cl < dec->num_class[ar]; cl++){
			for(co = 0; co < dec->num_kind_prd; co++){
				count++;
				sum += dec->num_nzcoef[cl][co][ar];				
			}
		}
	}
	ave = sum / count;
	printf("nzcoef_ave = %f\n", ave);
	
}

void print_dispersion(DECODER *dec)
{
	int cl, co, ar, k, count;
	long double sum, ave, s;
	char *name;
	name = "./LOG/dispersion.csv";
	FILE *fp;
	fp = fileopen(name, "w");
	sum = 0;
	ave = 0;
	s = 0;
	for(k = 0; k < dec->max_prd_order_all; k++){
		for(ar = 0; ar < dec->num_of_vgen; ar++){
			for(cl = 0; cl < dec->num_class[ar]; cl++){
				for(co = 0; co < dec->num_kind_prd; co++){
					count++;
					sum += dec->predictor[cl][co][ar][k];				
				}
			}
		}
		ave = sum / (long double)count;
		sum = 0;
		for(ar = 0; ar < dec->num_of_vgen; ar++){
			for(cl = 0; cl < dec->num_class[ar]; cl++){
				for(co = 0; co < dec->num_kind_prd; co++){
					sum += ((long double)dec->predictor[cl][co][ar][k] - ave) * ((long double)dec->predictor[cl][co][ar][k] - ave);				
				}
			}
		}
		s = sum / (long double)count;
		fprintf(fp, "s[%d] = %f\n", k, (double)s);
		sum = count = 0;
	}
	fclose(fp);
}

int main(int argc, char **argv)
{
    int i/*,k,cl,co,ar*/;
    IMAGE *img;
    DECODER *dec;
    char *infile, *outfile;
    FILE *fp;

    cpu_time();
printf("1\n");
    setbuf(stdout, 0);
    infile = outfile = NULL;
    for (i = 1; i < argc; i++) {
			if (infile == NULL) {
			    infile = argv[i];
			} else {
			    outfile = argv[i];
			}
    }
    if (infile == NULL || outfile == NULL) {
		  printf(BANNER"\n", 0.1 * VERSION);
			printf("usage: decmrp infile outfile\n");
			printf("infile:     Input file\n");
			printf("outfile:    Output file\n");
		        exit(0);
    }
    fp = fileopen(infile, "rb");
printf("******** Start decoding ********\n");
		dec = init_decoder(fp);//配列の動的確保
		init_ref_offset5_2(dec->max_prd_order_all, dec->height, dec->width, dec);
    if (dec->f_huffman == 0) {
			dec->rc = rc_init();//rc構造体の初期化
			rc_startdec(fp, dec->rc);
    }
    decode_class_vseg(fp, dec);
    decode_predictor2(fp, dec);
	print_dispersion(dec);
    decode_threshold(fp, dec);
    dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel,
									dec->pm_accuracy, dec->pm_idx, dec->sigma,
									dec->maxval + 1);
    img = decode_image(fp, dec);
    fclose(fp);
	print_nzcoef(dec);
    write_pgm(img, outfile);
    printf("cpu time: %.2f sec.\n", cpu_time());
    //print_predictor(dec->predictor, dec->max_prd_order_all, dec->num_class,
                    //dec->num_kind_prd, dec->max_coef, outfile);
    print_threshold(dec->num_of_vgen, dec->th, dec->num_group, dec->num_class,
		    dec->num_kind_prd, NULL, dec->pm_idx, outfile);
    print_class(dec->class, dec->height, dec->width, outfile);
    //print_color(dec->height, dec->width, dec->color, outfile);
//printf("　　　　　　　　　　　　　　　            --｡、　　　　       --｡、\n");
//printf("　　　　　　　　　　　　　　　        ／:::::::＼／: :: ::.:::＼\n");
//printf("　　　　　　　　　　       人人 ／ : : : : : : : : :.:.:: 　Ｘ＾Vヽ.\n");
//printf("　　　　　　　　   ┌―-こ　 ﾄ､   ヽ::.: :ｧ''￣: : :: .:.:　／\n");
//printf("　　　　　　　   r┘.:.:ﾄJ　＼　　' .::〃　: : :ｧ''￣ヽ:: /　　,､ 　　 )、\n");
//printf("　　　　　　　   | .:.:j.:   |:::|| ::  〃　::　　）:: ｛／　し､ 　　j入\n");
//printf("　　　　　　    |   ／:    jﾚ'_,r宀广{{广广冖=‐     ､)': ｌ: 仏ｨ'⌒ヽ\n");
//printf("　　　　　　　   |　　/  |  ..ア'    ￣ヾ小ｯ'￣ ￣`ヽ ＼.: . |:. '′　　厶\n");
//printf("　　　     _厶.o/　　| ::.／　 ,       　ｯ''￣ヽ　 ＼ ）:. |ﾟＯ:.  　　　 ｝\n");
//printf("　　　   ／　 ｏ〉   ,勹::〉′   /  /  ,ｲ!　　  ￣ヽ.   ヽ. ,尨ﾟO:. 　　　/\n");
//printf("　　　　/ ,.　o)′　　 儿′/  / ,  ' , ｲ爪　､　、 ､   ｌ 　ｌ: ,尨ミﾟO:. ..ィ\n");
//printf("  　　｛,ｲ　 /    ,r勹　 /  │ │ / /|│| │ │ | │　  │ 《冬ミﾟO　   ⌒ヽ.\n");
//printf("　　　　 辻 /    /  （　′l: |│ │ ｛ { 乂  |_│ │ |   │.VﾐOﾟ\n");
//printf(".　　　   》′   /  r勹 |  ﾚ'ﾌ仄下 ヾ（　 ヽ丁 仄卞ﾒ,ﾉ 　│:.: '｡ ｡　　   |\n");
//printf("　　　　 /　　　  └/　人  乂　[从ｧ'てうく　 　　   ｧ'てう:J　  |:.:.:｝ ﾄ､　　ト､ﾉ\n");
//printf("　　　 ,ｲ　.x(￣ヽ　Ｙ　　） 　刈 ｛{し ｨﾟ)｝'    '{し ｨﾟ)｝　 l＾Y:.:/ ,`ー′\n");
//printf(".  /|  ｛  )　　｝　ト-r仆!ﾊ  ゞ ﾆ°　　      ゞ ﾆﾟｨ′   小ﾉ::)ｲ /\n");
//printf("　　 ｌ 人　ヾ(__ノ 　人_｛从乂小、      r ｀┐　     彡ｲリ :勹′\n");
//printf("　　 レｰへ、　　_  x《 厂＼ヾ从介:｡.     ､_ﾉ　 　　 .ｨ仆.从彡'′\n");
//printf("　　 |  　ヾ, ￣　　〕｛ 　厂＼ヾ ミ辿＞;｡. __.. イ辿彡'′\n");
//printf("　　 |　　 　　 ￣＞宀ｰ へ、   ＼ｰこﾉ  └yｬ‐ｩy┘ 廴_r‐ｧ''￣￣｀\n");
//printf("　　 |      ／　　,｡　￣　＼    ヽ｢ ヽ（＾辷)＾）   　＼::;:  : :＼\n");
//printf("　　 l　   ／   ／: ::..]＼   │ │r介r'       V ::::     、\n");
//printf("    ヽ ／   ／:  : : :丨  ＼ |　｢ >'￣)      [::: :　　   |\n");
//printf("　　　　　　    （:: : : :.::｝　　 ﾞ |  ﾚ′/￣ヾ    ﾄ､厂::: :  　　八\n");
//printf("          　）　: : : : :〈_rへ 人　  ′  /￣rv―-＜:.:.: :  : ）\n");
//printf("　　　　　   　（: :   : : ::,' : └小、　　 ′ /7　rｧ  ＼:.: ::,.イ\n");
//printf("　　　　　　　  ヾ　 ＼_　:_人:: ＜￣￣ _介x　  〈/./  /  ,V￣￣`'く 八\n");
//printf("　　　      ,｀＞'⌒7＾ 7＾入:  `ｧ'￣    `ー 　ﾄ-ｲ　/  / 　｝ー　　　 〉 ﾉ\n");
//printf(".　　　   　〃:: : :′ ﾂ′ ` (    ,     (￣`''ﾍ ノ　x《　　   `く （\n");
//printf("　　　　　   `Ｚﾆイ: :x《    　　` Ｚ＿,.ｨ'   `7冖'今く　　＼ 　　　　 V\n");
//printf("　　　　　　　　    ヾ／ ’､     / |    |   /        ＼　　＼　 ┌┘\n");
//printf("   ＿＿　＿　＿   ＿　》｡   ′　　 ＼　  人_人　   人     ＼　　＼|\n");
//printf("   : :  ::.:.:::￣ ヽ,　　　　　／ l:  ::`)イ|  　ヾ,人ノ ＼　　＼\n");
///printf("   : : : -‐　 : : :￣ ＼　　／: :′  : : ｛.八        /ｰ ＼　　＼\n");
//printf("  　_　   -           　｀7´ : :/ 　　 : 、: :  ＼  ／ :    ＼　　＼　‐-　　_\n");
//printf("　　　　               : ′: : :′ : : /　　 : ＼: ｀て: :      ＼　　 ＼::::\n");
//printf("　　　　　　　　　     : : :／: : :/   : : /    : ＼ : : :　　　    ＼　　 ＼::\n");
    return(0);
}
