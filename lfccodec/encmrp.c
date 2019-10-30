#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "mrp.h"

//POINT_REF dyx3[NUM_KIND_PRD][NUM_OF_VGEN][MAX_PRD_ORDER_ALL];
//extern POINT_REF dyx[NUM_OF_VGEN][126];
extern POINT_REF dyx2[];
extern POINT lens_center[];
extern int num_ref[];
//extern POINT_REF ref_center[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];
void set_prd_pels(ENCODER *);
int mindy, mindx, maxdx;
int mindy2, mindx2, maxdx2;
//int roff_const[NUM_OF_VGEN][MAX_PRD_ORDER];
int ***roff_up;
int ***roff_left;
int ***roff_right;
int ***roff_up2;
int ***roff_left2;
int ***roff_right2;
int roff_const2[MAX_PRD_ORDER];
//int roff_const3[NUM_KIND_PRD][NUM_OF_VGEN][MAX_PRD_ORDER];
int dyx_y[MAX_PRD_ORDER];
int dyx_x[MAX_PRD_ORDER];
int dyx_p[MAX_PRD_ORDER];



double ****calc_entropy_of_conditional_probability(PMODEL ***pmodels, int num_group,
						  int num_pmodel, int pm_accuracy,
						  int maxval)
//条件確率のエントロピー計算？
{
    int i, j, k, l, gr, total;
		uint *cumfreq_p, *freq_p;
    PMODEL *pm_p, *pm;
    double *logfreq, logtotal, log2 = log(2.0);
    double **ptr1, *ptr2, ****c_prob, entropy;

    /* alloc 4D memory */
    c_prob = (double ****)alloc_2d_array(num_group, num_pmodel, sizeof(double **));
    ptr1 = (double **)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy)
			       * sizeof(double *));
    for (gr = 0; gr < num_group; gr++) {
			for (i = 0; i < num_pmodel; i++) {
	    	c_prob[gr][i] = ptr1;
	    	ptr1 += (1 << pm_accuracy);
			}
    }
    ptr2 = (double *)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy)
			      * (maxval + 1) * sizeof(double));
    for (gr = 0; gr < num_group; gr++) {
			for (i = 0; i < num_pmodel; i++) {
	    	for (j = 0; j < (1 << pm_accuracy); j++) {
					c_prob[gr][i][j] = ptr2;
					ptr2 += maxval + 1;
	    	}
			}
    }
    /* calc entropy of conditional probability */
    logfreq = (double *)alloc_mem(((maxval << 1) + 1) * sizeof(double));
    for (gr = 0; gr < num_group; gr++) {
			for (i = 0; i < num_pmodel; i++) {
	    	pm_p = pmodels[gr][i];
	    	for (j = 0; j < (1 << pm_accuracy); j++) {
					pm = pm_p + j;
					freq_p = pm->freq;
					cumfreq_p = pm->cumfreq;
					for (k = 0; k < (maxval << 1) + 1; k++) {
		    		logfreq[k] = log(freq_p[k]);
					}
					for (k = 0; k < maxval + 1; k++) {
		    		total = (cumfreq_p[k + maxval + 1]) - (cumfreq_p[k]);
		    		logtotal = log(total);
		    		entropy = 0;
		    		for (l = k; l < k + maxval + 1; l++) {
							entropy += (freq_p[l]) * (logtotal - logfreq[l]);
		    		}
		    		entropy /= (total * log2);
		    		c_prob[gr][i][j][k] = entropy;
					}
	    	}
			}
    }
    return c_prob;
}

void calc_ratio_of_model_to_rate(ENCODER *enc)
{
    int y, x, prd, gr, e, base, frac, /*cl,*/ totpel, *gr_pel;
    double *entropy, *cost, *ratio, totent, totcost, totratio;
    PMODEL *pm;
    int calc_entropy_flag = 0;
    double ****entropy_of_conditional_probability;


    if (calc_entropy_flag == 0) {
			entropy_of_conditional_probability =
			  calc_entropy_of_conditional_probability(enc->pmodels, enc->num_group,
								  enc->num_pmodel,enc->pm_accuracy,enc->maxval);
			calc_entropy_flag = 1;
    }
    entropy = (double *)alloc_mem(enc->num_group * sizeof(double));
    cost = (double *)alloc_mem(enc->num_group * sizeof(double));
    ratio = (double *)alloc_mem(enc->num_group * sizeof(double));
    gr_pel = (int *)alloc_mem(enc->num_group * sizeof(int));
    for (gr = 0; gr < enc->num_group; gr++) {
			entropy[gr] = 0;//グループ毎のモデル関数のエントロピー？
			cost[gr] = 0;//グループ毎の実際の符号量？
			gr_pel[gr] = 0;//グループ毎の画素数
    }
    for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
			    //cl = enc->class[y][x];
			    gr = enc->group[y][x];
			    e = enc->encval[y][x];//eには輝度値が入るが予測値で山を立てるので、予測誤差を符号化していることと数学的には同義
			    prd = enc->prd[y][x];
			    if (prd < 0) prd = 0;
			    else if (prd > enc->maxprd) prd = enc->maxprd;
			    base = enc->bconv[prd];
			    frac = enc->fconv[prd];
			    pm = enc->pmlist[gr] + frac;
			    cost[gr] += pm->cost[base + e] + pm->subcost[base];//符号量を計算
			    entropy[gr] += entropy_of_conditional_probability[gr][pm->id][frac][base];
			    gr_pel[gr]++;
			}
    }
    /* calc ratio */
    totpel = 0;
    totcost = totent = 0;
    for (gr = 0; gr < enc->num_group; gr++){
		//printf("entropy[%d] = %f\n", gr, entropy[gr]);
		//printf("cost[%d] = %f\n", gr, cost[gr]);
			if (entropy[gr] != 0.0)
			  	ratio[gr] = 1.0 - fabs(entropy[gr] - cost[gr]) / entropy[gr];
			else
			  	ratio[gr] = 1.0;
			totent += entropy[gr];
			totcost += cost[gr];
			totpel += gr_pel[gr];
    }
    totratio = 1.0 - fabs(totent - totcost) / totent;
    /* print */
    printf("******* differences between entropy and rate *******\n");
    printf("(gr)  [shape]\tentropy\t\t|\trate\t\t|fitness\t|pel\n");
    printf("------------------------------------------------------------------------------\n");
    for (gr = 0; gr < enc->num_group; gr++) {
			printf("(%2d)  [%.1f]\t%10.2f\t|\t%10.2f\t|   %.3f\t|%10d\n", gr,
			       0.2 * (enc->pmlist[gr]->id + 1), entropy[gr],
			       cost[gr], ratio[gr], gr_pel[gr]);
    }
    printf("------------------------------------------------------------------------------\n");
    printf("all         \t%10.2f\t|\t%10.2f\t|   %.3f\t|%10d\n", totent,
	   totcost, totratio, totpel);
    free(entropy);
    free(cost);
    free(ratio);
    free(gr_pel);
}

IMAGE *read_pgm(char *filename)
{
    int i, j, width, height, maxval;
    char tmp[256];
    IMAGE *img;
    FILE *fp;

    fp = fileopen(filename, "rb");
    fgets(tmp, 256, fp);
    if (tmp[0] != 'P' || tmp[1] != '5') {
	fprintf(stderr, "Not a PGM file!\n");
	exit(1);
    }
    while (*(fgets(tmp, 256, fp)) == '#');
    sscanf(tmp, "%d %d", &width, &height);
    while (*(fgets(tmp, 256, fp)) == '#');
    sscanf(tmp, "%d", &maxval);
    if ((width % BASE_BSIZE) || (height % BASE_BSIZE)) {
	fprintf(stderr, "Image width and height must be multiples of %d!\n",
		BASE_BSIZE);
	exit(1);
    }
    if (maxval > 255) {
	fprintf(stderr, "Sorry, this version only supports 8bpp images!\n");
	exit(1);
    }
    img = alloc_image(width, height, maxval);
    for (i = 0; i < img->height; i++) {
	for (j = 0; j < img->width; j++) {
	    img->val[i][j] = (img_t)fgetc(fp);
	}
    }
    fclose(fp);
    return (img);
}



/*void init_ref_offset4(int prd_order, int width)
{	
	int k, dy, dx;
	mindy = 0;
	mindx = 0;
	maxdx = 0;
	for (k = 0; k < prd_order; k++) {
		dy = dyx_y[k] = dyx[k].y;
		dx = dyx_x[k] = dyx[k].x;
		dyx_p[k] = dyx[k].p;
		if (dy < mindy) mindy = dy;
		if (dx < mindx) mindx = dx;
		if (dx > maxdx) maxdx = dx;
		roff_const[k] = dy * width + dx;
	}
	//mindy--;
	//mindx--;
	//maxdx++;
	printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
}*/

/*void fill_in_dyx_p(POINT_REF ***d, int max_prd_order_all)
{
	int co, ar, k;
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < NUM_OF_VGEN; ar++){
			for(k = 0; k < max_prd_order_all; k++){
				if(abs(d[co][ar][k].y) % 2 == 0){
					if(abs(d[co][ar][k].x) % 2 == 0){
						d[co][ar][k].p = 0;
					}else{
						d[co][ar][k].p = 1;
					}
				}else{
					if(abs(d[co][ar][k].x) % 2 == 0){
						d[co][ar][k].p = 3;
					}else{
						d[co][ar][k].p = 2;
					}
				}
			}
		}
	}
}*/

/*int cmp( const void *p, const void *q ) {
    int i;
	double dist_xp=0.0, dist_yp=0.0, dist_xq=0.0, dist_yq=0.0;
	double dist_p=0.0, dist_q=0.0;
	i=0;
	if(((POINT_CORRELATION*)p)->coef > ((POINT_CORRELATION*)q)->coef){
		i=-1;
	}else if(((POINT_CORRELATION*)p)->coef < ((POINT_CORRELATION*)q)->coef){
		i=1;
	}else{
		dist_xp = (double)(lens_center[((POINT_CORRELATION*)p)->lens_n].x - ((POINT_CORRELATION*)p)->x);
		dist_xq = (double)(lens_center[((POINT_CORRELATION*)q)->lens_n].x - ((POINT_CORRELATION*)q)->x);
		dist_yp = (double)(lens_center[((POINT_CORRELATION*)p)->lens_n].y - ((POINT_CORRELATION*)p)->y);
		dist_yq = (double)(lens_center[((POINT_CORRELATION*)q)->lens_n].y - ((POINT_CORRELATION*)q)->y);
		dist_p = dist_xp * dist_xp + dist_yp * dist_yp;
		dist_q = dist_xq * dist_xq + dist_yq * dist_yq;
		if(dist_p > dist_q) i=1;
		else if(dist_p > dist_q) i=-1;
		else i=0;
	}
	return(i);
}*/

/*void make_auto_ref_offset(ENCODER *enc){
	printf("make_auto_ref_offset\n");
	int i, j, k, dy, dx, tlx, tly, brx, bry, mask_count=0, range_y,range_x, range_k, range;
	int ar, co;
	int num_lens, n, mask_mhd, *mask;
	int ****p_ref;
	POINT_REF *d_p[NUM_KIND_PRD][NUM_OF_VGEN], **d_pp[NUM_KIND_PRD];
	num_lens = NUM_REF_LENS;
	mask_mhd = LENS_MASK_MHD;//マスクの範囲を決めるマンハッタンディスタンス
	range_y = SEARCH_RANGE_Y;
	range_x = SEARCH_RANGE_X;
	range = range_y * range_x;
	tly = -range_y+1;
	tlx = -range_x/2;
	bry = 0+1;
	brx = range_x/2;
	i = j = 0;
	mask = (int *)alloc_mem(range * sizeof(int));
	for(range_k = 0; range_k < range; range_k++) mask[range_k] = -1;
	p_ref = (int ****)alloc_4d_array(NUM_KIND_PRD, NUM_OF_VGEN, range_y, range_x, sizeof(int));
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < NUM_OF_VGEN; ar++){
			for(i = 0; i < range_y;i++){
				for(j = 0; j < range_x; j++){
					p_ref[co][ar][i][j] = 0;
					if(i == range_y-1 && j >= range_x/2) p_ref[co][ar][i][j] = -1;
				}
			}
		}
	}
	
	POINT_CORRELATION *ref_correlation;
	ref_correlation = (POINT_CORRELATION *)alloc_mem(1000 * sizeof(POINT_CORRELATION));
	
	for(n = 0; n < num_lens; n++){
		for(dy = tly; dy < bry; dy++){
			for(dx = tlx; dx < brx; dx++){
				if(dy == 0 && dx >= 0) continue;
				if(abs(lens_center[n].y - dy) + abs(lens_center[n].x - dx) <= mask_mhd){
					range_k = (dy + range_y-1) * range_x + (dx + range_x/2);
					mask[range_k] = n;
				}
			}
		}
	}
	for(range_k=0; range_k < range; range_k++){
		printf("%2d ", mask[range_k]);
		if(range_k % range_x == range_x-1)printf("\n");
	}
	for(co = 0; co < NUM_KIND_PRD; co++){
		printf("co = %d\n",co);
		for(ar = 0; ar < NUM_OF_VGEN; ar++){
			printf("ar = %d\n",ar);
			mask_count=0;
			for(dy = tly; dy < bry; dy++){
				for(dx = tlx; dx < brx; dx++){
					if(dy == 0 && dx >= 0)continue;
					range_k = (dy + range_y-1) * range_x + (dx + range_x/2);
					if(mask[range_k] >= 0){
						ref_correlation[mask_count].y = dy;
						ref_correlation[mask_count].x = dx;
						ref_correlation[mask_count].lens_n = mask[range_k];
						ref_correlation[mask_count].coef = enc->correlation_coef[co][ar][range_k];
						//printf("dy, dx, coef[%d] = %lf, lens_n=%d\n", dy, dx, mask_count, ref_correlation[mask_count].coef,ref_correlation[mask_count].lens_n);
						mask_count++;
						//printf("dy ,dx = %d,%d\n", dy, dx);
					}
				}
			}

			qsort(ref_correlation, mask_count, sizeof(POINT_CORRELATION), cmp);
			for(k = 0; k < MAX_PRD_ORDER_ALL; k++){	
				printf("dy, dx, ref_correlation[%d].coef = %d,%d,%lf\n", 
				k, ref_correlation[k].y, ref_correlation[k].x, ref_correlation[k].coef);
				dyx3[co][ar][k].y = ref_correlation[k].y;
				dyx3[co][ar][k].x = ref_correlation[k].x;
			}
			//}
			for(k = 0; k < MAX_PRD_ORDER_ALL; k++){
				printf("dyx3 (dy,dx) = (%d,%d)\n", dyx3[co][ar][k].y, dyx3[co][ar][k].x);
			}
			d_p[co][ar] = dyx3[co][ar];
			d_pp[co] = d_p[co];
			for(k = 0; k < MAX_PRD_ORDER_ALL; k++){
				p_ref[co][ar][range_y-1+dyx3[co][ar][k].y][range_x/2+dyx3[co][ar][k].x] = 1;
			}
		}
	}
	printf("check0\n");
	fill_in_dyx_p(d_pp, MAX_PRD_ORDER_ALL);
	printf("check1\n");
	FILE *fp;
	fp = fileopen("./LOG/print_ref.csv", "wb");

	for(ar = 0; ar < NUM_OF_VGEN; ar++){
		fprintf(fp, "ar = %d\n", ar);
		for(co = 0; co < NUM_KIND_PRD; co++){
			fprintf(fp, "co = %d\n", co);
			for(i = 0; i < range_y; i++){
				for(j = 0; j < range_x; j++){
					fprintf(fp, "%d,", p_ref[co][ar][i][j]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	printf("check2\n");
	fclose(fp);
	free(ref_correlation);
	free(mask);
}*/



/*void init_ref_offset5(POINT_REF **d_pp, int *prd_order, int height, int width, ENCODER *enc)
{	
	printf("init_ref_offset5->");
	int k, dy, dx, y, x, flag, temp_k, temp_dy, temp_dx, yy, xx, ar;
	mindy = 0;
	mindx = 0;
	maxdx = 0;
	for(ar = 0; ar < NUM_OF_VGEN; ar++){
		for (k = 0; k < prd_order[ar]; k++) {
			dy = d_pp[ar][k].y;
			dx = d_pp[ar][k].x;
			if(dy > 0) printf("you must change dy dyx[%d][%d] in common.c\n", ar, k);
			if(dy == 0 && dx >= 0) printf("you must change dy or dx dyx[%d][%d] in common.c\n", ar, k);
			if (dy < mindy) mindy = dy;
			if (dx < mindx) mindx = dx;
			if (dx > maxdx) maxdx = dx;
			roff_const[ar][k] = dy * width + dx;
		}
	}
	
	roff_up = (int ***)alloc_3d_array(abs(mindy), width, enc->max_prd_order_all, sizeof(int));
	roff_left = (int ***)alloc_3d_array(height-abs(mindy), abs(mindx), enc->max_prd_order_all, sizeof(int));
	roff_right = (int ***)alloc_3d_array(height-abs(mindy), abs(maxdx), enc->max_prd_order_all, sizeof(int));
	for(k = 0; k < enc->max_prd_order_all; k++){
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
	
	//printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
	printf("roff_up->");
	for(y = 0; y < abs(mindy); y++){
		for(x = 0; x < width; x++){
			ar = enc->area[y][x];
			for(k = 0; k < prd_order[ar]; k++){
				if(y == 0 && x == 0){
					dx = 0;
					dy = height;
					roff_up[y][x][k] = dy * width + dx;
				}else{
					dy = d_pp[ar][k].y;
					dx = d_pp[ar][k].x;
					if(y + dy < 0 || x + dx < 0 || x + dx >= width){
						flag = 0;
						temp_k = k;
						temp_k--;
						if(color(y, x) == 1 || color(y, x) == 2){
							while(temp_k >= 0){
								temp_dy = d_pp[ar][temp_k].y;
								temp_dx = d_pp[ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
										roff_up[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
							}
						}else{
							while(temp_k >= 0){
								temp_dy = d_pp[ar][temp_k].y;
								temp_dx = d_pp[ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(d_pp[ar][k].p == 0 || d_pp[ar][k].p == 2){
										if(d_pp[ar][temp_k].p == 0 || d_pp[ar][temp_k].p == 2){
											roff_up[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}else{
										if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
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
				//printf("roff_up[%d][%d][%d] = %d\n", y, x, k, roff_up[y][x][k]);
			}
		}
	}
	printf("roff_left->");
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(mindx); x++){
			ar = enc->area[yy][x];
			for(k = 0; k < prd_order[ar]; k++){
				dy = d_pp[ar][k].y;
				dx = d_pp[ar][k].x;
				if(yy + dy < 0 || x + dx < 0 || x + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, x) == 1 || color(yy, x) == 2){
						while(temp_k >= 0){
							temp_dy = d_pp[ar][temp_k].y;
							temp_dx = d_pp[ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
									roff_left[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = d_pp[ar][temp_k].y;
							temp_dx = d_pp[ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(d_pp[ar][k].p == 0 || d_pp[ar][k].p == 2){
									if(d_pp[ar][temp_k].p == 0 || d_pp[ar][temp_k].p == 2){
										roff_left[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
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
	printf("roff_right->");
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(maxdx); x++){
			xx = x + (width-abs(maxdx));
			ar = enc->area[yy][xx];
			for(k = 0; k < prd_order[ar]; k++){
				dy = d_pp[ar][k].y;
				dx = d_pp[ar][k].x;
				if(yy + dy < 0 || xx + dx < 0 || xx + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, xx) == 1 || color(yy, xx) == 2){
						while(temp_k >= 0){
							temp_dy = d_pp[ar][temp_k].y;
							temp_dx = d_pp[ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
									roff_right[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = d_pp[ar][temp_k].y;
							temp_dx = d_pp[ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(d_pp[ar][k].p == 0 || d_pp[ar][k].p == 2){
									if(d_pp[ar][temp_k].p == 0 || d_pp[ar][temp_k].p == 2){
										roff_right[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(d_pp[ar][k].p == d_pp[ar][temp_k].p){
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
	printf("ok\n");
}*/

/*void init_ref_offset5_eachcolor(POINT_REF ***d_ppp, int *prd_order, int height, int width, ENCODER *enc)
{	
	printf("init_ref_offset5_eachcolor\n");
	int k, dy, dx, y, x, flag, temp_k, temp_dy, temp_dx, yy, xx, ar, co;
	mindy = 0;
	mindx = 0;
	maxdx = 0;
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < NUM_OF_VGEN; ar++){
			for (k = 0; k < prd_order[ar]; k++) {
				dy = d_ppp[co][ar][k].y;
				dx = d_ppp[co][ar][k].x;
				if(dy > 0) printf("you must change dy dyx[%d][%d] in common.c\n", ar, k);
				if(dy == 0 && dx >= 0) printf("you must change dy or dx dyx[%d][%d] in common.c\n", ar, k);
				if (dy < mindy) mindy = dy;
				if (dx < mindx) mindx = dx;
				if (dx > maxdx) maxdx = dx;
				roff_const3[co][ar][k] = dy * width + dx;
			}
		}
	}
	
	roff_up = (int ***)alloc_3d_array(abs(mindy), width, enc->max_prd_order_all, sizeof(int));
	roff_left = (int ***)alloc_3d_array(height-abs(mindy), abs(mindx), enc->max_prd_order_all, sizeof(int));
	roff_right = (int ***)alloc_3d_array(height-abs(mindy), abs(maxdx), enc->max_prd_order_all, sizeof(int));
	for(k = 0; k < enc->max_prd_order_all; k++){
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
	
	//printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
	printf("roff_up\n");
	for(y = 0; y < abs(mindy); y++){
		for(x = 0; x < width; x++){
			co = color(y,x);
			ar = enc->area[y][x];
			for(k = 0; k < prd_order[ar]; k++){
				if(y == 0 && x == 0){
					dx = 0;
					dy = height;
					roff_up[y][x][k] = dy * width + dx;
				}else{
					dy = d_ppp[co][ar][k].y;
					dx = d_ppp[co][ar][k].x;
					if(y + dy < 0 || x + dx < 0 || x + dx >= width){
						flag = 0;
						temp_k = k;
						temp_k--;
						if(color(y, x) == 1 || color(y, x) == 2){
							while(temp_k >= 0){
								temp_dy = d_ppp[co][ar][temp_k].y;
								temp_dx = d_ppp[co][ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
										roff_up[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
							}
						}else{
							while(temp_k >= 0){
								temp_dy = d_ppp[co][ar][temp_k].y;
								temp_dx = d_ppp[co][ar][temp_k].x;
								if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(d_ppp[co][ar][k].p == 0 || d_ppp[co][ar][k].p == 2){
										if(d_ppp[co][ar][temp_k].p == 0 || d_ppp[co][ar][temp_k].p == 2){
											roff_up[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}else{
										if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
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
				//printf("roff_up[%d][%d][%d] = %d\n", y, x, k, roff_up[y][x][k]);
			}
		}
	}
	printf("roff_left\n");
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(mindx); x++){
			co = color(yy, x);
			ar = enc->area[yy][x];
			for(k = 0; k < prd_order[ar]; k++){
				dy = d_ppp[co][ar][k].y;
				dx = d_ppp[co][ar][k].x;
				if(yy + dy < 0 || x + dx < 0 || x + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, x) == 1 || color(yy, x) == 2){
						while(temp_k >= 0){
							temp_dy = d_ppp[co][ar][temp_k].y;
							temp_dx = d_ppp[co][ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
									roff_left[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = d_ppp[co][ar][temp_k].y;
							temp_dx = d_ppp[co][ar][temp_k].x;
							if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
								if(d_ppp[co][ar][k].p == 0 || d_ppp[co][ar][k].p == 2){
									if(d_ppp[co][ar][temp_k].p == 0 || d_ppp[co][ar][temp_k].p == 2){
										roff_left[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
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
			co = color(yy,xx);
			ar = enc->area[yy][xx];
			for(k = 0; k < prd_order[ar]; k++){
				dy = d_ppp[co][ar][k].y;
				dx = d_ppp[co][ar][k].x;
				if(yy + dy < 0 || xx + dx < 0 || xx + dx >= width){
					flag = 0;
					temp_k = k;
					temp_k--;
					if(color(yy, xx) == 1 || color(yy, xx) == 2){
						while(temp_k >= 0){
							temp_dy = d_ppp[co][ar][temp_k].y;
							temp_dx = d_ppp[co][ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
									roff_right[y][x][k] = temp_dy * width + temp_dx;
									flag = 1;
									break;
								}
							}
							temp_k--;
						}
					}else{
						while(temp_k >= 0){
							temp_dy = d_ppp[co][ar][temp_k].y;
							temp_dx = d_ppp[co][ar][temp_k].x;
							if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
								if(d_ppp[co][ar][k].p == 0 || d_ppp[co][ar][k].p == 2){
									if(d_ppp[co][ar][temp_k].p == 0 || d_ppp[co][ar][temp_k].p == 2){
										roff_right[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}else{
									if(d_ppp[co][ar][k].p == d_ppp[co][ar][temp_k].p){
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
}*/

void init_ref_offset5_2(int prd_order, int height, int width, ENCODER *enc)
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

/*void init_ref_offset6(int prd_order, int height, int width, ENCODER *enc)
{	
	printf("init_ref_offset6\n");
	int k, dy, dx, y, x, flag, temp_k, temp_dy, temp_dx, yy, xx;
	int rx, ry, id, temp_rx, temp_ry, temp_id, temp_roff, max_roff;
	double dist, min_dist, dist_x, dist_y, error;
	int k_buf[prd_order], k_buf2[prd_order], id_buf[prd_order], i=0, count=0, count2=0, count3=0;
	mindy = 0;
	mindx = 0;
	maxdx = 0;
	error = 2.0;
	for (k = 0; k < prd_order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;
		if(dy > 0) printf("you must change dyx[] in common.c\n");
		if(dy == 0 && dx >= 0) printf("you must change dyx[] in common.c\n");
		if (dy < mindy) mindy = dy;
		if (dx < mindx) mindx = dx;
		if (dx > maxdx) maxdx = dx;
		roff_const[k] = dy * width + dx;
	}
	roff_up = (int ***)alloc_3d_array(abs(mindy), width, prd_order, sizeof(int));
	roff_left = (int ***)alloc_3d_array(height-abs(mindy), abs(mindx), prd_order, sizeof(int));
	roff_right = (int ***)alloc_3d_array(height-abs(mindy), abs(maxdx), prd_order, sizeof(int));
	//printf("mindy = %d, mindx = %d, maxdx = %d\n", mindy, mindx, maxdx);
	printf("roff_up\n");
	for(y = 0; y < abs(mindy); y++){
		for(x = 0; x < width; x++){
			for(k = 0; k < prd_order; k++){
				if(y == 0 && x == 0){
					dx = 0;
					dy = height;
					roff_up[y][x][k] = dy * width + dx;
				}else{
					dy = dyx[k].y;
					dx = dyx[k].x;
					//printf("y,x,k = %d,%d,%d\n", y,x,k);
					if(y + dy < 0 || x + dx < 0 || x + dx >= width){
						ry = dyx[k].ry;
						rx = dyx[k].rx;
						count = 0;
						for(temp_k = 0; temp_k < prd_order; temp_k++){
							if(temp_k == k) continue;
							temp_dy = dyx[temp_k].y;
							temp_dx = dyx[temp_k].x;
							temp_roff = temp_dy * width + temp_dx;
							dist_y = dyx[k].ry - dyx[temp_k].ry;
							dist_x = dyx[k].rx - dyx[temp_k].rx;
							dist = sqrt(dist_y * dist_y + dist_x * dist_x);
							if(y + temp_dy >= 0 &&  x + temp_dx >= 0 && x + temp_dx < width && temp_roff < 0){
								if(dist < error && dyx[k].p == dyx[temp_k].p){
									k_buf[count] = temp_k;
									id_buf[count] = dyx[temp_k].id;
									count++;
								}
							}
							
						}
						
						min_dist = DBL_MAX;
						if(count != 0){
							
							count2 = 0;
							for(i = 0; i < count; i++){
								dist_y = ref_center[dyx[k].id].y - ref_center[id_buf[i]].y;
								dist_x = ref_center[dyx[k].id].x - ref_center[id_buf[i]].x;
								dist = sqrt(dist_y * dist_y + dist_x * dist_x);
								if(dist < min_dist){
									min_dist = dist;
									count2 = 0;
									k_buf2[count2] = k_buf[i];
									count2++;
								}else if(dist == min_dist){
									k_buf2[count2++] = k_buf[i];
								}
							}
							
							if(count2 == 1){
								temp_k = k_buf2[0];
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
							}else{
								
								min_dist = DBL_MAX;
								count3 = 0;
								for(i = 0; i < count2; i++){
									dist_y = dy - dyx[k_buf2[i]].y;
									dist_x = dx - dyx[k_buf2[i]].x;
									dist = sqrt(dist_y * dist_y + dist_x * dist_x);
									if(dist < min_dist){
										min_dist = dist;
										count3 = 0;
										k_buf[count3++] = k_buf2[i];
									}else if(dist == min_dist){
										k_buf[count3++] = k_buf2[i];
									}
								}
								
								if(count3 == 1){
									temp_k = k_buf[0];
									temp_dy = dyx[temp_k].y;
									temp_dx = dyx[temp_k].x;
								}else{
									
									max_roff = -INT_MAX;
									for(i = 0; i < count3; i++){
										temp_roff = dyx[k_buf[i]].y * width + dyx[k_buf[i]].x;
										if(temp_roff > max_roff){
											max_roff = temp_roff;
											temp_k = k_buf[i];
										}
									}
									temp_dy = dyx[temp_k].y;
									temp_dx = dyx[temp_k].x;
								}
							}
						}else{
							//if(y == 0 && x == 1 && k == 6) printf("check\n");
							flag = 0;
							temp_k = k;
							temp_k--;
					
								while(temp_k >= 0){
									temp_dy = dyx[temp_k].y;
									temp_dx = dyx[temp_k].x;
									if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
										if(dyx[k].p == dyx[temp_k].p){
											roff_up[y][x][k] = temp_dy * width + temp_dx;
											flag = 1;
											break;
										}
									}
									temp_k--;
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
						}
						
						roff_up[y][x][k] = temp_dy * width + temp_dx;
					}else{
						roff_up[y][x][k] = dy * width + dx;
					}
				}
				//printf("roff_up[%d][%d][%d] = %d\n", y, x, k, roff_up[y][x][k]);
			}
		}
	}
	
	printf("roff_left\n");
	//exit(1);
	for(y = 0; y < height-abs(mindy); y++){
		yy = y + abs(mindy);
		for(x = 0; x < abs(mindx); x++){
			for(k = 0; k < prd_order; k++){
				dy = dyx[k].y;
				dx = dyx[k].x;
				//printf("y,x,k = %d,%d,%d\n", y,x,k);
				if(yy + dy < 0 || x + dx < 0 || x + dx >= width){
					ry = dyx[k].ry;
					rx = dyx[k].rx;
					count = 0;
					for(temp_k = 0; temp_k < prd_order; temp_k++){
						if(temp_k == k) continue;
						temp_dy = dyx[temp_k].y;
						temp_dx = dyx[temp_k].x;
						temp_roff = temp_dy * width + temp_dx;
						dist_y = dyx[k].ry - dyx[temp_k].ry;
						dist_x = dyx[k].rx - dyx[temp_k].rx;
						dist = sqrt(dist_y * dist_y + dist_x * dist_x);
						if(yy + temp_dy >= 0 &&  x + temp_dx >= 0 && x + temp_dx < width && temp_roff < 0){
							if(dist < error && dyx[k].p == dyx[temp_k].p){
								k_buf[count] = temp_k;
								id_buf[count] = dyx[temp_k].id;
								count++;
							}
						}
						
					}
					
					min_dist = DBL_MAX;
					if(count != 0){
						
						count2 = 0;
						for(i = 0; i < count; i++){
							dist_y = ref_center[dyx[k].id].y - ref_center[id_buf[i]].y;
							dist_x = ref_center[dyx[k].id].x - ref_center[id_buf[i]].x;
							dist = sqrt(dist_y * dist_y + dist_x * dist_x);
							if(dist < min_dist){
								min_dist = dist;
								count2 = 0;
								k_buf2[count2] = k_buf[i];
								count2++;
							}else if(dist == min_dist){
								k_buf2[count2++] = k_buf[i];
							}
						}
						
						if(count2 == 1){
							temp_k = k_buf2[0];
							temp_dy = dyx[temp_k].y;
							temp_dx = dyx[temp_k].x;
						}else{
							
							min_dist = DBL_MAX;
							count3 = 0;
							for(i = 0; i < count2; i++){
								dist_y = dy - dyx[k_buf2[i]].y;
								dist_x = dx - dyx[k_buf2[i]].x;
								dist = sqrt(dist_y * dist_y + dist_x * dist_x);
								if(dist < min_dist){
									min_dist = dist;
									count3 = 0;
									k_buf[count3++] = k_buf2[i];
								}else if(dist == min_dist){
									k_buf[count3++] = k_buf2[i];
								}
							}
							
							if(count3 == 1){
								temp_k = k_buf[0];
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
							}else{
								
								max_roff = -INT_MAX;
								for(i = 0; i < count3; i++){
									temp_roff = dyx[k_buf[i]].y * width + dyx[k_buf[i]].x;
									if(temp_roff > max_roff){
										max_roff = temp_roff;
										temp_k = k_buf[i];
									}
								}
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
							}
						}
					}else{
						//if(y == 0 && x == 1 && k == 6) printf("check\n");
						flag = 0;
						temp_k = k;
						temp_k--;
					
							while(temp_k >= 0){
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
								if(yy + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
									if(dyx[k].p == dyx[temp_k].p){
										roff_up[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
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
			for(k = 0; k < prd_order; k++){
				dy = dyx[k].y;
				dx = dyx[k].x;
				//printf("y,x,k = %d,%d,%d\n", y,x,k);
				if(yy + dy < 0 || xx + dx < 0 || xx + dx >= width){
					ry = dyx[k].ry;
					rx = dyx[k].rx;
					count = 0;
					for(temp_k = 0; temp_k < prd_order; temp_k++){
						if(temp_k == k) continue;
						temp_dy = dyx[temp_k].y;
						temp_dx = dyx[temp_k].x;
						temp_roff = temp_dy * width + temp_dx;
						dist_y = dyx[k].ry - dyx[temp_k].ry;
						dist_x = dyx[k].rx - dyx[temp_k].rx;
						dist = sqrt(dist_y * dist_y + dist_x * dist_x);
						if(yy + temp_dy >= 0 &&  x + temp_dx >= 0 && x + temp_dx < width && temp_roff < 0){
							if(dist < error && dyx[k].p == dyx[temp_k].p){
								k_buf[count] = temp_k;
								id_buf[count] = dyx[temp_k].id;
								count++;
							}
						}
						
					}
					
					min_dist = DBL_MAX;
					if(count != 0){
						
						count2 = 0;
						for(i = 0; i < count; i++){
							dist_y = ref_center[dyx[k].id].y - ref_center[id_buf[i]].y;
							dist_x = ref_center[dyx[k].id].x - ref_center[id_buf[i]].x;
							dist = sqrt(dist_y * dist_y + dist_x * dist_x);
							if(dist < min_dist){
								min_dist = dist;
								count2 = 0;
								k_buf2[count2] = k_buf[i];
								count2++;
							}else if(dist == min_dist){
								k_buf2[count2++] = k_buf[i];
							}
						}
						
						if(count2 == 1){
							temp_k = k_buf2[0];
							temp_dy = dyx[temp_k].y;
							temp_dx = dyx[temp_k].x;
						}else{
							
							min_dist = DBL_MAX;
							count3 = 0;
							for(i = 0; i < count2; i++){
								dist_y = dy - dyx[k_buf2[i]].y;
								dist_x = dx - dyx[k_buf2[i]].x;
								dist = sqrt(dist_y * dist_y + dist_x * dist_x);
								if(dist < min_dist){
									min_dist = dist;
									count3 = 0;
									k_buf[count3++] = k_buf2[i];
								}else if(dist == min_dist){
									k_buf[count3++] = k_buf2[i];
								}
							}
							
							if(count3 == 1){
								temp_k = k_buf[0];
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
							}else{
								
								max_roff = -INT_MAX;
								for(i = 0; i < count3; i++){
									temp_roff = dyx[k_buf[i]].y * width + dyx[k_buf[i]].x;
									if(temp_roff > max_roff){
										max_roff = temp_roff;
										temp_k = k_buf[i];
									}
								}
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
							}
						}
					}else{
						//if(y == 0 && x == 1 && k == 6) printf("check\n");
						flag = 0;
						temp_k = k;
						temp_k--;
						
							while(temp_k >= 0){
								temp_dy = dyx[temp_k].y;
								temp_dx = dyx[temp_k].x;
								if(yy + temp_dy >= 0 && xx + temp_dx >= 0 && xx + temp_dx < width){
									if(dyx[k].p == dyx[temp_k].p){
										roff_up[y][x][k] = temp_dy * width + temp_dx;
										flag = 1;
										break;
									}
								}
								temp_k--;
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
					}
					roff_right[y][x][k] = temp_dy * width + temp_dx;
				}else{
					roff_right[y][x][k] = dy * width + dx;
				}
			}
		}
	}
}*/

/*int ref_offset4(int height, int width, int y, int x, int k, int co)
{
	int dy, dx, temp_dy, temp_dx, temp_k;
	
	if(y + mindy >= 0 && x + mindx >= 0 && x + maxdx < width){
		return(roff_const[k]);
	}else{
		if(y == 0 && x == 0){
			dx = 0;
			dy = height;
			return(dy * width + dx);
		}
		dy = dyx_y[k];
		dx = dyx_x[k];
		if(y + dy < 0 || x + dx < 0 || x + dx >= width){
			temp_k = k;
			temp_k--;
			switch(co){
				case 1:
				case 2:
					while(temp_k >= 0){
						temp_dy = dyx_y[temp_k];
						temp_dx = dyx_x[temp_k];
						if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
							if(dyx_p[k] == dyx_p[temp_k]){
								return(temp_dy * width + temp_dx);
							}
						}
						temp_k--;
					}
					break;
				default:
					while(temp_k >= 0){
						temp_dy = dyx_y[temp_k];
						temp_dx = dyx_x[temp_k];
						if(y + temp_dy >= 0 && x + temp_dx >= 0 && x + temp_dx < width){
							switch(dyx_p[k]){
								case 0:
								case 2:
									switch(dyx_p[temp_k]){
										case 0:
										case 2:
											return(temp_dy * width + temp_dx);
									}
									break;
								default:
									if(dyx_p[k] == dyx_p[temp_k]){
										return(temp_dy * width + temp_dx);
									}
									break;
							}
						}
						temp_k--;
					}
					break;
			}
			switch(x){
				case 0:
					dy = -1;
					dx = 0;
					break;
				default:
					dy = 0;
					dx = -1;
					break;
			}
		}
	}
	return(dy * width + dx);
}*/

/*inline int ref_offset5(int height, int width, int y, int x, int ar, int k)
{
	//int roff = 0;
	//if(k >=  MAX_PRD_ORDER) {printf("over\n"); exit(1);}
	//int ar=0;
	if(y + mindy >= 0){
		if(x + mindx >= 0){
			if(x + maxdx < width){
				return(roff_const[ar][k]);
				//roff = roff_const[k];
			}else{
				return(roff_right[y+mindy][x-width+maxdx][k]);
				//roff = roff_right[y+mindy][x-width+maxdx][k];
			}
		}else{
			return(roff_left[y+mindy][x][k]);
			//roff = roff_left[y+mindy][x][k];
		}
	}else{
		return(roff_up[y][x][k]);
		//roff = roff_up[y][x][k];
	}
	
		
}*/

/*inline int ref_offset5_eachcolor(int height, int width, int y, int x, int co, int ar, int k)
{
	//int roff = 0;
	//if(k >=  MAX_PRD_ORDER) {printf("over\n"); exit(1);}
	//int ar=0;
	if(y + mindy >= 0){
		if(x + mindx >= 0){
			if(x + maxdx < width){
				return(roff_const3[co][ar][k]);
				//roff = roff_const[k];
			}else{
				return(roff_right[y+mindy][x-width+maxdx][k]);
				//roff = roff_right[y+mindy][x-width+maxdx][k];
			}
		}else{
			return(roff_left[y+mindy][x][k]);
			//roff = roff_left[y+mindy][x][k];
		}
	}else{
		return(roff_up[y][x][k]);
		//roff = roff_up[y][x][k];
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

void enc_make_label(ENCODER *enc)//vseg
{		
	printf("make_label\n");
		FILE *fp1, *fp2;
		//int **label;
		//POINT *tl, *br;
		int x, y, height, width;
		x = 0;
		y = 0;
		int xg, yg, area;
		char *area_p;
		height = enc->height;
		width = enc->width;
		enc->max_geny = 0;
		enc->max_genx = 0;
		fp1 = fileopen("outfile_label.txt", "rb");
		//fp1 = fileopen("outfile_label.txt", "rb");
		fp2 = fileopen("outfile_range.txt", "rb");
		fscanf(fp1, "%d\n", &enc->num_of_vgen);
		fscanf(fp1, "%d\n", &enc->num_of_vgenc);
		printf("vgen,vgenc = %d,%d\n", enc->num_of_vgen, enc->num_of_vgenc);
		enc->num_of_geny = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
		enc->num_of_genx = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
		for(area = 0; area < enc->num_of_vgen; area++){
			fscanf(fp1, "%d\n",&enc->num_of_geny[area]);
			fscanf(fp1, "%d\n",&enc->num_of_genx[area]);
			printf("geny,genx = %d,%d\n",enc->num_of_geny[area], enc->num_of_genx[area]);
			if(enc->max_geny < enc->num_of_geny[area]) enc->max_geny = enc->num_of_geny[area];
			if(enc->max_genx < enc->num_of_genx[area]) enc->max_genx = enc->num_of_genx[area];
		}
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				fscanf(fp1, "%d\n", &enc->label[y][x].x);
				fscanf(fp1, "%d\n", &enc->label[y][x].y);
				fscanf(fp1, "%d\n", &enc->label[y][x].area);
			}
		}
		enc->gen_tl = (POINT ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(POINT));//vseg
		enc->gen_br = (POINT ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(POINT));//vseg
		for (area = 0; area < enc->num_of_vgen; area++) {
			for (yg = 0; yg < enc->num_of_geny[area]; yg++) { //vseg
				for (xg = 0; xg < enc->num_of_genx[area]; xg++) {
						enc->gen_tl[area][yg][xg].x = width-1;
						enc->gen_tl[area][yg][xg].y = height-1;
						enc->gen_tl[area][yg][xg].area = area;
						enc->gen_br[area][yg][xg].x = 0;
						enc->gen_br[area][yg][xg].y = 0;
						enc->gen_br[area][yg][xg].area = area;
				}
			}
		}
		for(area = 0; area < enc->num_of_vgen; area++){
			for(yg = 0; yg < enc->num_of_geny[area]; yg++){
				for(xg = 0; xg < enc->num_of_genx[area]; xg++){
					fscanf(fp2, "%d\n", &enc->gen_tl[area][yg][xg].x);
					fscanf(fp2, "%d\n", &enc->gen_tl[area][yg][xg].y);
					fscanf(fp2, "%d\n", &enc->gen_br[area][yg][xg].x);
					fscanf(fp2, "%d\n", &enc->gen_br[area][yg][xg].y);
					if(enc->gen_tl[area][yg][xg].x < 0){enc->gen_tl[area][yg][xg].x = 0;
					}else if(enc->gen_tl[area][yg][xg].x > width - 1){enc->gen_tl[area][yg][xg].x = width - 1;}
					if(enc->gen_tl[area][yg][xg].y < 0) {enc->gen_tl[area][yg][xg].y = 0;
					}else if(enc->gen_tl[area][yg][xg].y > height - 1){enc->gen_tl[area][yg][xg].y = height - 1;}
					if(enc->gen_br[area][yg][xg].x < 0) {enc->gen_br[area][yg][xg].x = 0;
					}else if(enc->gen_br[area][yg][xg].x > width - 1){enc->gen_br[area][yg][xg].x = width - 1;}
					if(enc->gen_br[area][yg][xg].y < 0) {enc->gen_br[area][yg][xg].y = 0;
					}else if(enc->gen_br[area][yg][xg].y > height - 1){enc->gen_br[area][yg][xg].y = height - 1;}
				}
			}
		}
		
		fclose(fp1);
		fclose(fp2);
	
    for(y = 0; y < height; y++){
		area_p = enc->area[y];
		for(x = 0; x < width; x++){
			area_p[x] = enc->label[y][x].area;
			//printf("%d\n", enc->area[y][x]);
		}
    }
		return;
}


void set_cost_model(ENCODER *enc, int f_mmse)//
{
    int gr, i, j, k;
    double a, b, c, var;
    int gauss_index = 9;        // for enc->pmodels
    PMODEL *pm;

    for (i = 0; i <= enc->maxval; i++) {
			for (j = 0; j <= (enc->maxval << 1); j++) { //
		    k = (i << 1) - j - 1;
		    if (k < 0) k = -(k + 1);
		    enc->econv[i][j] = k;
		}
    }
    enc->encval = enc->err;
    for (gr = 0; gr < enc->num_group; gr++) {//enc->num_group
			var = enc->sigma[gr] * enc->sigma[gr];
			if (f_mmse) {
		    a = 0;
		    b = 1.0;
			} else {
		    a = 0.5 * log(2 * M_PI * var) / log(2.0);
		    b = 1.0 / (2.0 * log(2.0) * var);
			}
			enc->pmlist[gr] = pm = enc->pmodels[gr][gauss_index];
			for (k = 0; k <= pm->size; k++) {
		    c = (double)k * 0.5 + 0.25;
		    pm->cost[k] = a + b * c * c;
			}
			pm->subcost[0] = 0.0;
    }
    for (k = 0; k <= enc->maxprd; k++) {
			enc->bconv[k] = 0;
			enc->fconv[k] = 0;
    }
    return;
}

ENCODER *init_encoder(IMAGE *img, int num_group,int prd_order,
                      int coef_precision, int f_huffman, int quadtree_depth,
                      int num_pmodel, int pm_accuracy, int num_kind_prd)
{
	printf("init_encoder\n");
    ENCODER *enc;
    int x, y, i, j, k, ar, co;
	int xg, yg, area;
    double c;
	int height, width, max_class, max_prd_order_all, *max_prd_order;

    enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
    height = enc->height = img->height;
    width = enc->width = img->width;
    enc->maxval = img->maxval;
    //enc->num_class = num_class;
    enc->num_group = num_group;
    enc->num_kind_prd = num_kind_prd;
	enc->num_center_class = 79;
	enc->label = (POINT **)alloc_2d_array(height, width, sizeof(POINT));//vseg
    enc->area = (char **)alloc_2d_array(height, width, sizeof(char));
	for(i = 0; i < height; i++){
      for(j = 0; j < width; j++){
        enc->label[i][j].area = 0;
		enc->label[i][j].y = 0;
		enc->label[i][j].x = 0;
        enc->area[i][j] = 0;
      }
    }
	enc_make_label(enc);//vseg
	max_prd_order = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
    //enc->base_prd_order = prd_order;
    max_prd_order_all = enc->max_prd_order_all = MAX_PRD_ORDER_ALL;
    //max_prd_order_all = enc->max_prd_order_all = 36;
	//////////////////////////////////////////////
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		max_prd_order[ar]=MAX_PRD_ORDER;
	}
	enc->max_num_pixel = 0;

	//////////////////////////////////////////////
    enc->coef_precision = coef_precision;
    enc->max_coef = (2 << coef_precision);
    enc->f_huffman = f_huffman;
    enc->num_pmodel = num_pmodel;
    enc->pm_accuracy = pm_accuracy;
    enc->maxprd = enc->maxval << enc->coef_precision;
	enc->num_class = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	enc->max_prd_order = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	enc->base_prd_order = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	for(ar = 0; ar < enc->num_of_vgen; ar++) {
		if(ar < enc->num_of_vgenc){
			enc->num_class[ar] = NUM_CLASS_CENTER;
		}else{
			enc->num_class[ar] = NUM_CLASS_EDGE;
		}
		enc->max_prd_order[ar] = max_prd_order[ar];
		enc->base_prd_order[ar] = BASE_PRD_ORDER;
	}
	max_class = enc->max_class = MAX_CLASS;
    enc->predictor = (int ****)alloc_4d_array(max_class, num_kind_prd, enc->num_of_vgen, max_prd_order_all, sizeof(int));
	//init_4d_array(enc->predictor, max_class, num_kind_prd, enc->num_of_vgen, max_prd_order_all, 0);
	enc->predict_out = (int **)alloc_2d_array(num_kind_prd, max_prd_order_all, sizeof(int));
	init_2d_array(enc->predict_out, num_kind_prd, max_prd_order_all, 0);
    enc->num_nzcoef = (int ***)alloc_3d_array(max_class, num_kind_prd, enc->num_of_vgen, sizeof(int));
	init_3d_array(enc->num_nzcoef, max_class, num_kind_prd, enc->num_of_vgen, max_prd_order_all);
	
	
	enc->num_pixel = (int ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(int));//vseg
	enc->gen_class = (int ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(int));//vseg
	enc->num_segment = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	init_array(enc->num_segment, enc->num_of_vgen, 0);
	for (area = 0; area < enc->num_of_vgen; area++) {
		for (yg = 0; yg < enc->max_geny; yg++) { //vseg
			for (xg = 0; xg < enc->max_genx; xg++) {
					enc->num_pixel[area][yg][xg] = 0;
					enc->gen_class[area][yg][xg] = 0;
			}
		}
	}
	
    enc->nzconv = (int ****)alloc_4d_array(max_class, num_kind_prd, enc->num_of_vgen, max_prd_order_all, sizeof(int));
	enc->num_search = (int ***)alloc_3d_array(max_class, num_kind_prd, enc->num_of_vgen, sizeof(int));
	init_3d_array(enc->num_search, max_class, num_kind_prd, enc->num_of_vgen, 30);
    enc->th = (int ****)alloc_4d_array(max_class, num_kind_prd, enc->num_of_vgen ,num_group, sizeof(int));
	init_4d_array(enc->th, max_class, num_kind_prd, enc->num_of_vgen, num_group, 0);
    for (i = 0; i < max_class; i++) {
		for (j = 0; j < num_kind_prd; j++) {
		    for(ar = 0; ar < enc->num_of_vgen; ar++){
				for (k = 0; k < max_prd_order_all; k++) {
					enc->nzconv[i][j][ar][k] = k;
					enc->predictor[i][j][ar][k] = 0;
				}
				enc->th[i][j][ar][enc->num_group - 1] = MAX_UPARA + 1;
			}
		}
    }
	enc->ord2mhd = (int *)alloc_mem(max_prd_order_all * sizeof(int));
    for (i = 0; i < max_prd_order_all; i++) {
        enc->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
    }
    enc->upara = (int **)alloc_2d_array(height, width, sizeof(int));
	init_2d_array(enc->upara, height, width, 0);
    enc->prd = (int **)alloc_2d_array(height, width, sizeof(int));
	init_2d_array(enc->prd, height, width, 0);
    enc->org = (int **)alloc_2d_array(height+1, width, sizeof(int));
	init_2d_array(enc->org, height+1, width, 0);
    enc->err = (int **)alloc_2d_array(height+1, width, sizeof(int));
	init_2d_array(enc->err, height+1, width, 0);
    enc->ctx_weight = init_ctx_weight();
    enc->class = (int **)alloc_2d_array(height, width, sizeof(int));
	init_2d_array(enc->class, height, width, 0);
    enc->group = (char **)alloc_2d_array(height, width,sizeof(char));
    for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			enc->group[y][x] = 6;
			enc->org[y][x] = (int)img->val[y][x];
		}
    }
    enc->org[height][0] = (enc->maxval + 1) >> 1;
    enc->err[height][0] = (enc->maxval + 1) >> 2;
    enc->uquant = (int ****)alloc_4d_array(max_class, num_kind_prd, enc->num_of_vgen, MAX_UPARA + 1, sizeof(int));
	init_4d_array(enc->uquant, max_class, num_kind_prd, enc->num_of_vgen, MAX_UPARA+1, 6/*enc->num_group - 1*/);
    enc->econv = (int **)alloc_2d_array(enc->maxval+1, (enc->maxval<<1)+1, sizeof(int));
	init_2d_array(enc->econv, enc->maxval+1, (enc->maxval<<1)+1, 0);
    enc->bconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
    enc->fconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
    enc->pmlist = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));
    enc->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
    enc->spm.cumfreq = &(enc->spm.freq[MAX_SYMBOL]);
    if (enc->f_huffman == 1) {
			enc->sigma = sigma_h;
    } else {
			enc->sigma = sigma_a;
    }
    enc->mtfbuf = (int *)alloc_mem(max_class * sizeof(int));
	init_array(enc->mtfbuf, max_class, 0);
	enc->zero_m = (int *)alloc_mem(max_prd_order_all * sizeof(int));
	init_array(enc->zero_m, max_prd_order_all, NUM_ZMODEL >> 1);
	enc->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for(i = 0; i < NUM_ZMODEL; i++){
		enc->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
    enc->coef_m = (int *)alloc_mem(max_prd_order_all * sizeof(int));
	init_array(enc->coef_m, max_prd_order_all, 0);
    enc->coef_cost = (cost_t **)alloc_2d_array(16, enc->max_coef + 1, sizeof(cost_t));
    enc->coef_cost2 = (cost_t ***)alloc_3d_array(NUM_ZMODEL, 16, enc->max_coef + 1, sizeof(cost_t));
	
	for (i = 0; i < 16; i++) {
#ifdef OPT_SIDEINFO
			if (enc->f_huffman == 1) {
			    for (j = 0; j <= enc->max_coef; j++) {
						enc->coef_cost[i][j] = ((j >> i) + i + 1);
						if (j > 0) enc->coef_cost[i][j] += 1.0;
			    }
			} else {
			    double p;
			    set_spmodel(&enc->spm, enc->max_coef + 1, i);
			    p = log(enc->spm.cumfreq[enc->max_coef + 1]);
			    for (j = 0; j <= enc->max_coef; j++) {
						enc->coef_cost[i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
						if (j > 0) enc->coef_cost[i][j] += 1.0;
			    }
			}
#else
			for (j = 0; j <= enc->max_coef; j++) {
			    enc->coef_cost[i][j] = 0;
			}
#endif
    }
	double p, zero, nonz;
	uint cumb = 0;
	for (k = 0; k < NUM_ZMODEL; k++) {
		for (i = 0; i < 16; i++) {
#ifdef OPT_SIDEINFO	
			nonz = log((double)TOT_ZEROFR / (double)enc->zero_fr[k]) / log(2.0);
			zero = log((double)TOT_ZEROFR / ((double)TOT_ZEROFR - (double)enc->zero_fr[k])) / log(2.0);
			set_spmodel(&enc->spm, enc->max_coef + 1, i);
			cumb = enc->spm.freq[0];
			p = log(enc->spm.cumfreq[enc->max_coef + 1] - cumb);
			for (j = 1; j <= enc->max_coef; j++) {
				enc->coef_cost2[k][i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
				enc->coef_cost2[k][i][j] += 1.0 + nonz;
			}
			enc->coef_cost2[k][i][0] = zero;
#else
			for (j = 0; j <= enc->max_coef; j++) {
				enc->coef_cost2[k][i][j] = 0;
			}
#endif
		}
	}
    enc->th_cost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
    for (i = 0; i < MAX_UPARA + 2; i++) {
			enc->th_cost[i] = 0;
    }
    enc->class_cost = (cost_t **)alloc_2d_array(enc->num_of_vgen, max_class, sizeof(cost_t));//vseg
    c = log((double)max_class) / log(2.0);
		for(j = 0; j < enc->num_of_vgen; j++){
    	for (i = 0; i < enc->num_class[j]; i++) {

				enc->class_cost[j][i] = c;//vseg
			}
    }
    enc->err_cost = NULL;
    enc->cl_hist = (int *)alloc_mem(max_class * sizeof(int));
	
	
	
	//*********符号量を計算
	enc->area_co_err = (int **)alloc_2d_array(enc->num_of_vgen, NUM_KIND_PRD, sizeof(int));
	enc->area_pel = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	enc->area_class_info = (double *)alloc_mem(enc->num_of_vgen * sizeof(double));
	enc->area_predictor_info = (double *)alloc_mem(enc->num_of_vgen * sizeof(double));
	enc->area_th_info = (double*)alloc_mem(enc->num_of_vgen * sizeof(double));
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		for(co = 0; co < NUM_KIND_PRD;co++){
			enc->area_co_err[ar][co] = 0;
		}
		enc->area_class_info[ar] = 0;
		enc->area_predictor_info[ar] = 0;
		enc->area_th_info[ar] = 0;
	}

	init_array(enc->area_pel, enc->num_of_vgen, 0);
	
    return (enc);
}


//vseg
void init_class_vseg(ENCODER *enc)
{//vseg
	printf("init_class->");
	int x, y, i, j, v, k, cl;
	int xg;
	int yg;
	int area;
	int tmp_n;
	//POINT tmp_p;
	double *var;
	double tmp_d, sum;
	int num_pix;
	int *num_block;
	int tly, tlx, bry, brx;
	num_block = (int *)alloc_mem((enc->max_geny) * (enc->max_genx) * sizeof(int));
	//var = (double *)alloc_mem((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) * sizeof(double));
	var = (double *)alloc_mem((enc->max_geny) * (enc->max_genx) * sizeof(double));
	//ptr_d = (double **)alloc_mem((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) * sizeof(double *));

		/*area = 0;
		i = 0;
		enc->num_segment[area] = 0;//マイクロレンズ内の相対的同位置の領域数
		for (yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++) {
			for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
				k = yg * (NUM_OF_GENx +RUN_OVER) + xg;
				if(k >= ((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER))) printf("ERROR\n");
				var[k] = 0;
				sum = 0;
				num_pix = 0;
				enc->num_pixel[area][yg][xg] = 0;
				tly = enc->gen_tl[area][yg][xg].y;
				tlx = enc->gen_tl[area][yg][xg].x;
				bry = enc->gen_br[area][yg][xg].y + 1;
				brx = enc->gen_br[area][yg][xg].x + 1;
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
							if (color(y, x) == 1 || color(y, x) == 2) {
										v = enc->org[y][x];
										sum += v;
										var[k] += v * v;
										num_pix++;//領域内の緑の画素数
							}
							enc->num_pixel[area][yg][xg]++;//領域内の画素数
						}
					}
				}
				if(enc->num_pixel[area][yg][xg] == 0){
					var[k] = -1;
					//ptr_d[k] = &(var[k]);
					num_block[k] = -1;
					enc->gen_class[area][yg][xg] = -1;
				}else{
					enc->num_segment[area]++;
					if(num_pix != 0){
						var[k] = var[k] / (double)num_pix;//var[k]に輝度値の2乗平均を記録しなおす。
						var[k] -= sum * sum / (double)(num_pix * num_pix);//var[k]からサンプルの平均を引くのでvar[k]には分散を格納	
					}else{
						var[k] = 0;
					}
					//ptr_d[k] = &(var[k]);
					num_block[k] = i;//形状が同じ領域の数を数える（領域内の画素数が0のものを除く）
					i++;
				}
			}
		}

		for (i = (NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) - 1; i > 0; i--) {//分散の大きさに従ってnum_blockを並び替える
			for (j = 0; j < i; j++) {
				//while(*ptr_d[j] < 0 && j < i - 1) j++;
				if(j >= ((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER))) printf("ERROR2\n");
				while(var[j] < 0 && j < i - 1) j++;//分散が1だったら飛ばす
				k = j + 1;
				//while(*ptr_d[k] < 0 && k < i ) k++;
				if(k >= ((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER))) printf("ERROR3\n");
				while(var[k] < 0 && k < i ) k++;
				//if(*ptr_d[k] < 0 || k > i + 1){
				if(var[k] < 0 || k > i + 1){
					continue;
				//}else if(*ptr_d[j] > *ptr_d[k]) {
				}else if(var[j] > var[k]) {
					
					tmp_d = var[j];
					var[j] = var[k];
					var[k] = tmp_d;
					tmp_n = num_block[j];
					num_block[j] = num_block[k];
					num_block[k] = tmp_n;
				}
			}
		}

		for (yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
			for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
				if(enc->num_pixel[area][yg][xg] == 0){
					continue;
				}else{
					k = yg * (NUM_OF_GENx +RUN_OVER) + xg;
					cl = ((num_block[k]) * enc->num_class[area]) / enc->num_segment[area];
					//if(cl >= 80) printf("cl = %d\n", cl);
					tly = enc->gen_tl[area][yg][xg].y;
					tlx = enc->gen_tl[area][yg][xg].x;
					bry = enc->gen_br[area][yg][xg].y + 1;
					brx = enc->gen_br[area][yg][xg].x + 1;
					for (y = tly; y < bry; y++) {
						for (x = tlx; x < brx; x++) {
							if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
								enc->class[y][x] = cl;
							}
						}
					}
					enc->gen_class[area][yg][xg] = cl;
				}
			}
		}
		*/
		for(area = 0; area < enc->num_of_vgen; area++){
			i = 0;
			enc->num_segment[area] = 0;

			for (yg = 0; yg < enc->num_of_geny[area]; yg++) {
				for (xg = 0; xg < enc->num_of_genx[area]; xg++) {
					tly = enc->gen_tl[area][yg][xg].y;
					tlx = enc->gen_tl[area][yg][xg].x;
					bry = enc->gen_br[area][yg][xg].y + 1;
					brx = enc->gen_br[area][yg][xg].x + 1;
					k = yg * enc->num_of_genx[area] + xg;
					var[k] = sum = 0;
					num_pix = 0;
					enc->num_pixel[area][yg][xg] = 0;
					for (y = tly; y < bry; y++) {
							for (x = tlx; x < brx; x++) {
								if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){//label
									if (color(y, x) == 1 || color(y, x) == 2) {
												v = enc->org[y][x];
												sum += v;
												var[k] += v * v;
												num_pix++;
									}
									enc->num_pixel[area][yg][xg]++;
								}
							}
					}
					if(enc->num_pixel[area][yg][xg] == 0){
						var[k] = -1;
						//ptr_d[k] = &(var[k]);
						num_block[k] = -1;
						enc->gen_class[area][yg][xg] = -1;
					}else{
						enc->num_segment[area]++;
						if(num_pix != 0){
							var[k] = var[k] / (double)num_pix;
							var[k] -= sum * sum / (double)(num_pix * num_pix);
						}else{
							var[k] = 0;
						}
						//ptr_d[k] = &(var[k]);
						num_block[k] = i;
						i++;
					}
				}
			}
			for (i = enc->num_of_geny[area] * enc->num_of_genx[area] - 1; i > 0; i--) {
				for (j = 0; j < i; j++) {
					//while(*ptr_d[j] < 0 && j < i - 1) j++;
					while(var[j] < 0 && j < i - 1) j++;
					k = j + 1;
					//while(*ptr_d[k] < 0 && k < i ) k++;
					while(var[k] < 0 && k < i ) k++;
					//if(*ptr_d[k] < 0 || k > i + 1){
					if(var[k] < 0 || k > i + 1){
						continue;
					//}else if(*ptr_d[j] > *ptr_d[k]) {
					}else if(var[j] > var[k]) {
						
						tmp_d = var[j]; //var
						var[j] = var[k];
						var[k] = tmp_d;
						tmp_n = num_block[j];
						num_block[j] = num_block[k];
						num_block[k] = tmp_n;
					}
				}
			}

			for (yg = 0; yg < enc->num_of_geny[area]; yg++) {
				for (xg = 0; xg < enc->num_of_genx[area]; xg++) {
					if(enc->num_pixel[area][yg][xg] == 0){
						continue;
					}else{
						k = yg * enc->num_of_genx[area] + xg;
						cl = ((num_block[k]) * enc->num_class[area]) / enc->num_segment[area];
						tly = enc->gen_tl[area][yg][xg].y;
						tlx = enc->gen_tl[area][yg][xg].x;
						bry = enc->gen_br[area][yg][xg].y + 1;
						brx = enc->gen_br[area][yg][xg].x + 1;
						for (y = tly; y < bry; y++) {
								for (x = tlx; x < brx; x++) {
									if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
										enc->class[y][x] = cl;
									}
								}
						}
						enc->gen_class[area][yg][xg] = cl;
					}
				}
			}
		}
		/*area = 0;
		for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg ++){
			for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg ++){
				if(enc->max_num_pixel < enc->num_pixel[area][yg][xg]) enc->max_num_pixel = enc->num_pixel[area][yg][xg];
			}
		}*/
		for(area = 0; area < enc->num_of_vgen; area++){
			for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
				for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
					if(enc->max_num_pixel < enc->num_pixel[area][yg][xg]) enc->max_num_pixel = enc->num_pixel[area][yg][xg];
				}
			}
		}

	free(var);
	//free(ptr_d);
	free(num_block);
	printf("ok\n");
}

/*void init_class_vseg(ENCODER *enc)
{//vseg
	printf("init_class\n");
	int x, y, i, j, v, k, cl;
	int xg;
	int yg;
	int area;
	int *tmp_n;
	//POINT tmp_p;
	double *var, **ptr_d, *tmp_d, sum;
	int num_pix;
	int **num_block;
	int tly, tlx, bry, brx;
	num_block = (int **)alloc_2d_array((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER), 1, sizeof(int));
	var = (double *)alloc_mem((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) * sizeof(double));
	ptr_d = (double **)alloc_mem((NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) * sizeof(double *));

		area = 0;
		i = 0;
		enc->num_segment[area] = 0;
		for (yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++) {
			for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
				k = yg * (NUM_OF_GENx +RUN_OVER) + xg;
				var[k] = sum = 0;
				num_pix = 0;
				enc->num_pixel[area][yg][xg] = 0;
				tly = enc->gen_tl[area][yg][xg].y;
				tlx = enc->gen_tl[area][yg][xg].x;
				bry = enc->gen_br[area][yg][xg].y + 1;
				brx = enc->gen_br[area][yg][xg].x + 1;
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
							if (color(y, x) == 1 || color(y, x) == 2) {
										v = enc->org[y][x];
										sum += v;
										var[k] += v * v;
										num_pix++;
							}
							enc->num_pixel[area][yg][xg]++;
						}
					}
				}
				if(enc->num_pixel[area][yg][xg] == 0){
					var[k] = -1;
					ptr_d[k] = &(var[k]);
					*num_block[k] = -1;
					enc->gen_class[area][yg][xg] = -1;
				}else{
					enc->num_segment[area]++;
					var[k] = var[k] / (double)num_pix;
					var[k] -= sum * sum / (double)(num_pix * num_pix);
					ptr_d[k] = &(var[k]);
					*num_block[k] = i;
					i++;
				}
			}
		}

		for (i = (NUM_OF_GENy +RUN_OVER) * (NUM_OF_GENx +RUN_OVER) - 1; i > 0; i--) {
			for (j = 0; j < i; j++) {
				while(*ptr_d[j] < 0 && j < i - 1) j++;
				k = j + 1;
				while(*ptr_d[k] < 0 && k < i ) k++;
				if(*ptr_d[k] < 0 || k > i + 1){
					continue;
				}else if(*ptr_d[j] > *ptr_d[k]) {
					tmp_d = ptr_d[j];
					ptr_d[j] = ptr_d[k];
					ptr_d[k] = tmp_d;
					tmp_n = num_block[j];
					num_block[j] = num_block[k];
					num_block[k] = tmp_n;
				}
			}
		}

		for (yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
			for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
				if(enc->num_pixel[area][yg][xg] == 0){
					continue;
				}else{
					k = yg * (NUM_OF_GENx +RUN_OVER) + xg;
					cl = ((*num_block[k]) * enc->num_class[area]) / enc->num_segment[area];
					tly = enc->gen_tl[area][yg][xg].y;
					tlx = enc->gen_tl[area][yg][xg].x;
					bry = enc->gen_br[area][yg][xg].y + 1;
					brx = enc->gen_br[area][yg][xg].x + 1;
					for (y = tly; y < bry; y++) {
						for (x = tlx; x < brx; x++) {
							if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
								enc->class[y][x] = cl;
							}
						}
					}
					enc->gen_class[area][yg][xg] = cl;
				}
			}
		}

		for(area = 1; area < enc->num_of_vgen; area++){
			i = 0;
			enc->num_segment[area] = 0;

			for (yg = 0; yg < NUM_OF_GEN_ey; yg++) {
				for (xg = 0; xg < NUM_OF_GEN_ex; xg++) {
					tly = enc->gen_tl[area][yg][xg].y;
					tlx = enc->gen_tl[area][yg][xg].x;
					bry = enc->gen_br[area][yg][xg].y + 1;
					brx = enc->gen_br[area][yg][xg].x + 1;
					k = yg * NUM_OF_GEN_ex + xg;
					var[k] = sum = 0;
					num_pix = 0;
					enc->num_pixel[area][yg][xg] = 0;
					for (y = tly; y < bry; y++) {
							for (x = tlx; x < brx; x++) {
								if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){//label
									if (color(y, x) == 1 || color(y, x) == 2) {
												v = enc->org[y][x];
												sum += v;
												var[k] += v * v;
												num_pix++;
									}
									enc->num_pixel[area][yg][xg]++;
								}
							}
					}
					if(enc->num_pixel[area][yg][xg] == 0){
						var[k] = -1;
						ptr_d[k] = &(var[k]);
						*num_block[k] = -1;
						enc->gen_class[area][yg][xg] = -1;
					}else{
						enc->num_segment[area]++;
						var[k] = var[k] / (double)num_pix;
						var[k] -= sum * sum / (double)(num_pix * num_pix);
						ptr_d[k] = &(var[k]);
						*num_block[k] = i;
						i++;
					}
				}
			}
			for (i = NUM_OF_GEN_ey * NUM_OF_GEN_ex - 1; i > 0; i--) {
				for (j = 0; j < i; j++) {
					while(*ptr_d[j] < 0 && j < i - 1) j++;
					k = j + 1;
					while(*ptr_d[k] < 0 && k < i ) k++;
					if(*ptr_d[k] < 0 || k > i + 1){
						continue;
					}else if(*ptr_d[j] > *ptr_d[k]) {
						tmp_d = ptr_d[j]; //var
						ptr_d[j] = ptr_d[k];
						ptr_d[k] = tmp_d;
						tmp_n = num_block[j];
						num_block[j] = num_block[k];
						num_block[k] = tmp_n;
					}
				}
			}

			for (yg = 0; yg < NUM_OF_GEN_ey; yg++) {
				for (xg = 0; xg < NUM_OF_GEN_ex; xg++) {
					if(enc->num_pixel[area][yg][xg] == 0){
						continue;
					}else{
						k = yg * NUM_OF_GEN_ex + xg;
						cl = ((*num_block[k]) * enc->num_class[area]) / enc->num_segment[area];
						tly = enc->gen_tl[area][yg][xg].y;
						tlx = enc->gen_tl[area][yg][xg].x;
						bry = enc->gen_br[area][yg][xg].y + 1;
						brx = enc->gen_br[area][yg][xg].x + 1;
						for (y = tly; y < bry; y++) {
								for (x = tlx; x < brx; x++) {
									if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
										enc->class[y][x] = cl;
									}
								}
						}
						enc->gen_class[area][yg][xg] = cl;
					}
				}
			}
		}


	free(var);
	free(ptr_d);
	free(num_block);
}*/



/*double fconv(int k, int maxprd, int shift, int mask){
	int i;
	i = (maxprd - k + (1 << shift) / 2) >> shift;
	return(i & mask);
}*/

void set_cost_rate(ENCODER *enc)
{
    int gr, k, i, j, mask, shift, num_spm;
    double a, c;
    PMODEL *pm;
	FILE *fp;
	fp = fileopen("./LOG/econv1.csv", "w");
    if (enc->pm_accuracy < 0) {
		printf("***pm_accuracy < 0***\n");
			for (i = 0; i <= enc->maxval; i++) {
			    for (j = 0; j <= (enc->maxval << 1); j++) {
						k = (j + 1) >> 1;
						enc->econv[i][j] = e2E(i - k, k, j&1, enc->maxval);
						fprintf(fp, "econv1[%d][%d] = %d\n", i, j, enc->econv[i][j]);
			    }
			}
    }
	fclose(fp);
    if (enc->pm_accuracy < 0) {
			num_spm = 1;
    } else {
			enc->encval = enc->org;
			mask = (1 << enc->pm_accuracy) - 1;
			shift = enc->coef_precision - enc->pm_accuracy;
			for (k = 0; k <= enc->maxprd; k++) {
			    i = (enc->maxprd - k + (1 << shift) / 2) >> shift;
			    enc->fconv[k] = (i & mask);
			    enc->bconv[k] = (i >> enc->pm_accuracy);
					//printf("enc->fconv[%d] = %d, enc->bconv[%d] = %d\n", k, enc->fconv[k], k, enc->bconv[k]);
			}
			num_spm = 1 << enc->pm_accuracy;
    }
    a = 1.0 / log(2.0);
    for (gr = 0; gr < enc->num_group; gr++) {
			for (i = 0; i < enc->num_pmodel; i++) {
			    pm = enc->pmodels[gr][i];
			    if (enc->f_huffman == 1) {
						for (k = 0; k < pm->size; k++) {
				        pm->cost[k] = enc->vlcs[gr][pm->id].len[k];
						}
						pm->subcost[0] = 0.0;
			    } else if (enc->pm_accuracy < 0) {
						for (k = 0; k < pm->size; k++) {
						    pm->cost[k] = -a * log(pm->freq[k]);
						}
				    c = pm->cumfreq[enc->maxval + 1];
						pm->subcost[0] = a * log(c);
			    } else {
						for (j = 0; j < num_spm; j++) {
						    for (k = 0; k < pm->size; k++) {
									pm->cost[k] = -a * log(pm->freq[k]);
						    }
						    for (k = 0; k <= enc->maxval; k++) {
									c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];
									pm->subcost[k] = a * log(c);
						    }
						    pm++;
						}
			    }
			}
    }
}

/*void predict_region2(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
    int x, y, k, l, cl, co, ar, prd, org;
    int *coef_p, *nzc_p;
    int *prd_p;
	int *roff_p, **roff_pp;
	int *err_p, *org_p;
    int *class_p;
	char *area_p;
    
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		area_p = &enc->area[y][tlx];
		org_p = &enc->org[y][tlx];
		roff_pp = &enc->roff[y][tlx];
		err_p = &enc->err[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			co = color(y, x);
			ar = *area_p;
			area_p = area_p + 1;
			roff_p = *roff_pp++;
			coef_p = enc->predictor[cl][co][ar];
			nzc_p = enc->nzconv[cl][co][ar];
			prd = 0;
			for (k = 0; k < enc->num_nzcoef[cl][co][ar]; k++) {
				l = nzc_p[k];
				prd += org_p[roff_p[l]] * (coef_p[l]);//prd += org_p[*roff_p++] * (*coef_p++);
			}
			org = *org_p++;
			*prd_p++ = prd;

			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;
			prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
		}
	}
}*/

void predict_region(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
    int x, y, k, l, cl, co, ar, prd, org, nzcoef/*, count = 0*/;
	int roff = 0;
    int *coef_p, *nzc_p;
    int *prd_p;
	int *err_p, *org_p;
    int *class_p;
	char *area_p;
    int width = enc->width;
	int height = enc->height;
	int maxprd, coef_precision;
	maxprd = enc->maxprd;
	coef_precision = enc->coef_precision;
	prd = org = nzcoef = 0;
	for (y = tly; y < bry; y++) {
		class_p = enc->class[y];
		area_p = enc->area[y];
		org_p = enc->org[y];
		err_p = enc->err[y];
		prd_p = enc->prd[y];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			co = color(y, x);
			ar = (int)*area_p++;
			coef_p = enc->predictor[cl][co][ar];
			//if(*coef_p == 128 || *coef_p == -128) count++;
			nzc_p = enc->nzconv[cl][co][ar];
			nzcoef = enc->num_nzcoef[cl][co][ar];
			prd = 0;
			for (k = 0; k < nzcoef; k++) {
				l = nzc_p[k];
				//roff = ref_offset4(height, width, y, x, l, co);
				roff = ref_offset5_2(height, width, y, x, l);
				prd += org_p[roff] * (coef_p[l]);
			}
			org = *org_p++;
			*prd_p++ = prd;
			if (prd < 0) prd = 0;
			else if (prd > maxprd) prd = maxprd;
			prd >>= (coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
		}
	}
	//printf("coef max count = %d\n", count);
}

inline int calc_uenc(ENCODER *enc, int y, int x, int co)
{
    int u, k, *err_p, *wt_p, height, width;
	int roff = 0;
    err_p = &enc->err[y][x];
	height = enc->height;
	width = enc->width;
    wt_p = enc->ctx_weight;
	//ar = enc->area[y][x];
    u = 0;
    for (k =0; k < NUM_UPELS; k++){
		//roff = ref_offset4(height, width, y, x, k, co);
		//roff = ref_offset5_eachcolor(height, width, y, x, co, ar, k);
		roff = ref_offset5_2(height, width, y, x, k);
		u += err_p[roff] * (*wt_p++);
    }
    u >>= 6;
    if (u > MAX_UPARA) u = MAX_UPARA;
    return (u);
}

cost_t calc_cost(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
    cost_t cost;
    int x, y, u, cl, co, ar, gr, prd, e, base, frac, maxprd;
    int *upara_p, *prd_p, *encval_p;
    int *class_p;
	char *group_p, *area_p;
	//img_t *bconv, *fconv;
    PMODEL *pm;
	maxprd = enc->maxprd;
	base = frac = prd = e = 0;
	/*bconv = (img_t *)alloc_mem((maxprd+1) * sizeof(img_t));
	fconv = (img_t *)alloc_mem((maxprd+1) * sizeof(img_t));
	for(prd = 0; prd <= maxprd; prd++){
		bconv[prd] = enc->bconv[prd];
		fconv[prd] = enc->fconv[prd];
	}*/
		if (tly < 0) tly = 0;
		if (bry > enc->height) bry = enc->height;
		if (tlx < 0) tlx = 0;
		if (brx > enc->width) brx = enc->width;
		cost = 0.0;
		for (y = tly; y < bry; y++) {
			class_p = enc->class[y];
			area_p = enc->area[y];
			group_p = enc->group[y];
			upara_p = enc->upara[y];
			encval_p = enc->encval[y];
			prd_p = enc->prd[y];
			for (x = tlx; x < brx; x++) {
				cl = *class_p++;
				co = color(y, x);
				ar = *area_p++;
				*upara_p++ = u = calc_uenc(enc, y, x, co);
				*group_p++ = gr = enc->uquant[cl][co][ar][u];
				e = *encval_p++;
				prd = *prd_p++;
				if (prd < 0) prd = 0;
				else if (prd > maxprd) prd = maxprd;
				base = enc->bconv[prd];
				frac = enc->fconv[prd];
				pm = enc->pmlist[gr] + frac;
				cost += pm->cost[base + e] + pm->subcost[base];
			}
		}
		/*free(bconv);
		free(fconv);*/
    return (cost);
}

cost_t calc_cost_vseg(ENCODER *enc, int tly, int tlx, int bry, int brx, int xg, int yg, int area, int cl)
{
	cost_t cost;
	int x, y, u, co, gr, prd, e, base, frac;
	PMODEL *pm;
	cost = 0;
	for (y = tly; y < bry; y++) {
			for (x = tlx; x < brx; x++) {
				if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
					co = color(y, x);
					//ar = enc->area[y][x];
					enc->upara[y][x] = u = calc_uenc(enc, y, x, co);//uparameter?
					enc->group[y][x] = gr = enc->uquant[cl][co][area][u];
					e = enc->encval[y][x];
					prd = enc->prd[y][x];
					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;
					base = enc->bconv[prd];
					frac = enc->fconv[prd];
					pm = enc->pmlist[gr] + frac;
					cost += pm->cost[base + e] + pm->subcost[base];
				}
			}
	}
	return (cost);

}

cost_t design_predictor(ENCODER *enc, int f_mmse)
{
    double w, e, d, pivot, wg, l;
	double *mat_p;
    int x, y, i, j, k, cl, co, gr, ar, pivpos, height, width, base_prd_order[enc->num_of_vgen], max_prd_order_all;
	int num_kind_prd, num_class[enc->num_of_vgen];
	char *area_p;
    int *org_p, *coef_p, *class_p, *color_p;
	int coef_precision, max_coef;
	int **org_buf, **color_buf;
	height = enc->height;
	width = enc->width;
	max_prd_order_all = enc->max_prd_order_all;
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		base_prd_order[ar] = enc->base_prd_order[ar];
	}
	double mat[max_prd_order_all][max_prd_order_all+1];
	double weight[enc->num_group];
	int index[max_prd_order_all];
	int org_roff[max_prd_order_all];
	coef_precision = enc->coef_precision;
	max_coef = enc->max_coef;
	num_kind_prd = enc->num_kind_prd;
	pivpos = w = e = d = pivot = wg = l = 0;
	for(i = 0; i < enc->num_of_vgen; i++){
		num_class[i] = enc->num_class[i];
	}
	//printf("cpu_time = %f\n", cpu_time());
	//printf("f_mmse = %d\n", f_mmse);
	org_buf = (int **)alloc_2d_array(height+1, width, sizeof(int));
	color_buf = (int **)alloc_2d_array(height, width, sizeof(int));
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			color_buf[y][x] = color(y, x);
			//if(enc->class[y][x] >= 80 || enc->class[y][x] < 0) printf("ERROR, y = %d, x = %d, class = %d\n", y, x, enc->class[y][x]);
		}
	}
	for(y = 0; y < height+1; y++){
		for(x = 0; x < width; x++){
			//if(enc->org[y][x] > 255 || 0 > enc->org[y][x]) printf("ERROR\n");
			org_buf[y][x] = enc->org[y][x];
		}
	}
		for (gr = 0; gr < enc->num_group; gr++) {
			if (f_mmse) {
				weight[gr] = 1.0;
			} else {
				weight[gr] = 1.0 / (double)((double)enc->sigma[gr] * (double)enc->sigma[gr]);//
			}
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for (cl = 0; cl < num_class[ar]; cl++) {
				for (co = 0; co < num_kind_prd; co++) {
					for (i = 0; i < max_prd_order_all; i++) {
						index[i] = 0;
						mat_p = mat[i];
						for (j = 0; j <= max_prd_order_all; j++) {
							mat_p[j] = 0.0;
						}
					}

						for (y = 0; y < height; y++) {
							class_p = enc->class[y];
							area_p = enc->area[y];
							color_p = color_buf[y];
							for (x = 0; x < width; x++) {
								if(class_p[x] == cl && (int)area_p[x] == ar && color_p[x] == co) {
									org_p = &org_buf[y][x];
									gr = (int)enc->group[y][x];
									wg = weight[gr];
									for(i = 0; i < base_prd_order[ar]; i++){
										org_roff[i] = org_p[ref_offset5_2(height, width, y, x, i)];
									}
									for (i = 0; i < base_prd_order[ar]; i++) {
										//roff = ref_offset5_eachcolor(height, width, y, x, i);
										w = wg * org_roff[i];
										mat_p = mat[i];
										for (j = i; j < base_prd_order[ar]; j++) {
											//roff = ref_offset5_eachcolor(height, width, y, x, j);
											mat_p[j] += w * (double)org_roff[j];
										}
										mat_p[base_prd_order[ar]] += w * (double)org_p[0];
									}
								}
								//class_p++;
								//area_p++;
								//color_p++;
								//org_p++;
							}
						}
					

					
					for (i = 0; i < base_prd_order[ar]; i++) {
						index[i] = i;
						for (j = 0; j < i; j++) {
							mat[i][j] = mat[j][i];
							//if(mat[i][j] < 0) printf("mat[%d][%d] < 0\n", i, j);
						}
					}
					
					for (i = 0; i < base_prd_order[ar]; i++) {
						pivpos = i;
						pivot = fabs(mat[index[i]][i]);
						for (k = i + 1; k < base_prd_order[ar]; k++) {
							if (fabs(mat[index[k]][i]) > pivot) {
								pivot = fabs(mat[index[k]][i]);
								pivpos = k;
							}
						}

						k = index[i];
						index[i] = index[pivpos];
						index[pivpos] = k;
						if (pivot > 1E-10) {
							d = mat[index[i]][i];
							for (j = i; j <= base_prd_order[ar]; j++){
								//if(mat[index[i]][j] == 0) printf("err\n");
								mat[index[i]][j] /= d;
							}
							for (k = 0; k < base_prd_order[ar]; k++) {
								if (k == i) continue;
								d = mat[index[k]][i];
								for (j = i; j <= base_prd_order[ar]; j++) {
									mat[index[k]][j] -= d * mat[index[i]][j];
								}	
							}
						}
					}
					
					
					
					w = (1 << coef_precision);
					e = 0.0;
					coef_p = enc->predictor[cl][co][ar];
					for (i = 0; i < base_prd_order[ar]; i++) {
						if (fabs(mat[index[i]][i]) > 1E-10) {
							d = mat[index[i]][base_prd_order[ar]] * w;
						} else {
							d = 0.0;
						}
						k = (int)d;
						if (k > d) k--;
						if (k < -max_coef) {
							k = d = -max_coef;
						} else if (k > max_coef) {
							k = d = max_coef;
						}
						coef_p[i] = k;
						d -= k;
						e += d;
						mat[index[i]][base_prd_order[ar]] = d;
					}
	/* minimize mean rounding errors */
					k = (int)(e + 0.5);
					for (;k > 0; k--) {
						d = 0.0;
						for (j = i = 0; i < base_prd_order[ar]; i++) {
							if (mat[index[i]][base_prd_order[ar]] > d) {
								d = mat[index[i]][base_prd_order[ar]];
								j = i;
							}
						}
						if (coef_p[j] < max_coef) coef_p[j]++;
						mat[index[j]][base_prd_order[ar]] = 0.0;
					}
				}
			}
		}
	free(org_buf);
	free(color_buf);
    set_prd_pels(enc);
    predict_region(enc, 0, 0, height, width);
	//printf("cpu_time = %f\n", cpu_time());
	return (calc_cost(enc, 0, 0, height, width));
}


cost_t optimize_group(ENCODER *enc)	
{
    cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
    int x, y, th1, th0, k, u, cl, co, ar, gr, prd, e, base, frac;
    int **trellis;
    PMODEL *pm, **pm_p;
	int height, width, maxprd;
	height = enc->height;
	width = enc->width;
	maxprd = enc->maxprd;
	int *class_p; 
	char *area_p;
    trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(int));
    dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
    cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));
    thc_p = enc->th_cost;
	for(gr = 0; gr < enc->num_group; gr++){
		for(u = 0; u < MAX_UPARA + 2; u++){
			trellis[gr][u] = 0;
			if(gr == 0){
				dpcost[u] = 0;
			}
		}
	}
    for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
    /* Dynamic programming */
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		for (cl = 0; cl < enc->num_class[ar]; cl++) {
			for (co = 0; co < enc->num_kind_prd; co++) {
				for (gr = 0; gr < enc->num_group; gr++) {
					cbuf_p = cbuf[gr];
					for (u = 0; u < MAX_UPARA + 2; u++) {
						cbuf_p[u] = 0;
					}
				}
				for (y = 0; y < height; y++) {
					class_p = enc->class[y];
					area_p = enc->area[y];
					for (x = 0; x < width; x++) {
						if (class_p[x] == cl && color(y, x) == co && area_p[x] == ar){
							u = enc->upara[y][x] + 1;
							e = enc->encval[y][x];
							prd = enc->prd[y][x];
							if (prd < 0) prd = 0;
							else if (prd > maxprd) prd = maxprd;
							base = enc->bconv[prd];
							frac = enc->fconv[prd];
							pm_p = enc->pmlist;
							for (gr = 0; gr < enc->num_group; gr++) {
								pm = (*pm_p++) + frac;
								cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
							}
						}
					}
				}
				for (gr = 0; gr < enc->num_group; gr++) {
				    cbuf_p = cbuf[gr];
				    for (u = 1; u < MAX_UPARA + 2; u++) {
						cbuf_p[u] += cbuf_p[u - 1];
				    }
				}
				cbuf_p = cbuf[0];
				for (u = 0; u < MAX_UPARA + 2; u++) {
				    dpcost[u] = cbuf_p[u] + thc_p[u];//dpcost:Tn(Thm(n))
				}
				for (gr = 1; gr < enc->num_group - 1; gr++) {
					cbuf_p = cbuf[gr];
					/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
					for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
						th0 = th1;
						min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
						for (k = 0; k < th1; k++) {
							cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
							if (cost < min_cost) {
								min_cost = cost;
								th0 = k;
							}
						}
						dpcost[th1] = min_cost + cbuf_p[th1];
						trellis[gr][th1] = th0;
					}
				}
				cbuf_p = cbuf[gr];
				/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
				th1 = MAX_UPARA + 1;
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1];
				for (k = 0; k < th1; k++) {
				    cost = dpcost[k] - cbuf_p[k];
				    if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
				    }
				}
				trellis[gr][th1] = th0;
				for (gr = enc->num_group - 1; gr > 0; gr--) {
				    th1 = trellis[gr][th1];
				    enc->th[cl][co][ar][gr - 1] = th1;
				}
			}
		}
	}
    /* set context quantizer */
	for(ar = 0; ar < enc->num_of_vgen; ar++){
    	for (cl = 0; cl < enc->num_class[ar]; cl++) {
			for (co = 0; co < enc->num_kind_prd; co++) {
				u = 0;
				for (gr = 0; gr < enc->num_group; gr++) {
					for (; u < enc->th[cl][co][ar][gr]; u++) {
						enc->uquant[cl][co][ar][u] = gr;
					}
				}
			}
    	}
	}
    /* renew groups */
    cost = 0;
    pm_p = enc->pmlist;
    for (y = 0; y < height; y++) {
		class_p = enc->class[y];
		area_p = enc->area[y];
		for (x = 0; x < width; x++) {
			cl = *class_p++;
			co = color(y, x);
			ar = *area_p++;
			u = enc->upara[y][x];
			enc->group[y][x] = gr = enc->uquant[cl][co][ar][u];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];
			if (prd < 0) prd = 0;
			else if (prd > maxprd) prd = maxprd;
			base = enc->bconv[prd];
			pm = pm_p[gr] + enc->fconv[prd];
			cost += pm->cost[base + e] + pm->subcost[base];
		}
    }
	
	
	
    /* optimize probability models */
    if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
		    cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel, sizeof(cost_t));
		}
		for (gr = 0; gr < enc->num_group; gr++) {
	    	for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
	    	}
		}
		for (y = 0; y < height; y++) {
	    	for (x = 0; x < width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				prd = enc->prd[y][x];
				if (prd < 0) prd = 0;
				else if (prd > enc->maxprd) prd = enc->maxprd;
				base = enc->bconv[prd];
				frac = enc->fconv[prd];
				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
	    	}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
		    pm = enc->pmodels[gr][0];
		    cost = cbuf[gr][0];
		    for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
	    	}
	    	pm_p[gr] = pm;
		}
		cost = 0.0;
		for (gr = 0; gr < enc->num_group; gr++) {
			cost += cbuf[gr][pm_p[gr]->id];
		}
    }
    free(cbuf);
    free(dpcost);
    free(trellis);
    return (cost);
}

inline void l_set_prdbuf_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area, int del_cl)
{
	int x, y, tlx, tly, brx, bry, cl, co, ar, k, l, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int org, *org_p, *nzc_p, roff = 0;
	int start_cl, end_cl, height, width;
		start_cl = 0;
		end_cl = enc->num_class[area];
	height = enc->height;
	width = enc->width;
	tlx = enc->gen_tl[area][yg][xg].x;
	tly = enc->gen_tl[area][yg][xg].y;
	brx = enc->gen_br[area][yg][xg].x + 1;
	bry = enc->gen_br[area][yg][xg].y + 1;
		for(cl = start_cl ;cl < end_cl; cl++){
			if(cl == del_cl) continue;
			prdbuf_p = prdbuf[cl];
			errbuf_p = errbuf[cl];
			for(y = tly; y < bry; y++){
				for(x = tlx; x < brx; x++){
					if(enc->area[y][x] == area && enc->label[y][x].y == yg && enc->label[y][x].x == xg){
						org_p = &enc->org[y][x];
						co = color(y, x);
						ar = area;
						nzc_p = enc->nzconv[cl][co][ar];
						if (cl == enc->class[y][x]) {
							*prdbuf_p++ = enc->prd[y][x];
							*errbuf_p++ = enc->err[y][x];
						} else {
							coef_p = enc->predictor[cl][co][ar];
							prd = 0;
							for (k = 0; k < enc->num_nzcoef[cl][co][ar]; k++) {
								l = nzc_p[k];
								//roff = ref_offset4(height, width, y, x, l, co);
								roff = ref_offset5_2(height, width, y, x, l);
								prd += org_p[roff] * (coef_p[l]);
							}
							org = *org_p;
							*prdbuf_p++ = prd;
							if (prd < 0) prd = 0;
							else if (prd > enc->maxprd) prd = enc->maxprd;
							prd >>= (enc->coef_precision - 1);
							*errbuf_p++ = enc->econv[org][prd];
						}
					}
				}
			}
		}

}
inline void set_prdbuf_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area)
{
	int x, y, tlx, tly, brx, bry, cl, co, k, l, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int org, *org_p, *nzc_p, roff = 0;
	int nzcoef, height, width, coef_precision, maxprd;
	height = enc->height;
	width = enc->width;
	coef_precision = enc->coef_precision;
	maxprd = enc->maxprd;
	int start_cl, end_cl;
	start_cl = 0;
	end_cl = enc->num_class[area];
	tlx = enc->gen_tl[area][yg][xg].x;
	tly = enc->gen_tl[area][yg][xg].y;
	brx = enc->gen_br[area][yg][xg].x + 1;
	bry = enc->gen_br[area][yg][xg].y + 1;
	
	for(cl = start_cl ;cl < end_cl; cl++){
		prdbuf_p = prdbuf[cl];
		errbuf_p = errbuf[cl];
		for(y = tly; y < bry; y++){
			for(x = tlx; x < brx; x++){
				if(enc->area[y][x] == area && enc->label[y][x].y == yg && enc->label[y][x].x == xg){
					if (cl == enc->class[y][x]) {
						*prdbuf_p++ = enc->prd[y][x];
						*errbuf_p++ = enc->err[y][x];
					} else {
						org_p = &enc->org[y][x];
						co = color(y, x);
						nzc_p = enc->nzconv[cl][co][area];
						coef_p = enc->predictor[cl][co][area];
						nzcoef = enc->num_nzcoef[cl][co][area];
						prd = 0;
						for (k = 0; k < nzcoef; k++) {
							l = nzc_p[k];
							//roff = ref_offset4(height, width, y, x, l, co);
							roff = ref_offset5_2(height, width, y, x, l);
							prd += org_p[roff] * (coef_p[l]);
						}
						org = *org_p;
						*prdbuf_p++ = prd;
						if (prd < 0) prd = 0;
						else if (prd > maxprd) prd = maxprd;
						prd >>= (coef_precision - 1);
						*errbuf_p++ = enc->econv[org][prd];
					}
				}

			}
		}
	}

}
inline int l_find_class_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area, cost_t *err_cost, int del_cl)
{
	int x, y, tlx, tly, brx, bry;
	int min_cl, cl;
	cost_t cost, min_cost;
	int *prdbuf_p, *errbuf_p;
	int start_cl = 0, end_cl = 0;

		start_cl = 0;
		end_cl = enc->num_class[area];

	tlx = enc->gen_tl[area][yg][xg].x;
	tly = enc->gen_tl[area][yg][xg].y;
	brx = enc->gen_br[area][yg][xg].x + 1;
	bry = enc->gen_br[area][yg][xg].y + 1;
	min_cost = DBL_MAX;
	min_cl = 0;


	for (cl = start_cl; cl < end_cl; cl++) {
		if(cl == del_cl) continue;
		prdbuf_p = prdbuf[cl];
		errbuf_p = errbuf[cl];
		for (y = tly; y < bry; y++) {
			for (x = tlx; x < brx; x++) {
				if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
					enc->err[y][x] = *errbuf_p++;
					enc->prd[y][x] = *prdbuf_p++;
					//enc->class[y][x]= cl;
				}
			}
		}
		cost = calc_cost_vseg(enc, tly, tlx, bry, brx, xg, yg, area, cl);
		err_cost[cl] = cost;
		//cost += enc->class_cost[area][enc->mtfbuf[cl]];
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
			enc->gen_class[area][yg][xg] = cl;
		}
	}

	prdbuf_p = prdbuf[min_cl];
	errbuf_p = errbuf[min_cl];
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			if(enc->label[y][x].x == xg && enc->label[y][x].y == yg && enc->label[y][x].area == area){
				enc->class[y][x] = min_cl;
				enc->prd[y][x] = *prdbuf_p++;
				enc->err[y][x] = *errbuf_p++;
			}
		}
	}
	return (min_cl);
}
inline int find_class_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area, cost_t *err_cost)
{
	int x, y, tlx, tly, brx, bry;
	int min_cl, cl;
	cost_t cost, min_cost;
	int *prdbuf_p, *errbuf_p;
	int start_cl = 0, end_cl = 0;
	start_cl = 0;
	end_cl = enc->num_class[area];
	tlx = enc->gen_tl[area][yg][xg].x;
	tly = enc->gen_tl[area][yg][xg].y;
	brx = enc->gen_br[area][yg][xg].x + 1;
	bry = enc->gen_br[area][yg][xg].y + 1;
	min_cost = DBL_MAX;
	min_cl = 0;
	for (cl = start_cl; cl < end_cl; cl++) {
		prdbuf_p = prdbuf[cl];
		errbuf_p = errbuf[cl];
		for (y = tly; y < bry; y++) {
				for (x = tlx; x < brx; x++) {
					if(enc->area[y][x] == area && enc->label[y][x].x == xg && enc->label[y][x].y == yg){
						enc->err[y][x] = *errbuf_p++;
						enc->prd[y][x] = *prdbuf_p++;
						//enc->class[y][x] = cl;
					}
				}
		}
		cost = calc_cost_vseg(enc, tly, tlx, bry, brx, xg, yg, area, cl);
		err_cost[cl] = cost;
		if(enc->optimize_loop == 2){
			cost += enc->class_cost[area][enc->mtfbuf[cl]];
		}
		if (cost < min_cost) {
				min_cost = cost;
				min_cl = cl;
				enc->gen_class[area][yg][xg] = cl;
		}
	}

	prdbuf_p = prdbuf[min_cl];
	errbuf_p = errbuf[min_cl];
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			if(enc->area[y][x] == area && enc->label[y][x].x == xg && enc->label[y][x].y == yg){
				enc->class[y][x] = min_cl;
				enc->prd[y][x] = *prdbuf_p++;
				enc->err[y][x] = *errbuf_p++;
			}
		}
	}
	return (min_cl);
}
inline void l_vbs_class_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area, int *blk, int del_cl)
{
	int i;
	cost_t *err_cost;
	err_cost = (cost_t *)alloc_mem(enc->num_class[area] * sizeof(cost_t));
	for(i = 0; i < enc->num_class[area]; i++){
		err_cost[i] = 0.0;
	}
	mtf_classlabel_vseg(enc->max_genx, enc->gen_class, enc->mtfbuf, xg, yg, area, enc->num_class);
	l_find_class_vseg(enc, prdbuf, errbuf, xg, yg, area, err_cost, del_cl);//vseg
	free(err_cost);
}

inline void vbs_class_vseg(ENCODER *enc, int **prdbuf, int **errbuf, int xg, int yg, int area, int *blk)
{
	int i;
	cost_t *err_cost;
	err_cost = (cost_t *)alloc_mem(enc->num_class[area] * sizeof(cost_t));
	for(i = 0; i < enc->num_class[area]; i++){
		err_cost[i] = 0.0;
	}
	mtf_classlabel_vseg(enc->max_genx,enc->gen_class, enc->mtfbuf, xg, yg, area, enc->num_class);
	

	find_class_vseg(enc, prdbuf, errbuf, xg, yg, area, err_cost);//vseg


	free(err_cost);
}


void count_cl(ENCODER *enc)
{
    int cl, y, x;

    for (cl = 0; cl < enc->max_class; cl++) enc->cl_hist[cl] = 0;
    for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
			    enc->cl_hist[(int)enc->class[y][x]]++;
			}
    }
}
inline void l_optimize_class_vseg(ENCODER *enc, int del_cl, int ar)
{
    int i, blk = 0;
    int **prdbuf, **errbuf;
	int xg, yg, area;
	int j, k;
    for (i = 0; i < enc->max_class; i++) {
        enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
    }
    prdbuf =(int **)alloc_2d_array(enc->max_class, enc->max_num_pixel,sizeof(int));
    errbuf =(int **)alloc_2d_array(enc->max_class, enc->max_num_pixel,sizeof(int));
	for(j = 0; j < enc->max_class; j++){
		for(k = 0; k < enc->max_num_pixel; k++){
			prdbuf[j][k] = 0;
			errbuf[j][k] = 0;
		}
	}

	for (i = 0; i < enc->max_class; i++) { //mtf
	enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
	}

	//if(ar == 0){
		/*area = 0;
		for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg ++){
			for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg ++){
				if(enc->num_pixel[area][yg][xg] == 0) continue;
				if(enc->gen_class[area][yg][xg] != del_cl) continue;
				l_set_prdbuf_vseg(enc, prdbuf, errbuf, xg, yg, area, del_cl);
				l_vbs_class_vseg(enc, prdbuf, errbuf, xg, yg, area, &blk, del_cl);
				blk++;
			}
		}
		for (i = 0; i < enc->max_class; i++) {
			enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
		}*/
	//}else{
		area = ar;
		for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
			for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
				if(enc->num_pixel[area][yg][xg] == 0) continue;
				if(enc->gen_class[area][yg][xg] != del_cl) continue;
				l_set_prdbuf_vseg(enc, prdbuf, errbuf, xg, yg, area, del_cl);
				l_vbs_class_vseg(enc, prdbuf, errbuf, xg, yg, area, &blk, del_cl);
				blk++;
			}
		}
		for (i = 0; i < enc->max_class; i++) {
		enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
		}
	//}
    free(errbuf);
    free(prdbuf);
    count_cl(enc);
}
//vseg
cost_t optimize_class_vseg(ENCODER *enc)
{
    int i, blk = 0;
    int **prdbuf, **errbuf;
	int xg, yg, area;
	int j, k;
    for (i = 0; i < enc->max_class; i++) {
        enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
    }
    prdbuf =(int **)alloc_2d_array(enc->max_class, enc->max_num_pixel, sizeof(int));
    errbuf =(int **)alloc_2d_array(enc->max_class, enc->max_num_pixel, sizeof(int));
	for(j = 0; j < enc->max_class; j++){
		for(k = 0; k < enc->max_num_pixel; k++){
			prdbuf[j][k] = 0;
			errbuf[j][k] = 0;
		}
	}
	for (i = 0; i < enc->max_class; i++) { //mtf
		enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
	}
	/*area = 0;
	for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg ++){
		for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg ++){
			if(enc->num_pixel[area][yg][xg] == 0) continue;
			set_prdbuf_vseg(enc, prdbuf, errbuf, xg, yg, area);
			vbs_class_vseg(enc, prdbuf, errbuf, xg, yg, area, &blk);
			blk++;
		}
	}
	for (i = 0; i < enc->max_class; i++) {
	  enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
	}*/
	for(area = 0; area < enc->num_of_vgen; area++){
		for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
			for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
				if(enc->num_pixel[area][yg][xg] == 0) continue;
				set_prdbuf_vseg(enc, prdbuf, errbuf, xg, yg, area);
				vbs_class_vseg(enc, prdbuf, errbuf, xg, yg, area, &blk);
				blk++;
			}
		}
		for (i = 0; i < enc->max_class; i++) {
			enc->mtfbuf[i] = i; //mtfbuf[]={0,1,2,3...}
		}
	}
    free(errbuf);
    free(prdbuf);
    count_cl(enc);
    return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

inline void set_prd_pels(ENCODER *enc)
{
    int cl, co, ar, k, i;
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for (cl = 0; cl < enc->num_class[ar]; cl++) {
				for (co = 0; co < enc->num_kind_prd; co++) {
					k = 0;
					for (i = 0; i < enc->max_prd_order[ar]; i++) {
						if (enc->predictor[cl][co][ar][i] != 0) {
							enc->nzconv[cl][co][ar][k] = i;
							k++;
						}
					}
					for (i = k; i < enc->max_prd_order[ar]; i++) {
						enc->nzconv[cl][co][ar][i] = 1E5;
					}
					enc->num_nzcoef[cl][co][ar] = k;
				}
			}
		}
}

inline void optimize_coef2(ENCODER *enc, int cl, int co, int ar, int pos, int *num_eff)
{//pos:ランダムに選ばれた次数
#define S_RANGE 2 //係数の修正の範囲。2だと+-1の範囲
#define SUBS_RANGE 2 //修正の範囲。（一つの係数のみを修正する場合の）
    cost_t *cbuf, *cbuf_p, *cbuf_p2, c_base, *coef_cost_p1, *coef_cost_p2;
    int i, j, k, l, x, y, df1, df2, df_f, *org_p, base;
    int coef_abs_pos, coef_abs, coef_pos;
    int num_swap_search, *swap_search;
    int prd, /*shift,*/ maxprd, *coef_p,/* *nzc_p,*/ prd_c;
    int *class_p, onoff, class, max_prd_order, max_coef, height, width;
	char *area_p, area;
    img_t *bconv_p, *fconv_p;
    PMODEL *pm, *pm_p;
    int diffy, diffx, *diff, *diff_p, *diff_p2, sgn;
	max_prd_order = enc->max_prd_order[ar];
	max_coef = enc->max_coef;
	height = enc->height;
	width = enc->width;
	//printf("***optimize_coef2***\n");
    num_swap_search = max_prd_order - enc->num_nzcoef[cl][co][ar];
	//cost_t cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search];
    cbuf = (cost_t *)alloc_mem((SUBS_RANGE + (max_prd_order * S_RANGE) + num_swap_search) * sizeof(cost_t));
	//int diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)+ num_swap_search]; 
    diff = (int *) alloc_mem ((SUBS_RANGE + (max_prd_order * S_RANGE)+ num_swap_search) * sizeof(int));
    //交換の対象となる係数の位置のテーブル。
	//int swap_search[num_swap_search];
    swap_search = (int *)alloc_mem(num_swap_search * sizeof(int));
	for(i = 0; i < num_swap_search; i++){
		swap_search[i] = 0;
	}
	for(i = 0; i < SUBS_RANGE + (max_prd_order * S_RANGE) + num_swap_search; i++){
		cbuf[i] = 0;
		diff[i] = 0;
	}
    //nzc_p = enc->nzconv[cl][co][ar];
    coef_p = enc->predictor[cl][co][ar];
    coef_pos = coef_p[pos];//coef_pos:ランダムに選ばれた係数
    coef_abs_pos = (coef_pos < 0) ? -coef_pos : coef_pos;//coef_abs_pos:ランダムに選ばれた係数の絶対値
    cbuf_p = cbuf;
    cbuf_p2 = &cbuf[SUBS_RANGE + (max_prd_order * S_RANGE)];//cbuf_p2:cbuf[num_swap_search]から指す
    diff_p = diff;
    diff_p2 = &diff[SUBS_RANGE + (max_prd_order * S_RANGE)];//diff_p2:diff_p2[num_sawp_search]から指す
    //i = enc->ord2mhd[pos];//ord2mhd:マンハッタンディスタンス？
    coef_cost_p1 = enc->coef_cost2[enc->zero_m[pos]][enc->coef_m[pos]];
    //printf("check10\n");
	/* 係数の分布が指数分布になることに基づいて付加情報を推定する？　*/
	/* ランダムに選んだ係数を単体で修正する場合　*/
	for (i = 0; i < (SUBS_RANGE >> 1); i++) {
		df1 = coef_pos - (i + 1);//ランダムに選んだ係数に微小変化を加える
		
		for (j = 0; j < 2; j++) {
			y = df1;
			sgn = 0;
			
			if (y < 0) {//係数が負の時はsgnに1を格納
				y = -y;
				sgn = 1;
			}
			if (y > max_coef) y = max_coef;
			*cbuf_p++ = coef_cost_p1[y];//修正した係数のコストを保存？
			*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;//微小変化分を保存？
			df1 += (i + 1) << 1;
		}
    }
	
	/* ランダムに選んだ係数と対になる係数をもう一つ選び、双方修正する場合　*/
    k = 0;
	//printf("check20\n");
    for (i = 0; i < max_prd_order; i++) {//交換する係数を探索するループ？
		onoff = (coef_p[i] == 0) ? 1 : 0;//coef_posが0かどうかでフラグを立てとく
		//j = enc->ord2mhd[i];//修正した次数のマンハッタンディスタンスを格納
		coef_cost_p2 = enc->coef_cost2[enc->zero_m[i]][enc->coef_m[i]];
		df2 = coef_p[i];
        coef_abs = (df2 < 0) ? -df2 : df2;
		c_base = coef_cost_p2[coef_abs];
		for (l = 0; l < (S_RANGE >> 1); l++) {
			df1 = coef_pos - (l + 1);
			df2 = coef_p[i] + (l + 1);
			
			for (j = 0; j < 2; j++) { //修正の+-がj=0の時と1の時で切り替わる
				y = df1;//posの係数の修正後がy
				x = df2;//上記と対になる係数の修正後がx
				sgn = 0;
				if (y < 0) {
					y = -y;
					sgn = 1;//sgnにyの修正後の符号を格納
				}
				diffy = y - max_coef;
				if (diffy < 0) diffy = 0;//y>max_coefのとき、はみ出た分がdiffyに値が入る
				if (y > max_coef) y = max_coef;//はみ出たらmax_coefに修正
				
				if (x < 0) x = -x;//xはかならずプラス？
				diffx = x - max_coef;
				if (diffx < 0) diffx = 0;//|x| > max_coefのとき、はみ出た分がdiffxに値が入る
				if (x > max_coef) x = max_coef;//はみ出たらmax_coefに修正
				if (diffy > 0 || diffx > 0) {//yかxがはみ出てた時
					if (diffy > diffx) {//はみ出てる度がyの方が多いとき
						x += (j) ? diffy - diffx : diffx - diffy;//修正が
						if (x < 0) {
							x = -x;
						}
					} else {
						y += (j) ? diffy - diffx : diffx - diffy;
						if (y < 0) {
							y = -y;
							sgn = (sgn) ? 0 : 1;
						}
					}
				}//この処理いみわからんちん
				
				*cbuf_p++ = coef_cost_p1[y] + coef_cost_p2[x] - c_base;//付加情報の増減だけ保存？
				*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;//yの微小変化分を保存
				df1 += (l + 1) << 1;
				df2 -= (l + 1) << 1;
			}
		}
		if (onoff == 1) {//coef_pos　== 0のとき
			*cbuf_p2++ = coef_cost_p1[0] + coef_cost_p2[coef_abs_pos] - coef_cost_p2[0];//p1がpos、p2が対になるやつ？
			*diff_p2++ = -coef_pos;
			swap_search[k++] = i;//swap_searchに対になるやつの次数が入る
		}
    }
    for (i = 0; i < S_RANGE; i++) {
		cbuf[SUBS_RANGE + pos * S_RANGE + i] = coef_cost_p1[coef_abs_pos];//cbufの配列の仕組みがいまいちわからん
		diff[SUBS_RANGE + pos * S_RANGE + i] = 0;
    }
    // ここまでは付加情報の計算
	
	//printf("check30\n");
    bconv_p = enc->bconv;
    fconv_p = enc->fconv;
    maxprd = enc->maxprd;
    //shift = enc->coef_precision - 1;
	class_p = enc->class[0];
	area_p = enc->area[0];
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            class = *class_p++;
			area = *area_p++;
			if (cl != class) continue;
			if (ar != (int)area) continue;
            if (co != color(y, x)) continue;
			prd = enc->prd[y][x];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[ref_offset5_2(height, width, y, x, pos)];
			cbuf_p = cbuf;
			cbuf_p2 = &cbuf[SUBS_RANGE + (max_prd_order * S_RANGE)];
			diff_p = diff;
			diff_p2 = &diff[SUBS_RANGE + (max_prd_order * S_RANGE)];
			df_f = df1;
			//係数を1つしか修正しない
			for (i = 0; i < (SUBS_RANGE >> 1); i++) {
				for (j = 0; j < 2; j++) { //修正の+-がj=0の時と1の時で切り替わる
					prd_c = prd + df_f * (*diff_p++);//prd_cは修正後の新しい予測値
					if (prd_c < 0) prd_c = 0;
					else if (prd_c > maxprd) prd_c = maxprd;//はみ出たら直す
					base = bconv_p[prd_c];
					pm = pm_p + fconv_p[prd_c];
					*cbuf_p++ += pm->cost[*org_p + base] + pm->subcost[base];//cbufに予測誤差符号量を保存
				}
			}
			//修正と修正の相手を探索
			for (i = 0; i < max_prd_order; i++) {
				df2 = org_p[ref_offset5_2(height, width, y, x, i)];
				df_f = df1 - df2;
				for (l = 0; l < (S_RANGE >> 1); l++) {
					for (j = 0; j < 2; j++) { //修正の+-がj=0の時と1の時で切り替わる
						prd_c = prd + df_f * (*diff_p++);
                        if (prd_c < 0) prd_c = 0;
                        else if (prd_c > maxprd) prd_c = maxprd;
                        base = bconv_p[prd_c];
                        pm = pm_p + fconv_p[prd_c];
                        *cbuf_p++ += pm->cost[*org_p + base] + pm->subcost[base];
					}
				}
			}
			//交換
			for (j = 0; j < num_swap_search; j++) {
				k = swap_search[j];
				prd_c = prd + (df1 - org_p[ref_offset5_2(height, width, y, x, k)]) * (*diff_p2++);
				if (prd_c < 0) prd_c = 0;
				else if (prd_c > maxprd) prd_c = maxprd;
				base = bconv_p[prd_c];
				pm = pm_p + fconv_p[prd_c];
				*cbuf_p2++ += pm->cost[*org_p + base] + pm->subcost[base];
			}
		}
    }
	//printf("check40\n");
    j = SUBS_RANGE + pos * S_RANGE;
    for (i = 0; i < SUBS_RANGE + (max_prd_order * S_RANGE) + num_swap_search; i++){//ここで総コストを計算？
		if (cbuf[i] < cbuf[j]) {
			j = i;
		}
    }
    
    if (j == SUBS_RANGE + pos * S_RANGE) {
		free(swap_search);
		cbuf_p = NULL;
		free(cbuf);
		//swap_search = NULL;
		diff_p = NULL;
		diff_p2 = NULL;
		free(diff);
		//printf("check41\n");
		return;
    }
    if (j < SUBS_RANGE) {//jがcbufの配列のうちSUBS_RANGEの部分なら
		class_p = enc->class[0];
		area_p = enc->area[0];
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
                class = *class_p++;
				area = *area_p++;
				if (cl == class && ar == (int)area && co == color(y, x)) {
					org_p = &enc->org[y][x];
					enc->prd[y][x] += org_p[ref_offset5_2(height, width, y, x, pos)] * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
    } else {
		if (j < SUBS_RANGE + (max_prd_order * S_RANGE)) {//jがcbufのswap_search以前なら
			i = j - SUBS_RANGE;
			i /= S_RANGE;
		} else {//jがcbufのswap_serch以降なら
			i = j - SUBS_RANGE - (max_prd_order * S_RANGE);
			i = swap_search[i];
		}
		class_p = enc->class[0];
		area_p = enc->area[0];
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
                class = *class_p++;
				area = *area_p++;
				if (cl == class && ar == (int)area && co == color(y, x)) {
					org_p = &enc->org[y][x];
					enc->prd[y][x] += (org_p[ref_offset5_2(height, width, y, x, pos)]
					- org_p[ref_offset5_2(height, width, y, x, i)]) * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
		coef_p[i] -= diff[j];
    }
	
    free(swap_search);
	diff_p = NULL;
	diff_p2 = NULL;
    free(diff);
	cbuf_p = NULL;
	free(cbuf);
	if (diff[j] != 0) (*num_eff)++;
	//printf("fin\n");
}

cost_t optimize_predictor2(ENCODER *enc)//予測器の最適化
{
    int cl, co, ar, pos, k, /*num_nzc,*/ num_eff;
	int height, width;
	height = enc->height;
	width = enc->width;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif
//num_searchは探索回数？
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		for (cl = 0; cl < enc->num_class[ar]; cl++) {
			for (co = 0; co < enc->num_kind_prd; co++){
				//num_nzc = enc->num_nzcoef[cl][co][ar];
				num_eff = 0;
				if (enc->cl_hist[cl] == 0) continue;
				for (k = 0; k < enc->num_search[cl][co][ar]; k++) {//num_searchを変更する
					if (enc->num_nzcoef[cl][co][ar] == 0) continue;
					pos = (int)(((double)rand() * enc->num_nzcoef[cl][co][ar]) / (RAND_MAX+1.0));
					pos = enc->nzconv[cl][co][ar][pos];
					optimize_coef2(enc, cl, co, ar, pos, &num_eff);//予測器係数の修正
					set_prd_pels(enc);
				}
				enc->num_search[cl][co][ar] = num_eff + 3;//num_searchを変更する
			}
		}
	}
	set_prd_pels(enc);
    predict_region(enc, 0, 0, height, width);//予測値と予測誤差を各画素ごとに保存
    return (calc_cost(enc, 0, 0, height, width));
}

void optimize_coef(ENCODER *enc, int cl, int co, int ar, int pos1, int pos2)
{
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
    cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
    float *pmcost_p;
    int i, j, k, x, y, df1, df2, base;
    int prd, prd_f, shift, maxprd, max_coef;
    int *coef_p, *econv_p, *org_p;
    int *class_p, class;
		char *area_p, area;
    img_t *bconv_p, *fconv_p;
    PMODEL *pm, *pm_p;
	int height, width;
	height = enc->height;
	width = enc->width;
	max_coef = enc->max_coef;
    cbuf_p = cbuf;
    coef_p = enc->predictor[cl][co][ar];
    k = 0;
    for (i = 0; i < SEARCH_RANGE; i++) {
    	y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);
      if (y < 0) y = -y;
      if (y > max_coef) y = max_coef;
			for (j = 0; j < SUBSEARCH_RANGE; j++) {
        x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1)) - (j - (SUBSEARCH_RANGE >> 1));
      	if (x < 0) x = -x;
        if (x > max_coef) x = max_coef;
        cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y] + enc->coef_cost[enc->coef_m[pos2]][x];
			}
    }
		
    bconv_p = enc->bconv;
    fconv_p = enc->fconv;
    maxprd = enc->maxprd;
    shift = enc->coef_precision - 1;
	class_p = enc->class[0];
	area_p = enc->area[0];
    for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
        class = *class_p++;
		area = *area_p++;
        if (cl != class) continue;
		if (ar != area) continue;
        if (co != color(y, x)) continue;
        prd = enc->prd[y][x];
        org_p = &enc->org[y][x];
        pm_p = enc->pmlist[(int)enc->group[y][x]];
        df1 = org_p[ref_offset5_2(height, width, y, x, pos1)];
        df2 = org_p[ref_offset5_2(height, width, y, x, pos2)];
        prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1) + df2 * (SUBSEARCH_RANGE >> 1);
        cbuf_p = cbuf;
	    	if (enc->pm_accuracy < 0) {
          econv_p = enc->econv[*org_p];
          pmcost_p = pm_p->cost;

					for (i = 0; i < SEARCH_RANGE; i++) {
		    		for (j = 0; j < SUBSEARCH_RANGE; j++) {
              prd = prd_f;
              if (prd < 0) prd = 0;
              else if (prd > maxprd) prd = maxprd;
              (*cbuf_p++) += pmcost_p[econv_p[prd >> shift]];
              prd_f -= df2;
		    		}
            prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
					}
	    	} else {
					for (i = 0; i < SEARCH_RANGE; i++) {
		    		for (j = 0; j < SUBSEARCH_RANGE; j++) {
              prd = prd_f;
              if (prd < 0) prd = 0;
              else if (prd > maxprd) prd = maxprd;
              base = bconv_p[prd];
              pm = pm_p + fconv_p[prd];
              (*cbuf_p++) += pm->cost[*org_p + base] + pm->subcost[base];
              prd_f -= df2;
		    		}
          	prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
					}
	    	}
			}
    }
    cbuf_p = cbuf;
    j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;
    for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
			if (cbuf_p[i] < cbuf_p[j]) {
        j = i;
			}
    }
    i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
    j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
    y = coef_p[pos1] + i;
    x = coef_p[pos2] - i - j;
    if (y < -max_coef) y = -max_coef;
    else if (y > max_coef) y = max_coef;
    if (x < -max_coef) x = -max_coef;
    else if (x > max_coef) x = max_coef;
    i = y - coef_p[pos1];
    j = x - coef_p[pos2];
    if (i != 0 || j != 0) {
	class_p = enc->class[0];
	area_p = enc->area[0];
	for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
        	class = *class_p++;
			area = *area_p++;
			if (cl == class && ar == area && co == color(y, x)) {
				org_p = &enc->org[y][x];
				enc->prd[y][x] += org_p[ref_offset5_2(height, width, y, x, pos1)]
				* i + org_p[ref_offset5_2(height, width, y, x, pos2)] * j;
			}
	    }
	}
      coef_p[pos1] += i;
      coef_p[pos2] += j;
    }
}

cost_t optimize_predictor(ENCODER *enc)
{
    int cl, co, ar, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif
	for(ar = 0; ar < enc->num_of_vgen; ar++){
    	for (cl = 0; cl < enc->num_class[ar]; cl++) {
			for (co = 0; co < enc->num_kind_prd; co++) {
				for (k = 0; k < enc->max_prd_order[ar]; k++) {
	        		retry:
					pos1 = (int)(((double)rand() * enc->max_prd_order[ar]) / (RAND_MAX+1.0));
					pos2 = (int)(((double)rand() * enc->max_prd_order[ar]) / (RAND_MAX+1.0));
					if (pos1 == pos2) goto retry;
					optimize_coef(enc, cl, co, ar, pos1, pos2);

	    		}
			}
		}
    }
	set_prd_pels(enc);
    predict_region(enc, 0, 0, enc->height, enc->width);
    return (calc_cost(enc, 0, 0, enc->height, enc->width));
}



int putbits(FILE *fp, int n, uint x)
{
    static int bitpos = 8;
    static uint bitbuf = 0;
    int bits;

    bits = n;
    if (bits <= 0) return (0);
    while (n >= bitpos) {
        n -= bitpos;
				if (n < 32) {
			            bitbuf |= ((x >> n) & (0xff >> (8 - bitpos)));
				}
        putc(bitbuf, fp);
        bitbuf = 0;
        bitpos = 8;
    }
    bitpos -= n;
    bitbuf |= ((x & (0xff >> (8 - n))) << bitpos);
    return (bits);
}

void remove_emptyclass_vseg(ENCODER *enc)
{
    int cl, co, i, k, x, y, max;
		int xg, yg, area, ar;
		int **buf;

		buf = (int **)alloc_2d_array(enc->num_of_vgen, MAX_CLASS, sizeof(int));
		printf("remove_emptyclass.\n");
		for(ar = 0; ar < enc->num_of_vgen; ar++){
    	for (cl = 0; cl < MAX_CLASS; cl++) {
				buf[ar][cl] = 0;
			}
    }
		/*area = 0;
		for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg ++){
			for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg ++){
				cl = (int)enc->gen_class[area][yg][xg];
				if(cl < 0) continue;
				buf[area][cl]++;
			}
		}*/
		for(area = 0; area < enc->num_of_vgen; area++){
			for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
				for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
					cl = (int)enc->gen_class[area][yg][xg];
					if(cl < 0) continue;
					buf[area][cl]++;
				}
			}
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
    	for (i = cl = 0; i < enc->num_class[ar]; i++) {
				//printf("buf[%d][%d] = %d\n", ar, i, buf[ar][i]);
				if (buf[ar][i] == 0) {
			  	buf[ar][i] = -1;
					//if(i < k) enc->num_center_class--;
				} else {
			  	buf[ar][i] = cl++;
				}
			}
    }
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(cl = 0; cl < enc->num_class[ar]; cl++){
    		if (buf[ar][cl] < 0){
					goto change;
	  		}else{

				}	/* no empty class */
			}
		}
		printf("num of class is not changed.\n");
		return;


		change:
		printf("num of class is changed.\n");
		/*area = 0;
		for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg ++){
			for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg ++){
				i = (int)enc->gen_class[area][yg][xg];
				if(i < 0) continue;
				enc->gen_class[area][yg][xg] = (int)buf[area][i];
			}
		}*/
		for(area = 0; area < enc->num_of_vgen; area++){
			for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
				for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
					i = (int)enc->gen_class[area][yg][xg];
					if(i < 0) continue;
					enc->gen_class[area][yg][xg] = (int)buf[area][i];
				}
			}
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
    	for (i = cl = 0; i < enc->num_class[ar]; i++) {
      	if (buf[ar][i] < 0) continue;
      	if (cl != i) {
			    for (co = 0; co < enc->num_kind_prd; co++) {
		        for (k = 0; k < enc->max_prd_order[ar]; k++) {
							enc->predictor[cl][co][ar][k] = enc->predictor[i][co][ar][k];
						}
						for (k = 0; k < enc->num_group - 1; k++) {
							enc->th[cl][co][ar][k] = enc->th[i][co][ar][k];
						}
					}
				}
				cl++;
			}
			enc->num_class[ar] = cl;
    }
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
			  yg = enc->label[y][x].y;
				xg = enc->label[y][x].x;
				ar = enc->label[y][x].area;
				i = (int)enc->gen_class[ar][yg][xg];

				enc->class[y][x] = i;

			}
    }
		max = 0;
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			if(max < enc->num_class[ar]){
				max = enc->num_class[ar];
			}
		}
		enc->max_class = max;
  
	set_prd_pels(enc);
	free(buf);
}


int write_header(ENCODER *enc, FILE *fp)
{
	printf("write_header...");
    int bits, ar;
    bits = putbits(fp, 16, MAGIC_NUMBER);
    bits += putbits(fp, 8, VERSION);
    bits += putbits(fp, 16, enc->width);
    bits += putbits(fp, 16, enc->height);
    bits += putbits(fp, 16, enc->maxval);
    bits += putbits(fp, 4, 1);	/* number of components (1 = monochrome) */
	for(ar = 0; ar <enc->num_of_vgen; ar++){
		bits += putbits(fp, 8, enc->num_class[ar]);
	}
    bits += putbits(fp, 6, enc->num_group);
	for(ar = 0; ar <enc->num_of_vgen; ar++){
		bits += putbits(fp, 8, enc->max_prd_order[ar]);
	}
    bits += putbits(fp, 3, enc->num_kind_prd);
    bits += putbits(fp, 6, enc->num_pmodel - 1);
    bits += putbits(fp, 4, enc->coef_precision - 1);
    bits += putbits(fp, 3, enc->pm_accuracy + 1);
    bits += putbits(fp, 1, enc->f_huffman);
    bits += putbits(fp, 1, (enc->quadtree_depth < 0)? 0 : 1);
	printf("ok\n");
    return (bits);
}

int encode_golomb(FILE *fp, int m, int v)
{
    int bits, p;

    bits = p = (v >> m) + 1;
    while (p > 32) {
	putbits(fp, 32, 0);
	p -= 32;
    }
    putbits(fp, p, 1);	/* prefix code */
    putbits(fp, m, v);
    return (bits + m);
}
int encode_class_vseg(FILE *fp, ENCODER *enc, int *qt_cost, int flag)
{
	//printf("encode_class->");
    int i, j, l, cl, ctx, bits, /*buff,*/***index;
    uint **hist;
    cost_t cost;
		int xg, yg, area;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif

    hist = (uint **)alloc_2d_array(enc->num_of_vgen, enc->max_class, sizeof(uint));
    index = (int ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(int));

    for (i = 0; i < enc->max_class; i++) {
		for(l = 0; l < enc->num_of_vgen; l++){
			hist[l][i] = 0;
		}
			enc->mtfbuf[i] = i;
    }
		
			for(area = 0; area < enc->num_of_vgen; area++){			
				for(yg = 0; yg < enc->num_of_geny[area]; yg ++){
					for(xg = 0; xg < enc->num_of_genx[area]; xg ++){
						cl = enc->gen_class[area][yg][xg];
						if(cl < 0) continue;
						mtf_classlabel_vseg(enc->max_genx, enc->gen_class, enc->mtfbuf, xg, yg, area,
							 									enc->num_class);
						i = enc->mtfbuf[cl];
						index[area][yg][xg] = i;
						hist[area][i]++;
					}
				}
				for (i = 0; i < enc->max_class; i++) {
					enc->mtfbuf[i] = i;
				}
			}
			//printf("3,");
    	bits = 0;
    			/********** Arithmetic (???δづつ算ﾃδづつ??)***********/
			PMODEL *pm;
			double p[enc->num_of_vgen], c[enc->num_of_vgen];//vseg
			//int qtree_code[QUADTREE_DEPTH << 2];
			int mtf_code[enc->num_of_vgen][enc->max_class];
			//cost_t qtflag_cost[QUADTREE_DEPTH << 3];
			cost_t **class_cost;

			class_cost = (cost_t **)alloc_2d_array(enc->num_of_vgen, enc->max_class, sizeof(cost_t));
			/* quantization of log-transformed probability */
			for(l = 0; l < enc->num_of_vgen; l++){
				c[l] = 0.0;//vseg
		  }
			for(l = 0; l < enc->num_of_vgen; l++){
				for (i = 0; i < enc->num_class[l]; i++){
			    c[l] += (double)hist[l][i];
				}
			}
			//printf("4,");
			for(l = 0; l < enc->num_of_vgen; l++){
				for (i = 0; i < enc->num_class[l]; i++) {
					//printf("area%d hist %d\n", l,hist[l][i]);
				    p[l] = (double)hist[l][i] / c[l];
				    if (p[l] > 0.0) {
							mtf_code[l][i] = -log(p[l]) / log(2.0)
							    							* (PMCLASS_LEVEL / PMCLASS_MAX);
							if (mtf_code[l][i] >= PMCLASS_LEVEL) {
							    mtf_code[l][i] = PMCLASS_LEVEL - 1;
							}
				    } else {
							mtf_code[l][i] = PMCLASS_LEVEL - 1;
				    }
				    p[l] = exp(-log(2.0) * ((double)mtf_code[l][i] + 0.5)
					    										* PMCLASS_MAX / PMCLASS_LEVEL);
				    class_cost[l][i] = -log(p[l]) / log(2.0);
				    hist[l][i] = p[l] * (1 << 10);
				    if (hist[l][i] <= 0) hist[l][i] = 1;
				}
			}
			//printf("ip");
			//printf("5,");
			if (fp == NULL) {
			    cost = 0.0;
					
					for(l = 0; l < enc->num_of_vgen; l++){
						enc->area_class_info[l] = 0;
						for(yg = 0; yg < enc->num_of_geny[l]; yg ++){
							for(xg = 0; xg < enc->num_of_genx[l]; xg ++){
								if(enc->gen_class[l][yg][xg] < 0) continue;
								i = index[l][yg][xg];
								cost += class_cost[l][i];
								enc->area_class_info[l] += class_cost[l][i];
							}
						}
						//printf("ip");
					
				    if (flag == 1) {
			              //  for (i = 0; i < (QUADTREE_DEPTH << 3); i++)
			              //    					enc->qtflag_cost[i] = qtflag_cost[i];
			                for (i = 0; i < enc->num_class[l]; i++)
			                  				enc->class_cost[l][i] = class_cost[l][i];
				    }
				  }
					//printf("8,");

			    bits = (int)cost;
			} else {	/* actually encode */
			printf("9,");
			    PMODEL cpm[1];
			    /* additional info. */
			    pm = &enc->spm;
					
					for(l = 0; l < enc->num_of_vgen; l++){
						enc->area_class_info[l] = 0;
						cpm->freq = (uint *)alloc_mem((enc->num_class[l] * 2 + 1) * sizeof(uint));
					  cpm->size = enc->num_class[l];
				    set_spmodel(pm, PMCLASS_LEVEL, -1);
						cpm->cumfreq = &(cpm->freq[cpm->size]);
						cpm->cumfreq[0] = 0;
				    for (i = 0; i < enc->num_class[l]; i++) {
							j = mtf_code[l][i];
							rc_encode(fp, enc->rc, pm->cumfreq[j], pm->freq[j],
								  pm->cumfreq[pm->size]);
							if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
							    for (j = 0; j < pm->size; j++) {
										if (j < mtf_code[l][i]) {
										    pm->freq[j] /= 2;
										} else {
										    pm->freq[j] *= 2;
										}
										if (pm->freq[j] <= 0) pm->freq[j] = 1;
										pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
							    }
							}
				    }
					
				    /* set prob. models */
						//printf("rc->code %llu\n", enc->rc->code);

				    for (i = 0; i < enc->num_class[l]; i++) {
							cpm->freq[i] = hist[l][i];
							cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
							//printf("hist %d\n", hist[l][i]);
							//printf("cpm->cumfreq[i + 1] %d\n", cpm->cumfreq[i + 1]);
				    }
						//printf("cpm->size %d\n", cpm->size);
					
			      //buff = 0;
						for (yg = 0; yg < enc->num_of_geny[l]; yg++) {
							for (xg = 0; xg < enc->num_of_genx[l]; xg++) {
								if(enc->gen_class[l][yg][xg] < 0) continue;
								i = index[l][yg][xg];
								if (i < 0) {
										i = -(i + 1);
										ctx = i & (~1);
										rc_encode(fp, enc->rc,
														pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i],
														pm->cumfreq[ctx + 2] - pm->cumfreq[ctx]);
																//	*qt_cost += enc->rc->code - buff;
								} else {
										rc_encode(fp, enc->rc, cpm->cumfreq[i], cpm->freq[i],
												cpm->cumfreq[cpm->size]);
								}
													//buff = enc->rc->code;
							}
						}
						
				    bits += enc->rc->code;
					enc->area_class_info[l] += enc->rc->code;
						//printf("rc->code %llu\n", enc->rc->code);
				    enc->rc->code = 0;
						free(cpm->freq);
					}
					printf("13,");

		  }

		free(class_cost);
    free(index);
    free(hist);
	//printf("ok\n");
    return (bits);
}




int encode_predictor2(FILE *fp, ENCODER *enc, int flag)
{//nonzerocoef is pre-encoded.
    int cl, co, ar, coef, sgn, k, m, min_m, bits, buff = 0.0;
    cost_t cost, min_cost, t_cost;
	uint cumb;
	int zrfreq, nzfreq;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    t_cost = 0.0;
	for(ar = 0; ar < enc->num_of_vgen;ar++) enc->area_predictor_info[ar] = 0;
    for (k = 0; k < enc->max_prd_order_all; k++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
	    	cost = 0.0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						if(k < enc->max_prd_order[ar]){
							coef = enc->predictor[cl][co][ar][k];
							if (coef < 0) coef = -coef;
							cost += enc->coef_cost2[enc->zero_m[k]][m][coef];//変更する
							//enc->area_predictor_info[ar] += enc->coef_cost2[enc->zero_m[k]][m][coef];
						}
					}
				}
	    	}
			if (cost < min_cost) {
					min_cost = cost;
					min_m = m;
			}
		}
		//t_cost += min_cost;//この文必要ないの？
		if (flag) enc->coef_m[k] = min_m;
    }
	for (k = 0; k < enc->max_prd_order_all; k++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < NUM_ZMODEL; m++) {
	    	cost = 0.0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				//enc->area_predictor_info[ar] = 0;
				for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						if(k < enc->max_prd_order[ar]){
							coef = enc->predictor[cl][co][ar][k];
							if (coef < 0) coef = -coef;
							cost += enc->coef_cost2[m][enc->coef_m[k]][coef];//変更する
							//enc->area_predictor_info[ar] += enc->coef_cost2[m][enc->coef_m[k]][coef];
						}
					}
				}
	    	}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		t_cost += min_cost;
		if (flag) enc->zero_m[k] = min_m;
    }
    bits = (int)t_cost;//符号量見積の計算はここまで？
    
	
	if (fp != NULL) {//実際に符号化
		bits = 0;
		/* Arithmetic */
	    	PMODEL *pm;
	    	pm = &enc->spm;
			for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_predictor_info[ar] = 0;
	    	for (k = 0; k < enc->max_prd_order_all; k++) {
				//printf("zero_m[%d] = %d\n",k ,enc->zero_m[k]);
				rc_encode(fp, enc->rc, enc->zero_m[k], 1, NUM_ZMODEL);
				rc_encode(fp, enc->rc, enc->coef_m[k], 1, 8);
				nzfreq = enc->zero_fr[enc->zero_m[k]];
				zrfreq = TOT_ZEROFR - nzfreq;
				set_spmodel(pm, enc->max_coef + 1, enc->coef_m[k]);
				cumb = pm->freq[0];//cumbにラプラス関数に0を入力したときの出力を格納
				for(ar = 0; ar < enc->num_of_vgen; ar++){
					
					for (cl = 0; cl < enc->num_class[ar]; cl++) {
						for (co = 0; co < enc->num_kind_prd; co++) {
							if(k < enc->max_prd_order[ar]){
								coef = enc->predictor[cl][co][ar][k];
								if(coef == 0){
									rc_encode(fp, enc->rc, 0, zrfreq, TOT_ZEROFR);
								}else{
									rc_encode(fp, enc->rc, zrfreq, nzfreq, TOT_ZEROFR);//
									sgn = (coef < 0) ? 1 : 0;
									if (coef < 0) coef = -coef;
									rc_encode(fp, enc->rc, pm->cumfreq[coef] - cumb,  pm->freq[coef],
									pm->cumfreq[pm->size] - cumb);//cumbを追加済み//累積分布と最大値からcumbを差し引く
									rc_encode(fp, enc->rc, sgn, 1, 2);
								}
							}
						}
					}
					enc->area_predictor_info[ar] += enc->rc->code - buff;
					buff = enc->rc->code;
				}
			}
			bits = (int)enc->rc->code;
			//temp = bits;
			//for(ar = 0; ar < enc->num_of_vgen; ar++) temp -= enc->area_predictor_info[ar];
			//for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_predictor_info[ar] += bits / enc->num_of_vgen; 
	    	enc->rc->code = 0;
    }
    return (bits);
}


int encode_predictor(FILE *fp, ENCODER *enc, int flag)
{
    int cl, co, ar, coef, sgn, k, m, min_m, bits, buff=0;
    cost_t cost, min_cost, t_cost;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    t_cost = 0.0;
	for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_predictor_info[ar] = 0;
    for (k = 0; k < enc->max_prd_order_all; k++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
	    	cost = 0.0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						coef = enc->predictor[cl][co][ar][k];
						if (coef < 0) coef = -coef;
		    			cost += enc->coef_cost[m][coef];
						//enc->area_predictor_info[ar] += enc->coef_cost[m][coef];
					}
				}
	    	}
			if (cost < min_cost) {
					min_cost = cost;
					min_m = m;
			}
		}
		t_cost += min_cost;
		if (flag) enc->coef_m[k] = min_m;
    }
    bits = t_cost;
    if (fp != NULL) {
		for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_predictor_info[ar] = 0;
		bits = 0;
		if (enc->f_huffman == 1) {      /* Huffman */
	    	for (k = 0; k < enc->max_prd_order_all; k++) {
				bits += putbits(fp, 4, enc->coef_m[k]);
				for(ar = 0; ar < enc->num_of_vgen; ar++){
					for (cl = 0; cl < enc->num_class[ar]; cl++) {
		    			for (co = 0; co < enc->num_kind_prd; co++) {
							coef = enc->predictor[cl][co][ar][k];
							sgn = (coef < 0)? 1 : 0;
							if (coef < 0) coef = -coef;
								bits += encode_golomb(fp, enc->coef_m[k], coef);
							if (coef != 0) {
								bits += putbits(fp, 1, sgn);
							}
						}
		    		}
				}
	    	}
		} else {               /* Arithmetic */
	    	PMODEL *pm;
	    	pm = &enc->spm;
	    	for (k = 0; k < enc->max_prd_order_all; k++) {
				set_spmodel(pm, enc->max_coef + 1, enc->coef_m[k]);
				rc_encode(fp, enc->rc, enc->coef_m[k], 1, 16);
				for(ar = 0; ar < enc->num_of_vgen; ar++){
					
					for (cl = 0; cl < enc->num_class[ar]; cl++) {
						for (co = 0; co < enc->num_kind_prd; co++) {
							coef = enc->predictor[cl][co][ar][k];
							sgn = (coef < 0)? 1 : 0;
							if (coef < 0) coef = -coef;
							rc_encode(fp, enc->rc, pm->cumfreq[coef],  pm->freq[coef],
							pm->cumfreq[pm->size]);
							if (coef > 0) {
								rc_encode(fp, enc->rc, sgn, 1, 2);
							}
						}
					}
					enc->area_predictor_info[ar] += enc->rc->code -buff;
					buff = enc->rc->code;
				}
			}
			bits = enc->rc->code;
			//temp = bits;
			//for(ar = 0; ar < enc->num_of_vgen; ar++) temp -= enc->area_predictor_info[ar];
			//for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_predictor_info[ar] += bits / enc->num_of_vgen; 
	    	enc->rc->code = 0;
		}
    }
    return (bits);
}


int encode_threshold(FILE *fp, ENCODER *enc, int flag)
{
    int cl, co, ar, gr, i, k, m, min_m, bits, buff=0;
    cost_t cost, min_cost;
    PMODEL *pm;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    if (enc->f_huffman == 1) {	/* Huffman */
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
	    	bits = 0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
	    		for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						k = 0;
		    			for (gr = 1; gr < enc->num_group; gr++) {
							i = enc->th[cl][co][ar][gr - 1] - k;
							bits++;
							if (i > 0) {
								bits += ((i - 1) >> m) + m + 1;
							}
							k += i;
							if (k > MAX_UPARA) break;
		    			}
					}
	  			}
			}
	    	if ((cost = bits) < min_cost) {
				min_cost = cost;
				min_m = m;
	    	}
		}
		if (fp == NULL) {
	    	enc->th_cost[0] = 1.0;
	    	for (i = 1; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] =((i - 1) >> min_m) + min_m + 1 + 1;
	    	}
	    	bits = min_cost;
		} else {
	    	bits = putbits(fp, 4, min_m);
			for(ar = 0; enc->num_of_vgen; ar++){
	    		for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						k = 0;
		    			for (gr = 1; gr < enc->num_group; gr++) {
							i = enc->th[cl][co][ar][gr - 1] - k;
							if (i == 0) {
			    				bits += putbits(fp, 1, 0);
							} else {
			    				bits += putbits(fp, 1, 1);
			    				bits += encode_golomb(fp, min_m, i - 1);
							}
							k += i;
							if (k > MAX_UPARA) break;
		    			}
					}
	    		}
			}
	    	if (enc->num_pmodel > 1) {
				for (k = 1; (1 << k) < enc->num_pmodel; k++);
					for (gr = 0; gr < enc->num_group; gr++) {
						pm = enc->pmlist[gr];
						bits += putbits(fp, k, pm->id);
					}
				}
			}
    } else {			/* Arithmetic */
		double p;
		for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_th_info[ar] = 0;
		pm = &enc->spm;
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
			set_spmodel(pm, MAX_UPARA + 2, m);
			cost = 0.0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						k = 0;
						for (gr = 1; gr < enc->num_group; gr++) {
							i = enc->th[cl][co][ar][gr - 1] - k;
							p = (double)pm->freq[i] / (pm->cumfreq[pm->size - k]);
							cost += -log(p);
							k += i;
							if (k > MAX_UPARA) break;
						}
					}
				}
			}
			cost /= log(2.0);
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		set_spmodel(pm, MAX_UPARA + 2, min_m);
		p = log(pm->cumfreq[MAX_UPARA + 2]);
		if (fp == NULL) {
			if (flag == 1){
				for (i = 0; i < MAX_UPARA + 2; i++) {
					enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);//enc->th_cost
				}
			}
			bits = min_cost;
		} else {
			buff = enc->rc->code;
			rc_encode(fp, enc->rc, min_m, 1, 16);
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				for (cl = 0; cl < enc->num_class[ar]; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						k = 0;
						for (gr = 1; gr < enc->num_group; gr++) {
							i = enc->th[cl][co][ar][gr - 1] - k;
							rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
							pm->cumfreq[pm->size - k]);
							k += i;
							if (k > MAX_UPARA) break;
						}
					}
				}
				enc->area_th_info[ar] += enc->rc->code - buff;
				buff = enc->rc->code;
			}
			if (enc->num_pmodel > 1) {//cの符号化かな？
				for (gr = 0; gr < enc->num_group; gr++) {
					pm = enc->pmlist[gr];
					rc_encode(fp, enc->rc, pm->id, 1, enc->num_pmodel);
				}
			}
			enc->area_th_info[0] += enc->rc->code - buff;//しょうがないので0に足す。
			bits = enc->rc->code;
			//temp = bits;
			//for(ar = 0; ar < enc->num_of_vgen; ar++) temp -= enc->area_th_info[ar];
			//for(ar = 0; ar < enc->num_of_vgen; ar++) enc->area_th_info[ar] += bits / enc->num_of_vgen; 
			enc->rc->code = 0;
		}
	}
    return (bits);
}

int encode_image(FILE *fp, ENCODER *enc)
{
    int x, y, e, co, ar, prd, base, bits, gr, cumbase, buff;
    PMODEL *pm;

    bits = buff = 0;

			/* Arithmetic */
			for (y = 0; y < enc->height; y++) {
			    for (x = 0; x < enc->width; x++) {
					gr = enc->group[y][x];
					co = color(y, x);
					ar = enc->area[y][x];
					prd = enc->prd[y][x];
					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;
					e = enc->encval[y][x];
					base = enc->bconv[prd];
					pm = enc->pmlist[gr] + enc->fconv[prd];
					cumbase = pm->cumfreq[base];
					rc_encode(fp, enc->rc,
					pm->cumfreq[base + e] - cumbase,
					pm->freq[base + e],
					pm->cumfreq[base + enc->maxval + 1] - cumbase);
					enc->area_co_err[ar][co] += enc->rc->code - buff;
					buff = enc->rc->code;
			    }
			}
			rc_finishenc(fp, enc->rc);
			
			enc->area_co_err[0][0] += enc->rc->code - buff;//しょうがないので00に足す。
			bits += enc->rc->code;

    return (bits);
}
inline cost_t opcl_vseg(ENCODER *enc, int k, int sub_cl, int ar, int restore,
	int **class, int ***gen_class)
{
    int x, y, i, j, l, m, **th_s, **prd_s;
    int **uq_s, *cl_p, *cls_p, *gcl_p, *gcls_p, **save_prd, **save_err;
		char *area_p;
		cost_t cost, min_cost;
		int xg, yg, flag=0;
		if(restore == 1) flag = 0;
		else flag = 1;
		min_cost = DBL_MAX;
    th_s = (int **)alloc_2d_array(enc->num_kind_prd, enc->num_group, sizeof(int));
    prd_s = (int **)alloc_2d_array(enc->num_kind_prd, enc->max_prd_order_all, sizeof(int));
    uq_s = (int **)alloc_2d_array(enc->num_kind_prd, (MAX_UPARA + 1), sizeof(int));
	save_prd = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	save_err = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	for (l = 0; l < enc->num_kind_prd; l++) {
		for (i = 0; i < enc->num_group; i++) {
			th_s[l][i] = 0;
		}
		for (i = 0; i < enc->max_prd_order_all; i++) {
			prd_s[l][i] = 0;
		}
		for (i = 0; i < (MAX_UPARA + 1); i++) {
			uq_s[l][i] = 0;
		}
    }
	for(y = 0; y < enc->height; y++){
		for(x = 0; x < enc->width; x++){
			save_prd[y][x]=enc->prd[y][x];
			save_err[y][x]=enc->err[y][x];
		}
	}
	l_optimize_class_vseg(enc, k, ar);
	enc->sub_cl = 0;
	for (l = 0; l < enc->num_kind_prd; l++) {
		for (i = 0; i < enc->num_group; i++) {
	    	th_s[l][i] = enc->th[k][l][ar][i];
		}
		for (i = 0; i < enc->max_prd_order_all; i++) {
	    	prd_s[l][i] = enc->predictor[k][l][ar][i];
		}
		for (i = 0; i < (MAX_UPARA + 1); i++) {
			uq_s[l][i] = enc->uquant[k][l][ar][i];
		}
    }
	enc->num_class[ar]--;
	for(m = k; m < enc->num_class[ar]; m++){
		for (l = 0; l < enc->num_kind_prd; l++) {
			for (j = 0; j < enc->num_group; j++) {
				enc->th[m][l][ar][j] = enc->th[m+1][l][ar][j];
			}
			for (j = 0; j < enc->max_prd_order_all; j++) {
				enc->predictor[m][l][ar][j] = enc->predictor[m+1][l][ar][j];
			}
			for (j = 0; j < MAX_UPARA + 1; j++) {
				enc->uquant[m][l][ar][j] = enc->uquant[m+1][l][ar][j];
			}
		}
	}
	enc->num_class[ar]++;
    for (y = 0; y < enc->height; y++) {
		cl_p = enc->class[y];
		area_p = enc->area[y];
		for (x = 0; x < enc->width; x++) {
			if(*area_p == ar){
			if(*cl_p > k){
					*cl_p = *cl_p - 1;
				}
			}
			cl_p++;
			area_p++;
		}
    }
	//if(ar == 0){
		for(yg = 0; yg < enc->num_of_geny[ar]; yg++){
			gcl_p = enc->gen_class[ar][yg];
			for(xg = 0; xg < enc->num_of_genx[ar]; xg++){
				if(*gcl_p >= 0){
					if(*gcl_p > k){
						*gcl_p = *gcl_p - 1;
					}
				}
				gcl_p++;
			}
		}
	//}else{
		/*for(yg = 0; yg < NUM_OF_GEN_ey; yg++){
			gcl_p = enc->gen_class[ar][yg];
			for(xg = 0; xg < NUM_OF_GEN_ex; xg++){
				if(*gcl_p >= 0){
					if(*gcl_p > k){
						*gcl_p = *gcl_p - 1;
					}
				}
				gcl_p++;
			}
		}*/
	//}
	enc->num_class[ar]--;
	set_prd_pels(enc);
	//predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);
	cost += encode_class_vseg(NULL, enc, NULL, flag);
	cost += encode_predictor2(NULL, enc, 1);
	cost += encode_threshold(NULL, enc, 0);
	min_cost = cost;
		
	if (restore == 1){
		enc->num_class[ar]++;
		for(y = 0; y < enc->height; y++){
			for(x = 0; x < enc->width; x++){
				enc->prd[y][x]=save_prd[y][x];
				enc->err[y][x]=save_err[y][x];
			}
		}
		for (y = 0; y < enc->height; y ++) {
			area_p = enc->area[y];
			cl_p = enc->class[y];
			cls_p = class[y];
			for (x = 0; x < enc->width; x ++) {
				if(*area_p == ar){
					if (*cl_p != *cls_p) {
						*cl_p = *cls_p;
					}
				}
				area_p++;
				cl_p++;
				cls_p++;
			}
			
		}
		for (yg = 0; yg < enc->num_of_geny[ar]; yg++) {
			gcl_p = enc->gen_class[ar][yg];
			gcls_p = gen_class[ar][yg];
			for (xg = 0; xg < enc->num_of_genx[ar]; xg++) {
				if(*gcl_p != *gcls_p){
					*gcl_p = *gcls_p;
				}
				gcl_p++;
				gcls_p++;
			}
		}
		


		for(m = (enc->num_class[ar] - 1); m > k; m--){
			for (l = 0; l < enc->num_kind_prd; l++) {
				for (j = 0; j < enc->num_group; j++) {
					enc->th[m][l][ar][j] = enc->th[m-1][l][ar][j];
				}
				for (j = 0; j < enc->max_prd_order_all; j++) {
					enc->predictor[m][l][ar][j] = enc->predictor[m-1][l][ar][j];
				}
				for (j = 0; j < MAX_UPARA + 1; j++) {
					enc->uquant[m][l][ar][j] = enc->uquant[m-1][l][ar][j];
				}
			}
		}
		
		for (l = 0; l < enc->num_kind_prd; l++) {
			for (i = 0; i < enc->num_group; i++) {
				enc->th[k][l][ar][i] = th_s[l][i];
			}
			for (i = 0; i < enc->max_prd_order_all; i++) {
				enc->predictor[k][l][ar][i] = prd_s[l][i];
			}
			for (i = 0; i < (MAX_UPARA + 1); i++) {
				enc->uquant[k][l][ar][i] = uq_s[l][i];
			}
		}
		set_prd_pels(enc);
		//predict_region(enc, 0, 0, enc->height, enc->width);
		calc_cost(enc, 0, 0, enc->height, enc->width);
	}
	free(save_prd);
	free(save_err);
    free(th_s);
    free(prd_s);
    free(uq_s);
    return(min_cost);
}

cost_t auto_del_class_vseg(ENCODER *enc, cost_t pre_cost)
{
    int x, y, i, j, k, del_cl, sub, ar, max, /*count,*/ r;
    cost_t cost, min_cost, pre_save;
    int **class, ***gen_class;
		int xg, yg, area;
		int **r_cl;
		int ndel_f[enc->num_of_vgen], n[enc->num_of_vgen];
		r_cl = (int **)alloc_2d_array(enc->num_of_vgen, enc->max_class ,sizeof(int));
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			ndel_f[ar] = 0;
			n[ar] = 0;
		}
		pre_save = pre_cost;
		
		restart:
		max = 0;
		sub = 0;
		printf("\n********auto_del_class_vseg*********\n");
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			n[ar] = enc->num_class[ar] / 5;
		}
		class = (int **)alloc_2d_array(enc->height, enc->width,sizeof(int));
		gen_class = (int ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(int));
		
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
			    class[y][x] = enc->class[y][x];
			}
		}
		
		for(area = 0; area < enc->num_of_vgen; area++){
			for(yg = 0; yg < enc->num_of_geny[area]; yg++){
				for(xg = 0; xg < enc->num_of_genx[area]; xg++){
					gen_class[area][yg][xg] = enc->gen_class[area][yg][xg];
				}
			}
		}
    
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			r = (int)(((double)rand() * enc->num_class[ar]) / ((RAND_MAX+1.0) * n[ar]));
			for(i = 0; i < n[ar]; i++){
				r_cl[ar][i] = r + (enc->num_class[ar] / n[ar]) * i;
			}
			if(r_cl[ar][n[ar] - 1] >= enc->num_class[ar]){
				r_cl[ar][n[ar] - 1] = enc->num_class[ar] - 1;
			}
			
			min_cost = DBL_MAX;
			del_cl = 0;
			j = 0;
			for (k = 0; k < enc->num_class[ar]; k++) {
				
					if(k != r_cl[ar][j]){
						continue;
					}else{
						j++;
					}
				
				cost = opcl_vseg(enc, k, enc->sub_cl, ar, 1, class, gen_class);
				if (cost < min_cost) {
			    min_cost = cost;
			    del_cl = k;
					sub = enc->sub_cl;
				}
			}
			if (pre_cost > min_cost) {
				pre_cost = min_cost;
				printf("ar = %d, del_cl = %d\n", ar, del_cl);
				cost = opcl_vseg(enc, del_cl, sub, ar, 0, class, gen_class);
				if(cost != min_cost){
					printf("auto_del_class ERROR!!\n");
					printf("now cost = %d, min_cost = %d\n", (int)cost, (int)min_cost);
					exit(1);
				}
			} else {
				ndel_f[ar]++;
			}
		}
		max = 0;
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			if(max < enc->num_class[ar]){
				max = enc->num_class[ar];

			}
		}
		enc->max_class = max;
		predict_region(enc, 0, 0, enc->height, enc->width);
		cost = calc_cost(enc, 0, 0, enc->height, enc->width)
				+ encode_class_vseg(NULL, enc, NULL, 1) 
				+ encode_predictor2(NULL, enc, 1)
				+ encode_threshold(NULL, enc, 1);
		free(class);
		free(gen_class);
		if(cost < pre_save) {
			pre_save = cost;
				goto restart;
		}else{
			i = 0;
			for(ar = 0; ar < enc->num_of_vgen; ar++){
				i += ndel_f[ar];
			}
			if(i == enc->num_of_vgen){
				goto restart;
			}
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			printf("num_class[%d] = %d\n",ar, enc->num_class[ar]);
		}
		free(r_cl);
    return(cost);
}

/*void correlation_coef_check(ENCODER *enc, int **edge_label)
{
	printf("correlation_coef_check\n");
	int y, x, z, k, i, j, s,t, ar, co;
	double temp;
	int height, width;
	height = enc->height;
	width = enc->width;
	int range_y, range_x, range, area;
	int *org_p;
	int *ref_count, max_ref_count = -INT_MAX;
	POINT **ref_coordinate;
	char *area_p;
	range_y = SEARCH_RANGE_Y;//40必ず偶数で記述
	range_x = SEARCH_RANGE_X;//70必ず偶数で記述
	area = enc->num_of_vgen;
	range = range_y*range_x;
	double ***sum, ***ave, ***sum0, ***ave0, ***variance_x, ***variance_y, ***covariance, ***coef;
	int ***count;
	sum = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	ave = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	sum0 = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	ave0 = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	variance_x = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	variance_y = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	covariance = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	coef = (double ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(double));
	count = (int ***)alloc_3d_array(NUM_KIND_PRD, area, range, sizeof(int));
	//ref_count = (int *)alloc_mem(area*sizeof(int));
	//printf("check\n");
	i = -range_y+1;
	j = -range_x/2;
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(k = 0; k < range; k++){
				sum[co][ar][k] = ave[co][ar][k] = sum0[co][ar][k] = ave0[co][ar][k] 
				= variance_y[co][ar][k] = variance_x[co][ar][k] = covariance[co][ar][k] 
				= coef[co][ar][k] = 0.0;
				count[co][ar][k] = 0;
			}
		}
	}

	//printf("check2\n");
	for(y = 0; y < height; y++){
		//org_p = enc->org[y];
		area_p = enc->area[y];
		for(x = 0; x < width; x++){
			
			ar = (int)area_p[x];
			if(ar > 0 && edge_label[y][x] > 0) {
				co = color(y, x);
				for(k = 0; k < range; k++){
					s = i + (int)(k / range_x);//y軸相対座標
					t = j + (int)(k % range_x);//x軸相対座標k　% range_x = - jの時t=0
					if(s+y >= 0 && s+y < height){
						if(t+x >= 0 && t+x < width){
							sum[co][ar][k] += enc->org[y+s][x+t];
							sum0[co][ar][k] += enc->org[y][x];
							count[co][ar][k]++;
						}				
					}
				}
			}
			//org_p++;
		}
	}
	//printf("check3\n");
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(k = 0; k < range; k++){
				ave[co][ar][k] = sum[co][ar][k] / (double)count[co][ar][k];
				ave0[co][ar][k] = sum0[co][ar][k] / (double)count[co][ar][k];
				//printf("k = %d, ar = %d, count = %d, ave = %lf, ave0 = %lf\n", k, ar, count[ar][k], ave[ar][k], ave0[ar][k]);
			}
		}
	}
	for(y = 0; y < height; y++){
		org_p = enc->org[y];
		area_p = enc->area[y];
		for(x = 0; x < width; x++){
			ar = area_p[x];
			if(ar > 0 && edge_label[y][x] > 0){
				co = color(y, x);
				for(k = 0; k < range; k++){
					s = i + (int)(k / range_x);//y軸相対座標
					t = j + (int)(k % range_x);//x軸相対座標k　% range_x = - jの時t=0
					if(s+y >= 0 && s+y < height){
						if(t+x >= 0 && t+x < width){
							covariance[co][ar][k] += ((double)org_p[0] - ave0[co][ar][k]) 
													* ((double)enc->org[y+s][x+t] - ave[co][ar][k]);  
							variance_x[co][ar][k] += ((double)org_p[0] - ave0[co][ar][k]) 
													* ((double)org_p[0] - ave0[co][ar][k]);
							variance_y[co][ar][k] += ((double)enc->org[y+s][x+t] - ave[co][ar][k]) 
													* ((double)enc->org[y+s][x+t] - ave[co][ar][k]);
						}
					}
				}
			}
			org_p++;
		}
	}
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(k = 0; k < range; k++){
				covariance[co][ar][k] /= (double)count[co][ar][k];
				variance_y[co][ar][k] /= (double)count[co][ar][k];
				variance_x[co][ar][k] /= (double)count[co][ar][k];
			}
		}
	}
	for(co = 0; co < NUM_KIND_PRD; co++){
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			//ref_count[ar] = 0;
			for(k = 0; k < range; k++){
				coef[co][ar][k] = covariance[co][ar][k] / (sqrt(variance_x[co][ar][k] * variance_y[co][ar][k]));
				//temp = coef[co][ar][k] * 50;
				//coef[co][ar][k] = round(temp) / 50;
				//printf("coef[%d][%4d] = %lf, covariance[%d][%4d] = %lf\n", ar, k, coef[ar][k], ar, k, covariance[ar][k]);
				if(k >= range-range_x/2) coef[co][ar][k] = 0.0;
				
				if(coef[co][ar][k] > 1 || coef[co][ar][k] < -1) printf("yabai\n");
			}
		}
	}
	

	enc->correlation_coef = coef;
	printf("check file out\n");
	FILE *fp;
	char *name = "./LOG/correlation_coef_check_edgeoff.csv";
	fp = fileopen(name, "wb");
	
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		fprintf(fp, "area = %d\n", ar);
		for(co = 0; co < NUM_KIND_PRD; co++){
			fprintf(fp, "co = %d\n", co);
			for(k = 0; k < range; k++){
				fprintf(fp, "%lf, ", coef[co][ar][k]);
				if((k % range_x) == range_x-1) fprintf(fp, "\n");
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	free(sum);
	free(ave);
	free(sum0);
	free(ave0);
	free(variance_x);
	free(variance_y);
	free(covariance);
	//free(coef);
	free(count);
	//free(ref_count);
	//free(ref_coordinate);
}*/

void edge_label_analysis(ENCODER *enc, int **edge_label)
{
	int y, x, ar, hist[256]={0}, height, width, max_val=256;
	int s, hist_area[2][256];
	int count_hist=0,count_hist_area[2]={0};
	height = enc->height;
	width = enc->width;
	for(ar = 0; ar < 2; ar++){
		for(s = 0; s < max_val; s++){
			hist_area[ar][s]=0;
		}
	}
	
	
	FILE *fp;
	char *name = "./LOG/edge_hist.csv";
	fp = fileopen(name, "wb");
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			hist[enc->org[y][x]]++;
			count_hist++;
		}
	}
	for(s = 0; s < max_val; s++) {
		fprintf(fp, "%d", s);
		if(s != max_val-1) fprintf(fp, ",");
		else fprintf(fp, "\n");
	}
	fprintf(fp, "image_all\n");
	for(s = 0; s < max_val; s++) {
		fprintf(fp, "%f", (double)hist[s]/(double)count_hist);
		if(s != max_val-1) fprintf(fp, ",");
		else fprintf(fp, "\n");
		hist[s]=0;
	}
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			if(enc->area[y][x] == 0) ar = 0;
			else ar =1;
			hist_area[ar][enc->org[y][x]]++;
			count_hist_area[ar]++;
		}
	}
	fprintf(fp, "each area\n");
	for(ar = 0; ar < 2; ar++){
		fprintf(fp, "%d\n", ar);
		for(s = 0; s < max_val; s++) {
			fprintf(fp, "%f", (double)hist_area[ar][s]/(double)count_hist_area[ar]);
			if(s != max_val-1) fprintf(fp, ",");
			else fprintf(fp, "\n");
			hist_area[ar][s]=0;
		}
		count_hist_area[ar]=0;
	}
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			if(edge_label[y][x] > 0) ar = 1;
			else ar = 0;
				hist_area[ar][enc->org[y][x]]++;
				count_hist_area[ar]++;
		}
	}
	fprintf(fp, "each area and edge\n");
	for(ar = 0; ar < 2; ar++){
		fprintf(fp, "%d\n", ar);
		for(s = 0; s < max_val; s++) {
			fprintf(fp, "%f", (double)hist_area[ar][s]/(double)count_hist_area[ar]);
			if(s != max_val-1) fprintf(fp, ",");
			else fprintf(fp, "\n");
		}
	}
	fclose(fp);
	name ="./edge_labelout.pgm";
	fp = fileopen(name, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	for (y = 0; y < height; y++){
		for (x = 0; x < width; x++){
			if(edge_label[y][x] == 0){
				putc(0,fp);
			}else{
				putc(enc->org[y][x], fp);
			}
		}
    }
	fclose(fp);
	
}

/*void roff_check(ENCODER *enc)
{	
	printf("roff_check -> ");
	int y, x, k, roff, roff2;
	for(y = 0; y < enc->height; y++){
		for(x = 0; x < enc->width; x++){
			for(k = 0; k < enc->max_prd_order_all; k++){
				roff = ref_offset4(enc->height, enc->width, y, x, k, color(y, x));
				roff2 = ref_offset5_eachcolor(enc->height, enc->width, y, x, k);
				if(roff2 >= 0){
					printf("ERROR y = %d, x = %d, k = %d, ref_offset5_eachcolor = %d, ref_offset4 = %d\n", y, x, k, roff2, roff);
					//exit(1);
				}
			}
		}
	}
	printf("ok\n");
	//printf("roff = %d\n", ref_offset5_eachcolor(enc->height, enc->width, 54, 5, 124));
}*/

void coding_rate_par_area(ENCODER *enc)
{
	int co, ar, y, x, area_err[enc->num_of_vgen];
	int total_pel=0;
	double bits=0.0, coding_length[enc->num_of_vgen], coding_rate[enc->num_of_vgen];
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		coding_length[ar] = 0.0;
		coding_rate[ar] = 0.0;
		area_err[ar] = 0;
	}
	double err=0.0, class=0.0, predictor=0.0, th=0.0;
	for(y = 0; y < enc->height; y++){
		for(x = 0; x < enc->width; x++){
			enc->area_pel[(int)enc->area[y][x]]++;
		}
	}
	printf("coding_rate_par_area(without header info.)\n");
	total_pel = enc->height * enc->width;
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		for(co = 0; co < NUM_KIND_PRD; co++){
			area_err[ar] += enc->area_co_err[ar][co];
		}
		coding_length[ar] = (double)area_err[ar] + enc->area_class_info[ar]
							+ enc->area_predictor_info[ar] + enc->area_th_info[ar];
		coding_rate[ar] = coding_length[ar] / (double)enc->area_pel[ar];
		bits += coding_length[ar];
		err += area_err[ar];
		class += enc->area_class_info[ar];
		predictor += enc->area_predictor_info[ar];
		th += enc->area_th_info[ar];
		printf("--------------------area = %d--------------------\n",ar);
		printf("pixel\t:\t%10d pixel (%4.2f %%)\n", (int)enc->area_pel[ar], (double)enc->area_pel[ar]/(double)total_pel*100);
		printf("class info.\t:\t%10d bits (%4.2f %%)\n", (int)enc->area_class_info[ar], enc->area_class_info[ar]/coding_length[ar]*100);
		printf("predictor\t:\t%10d bits (%4.2f %%)\n", (int)enc->area_predictor_info[ar], enc->area_predictor_info[ar]/coding_length[ar]*100);
		printf("thresholds\t:\t%10d bits (%4.2f %%)\n",(int)enc->area_th_info[ar], enc->area_th_info[ar]/coding_length[ar]*100);	
		printf("G1\t\t pred. errors \t:%10d bits\n", enc->area_co_err[ar][0]);
		printf("G2\t\t pred. errors \t:%10d bits\n", enc->area_co_err[ar][3]);
		printf("R\t\t pred. errors \t:%10d bits\n", enc->area_co_err[ar][2]);
		printf("B\t\t pred. errors \t:%10d bits\n", enc->area_co_err[ar][1]);
		printf("total\t\t pred. errors \t:%10d bits (%4.2f %%)\n", area_err[ar], (double)area_err[ar]/coding_length[ar]*100);
		printf("coding rate\t%10.5f bits/pel\n", coding_rate[ar]);
	}
	printf("\n");
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		printf("num_class[%d]=%d\n",ar, enc->num_class[ar]);
	}
	/*printf("bits = %f\n", bits);
	printf("err = %f\n", err);
	printf("class = %f\n", class);
	printf("predictor = %f\n", predictor);
	printf("th = %f\n", th);*/
}

void print_coding_condition(ENCODER *enc)
{
	int ar, y, x;
	printf("BAYER_cfa\n");
	printf("color(y,x) mean ... \n1->B, 2->R, 0->G1, 3->G2\n");
	for(y = 0; y < 2; y++){
		for(x = 0; x < 2; x++){
			printf("%d,", color(y,x));
		}
		printf("\n");
	}
	printf("max_prd_order_all = %d\n", MAX_PRD_ORDER_ALL); 
	printf("coef_precision = %d\nnum_pmodel = %d\npm_accuracy = %d\n",COEF_PRECISION, NUM_PMODEL, PM_ACCURACY);
	printf("MAX_LOOP_NUM = %d\n", MAX_ITERATION);
	printf("\n-----------------------------------------------------------------\n");
	printf("\t|\t\tarea");
	printf("\n\t---------------------------------------------------------\n");
	printf("\t|");
	for(ar = 0; ar < enc->num_of_vgen; ar++) {
		printf("%3d\t|", ar);
	}
	printf("\n-----------------------------------------------------------------\n");
	printf("class\t|");
	for(ar = 0; ar < enc->num_of_vgen; ar++) {
		printf("%3d\t|", enc->num_class[ar]);
	}
	printf("\n-----------------------------------------------------------------\n");
	printf("base_K\t|");
	for(ar = 0; ar < enc->num_of_vgen; ar++) {
		printf("%3d\t|", enc->base_prd_order[ar]);
	}
	printf("\n-----------------------------------------------------------------\n");
	printf("max_K\t|");
	for(ar = 0; ar < enc->num_of_vgen; ar++) {
		printf("%3d\t|", enc->max_prd_order[ar]);
	}
	printf("\n-----------------------------------------------------------------\n");
}

int main(int argc, char **argv)
{
    cost_t cost, min_cost, side_cost, sc;
	cost = side_cost = sc = 0.0;
    int i, j, k, x, y, yg, xg,/*z, yg, xg,*/ cl, co, ar, /*area,*/ bits, ****prd_save, ****th_save;
		/*int xt, yt;*/
	int gr; 
	PMODEL **pmlist_save;
    int prd_cost, class_cost, th_cost, err_cost, qt_cost;
	int height, width;
	//POINT_REF *d_p[NUM_OF_VGEN];
    int **class_save, ***gen_class_save, save_max_class=MAX_CLASS, *save_num_class;
    double rate;
    IMAGE *img, *img2;
    ENCODER *enc;
    double elapse = 0.0;
    int f_mmse = 0;
    int f_optpred = 0;
    int f_huffman = 0;
    int quadtree_depth = QUADTREE_DEPTH;
	quadtree_depth = -1;
    int num_class=MAX_CLASS;
    int num_group = NUM_GROUP;
    //int prd_order = PRD_ORDER;
	int prd_order = BASE_PRD_ORDER;
    int coef_precision = COEF_PRECISION;
    int num_pmodel = NUM_PMODEL;
    int pm_accuracy = PM_ACCURACY;
    int max_iteration = MAX_ITERATION;
    int num_kind_prd = NUM_KIND_PRD;
	int **buf;
    char *infile, *outfile;
    FILE *fp;
    cpu_time();
    setbuf(stdout, 0);
    infile = outfile = NULL;
    for (i = 1; i < argc; i++) {
			if (argv[i][0] == '-') {
	    	switch (argv[i][1]) {
		case 'M':
		    num_class = atoi(argv[++i]);
		    if (num_class <= 0 || num_class > MAX_CLASS) {
					num_class = MAX_CLASS;
		    }
		    break;
		case 'K':
		    prd_order = atoi(argv[++i]);
		    	if (prd_order <= 0 || prd_order > 126) {
						prd_order = BASE_PRD_ORDER;
		    	}
		    	break;
		case 'P':
		    coef_precision = atoi(argv[++i]);
		    if (coef_precision <= 0 || coef_precision > 16) {
					coef_precision = COEF_PRECISION;
		    }
		    break;
		case 'V':
		    num_pmodel = atoi(argv[++i]);
		    if (num_pmodel <= 0 || num_pmodel > 64) {
					num_pmodel = NUM_PMODEL;
		    }
		    break;
		case 'A':
		    pm_accuracy = atoi(argv[++i]);
		    if (pm_accuracy < -1 || pm_accuracy > 6) {
					pm_accuracy = PM_ACCURACY;
		    }
		    break;
		case 'I':
		    max_iteration = atoi(argv[++i]);
		    if (max_iteration <= 0) {
					max_iteration = MAX_ITERATION;
		    }
		    break;
		case 'm':
		    f_mmse = 1;
		    break;
		case 'o':
		    f_optpred = 1;
		    break;
		case 'h':
		    f_huffman = 1;
		    break;
		case 'f':
		    quadtree_depth = -1;
		    break;
		default:
		    fprintf(stderr, "Unknown option: %s!\n", argv[i]);
		    exit (1);
	  		}
			} else {
	    	if (infile == NULL) {
					infile = argv[i];
	    	} else {
					outfile = argv[i];
	    	}
			}
    }
    if (f_huffman == 1) pm_accuracy = -1;
    if (pm_accuracy > coef_precision) pm_accuracy = coef_precision;
    if (infile == NULL || outfile == NULL) {
			printf(BANNER"\n", 0.1 * VERSION);
			printf("usage: encmrp [options] infile outfile\n");
			printf("options:\n");
			printf("    -M num  Number of predictors [%d]\n", num_class);
			printf("    -K num  Prediction order [%d]\n", prd_order);
			printf("    -P num  Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
			printf("    -V num  Number of probability models [%d]\n", num_pmodel);
			printf("    -A num  Accuracy of probability models [%d]\n", pm_accuracy);
			printf("    -I num  Maximum number of iterations [%d]\n", max_iteration);
			printf("    -m      Use MMSE predictors\n");
			printf("    -h      Use Huffman coding\n");
			printf("    -f      Fixed block-size for adaptive prediction\n");
			printf("    -o      Further optimization of predictors (experimental)\n");
			printf("infile:     Input file (must be in a raw PGM format)\n");
			printf("outfile:    Output file\n");
			exit(0);
    }
    img = read_pgm(infile);
	height = img->height;
	width = img->width;
	img2 = read_pgm("./outfile_edge_label.pgm");
	int **edge_label;
	edge_label = (int **)alloc_2d_array(height, width, sizeof(int));
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			edge_label[y][x] = (int)img2->val[y][x];
		}
	}
	free(img2->val);
	free(img2);
    fp = fileopen(outfile, "wb");
    printf("%s -> %s (%dx%d)\n", infile, outfile, img->width, img->height);
	enc = init_encoder(img, num_group, prd_order, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, num_kind_prd);
	print_coding_condition(enc);
	enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel, enc->pm_accuracy, NULL, enc->sigma, enc->maxval + 1);
    if (enc->f_huffman == 1) enc->vlcs = init_vlcs(enc->pmodels, enc->num_group, enc->num_pmodel);
	
    set_cost_model(enc, f_mmse);
    init_class_vseg(enc); //vseg
	//correlation_coef_check(enc, edge_label);
	/*edge_label_analysis(enc, edge_label);
	free(edge_label);
	exit(1);*/
	//make_auto_ref_offset(enc);
	//for(ar = 0; ar < enc->num_of_vgen; ar++) d_p[ar] = dyx[ar];
	//init_ref_offset5_eachcolor(d_pp, enc->max_prd_order, enc->height, enc->width, enc);
	//init_ref_offset5(d_p, enc->max_prd_order, enc->height, enc->width, enc);
	init_ref_offset5_2(enc->max_prd_order_all, height, width, enc);
	//char *outfil = "optimize_class_var";//init_class_vseg
    buf = (int **)alloc_2d_array(enc->num_of_vgen, MAX_CLASS, sizeof(int));
	for(ar = 0; ar < enc->num_of_vgen; ar++){
		for (cl = 0; cl < MAX_CLASS; cl++) {
			buf[ar][cl] = 0;
		}
    }
	save_num_class = (int *)alloc_mem(enc->num_of_vgen * sizeof(int));
	prd_save = (int ****)alloc_4d_array(enc->max_class, enc->num_kind_prd, enc->num_of_vgen, enc->max_prd_order_all, sizeof(int));
    th_save = (int ****)alloc_4d_array(enc->max_class, enc->num_kind_prd, enc->num_of_vgen, enc->num_group, sizeof(int));
    class_save = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
    gen_class_save = (int ***)alloc_3d_array(enc->num_of_vgen, enc->max_geny, enc->max_genx, sizeof(int));
	pmlist_save = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));
	/********************** 1st loop ****************************************/
	//	roff_check(enc);
	
	printf("[ 1st loop ]\n");
    enc->optimize_loop = 1;
    min_cost = DBL_MAX;
    for (i = j = 0; i < max_iteration; i++) {
		printf("[%2d] cost =", i);
		cost = design_predictor(enc, f_mmse);
		printf(" %d ->", (int)cost);
		cost = optimize_group(enc);
		printf(" %d ->", (int)cost);
		cost = optimize_class_vseg(enc);
		printf(" %d", (int)cost);
		if (cost < min_cost) {
		printf(" *\n");
		min_cost = cost;
		j = i;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				class_save[y][x] = enc->class[y][x];
			}
		}
		for (cl = 0; cl < enc->max_class; cl++) {
			for (co = 0; co <enc->num_kind_prd; co++) {
				for(ar = 0; ar < enc->num_of_vgen; ar++) {
					for (k = 0; k < enc->max_prd_order_all; k++) {
						prd_save[cl][co][ar][k] = enc->predictor[cl][co][ar][k];
					}
						for (k= 0; k < enc->num_group; k++) {
							th_save[cl][co][ar][k] = enc->th[cl][co][ar][k];
						}
					}
				}
			}
		} else {
			printf("\n");
		}
		if (i - j >= EXTRA_ITERATION) break;
		elapse += cpu_time();
  	}
  	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
	    	enc->class[y][x] = class_save[y][x];
		}
  	}
  	for (cl = 0; cl < enc->max_class; cl++) {
		for (co = 0; co < enc->num_kind_prd; co++) {
			for(ar = 0; ar < enc->num_of_vgen; ar++) {
				for (k = 0; k < enc->max_prd_order_all; k++) {
					enc->predictor[cl][co][ar][k] = prd_save[cl][co][ar][k];
	    		}
				i = 0;
				for (k = 0; k < enc->num_group; k++) {
					enc->th[cl][co][ar][k] = th_save[cl][co][ar][k];
					for (; i < enc->th[cl][co][ar][k]; i++) {
						enc->uquant[cl][co][ar][i] = k;
					}
				}
			}
		}
  	}
    set_cost_rate(enc);
    set_prd_pels(enc);
    predict_region(enc, 0, 0, height, width);
    cost = calc_cost(enc, 0, 0, height, width);
    /************************** 2nd loop ***********************************/
	printf("[ 2nd loop ]\n");
	enc->optimize_loop = 2;
	min_cost = INT_MAX;
	for (i = j = 0; i < max_iteration; i++) {
		printf("(%2d) cost =", i);
		if (f_optpred) {
			cost = optimize_predictor2(enc);
			printf(" %d", (int)cost);
		}
		side_cost = sc = (int)encode_predictor2(NULL, enc, 1);
		printf("[%.0f] ->", (double)sc);
		cost = optimize_group(enc);
		side_cost += sc = encode_threshold(NULL, enc, 1);
		printf(" %d [%d] ->", (int)cost, (int)sc);
		cost = optimize_class_vseg(enc);
		side_cost += sc = encode_class_vseg(NULL, enc, NULL, 1);
		printf(" %d [%d] (%.0f)", (int)cost, (int)sc, (double)side_cost);
		cost += side_cost;
		printf(" -> %d", (int)cost);
		//cost = auto_del_class_vseg(enc, cost);
		//printf(" -> %d", (int)cost);
		if (cost < min_cost) {
			//if (i==0) {
			printf(" *\n");
			min_cost = cost;
			j = i;
			if (f_optpred) {
				save_max_class = enc->max_class;
				for(ar = 0; ar < enc->num_of_vgen; ar++) save_num_class[ar] = enc->num_class[ar];
				for (y = 0; y < height; y++) {
					for (x = 0; x < width; x++) {
						class_save[y][x] = enc->class[y][x];
					}
				}
				for (cl = 0; cl < enc->max_class; cl++) {
					for (co = 0; co < enc->num_kind_prd; co++) {
						for(ar = 0; ar < enc->num_of_vgen; ar++){
							for (k= 0; k < enc->max_prd_order_all; k++) {
								prd_save[cl][co][ar][k] = enc->predictor[cl][co][ar][k];
							}
							for (k= 0; k < enc->num_group; k++) {
								th_save[cl][co][ar][k] = enc->th[cl][co][ar][k];
							}
						}
					}
				}
				for(ar = 0; ar < enc->num_of_vgen; ar++){
					for(yg = 0; yg < enc->num_of_geny[ar];yg++){
						for(xg = 0; xg < enc->num_of_genx[ar]; xg++){
							gen_class_save[ar][yg][xg] = enc->gen_class[ar][yg][xg];
						}
					}
				}
				for(gr = 0; gr < enc->num_group; gr++){
					pmlist_save[gr] = enc->pmlist[gr];
				}
			}
		} else {
			cost = auto_del_class_vseg(enc, cost);
			printf(" -> %d", (int)cost);
			if (cost < min_cost) {
				printf(" *\n");
				min_cost = cost;
				if (f_optpred) {
					save_max_class = enc->max_class;
					for(ar = 0; ar < enc->num_of_vgen; ar++) save_num_class[ar] = enc->num_class[ar];
					for (y = 0; y < height; y++) {
						for (x = 0; x < width; x++) {
							class_save[y][x] = enc->class[y][x];
						}
					}
					for (cl = 0; cl < enc->max_class; cl++) {
						for (co = 0; co < enc->num_kind_prd; co++) {
							for(ar = 0; ar < enc->num_of_vgen; ar++){
								for (k= 0; k < enc->max_prd_order_all; k++) {
									prd_save[cl][co][ar][k] = enc->predictor[cl][co][ar][k];
								}
								for (k= 0; k < enc->num_group; k++) {
									th_save[cl][co][ar][k] = enc->th[cl][co][ar][k];
								}
							}
						}
					}
					for(ar = 0; ar < enc->num_of_vgen; ar++){
						for(yg = 0; yg < enc->num_of_geny[ar];yg++){
							for(xg = 0; xg < enc->num_of_genx[ar]; xg++){
								gen_class_save[ar][yg][xg] = enc->gen_class[ar][yg][xg];
							}
						}
					}
					for(gr = 0; gr < enc->num_group; gr++){
						pmlist_save[gr] = enc->pmlist[gr];
					}
				}
			}else printf("\n");
		}
		if (f_optpred) {
			if (i - j >= EXTRA_ITERATION) break;
		} else {
			if (i > j) break;
		}
		elapse += cpu_time();
	}
	if (f_optpred) {
		enc->max_class = save_max_class;
		for(ar = 0; ar < enc->num_of_vgen; ar++) enc->num_class[ar] = save_num_class[ar];
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				enc->class[y][x] = class_save[y][x];
			}
		}
		for (cl = 0; cl < enc->max_class; cl++) {
			for (co = 0; co < enc->num_kind_prd; co++) {
				for(ar = 0; ar < enc->num_of_vgen; ar++){
					for (k= 0; k < enc->max_prd_order_all; k++) {
						enc->predictor[cl][co][ar][k] = prd_save[cl][co][ar][k];
					}
					i = 0;
					for (k= 0; k < enc->num_group; k++) {
						enc->th[cl][co][ar][k] = th_save[cl][co][ar][k];
						for (; i < enc->th[cl][co][ar][k]; i++) {
							enc->uquant[cl][co][ar][i] = k;
						}
					}
				}
			}
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(yg = 0; yg < enc->num_of_geny[ar];yg++){
				for(xg = 0; xg < enc->num_of_genx[ar]; xg++){
					enc->gen_class[ar][yg][xg] = gen_class_save[ar][yg][xg];
				}
			}
		}
		for(gr = 0; gr < enc->num_group; gr++){
			enc->pmlist[gr] = pmlist_save[gr];
		}
		set_prd_pels(enc);
		predict_region(enc, 0, 0, enc->height, enc->width);
		calc_cost(enc, 0, 0, enc->height, enc->width);
	}
  	remove_emptyclass_vseg(enc);
/****************************************************/
	char *name;//vseg
	FILE *fl;//vseg
		name = "./LOG/class&index_hist.csv";
		fl = fileopen(name, "wb");
		fprintf(fl, "num_segment\n");
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			fprintf(fl, "%d", ar);
			if(ar == enc->num_of_vgen-1) fprintf(fl, "\n");
			else fprintf(fl, ",");
		}
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			fprintf(fl, "%d", enc->num_segment[ar]);
			if(ar == enc->num_of_vgen-1) fprintf(fl, "\n");
			else fprintf(fl, ",");
		}

		fprintf(fl, "area,class,hist\n");
		
		for(i = 0; i < enc->num_of_vgen; i++){
			for (cl = 0; cl < enc->num_class[i]; cl++) {
				enc->mtfbuf[cl] = 0;
			}
			for(y = 0; y < enc->num_of_geny[i]; y++){
				for(x = 0; x < enc->num_of_genx[i]; x++){//mtf_buf
					cl = enc->gen_class[i][y][x];
					if(cl < 0)continue;
					enc->mtfbuf[cl]++;
				}
			}
			for(cl = 0; cl < enc->num_class[i]; cl++){
				fprintf(fl, "%d,%d,%d\n", i, cl, enc->mtfbuf[cl]);
			}
		}
		fprintf(fl, "area,index,hist\n");
		int hist[enc->max_class];

		for(i = 0; i < enc->num_of_vgen; i++){
			for (cl = 0; cl < enc->num_class[i]; cl++) {
				enc->mtfbuf[cl] = cl;
				hist[cl] = 0;
			}
			for(y = 0; y < enc->num_of_geny[i]; y++){
				for(x = 0; x < enc->num_of_genx[i]; x++){//mtf_buf
					cl = enc->gen_class[i][y][x];
					if(cl < 0)continue;
					mtf_classlabel_vseg(enc->max_genx, enc->gen_class, enc->mtfbuf, x, y, i, enc->num_class);
					j = enc->mtfbuf[cl];
					hist[j]++;
				}
			}
			for(cl = 0; cl < enc->num_class[i]; cl++){
				fprintf(fl, "%d,%d,%d\n", i, cl, hist[cl]);
			}
		}
		fclose(fl);
/**************************************************/


  	qt_cost = 0;
	printf("\n");
  	printf("******** Start encoding ********\n");
  	bits = k = write_header(enc, fp);
  	if (enc->f_huffman == 0) {
		enc->rc = rc_init();
		putbits(fp, 7, 0);	/* byte alignment for the rangecoder */
  	}
	printf("encode_class...");
    bits += class_cost = encode_class_vseg(fp, enc, &qt_cost, 1);
	printf("ok\n");
	printf("encode_predictor...");
	bits += prd_cost= encode_predictor2(fp, enc, 1);
	printf("ok\n");
	printf("encode_threshold...");
	bits += th_cost = encode_threshold(fp, enc, 1);
	printf("ok\n");
	printf("encode_image...");
	bits += err_cost = encode_image(fp, enc);
	printf("ok\n");
	fclose(fp);
	rate = (double)bits / (enc->height * enc->width);
	printf("header info.\t:%10d bits (%4.2f %%)\n", k, (double)k/(double)bits*100);
    printf("class info.\t:%10d bits (%4.2f %%)\n", class_cost - qt_cost, (double)(class_cost - qt_cost)/(double)bits*100);
    if (enc->quadtree_depth > 0) printf("quadtree info.\t:%10d bits\n", qt_cost);
    printf("predictors\t:%10d bits (%4.2f %%)\n", prd_cost, (double)prd_cost/(double)bits*100);    
    printf("thresholds\t:%10d bits (%4.2f %%)\n", th_cost, (double)th_cost/(double)bits*100);
    printf("pred. errors\t:%10d bits (%4.2f %%)\n", err_cost, (double)err_cost/(double)bits*100);
    printf("------------------------------\n");
    printf("total\t\t:%10d bits\n", bits);
    printf("coding rate\t:%10.5f b/p\n", rate);
    elapse += cpu_time();
    printf("cpu time: %.2f sec.\n", elapse);
	coding_rate_par_area(enc);
	calc_ratio_of_model_to_rate(enc);
	print_error(enc->err, enc->height, enc->width, outfile);
	print_prd(enc->prd, enc->height, enc->width, outfile);
    //print_predictor(enc->predictor, enc->max_prd_order_all, enc->num_class,
                    //enc->num_kind_prd, enc->max_coef, outfile);
	for(co = 0; co < enc->num_kind_prd; co++){
		for(ar = 0; ar < enc->num_of_vgen; ar++){
			for(cl = 0; cl < enc->num_class[ar]; cl++){
				for(k = 0; k < enc->max_prd_order[ar]; k++){
					enc->predict_out[co][k] += abs(enc->predictor[cl][co][ar][k]);
				}
			}
		}
	}
	//print_predict_out(enc->predict_out, enc->max_prd_order_all, enc->num_kind_prd, outfile);
	//print_predict_outall(enc->predict_out, enc->max_prd_order_all, enc->num_kind_prd, outfile);
	print_threshold(enc->num_of_vgen, enc->th, enc->num_group, enc->num_class,
		    						enc->num_kind_prd, enc->pmlist, NULL, outfile);
	print_class(enc->class, enc->height, enc->width, outfile);
    print_block_size(img->val, enc->height, enc->width, enc->quadtree_depth,
		    							enc->qtmap, outfile);
    //for(ar = 0; ar < enc->num_of_vgen; ar++){
    	//printf("class[%d] = %d\n", ar, enc->num_class[ar]);
	//}
	//printf("enc->max_class = %d\n", enc->max_class);
    return (0);
}
