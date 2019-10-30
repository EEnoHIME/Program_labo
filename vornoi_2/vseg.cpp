#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "mrp1.h"
#include <float.h>
#include "vseg.h"


FILE *fileopen(char *filename, char *mode)
{
    FILE *fp;
    fp = fopen(filename, mode);
    if (fp == NULL) {
        fprintf(stderr, "Can\'t open %s!\n", filename);
        exit(1);
    }
    return (fp);
}



void *alloc_mem(size_t size)
{
    void *ptr;
    if ((ptr = (void *)malloc(size)) == NULL) {
        fprintf(stderr, "Can\'t allocate memory (size = %u)!\n", (int)size);//lfc
        exit(1);
    }
    return (ptr);
}

void **alloc_2d_array(int height, int width, size_t size)
{
    void **mat;
    void *ptr;
    int k;

    mat = (void **)alloc_mem(sizeof(*mat) * height + height * width * size);
    ptr = (void *)(mat + height);
    for (k = 0; k < height; k++) {
			mat[k] =  ptr;
			ptr += width * size;
    }
    return (mat);
}

void ***alloc_3d_array(int height, int width, int depth, size_t size)
{
		void ***mat, **ptr2;
		void *ptr;
		int t, k, l;
		t = size * depth;

		mat = (void ***)alloc_mem((sizeof(*mat) + sizeof(**mat) * width + t * width) * height);
		ptr2 = (void **)(mat + height);
		ptr = (void *)(ptr2 + height * width);

		for(k = 0; k < height; k++){
			mat[k] = ptr2;
			for(l = 0; l < width; l++){
				ptr2[l] = ptr;
				ptr += t;
			}
			ptr2 += width;
		}
		return (mat);
}

IMAGE *alloc_image(int width, int height, int maxval)
{
		int y,x;
    IMAGE *img;
    img = (IMAGE *)alloc_mem(sizeof(IMAGE));
    img->width = width;
    img->height = height;
    img->maxval = maxval;
    img->val = (img_t **)alloc_2d_array(img->height, img->width,
                                        sizeof(img_t));
		for(y = 0; y < img->height; y++){
			for(x = 0; x < img->width; x++){
				img->val[y][x] = 0;
			}
		}
    return (img);
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


double cpu_time(void)
{
#include <time.h>
#ifndef HAVE_CLOCK
#  include <sys/times.h>
    struct tms t;
#endif
#ifndef CLK_TCK
#  define CLK_TCK 60
#endif
    static clock_t prev = 0;
    clock_t cur, dif;

#ifdef HAVE_CLOCK
    cur = clock();
#else
    times(&t);
    cur = t.tms_utime + t.tms_stime;
#endif
    if (cur > prev) {
	dif = cur - prev;
    } else {
	dif = (unsigned)cur - prev;
    }
    prev = cur;

#ifdef HAVE_CLOCK
    return ((double)dif / CLOCKS_PER_SEC);
#else
    return ((double)dif / CLK_TCK);
#endif
}

void print_class(char **class, int height, int width)
{
    int i, j;
    char name[100];
    FILE *fp;

    sprintf(name, "outfile_class.pgm");
    fp = fileopen(name, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
            putc(class[i][j]*4, fp);
	     }
    }
    fclose(fp);
    return;
}

void print_label(LABEL *label, int height, int width, int *num_of_vgenc,
int *num_of_vgen, int *num_of_geny, int *num_of_genx)
{
	printf("***print_label***\n");
    int i, j, ar;
    char name[100];
    FILE *fp;

    sprintf(name, "outfile_label.txt");
    fp = fileopen(name, "wb");
	fprintf(fp, "%d\n",*num_of_vgen);
	fprintf(fp, "%d\n",*num_of_vgenc);
	for(ar = 0; ar < *num_of_vgen; ar++){
		fprintf(fp, "%d\n",num_of_geny[ar]);
		fprintf(fp, "%d\n",num_of_genx[ar]);
	}
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
         fprintf(fp, "%d\n",(int)(label->value[i][j].x));
         fprintf(fp, "%d\n",(int)(label->value[i][j].y));
         fprintf(fp, "%d\n",(int)(label->value[i][j].area));
	     }
    }
    fclose(fp);
}

void print_range(POINT ***gen_tl, POINT ***gen_br, 
int *num_of_vgen, int *num_of_geny, int *num_of_genx, char *filename)
{
	
	printf("***print_range***\n");
    FILE *fp;
    int xg, yg, ar;

    fp = fileopen(filename,"wb");
	for(ar = 0; ar < *num_of_vgen; ar++){
		for(yg = 0; yg < num_of_geny[ar]; yg++){//探索
			for(xg = 0; xg < num_of_genx[ar]; xg++){
				fprintf(fp, "%d\n", (int)gen_tl[yg][xg][ar].x);
				fprintf(fp, "%d\n", (int)gen_tl[yg][xg][ar].y);
				fprintf(fp, "%d\n", (int)gen_br[yg][xg][ar].x);
				fprintf(fp, "%d\n", (int)gen_br[yg][xg][ar].y);
			}
		}
	}
	
    fclose(fp);
	

}

void print_lens(IMAGE *img, LABEL *label, char *outfile)
{
	printf("***print_lens***\n");
    int i, j;
    char name[100];
    FILE *fp;
	int **rorg, **gorg, **borg;
	rorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	gorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	borg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	
    sprintf(name, "./infofile/result/%s_lens.ppm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", img->width, img->height);
    for (i = 0; i < img->height; i++) {
			for (j = 0; j < img->width; j++) {
				
				rorg[i][j] = gorg[i][j] = borg[i][j] = (int)img->val[i][j];
					if((int)label->value[i][j].y  % 2 == 1){
						if((int)label->value[i][j].x % 2 == 1){
							borg[i][j] = 0;
						}else{
							gorg[i][j] = 0;
						}
					}else if((int)label->value[i][j].y % 2 == 0){
						if((int)label->value[i][j].x % 2 == 1){
							rorg[i][j] = 0;
						}else{
							rorg[i][j] = gorg[i][j] = 0;
						}
					}
					if(label->value[i][j].area == 6) gorg[i][j] = 255;
					if(label->value[i][j].area == 7) rorg[i][j] = 255;
		            putc(rorg[i][j], fp);//R
					putc(gorg[i][j], fp);//G
					putc(borg[i][j], fp);//B
					
			}
    }
    fclose(fp);
	free(rorg);
	free(gorg);
	free(borg);
    return;
}

void print_lens_center(IMAGE *img, VORONOI **voronoi, char *outfile)
{
	printf("***print_lens_center***\n");
    int i, j;
    char name[100];
    FILE *fp;
	int **rorg, **gorg, **borg, **flag, max;
	max = -10000000;
	int yg, xg;
	double tlx = TOP_LEFT_x, trx = TOP_RIGHT_x, blx = BOTTOM_LEFT_x, brx = BOTTOM_RIGHT_x;
	double tly = TOP_LEFT_y, try = TOP_RIGHT_y, bly = BOTTOM_LEFT_y, bry = BOTTOM_RIGHT_y;
	rorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	gorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	borg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	flag = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	for(i = 0; i < img->height; i++){
		for(j = 0; j< img->width; j++){
			flag[i][j] = 0;
			
		}
	}
	for(yg = 0; yg < NUM_OF_GENy+RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx+RUN_OVER; xg++){
			if((int)round(voronoi[yg][xg].gen.y) < 0) continue;
			if((int)round(voronoi[yg][xg].gen.y) >= (int)img->height) continue;
			if((int)round(voronoi[yg][xg].gen.x) < 0) continue;
			if((int)round(voronoi[yg][xg].gen.x) >= (int)img->width) continue;	
			flag[(int)round(voronoi[yg][xg].gen.y)][(int)round(voronoi[yg][xg].gen.x)]++;
		}
	}
	flag[(int)round(tly)][(int)round(tlx)]=2;
	flag[(int)round(try)][(int)round(trx)]=2;
	flag[(int)round(bly)][(int)round(blx)]=2;
	flag[(int)round(bry)][(int)round(brx)]=2;
    sprintf(name, "./infofile/result/%s_lens_center.ppm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", img->width, img->height);
    for (i = 0; i < img->height; i++) {
			for (j = 0; j < img->width; j++) {
				
				rorg[i][j] = gorg[i][j] = borg[i][j] = (int)img->val[i][j];
					if(max < (int)img->val[i][j]) max = (int)img->val[i][j];
					
					if((int)flag[i][j] != 0){
						rorg[i][j] = 0;
						gorg[i][j] = borg[i][j] = 255;
						if((int)flag[i][j] > 1) borg[i][j] = 0;
					}
		            putc(rorg[i][j], fp);//R
					putc(gorg[i][j], fp);//G
					putc(borg[i][j], fp);//B
					
			}
    }
	printf("max = %d\n", max);
    fclose(fp);
	free(rorg);
	free(gorg);
	free(borg);
	free(flag);
    return;
}

RBS **roughly_belongs_search(int height, int width, VORONOI **voronoi)
{
	int xg, yg, y, x, count;
	VORONOI *vp;
	RBS **rbs;
	int start_y, start_x, end_y, end_x, search_mask;
	search_mask = 6;
	rbs = (RBS **)alloc_2d_array(height, width, sizeof(RBS));
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			rbs[y][x].count = 0;
		}
	}
	for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
			vp = &voronoi[yg][xg];
			start_y = (int)round(vp->gen.y) - search_mask;
			start_x = (int)round(vp->gen.x) - search_mask;
			end_y = (int)round(vp->gen.y) + search_mask;
			end_x = (int)round(vp->gen.x) + search_mask;
			for(y = start_y; y <= end_y; y++){
				for(x = start_x; x <= end_x; x++){
					if(y < 0) continue;
					if(y >= height)continue;
					if(x < 0) continue;
					if(x >= width) continue;
					count = rbs[y][x].count;
					if(rbs[y][x].count >= 10){
						printf("err\n");
						exit(1);
					}
					rbs[y][x].id[count].x = (double)xg;
					rbs[y][x].id[count].y = (double)yg;
					rbs[y][x].count++;
				}
			}
		}
	}
	return(rbs);
}

inline double comp_dist_cp(int x, int y, double min_dist, double min_dist2,
	      VORONOI *vpt, VORONOI *vps)
{  
    double dx=0.0, dy=0.0;
    double dist2 = DBL_MAX;

    dx = DIFFPG((double)x, vpt->gen.x);
	if(dx < 0.0) printf("ERRdx");
    if (dx > min_dist) return (-1.0);
    dy = DIFFPG((double)y, vpt->gen.y);
	if(dy < 0.0) printf("ERRdy");
    if (dy > min_dist) return (-1.0);
    dist2 = (dx*dx + dy*dy);
	
	
    if (dist2 > min_dist2) return(-1.0);
    if (dist2 < min_dist2) return(dist2);
    if (vpt->gen.y > vps->gen.y) return(-1.0);
    if (vpt->gen.y < vps->gen.y) return(dist2);
    if (vpt->gen.x > vps->gen.x) return(-1.0);
    if (vpt->gen.x < vps->gen.x) return(dist2);
    
    return (dist2);
}

inline double comp_dist_cp2(int x, int y, double min_dist, double min_dist2,
	      double vptgeny, double vptgenx, double vpsgeny, double vpsgenx)
{  
    double dx=0.0, dy=0.0;
    double dist2 = DBL_MAX;

    dx = fabs((double)x - vptgenx);
	
    if (dx > min_dist) return (-1.0);
    dy = fabs((double)y - vptgeny);
	
    if (dy > min_dist) return (-1.0);
    dist2 = (dx*dx + dy*dy);
	
	
    if (dist2 > min_dist2) return(-1.0);
    else if (dist2 < min_dist2) return(dist2);
    else{ 
		if (vptgeny > vpsgeny) return(-1.0);
		if (vptgeny < vpsgeny) return(dist2);
		if (vptgenx > vpsgenx) return(-1.0);
		if (vptgenx < vpsgenx) return(dist2);
	}
    
    return (dist2);
}


void partition_center(VORONOI **voronoi, VORONOI *voronoi_center, LABEL *label)
{
	int x, y;
    int xg, yg, l;
    int height, width;
	height = label->height;
	width = label->width;
    double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp, *vpc, *ivpc;
    printf("c\n");
    printf("h,w = %d,%d\n",label->height, label->width);
	ivp = &voronoi[0][0];
	vp = &voronoi[0][0];
	vpc = &voronoi_center[0];
	ivpc = &voronoi_center[0];
	if(NUM_OF_VGENc >= 1.0){
		for (y = 0; y < height; y++) { //label範囲のすべての画素について
			for (x = 0; x < width; x++) {
				if(label->value[y][x].area == 100) continue;
				xg = (int)label->value[y][x].x;
				yg = (int)label->value[y][x].y;
				vp = &voronoi[yg][xg];
				min_dist2 = DBL_MAX;
				min_dist = DBL_MAX;
				/* searching */
				ivp = &voronoi[0][0];
				ivpc = &voronoi_center[0];
				for (l = 0;  l < NUM_OF_VGENc; l++) { //マイクロレンズ内の７つの母点でボロノイ分割
					vpc = &voronoi_center[l];
					dist2 = comp_dist_cp2(x, y, min_dist, min_dist2, 
					vp->gen.y+vpc->gen.y, vp->gen.x+vpc->gen.x, 
					ivp->gen.y+ivpc->gen.y, ivp->gen.x+ivpc->gen.x); //x,yと最も近い母点との距
					if (dist2 >= 0){
						ivp = vp;
						ivpc = vpc;
						min_dist2 = dist2;
						min_dist = (sqrt((double)min_dist2));
					}
				}
				label->value[y][x].area = ivpc->id.area; //labelのvalueをidに（最も近い母点のid）
			}
		}
	}
}
void partition_edge(VORONOI **voronoi, VORONOI *voronoi_edge, LABEL *label)
{
	int x, y;
    int xg, yg, l;
    int height, width;
	height = label->height;
	width = label->width;
    double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp, *vpc, *ivpc;
    printf("c\n");
    printf("h,w = %d,%d\n",label->height, label->width);
	ivp = &voronoi[0][0];
	vp = &voronoi[0][0];
	vpc = &voronoi_edge[0];
	ivpc = &voronoi_edge[0];
	if((NUM_OF_VGEN - NUM_OF_VGENc) >= 1.0){
		for (y = 0; y < height; y++) { //label範囲のすべての画素について
			for (x = 0; x < width; x++) {
				if(label->value[y][x].area == 100){
					xg = (int)label->value[y][x].x;
					yg = (int)label->value[y][x].y;
					vp = &voronoi[yg][xg];
					min_dist2 = DBL_MAX;
					min_dist = DBL_MAX;
					/* searching */
					ivp = &voronoi[0][0];
					ivpc = &voronoi_edge[0];
					for (l = 0;  l < NUM_OF_VGEN - NUM_OF_VGENc; l++) { //マイクロレンズ内の７つの母点でボロノイ分割
						vpc = &voronoi_edge[l];
						dist2 = comp_dist_cp2(x, y, min_dist, min_dist2, 
						vp->gen.y+vpc->gen.y, vp->gen.x+vpc->gen.x, 
						ivp->gen.y+ivpc->gen.y, ivp->gen.x+ivpc->gen.x); //x,yと最も近い母点との距
						if (dist2 >= 0){
							ivp = vp;
							ivpc = vpc;
							min_dist2 = dist2;
							min_dist = (sqrt((double)min_dist2));
						}
					}
					label->value[y][x].area = ivpc->id.area; //labelのvalueをidに（最も近い母点のid）
				}
			}
		}
	}
}

void fill_recur2(RBS **rbs, VORONOI **voronoi, LABEL *label)
{
	int x, y, height, width, rbsidx, rbsidy, i, count;
	double min_dist, min_dist2, dist2;
	VORONOI *vp, *ivp;
	height = label->height;
	width = label->width;
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			min_dist2 = DBL_MAX;
			min_dist = DBL_MAX;
			count = rbs[y][x].count;
			if(rbs[y][x].count == 1){
				label->value[y][x].x = rbs[y][x].id[0].x;
				label->value[y][x].y = rbs[y][x].id[0].y;
			}else{
				vp = &voronoi[0][0];
				ivp = vp;
				for(i = 0; i < count; i++){
					rbsidy = (int)rbs[y][x].id[i].y;
					rbsidx = (int)rbs[y][x].id[i].x;
					vp = &voronoi[rbsidy][rbsidx];
					dist2 = comp_dist_cp2(x, y, min_dist, min_dist2, 
					vp->gen.y, vp->gen.x, ivp->gen.y, ivp->gen.x);
					if(dist2 >= 0.0){
						ivp = vp;
						min_dist2 = dist2;
						min_dist = sqrt((double)min_dist2);
					}
				}
				label->value[y][x].x = ivp->id.x;
				label->value[y][x].y = ivp->id.y;
			}
			label->value[y][x].area = 0;
			ivp = &voronoi[(int)label->value[y][x].y][(int)label->value[y][x].x];
			ivp->tl.x = MIN(ivp->tl.x, (double)x);
			ivp->tl.y =	MIN(ivp->tl.y, (double)y);
            ivp->br.x = MAX(ivp->br.x, (double)x);
            ivp->br.y = MAX(ivp->br.y, (double)y);
		}
	}
}

void partition_vignetting(VORONOI **voronoi, LABEL *label, IMAGE *img)
{
	printf("partition_vignetting->");
	int y, x, yg, xg, i;
	int tlx, tly, bry, brx, height, width;
	double dy, dx, dist2, dist, mla_radius;
	double hist[2][256], prob[2][256];
	double log2 = log(2.0);
	int count[2]={0};
	int id, idmax, save_id=0;
	idmax=20;
	double cost, min_cost;
	
	//mla_radius = MLADIST * 0.371;//1:1
	
	VORONOI *vp;
	height = label->height;
	width = label->width;
	min_cost = DBL_MAX;
	for(id = 0; id < idmax; id++){
		for(y = 0; y < 2; y++){
			for(x = 0; x < 256; x++){
				hist[y][x]=prob[y][x]=0.0;
			}
			count[y] =0;
		}
		mla_radius = MLADIST * 0.3 + MLADIST * 0.01 * (double)id;
	
		cost = 0.0;
		for(yg = 0; yg < NUM_OF_GENy+RUN_OVER; yg++){
			for(xg = 0; xg < NUM_OF_GENx+RUN_OVER; xg++){
				vp = &voronoi[yg][xg];
				tly = vp->tl.y;
				tlx = vp->tl.x;
				bry = vp->br.y+1;
				brx = vp->br.x+1;
				for(y = tly; y < bry; y++){
					for(x = tlx; x < brx; x++){
						if(y < 0 || y >= height) continue;
						if(x < 0 || x >= width) continue;
						if(yg == label->value[y][x].y && xg == label->value[y][x].x){
							dy = y - vp->gen.y;
							dx = x - vp->gen.x;
							dist2 = pow(dy,2)+pow(dx,2);
							dist = sqrt(dist2);
							if(dist < mla_radius) label->value[y][x].area = 0;
							else label->value[y][x].area = 1;
						}
					}
				}
			}
		}
		for(y = 0; y < height; y++){
			for(x = 0; x < width; x++){
				hist[(int)label->value[y][x].area][(int)img->val[y][x]]++;
				count[(int)label->value[y][x].area]++;
			}
		}
		for(y = 0; y < 2; y++){
			for(x = 0; x < 256; x++){
				prob[y][x] = hist[y][x] / (double)count[y];
				if(prob[y][x] == 0.0) prob[y][x] = DBL_MIN;
			}
			for(i = 0; i < 256; i++){
				cost += -log(prob[y][i]) * hist[y][i];
			}
		}
		cost /= log2;
		printf("cost=%f\n", cost);
		if(cost < min_cost){
			min_cost = cost;
			save_id = id;
		}
	}
	
	mla_radius = MLADIST * 0.3 + MLADIST * 0.01 * (double)save_id;
	//mla_radius = MLADIST * 0.371;//1:1
	printf("optimize_mla_radius = %f\n", mla_radius);
	for(yg = 0; yg < NUM_OF_GENy+RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx+RUN_OVER; xg++){
			vp = &voronoi[yg][xg];
			tly = vp->tl.y;
			tlx = vp->tl.x;
			bry = vp->br.y+1;
			brx = vp->br.x+1;
			for(y = tly; y < bry; y++){
				for(x = tlx; x < brx; x++){
					if(y < 0 || y >= height) continue;
					if(x < 0 || x >= width) continue;
					if(yg == label->value[y][x].y && xg == label->value[y][x].x){
						dy = y - vp->gen.y;
						dx = x - vp->gen.x;
						dist2 = pow(dy,2)+pow(dx,2);
						dist = sqrt(dist2);
						if(dist < mla_radius) label->value[y][x].area = 0;
						else label->value[y][x].area = 100;
					}
				}
			}
		}
	}
	
	printf("ok\n");
}

void paritition_along_mla(VORONOI **voronoi, LABEL *label, RBS **rbs)
{	
	printf("paritition_along_mla->");
    int x, y;
    int yg, xg;
    VORONOI *vp;
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
		for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
			vp = &voronoi[yg][xg];
      		vp->tl.x = vp->br.x = vp->gen.x;
      		vp->tl.y = vp->br.y = vp->gen.y;
      		vp->num_pixel = 0;
		}
    }
	
    for (y = 0; y < label->height; y++) {
  		for (x = 0; x < label->width; x++) {
        label->value[y][x].x = -1;
        label->value[y][x].y = -1;
  		label->value[y][x].area = -1;
  		}
    }
	
	fill_recur2(rbs, voronoi, label);
	//fill_recur(voronoi, label);
		
	printf("ok\n");
    return;
}





void label_check(ENCODER *enc, LABEL *label)
{
  FILE *fp;
  int i, j;
  int x, y;

  fp = fileopen("./outfile_label.txt", "rb");

  for (i = 0; i < enc->height; i++) {
    for (j = 0; j < enc->width; j++) {
      fscanf(fp, "%d\n", &x);
      fscanf(fp, "%d\n", &y);
      if((int)label->value[i][j].x != x)printf("errx\n");
	  if((int)label->value[i][j].y != y)printf("erry\n");
    }
  }

  fclose(fp);
}

void make_edge_label(IMAGE *img, ENCODER *enc, int **edge_label, LABEL *label)
{
	int y, x, height, width;
	int i, j, tly, tlx, bry, brx, karnel_size;
	FILE *fp;
	char name[100];
	karnel_size = 1;
	height = enc->height;
	width = enc->width;
	printf("h.w = %d,%d\n", height, width);
	//for(y = 0; y < height; y++){
		//printf("%d\n",y);
		//tly = y - karnel_size;
		//bry = y + karnel_size;
		//for(x = 0; x < width; x++){
			/*tlx = x - karnel_size;
			brx = x + karnel_size;
			if(tly > 0){
				i = tly;
				j = x;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(tlx > 0){
				i = y;
				j = tlx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(tly > 0 && tlx > 0){
				i = tly;
				j = tlx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(bry < height){
				i = bry;
				j = x;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(brx < width){
				i = y;
				j = brx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(bry < height && brx < width){
				i = bry;
				j = brx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}
			if(tly > 0 && tlx > 0 && bry < height && brx < width){
				i = tly;
				j = brx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
				i = bry;
				j = tlx;
				if((save_labelx[i][j] != save_labelx[y][x]) 
				|| (save_labely[i][j] != save_labely[y][x])) edge_label[y][x]=100;
			}*/
			//edge_label[y][x] = label->value[y][x].area*100;
		//}
		
	//}
	
	sprintf(name, "./outfile_edge_label.pgm");
	fp = fileopen(name, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
            putc(edge_label[i][j], fp);
	     }
    }
    fclose(fp);
	sprintf(name, "./outfile_edge_label.ppm");
    fp = fileopen(name, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
		            
					if(edge_label[i][j] > 0){
						putc(0,fp);
						putc(0,fp);
						putc(255, fp);
					}else{
						putc((int)img->val[i][j], fp);//B
						putc((int)img->val[i][j], fp);//R
						putc((int)img->val[i][j], fp);//G
					}
					
					
			}
    }
	fclose(fp);
	
}

VORONOI **alloc_voronoi(LABEL *label)
{
	printf("alloc_voronoi->");
	printf("MLADIST = %lf", MLADIST);
	VORONOI **voronoi, *vp;
	int xg, yg, height, width;
	height = label->height;
	width = label->width;
	double tlx = TOP_LEFT_x, trx = TOP_RIGHT_x, blx = BOTTOM_LEFT_x, brx = BOTTOM_RIGHT_x;
	double tly = TOP_LEFT_y, try = TOP_RIGHT_y, bly = BOTTOM_LEFT_y, bry = BOTTOM_RIGHT_y;
	double r = RADIUS;
	double shift_y, shift_x;
	double sensor_offset_y, sensor_offset_x;
	double rd;
	double dblxg, dblyg;
	double cos1, sin1, cos2, cos3, sin3;
	double cent_dist_x, cent_dist_y;
	double rotated_bly, rotated_blx, rotated_bry, rotated_brx;
	double slope_a,slope_b,slope_c,slope_d;
	cent_dist_x = ((trx - tlx)/NUM_OF_GENx+(brx - blx)/NUM_OF_GENx)/ 2.0;
	cent_dist_y = ((bly - tly)/NUM_OF_GENy + (bry - try)/NUM_OF_GENy)/2.0;
	cos1 = ((trx - tlx)/(sqrt(pow(trx - tlx,2.0) + pow(try - tly, 2.0)))
		 + (brx - blx)/(sqrt(pow(brx - blx,2.0) + pow(bry - bly, 2.0))))/ 2.0;
	sin1 = ((tly - try)/(sqrt(pow(trx - tlx,2.0) + pow(try - tly, 2.0)))
		 + (bly - bry)/(sqrt(pow(trx - tlx,2.0) + pow(try - tly, 2.0))))/ 2.0;
	cos2 = ((bly - tly)/(sqrt(pow(blx - tlx,2.0) + pow(bly - tly, 2.0)))
		 + (bry - try)/(sqrt(pow(brx - trx,2.0) + pow(bry - try, 2.0))))/ 2.0;
	rotated_blx = (blx - tlx) * cos1 - (bly - tly) * sin1 + tlx;
	rotated_bly = (blx - tlx) * sin1 + (bly - tly) * cos1 + tly;
	rotated_brx = (brx - tlx) * cos1 - (bry - tly) * sin1 + tlx;
	rotated_bry = (brx - tlx) * sin1 + (bry - tly) * cos1 + tly;
	printf("rotated_blx,y = %f,%f\n",rotated_bly, rotated_blx);
	cos3 = ((rotated_bly - tly)/(sqrt(pow(rotated_blx - tlx, 2.0) + pow(rotated_bly - tly, 2.0)))
		  + (rotated_bry - try)/(sqrt(pow(rotated_brx - trx, 2.0) + pow(rotated_bry - try, 2.0))))/ 2.0;
	sin3 = ((rotated_blx - tlx)/(sqrt(pow(rotated_blx - tlx, 2.0) + pow(rotated_bly - tly, 2.0)))
		 + (rotated_brx - trx)/(sqrt(pow(rotated_brx - trx, 2.0) + pow(rotated_bry - try, 2.0))))/ 2.0;
	printf("x,y = %f,%f\n",cent_dist_x, cent_dist_y);
	printf("cos1,sin1,cos2,cos3 = %f,%f,%f,%f\n",cos1, sin1,cos2, cos3);
	rd = ROTATION;
	slope_a = (try - tly)/(trx - tlx);
	slope_b = (bly -  tly)/(blx - tlx);
	slope_c = (bry - bly)/(brx - blx);
	slope_d = (bry - try)/(brx - trx);
	printf("a,b,c,d = %f,%f,%f,%f\n", slope_a, slope_b, slope_c, slope_d);
	sensor_offset_y = SENSOR_OFFSET_y;
	sensor_offset_x = SENSOR_OFFSET_x;
	shift_y = (height-1.0)/2.0 + sensor_offset_y;
	shift_x = (width-1.0)/2.0 + sensor_offset_x;
	voronoi = (VORONOI **)alloc_2d_array(NUM_OF_GENy+RUN_OVER, NUM_OF_GENx+RUN_OVER, sizeof(VORONOI));

	for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
			dblxg = (double)xg;
			dblyg = (double)yg;
			if(yg % 2 == 0){
				vp = &voronoi[yg][xg];
				vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) + (dblxg * cent_dist_x/cos1) + sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2));
				vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
				vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
				vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
				vp->id.x = xg;
				vp->id.y = yg;
				vp->id.area = 0;
			}
			if(yg % 2 == 1){
				vp = &voronoi[yg][xg];
				vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) + sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2));
				vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
				vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
				vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
				vp->id.x = xg;
				vp->id.y = yg;
				vp->id.area = 0;
			}
		}
	}
	printf("ok\n");
	return(voronoi);
}


VORONOI *alloc_voronoi_center(void)
{
	printf("alloc_voronoi_center->");
	VORONOI *voronoi_center, *vp;
	int ar;
	double r = RADIUS, rd = 2*M_PI/NUM_OF_VGENc;
	double g_offset, offset;
	g_offset = -M_PI/2.0;
	offset = rd / 2.0;
	
	voronoi_center = (VORONOI *)alloc_mem(NUM_OF_VGENc *sizeof(VORONOI));
	if(NUM_OF_VGENc > 1.0){
		for(ar = 0; ar < NUM_OF_VGENc; ar++){
			vp = &voronoi_center[ar];
			vp->gen.x = r * cos(rd * ar + g_offset + offset);
			vp->gen.y = r * sin(rd * ar + g_offset + offset);
			vp->id.x = 0;
			vp->id.y = 0;
			vp->id.area = ar;
		}
	}else{
		vp = &voronoi_center[0];
			vp->gen.x = 0;
			vp->gen.y = 0;
			vp->id.x = 0;
			vp->id.y = 0;
			vp->id.area = 0;
	}
	printf("ok\n");
	return(voronoi_center);
}

VORONOI *alloc_voronoi_edge(void)
{
	printf("alloc_voronoi_edge->");
	VORONOI *voronoi_edge, *vp;
	int ar;
	double g_offset, offset;
	double r = RADIUS, rd = 2*M_PI/(NUM_OF_VGEN - NUM_OF_VGENc);
	g_offset = -M_PI/2.0;
	offset = rd / 2.0;
	
	voronoi_edge = (VORONOI *)alloc_mem((NUM_OF_VGEN - NUM_OF_VGENc) *sizeof(VORONOI));
	if((NUM_OF_VGEN - NUM_OF_VGENc) > 1.0){
		for(ar = 0; ar < NUM_OF_VGEN - NUM_OF_VGENc; ar++){
			vp = &voronoi_edge[ar];
			vp->gen.x = r * cos(rd * ar + g_offset + offset);
			vp->gen.y = r * sin(rd * ar + g_offset + offset);
			vp->id.x = 0;
			vp->id.y = 0;
			vp->id.area = ar+NUM_OF_VGENc;
		}
	}else{
		vp = &voronoi_edge[0];
			vp->gen.x = 0;
			vp->gen.y = 0;
			vp->id.x = 0;
			vp->id.y = 0;
			vp->id.area = (int)(NUM_OF_VGEN-1);
	}
	printf("ok\n");
	return(voronoi_edge);
}

LABEL *alloc_label(int height, int width, int init)
{
    LABEL *label;
	printf("alloc_label->");
    label = (LABEL *)alloc_mem(sizeof(LABEL));
    label->height = height;
	label->width = width;
    label->value = (POINT **)alloc_2d_array(height, width, sizeof(POINT));
	printf("ok\n");
    return (label);
}




void merge_center(int *num_of_geny, int *num_of_genx, LABEL *label, VORONOI **voronoi, POINT ***gen_tl, POINT ***gen_br)
{
	printf("merge_center->");
	int yg, xg, ar, y, x, height, width;
	double max_y, max_x;
	POINT cp;
	height = label->height;
	width = label->width;
	max_y = 0.0;
	max_x = 0.0;
	for(yg = 0; yg < NUM_OF_GENy+RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx+RUN_OVER; xg++){
			for(ar = 0; ar < NUM_OF_VGENc; ar++){
				gen_tl[yg][xg][ar].x = width-1;
				gen_tl[yg][xg][ar].y = height-1;
				gen_br[yg][xg][ar].x = 0;
				gen_br[yg][xg][ar].y = 0;
			}
		}
	}
	max_y = (int)((NUM_OF_GENy+RUN_OVER-1.0) / MERGE_CENTERy);
	max_x = (int)((NUM_OF_GENx+RUN_OVER-1.0) / MERGE_CENTERx);
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			if(label->value[y][x].area < NUM_OF_VGENc){
				label->value[y][x].x = (int)(label->value[y][x].x / MERGE_CENTERx);
				label->value[y][x].y = (int)(label->value[y][x].y / MERGE_CENTERy);
				cp = label->value[y][x];
				gen_tl[(int)cp.y][(int)cp.x][cp.area].x = 
							MIN(x, gen_tl[(int)cp.y][(int)cp.x][cp.area].x);
				gen_tl[(int)cp.y][(int)cp.x][cp.area].y = 
							MIN(y, gen_tl[(int)cp.y][(int)cp.x][cp.area].y);
				gen_br[(int)cp.y][(int)cp.x][cp.area].x = 
							MAX(x, gen_br[(int)cp.y][(int)cp.x][cp.area].x);
				gen_br[(int)cp.y][(int)cp.x][cp.area].y = 
							MAX(y, gen_br[(int)cp.y][(int)cp.x][cp.area].y);
			}
		}
	}
	for(ar = 0; ar < NUM_OF_VGENc; ar++){
		num_of_geny[ar] = max_y+1;
		num_of_genx[ar] = max_x+1;
		printf("geny,genx=%d,%d\n", num_of_geny[ar], num_of_genx[ar]);
	}
	printf("ok\n");
}

void merge_edge(int *num_of_geny, int *num_of_genx, LABEL *label, VORONOI **voronoi, POINT ***gen_tl, POINT ***gen_br)
{
	printf("merge_edge->");
	int yg, xg, ar, y, x, height, width;
	double max_y, max_x;
	POINT cp;
	height = label->height;
	width = label->width;
	max_y = 0.0;
	max_x = 0.0;
	for(yg = 0; yg < NUM_OF_GENy+RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx+RUN_OVER; xg++){
			for(ar = NUM_OF_VGENc; ar < NUM_OF_VGEN; ar++){
				gen_tl[yg][xg][ar].x = width-1;
				gen_tl[yg][xg][ar].y = height-1;
				gen_br[yg][xg][ar].x = 0;
				gen_br[yg][xg][ar].y = 0;
			}
		}
	}
	max_y = (int)((NUM_OF_GENy+RUN_OVER-1.0) / MERGE_EDGEy);
	max_x = (int)((NUM_OF_GENx+RUN_OVER-1.0) / MERGE_EDGEx);
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			if(label->value[y][x].area >= NUM_OF_VGENc && label->value[y][x].area < NUM_OF_VGEN){
				label->value[y][x].x = (int)(label->value[y][x].x / MERGE_EDGEx);
				label->value[y][x].y = (int)(label->value[y][x].y / MERGE_EDGEy);
				cp = label->value[y][x];
				gen_tl[(int)cp.y][(int)cp.x][cp.area].x = 
							MIN(x, gen_tl[(int)cp.y][(int)cp.x][cp.area].x);
				gen_tl[(int)cp.y][(int)cp.x][cp.area].y = 
							MIN(y, gen_tl[(int)cp.y][(int)cp.x][cp.area].y);
				gen_br[(int)cp.y][(int)cp.x][cp.area].x = 
							MAX(x, gen_br[(int)cp.y][(int)cp.x][cp.area].x);
				gen_br[(int)cp.y][(int)cp.x][cp.area].y = 
							MAX(y, gen_br[(int)cp.y][(int)cp.x][cp.area].y);
			}
		}
	}
	for(ar = NUM_OF_VGENc; ar < NUM_OF_VGEN; ar++){
		num_of_geny[ar] = max_y+1;
		num_of_genx[ar] = max_x+1;
		printf("geny,genx=%d,%d\n", num_of_geny[ar], num_of_genx[ar]);
	}
	printf("ok\n");
}

void merge_label(LABEL *label, int *num_of_vgenc, int *num_of_vgen)
{
	int y,x, height, width;
	height = label->height;
	width = label->width;
	//統合する領域を記述
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			if(label->value[y][x].area == 0) label->value[y][x].area = 0;
		}
	}
	//領域番号を修正（番号をつめる）
	//領域数を修正
}

int main(int argc, char **argv)
{
  ENCODER *enc;
  IMAGE *img;
  int height, width, *num_of_geny, *num_of_genx, num_of_vgenc, num_of_vgen;
  num_of_geny = (int *)alloc_mem(NUM_OF_VGEN * sizeof(int));
  num_of_genx = (int *)alloc_mem(NUM_OF_VGEN * sizeof(int));
  char *infile = NULL, *outfile = NULL;
  enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
  int **edge_label;
  int edge_pix=0, cent_pix=0;
  int x,y,ar;
  double elapse = 0.0;
  char *class_p;
  infile = argv[1];
  outfile = argv[2];
  printf("***infile -> %s***\n", infile);
  printf("***outfile_name -> %s***\n", outfile);
  for(ar = 0; ar < NUM_OF_VGEN; ar++){
	  num_of_geny[ar] = NUM_OF_GENy+RUN_OVER;
	  num_of_genx[ar] = NUM_OF_GENx+RUN_OVER;
  }
  num_of_vgen = NUM_OF_VGEN;
  num_of_vgenc = NUM_OF_VGENc;
  img = read_pgm(infile);
  height = enc->height = img->height;
  width = enc->width = img->width;
  enc->class = (char **)alloc_2d_array(height, width,
         sizeof(char));
  edge_label = (int **)alloc_2d_array(height, width, sizeof(int));
  for(y = 0; y < height; y++){
	  for(x = 0; x < width; x++){
		  edge_label[y][x] = 0;
	  }
  }
  VORONOI **voronoi,*voronoi_center, *voronoi_edge;
  LABEL *label;
  POINT ***gen_tl, ***gen_br;
  gen_tl = (POINT ***)alloc_3d_array(NUM_OF_GENy+RUN_OVER, NUM_OF_GENx+RUN_OVER, NUM_OF_VGEN, sizeof(POINT));
  gen_br = (POINT ***)alloc_3d_array(NUM_OF_GENy+RUN_OVER, NUM_OF_GENx+RUN_OVER, NUM_OF_VGEN, sizeof(POINT));
  RBS **rbs;
  label = alloc_label(height, width, 0);
  voronoi = alloc_voronoi(label);
  if(NUM_OF_VGENc >= 1.0){
	voronoi_center = alloc_voronoi_center();
  }
  if(NUM_OF_VGEN - NUM_OF_VGENc >= 1.0){
	voronoi_edge = alloc_voronoi_edge();
  }
  rbs = roughly_belongs_search(height, width, voronoi);
  //ここまで宣言と初期化
  
  
  paritition_along_mla(voronoi, label, rbs);
  free(rbs);
  //ここまでマイクロレンズに沿った分割
  if(NUM_OF_VGEN - NUM_OF_VGENc >= 1.0){
	partition_vignetting(voronoi, label, img);
  }else{
	  printf("周辺部は作られません");
  }
  //ここまでマイクロレンズ内に中心部と周辺部を形成
  for(y = 0;y < height; y++){
	class_p = enc->class[y];
    for(x =0; x < width; x++){
      class_p[x] = ((char)(((int)label->value[y][x].x % 32) + ((int)label->value[y][x].y % 32)))*label->value[y][x].area;
	  if(label->value[y][x].area == 0) cent_pix++;
	  else edge_pix++;
	  edge_label[y][x] = label->value[y][x].area;
	  //label->value[y][x].area = 0;
	}
  }
  printf("vignetting_ratio = %f\n", (double)edge_pix/(double)cent_pix);
  make_edge_label(img, enc, edge_label, label);
  //ここまで中心部と周辺部の画像の書き出し
  if(NUM_OF_VGENc >= 1.0){
	partition_center(voronoi, voronoi_center, label);
  }else{
	  printf("中心部は分割されません\n");
  }
  //ここまで中心部をNUM_OF_VGENc領域に分割
  if(NUM_OF_VGEN - NUM_OF_VGENc >= 1.0){
	partition_edge(voronoi, voronoi_edge, label);
  }else{
	  printf("周辺部は分割されません\n");
  }
  //ここまで周辺部をNUM_OF_VGEN-NUM_OF_VGENc領域に分割
  if(NUM_OF_VGENc >= 1.0){
	merge_center(num_of_geny, num_of_genx, 
				label, voronoi, gen_tl, gen_br);
  }else{
	  printf("中心部は統合されません\n");
  }
  //ここまで中心部の同一形状ブロック共有
  if(NUM_OF_VGEN - NUM_OF_VGENc >= 1.0){
	merge_edge(num_of_geny, num_of_genx, 
				label, voronoi, gen_tl, gen_br);
  }else{
	  printf("周辺部は統合されません\n");
  }
  //ここまで周辺部の同一形状ブロック共有
  merge_label(label, &num_of_vgenc, &num_of_vgen);
  //ここまで異なる形状同士で領域番号を統合
  
  print_label(label, enc->height, enc->width, &num_of_vgenc,
				&num_of_vgen, num_of_geny, num_of_genx);
  print_range(gen_tl, gen_br, &num_of_vgen, num_of_geny, 
				num_of_genx, "outfile_range.txt");
  print_class(enc->class, enc->height, enc->width);
  print_lens(img, label, outfile);
  print_lens_center(img, voronoi, outfile);
  free(label->value);
  free(label);
  free(voronoi);
  free(enc->class);
  free(enc);
  free(edge_label);
  return 0;
}