#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "mrp1.h"
#include <float.h>
#include "vseg.h"

void fill_label_grobal(VORONOI ***, LABEL *, RBS **);
//void fill_recur(VORONOI ***, LABEL *);
void fill_recur2(RBS **, VORONOI ***, LABEL *);
int set_range(VORONOI *, LABEL *, RANGE *, ENCODER *, int );
RANGE *alloc_range(int );
LABEL *alloc_label(int , int , int );
VORONOI ***alloc_voronoi(void);
long comp_dist(int , int , long , long , VORONOI *, VORONOI *);
inline double comp_dist_cp(int , int , double , double , VORONOI *, VORONOI *);
inline double comp_dist_cp2(int, int, double, double, double, double, double, double);


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
/*
void **alloc_2d_array(int height, int width, int size)
{
    void **mat;
    char *ptr;
    int k;

    mat = (void **)alloc_mem(sizeof(void *) * height + height * width * size);
    ptr = (char *)(mat + height);
    for (k = 0; k < height; k++) {
	mat[k] =  ptr;
	ptr += width * size;
    }
    return (mat);
}
void ***alloc_3d_array(int height, int width, int depth, int size)
{
    void ***mat;
    char *ptr;
    int k, l;

    mat = (void ***)alloc_2d_array(height, width, sizeof(void *));
    ptr = (char *)alloc_mem(height * width * depth * size);
    for (k = 0; k < height; k++) {
        for (l = 0; l < width; l++) {
            mat[k][l] = ptr;
            ptr += depth * size;
	}
    }
    return(mat);
}*/

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

void print_label(LABEL *label, int height, int width)
{
    int i, j;
    char name[100];
    FILE *fp;

    sprintf(name, "outfile_label.txt");
    fp = fileopen(name, "wb");
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
         fprintf(fp, "%d\n",(int)(label->value[i][j].x));
         fprintf(fp, "%d\n",(int)(label->value[i][j].y));
         fprintf(fp, "%d\n",(int)(label->value[i][j].area));
	     }
    }
    fclose(fp);
    return;
}

void print_range(VORONOI ***voronoi, int num_of_gen, char *filename)
{
    FILE *fp;
    VORONOI *vp;
    int xg, yg, l;

    fp = fileopen(filename,"wb");
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){//探索
      for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
          vp = &voronoi[yg][xg][0];
          fprintf(fp, "%d\n", (int)vp->tl.x);
          fprintf(fp, "%d\n", (int)vp->tl.y);
          fprintf(fp, "%d\n", (int)vp->br.x);
          fprintf(fp, "%d\n", (int)vp->br.y);
      }
    }
    fclose(fp);
    return;
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

void print_lens_center(IMAGE *img, VORONOI ***voronoi, char *outfile)
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
			if((int)round(voronoi[yg][xg][0].gen.y) < 0) continue;
			if((int)round(voronoi[yg][xg][0].gen.y) >= (int)img->height) continue;
			if((int)round(voronoi[yg][xg][0].gen.x) < 0) continue;
			if((int)round(voronoi[yg][xg][0].gen.x) >= (int)img->width) continue;	
			flag[(int)round(voronoi[yg][xg][0].gen.y)][(int)round(voronoi[yg][xg][0].gen.x)]++;
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

RBS **roughly_belongs_search(int height, int width, VORONOI ***voronoi)
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
			vp = &voronoi[yg][xg][0];
			start_y = (int)round(vp->gen.y) - search_mask;
			start_x = (int)round(vp->gen.x) - search_mask;
			end_y = (int)round(vp->gen.y) + search_mask;
			end_x = (int)round(vp->gen.x) + search_mask;
			/*if(start_y < 0) start_y = 0;
			if(start_x < 0) start_x = 0;
			if(end_y < 0) end_y = 0;
			if(end_x < 0) end_x = 0;
			if(start_y >= height) start_y = height -1;
			if(start_x >= width) start_x = width -1;
			if(end_y >= height) end_y = height -1;
			if(end_x >= width) end_x = width -1;*/
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


/*void fill_recur(VORONOI ***voronoi, LABEL *label)
{	
	printf("fill_recur->");
    int x, y;
    int xg, yg, l;
	int height, width;
	double tlx = TOP_LEFT_x, trx = TOP_RIGHT_x, blx = BOTTOM_LEFT_x, brx = BOTTOM_RIGHT_x;
	//double tly = TOP_LEFT_y, try = TOP_RIGHT_y, bly = BOTTOM_LEFT_y, bry = BOTTOM_RIGHT_y;
	double cent_dist_x = ((trx - tlx)/329.0+(brx - blx)/329.0)/ 2.0;
	double absolute_flag = cent_dist_x * 0.4; 
    double min_dist, min_dist2, dist2;
    VORONOI *vp;
    VORONOI *ivp;
	height = label->height;
	width = label->width;
    //printf("c\n");
    l = 0;
    
	printf("height, width = %d, %d\n", height, width);
    for (y = 0; y < height; y++) {
	printf("%d\n",(int)y);
  		for (x = 0; x < width; x++) {
  		
		//printf("%d ", x);
        min_dist2 = DBL_MAX;
        min_dist = DBL_MAX;
		
		vp = &voronoi[0][0][l];
		ivp = vp;
		//printf("num_yg, num_xg = %d,%d", NUM_OF_GENy, NUM_OF_GENx);
		
			for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
				if(min_dist <= absolute_flag) break;
				for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
					vp = &voronoi[yg][xg][l];
					if(fabs((double)y - vp->gen.y) > 10.0 || fabs((double)x - vp->gen.x) > 10.0){
						continue;
					}
					
					
					dist2 = comp_dist_cp2(x, y, min_dist, min_dist2, vp->gen.y, vp->gen.x, ivp->gen.y, ivp->gen.x); //
					//printf("E ");
					if (dist2 >= 0.0){
						ivp = vp;
						min_dist2 = dist2;
						min_dist = sqrt((double)min_dist2);
						if(min_dist <= absolute_flag) break;
					}
				}
			}
			
  			label->value[y][x].x = ivp->id.x;
			label->value[y][x].y = ivp->id.y;
			label->value[y][x].area = ivp->id.area;
  			ivp->num_pixel++;
			ivp->tl.x = MIN(ivp->tl.x, (double)x);
			ivp->tl.y =	MIN(ivp->tl.y, (double)y);
            ivp->br.x = MAX(ivp->br.x, (double)x);
            ivp->br.y = MAX(ivp->br.y, (double)y);
			//printf("E");
  		}
	
		
    }
}*/
void fill_recur2(RBS **rbs, VORONOI ***voronoi, LABEL *label)
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
				vp = &voronoi[0][0][0];
				ivp = vp;
				for(i = 0; i < count; i++){
					rbsidy = (int)rbs[y][x].id[i].y;
					rbsidx = (int)rbs[y][x].id[i].x;
					vp = &voronoi[rbsidy][rbsidx][0];
					dist2 = comp_dist_cp2(x, y, min_dist, min_dist2, vp->gen.y, vp->gen.x, ivp->gen.y, ivp->gen.x);
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
			ivp = &voronoi[(int)label->value[y][x].y][(int)label->value[y][x].x][0];
			ivp->tl.x = MIN(ivp->tl.x, (double)x);
			ivp->tl.y =	MIN(ivp->tl.y, (double)y);
            ivp->br.x = MAX(ivp->br.x, (double)x);
            ivp->br.y = MAX(ivp->br.y, (double)y);
		}
	}
}

void partition_vignetting(VORONOI ***voronoi, LABEL *label, IMAGE *img)
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
				vp = &voronoi[yg][xg][0];
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
			vp = &voronoi[yg][xg][0];
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
	
	printf("ok\n");
}

void fill_label_grobal(VORONOI ***voronoi, LABEL *label, RBS **rbs)
{	
	printf("fill_label_grobal->");
    int x, y, l;
    int yg, xg;
    VORONOI *vp;
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
      for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
        for(l = 0; l < 7; l++){
          vp = &voronoi[yg][xg][l];
      		vp->tl.x = vp->br.x = vp->gen.x;
      		vp->tl.y = vp->br.y = vp->gen.y;
      		vp->num_pixel = 0;
        }
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

int main(int argc, char **argv)
{
  ENCODER *enc;
  IMAGE *img;
  int height, width;
  char *infile = NULL, *outfile = NULL;
  enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
  int **edge_label;
  int edge_pix=0, cent_pix=0;
  int x,y;
  double elapse = 0.0;
  char *class_p;
  infile = argv[1];
  outfile = argv[2];
  printf("***infile -> %s***\n", infile);
  printf("***outfile_name -> %s***\n", outfile);
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
  VORONOI ***voronoi;
  LABEL *label;
  RBS **rbs;

  label = alloc_label(width,height,0);
  voronoi = alloc_voronoi();
  rbs = roughly_belongs_search(height, width, voronoi);
  

  fill_label_grobal(voronoi, label, rbs);
  free(rbs);
  partition_vignetting(voronoi, label, img);
  //label_check(enc, label);
  for(y = 0;y < height; y++){
	class_p = enc->class[y];
    for(x =0; x < width; x++){
      class_p[x] = ((char)(((int)label->value[y][x].x % 32) + ((int)label->value[y][x].y % 32)))*label->value[y][x].area;
	  if(label->value[y][x].area == 0) cent_pix++;
	  else edge_pix++;
	  edge_label[y][x] = label->value[y][x].area*100;
	  //label->value[y][x].area = 0;
	}
  }
  printf("vignetting_ratio = %f\n", (double)edge_pix/(double)cent_pix);
  make_edge_label(img, enc, edge_label, label);
  print_label(label, enc->height, enc->width);
  print_range(voronoi, NUM_OF_GEN, "outfile_range.txt");
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





VORONOI ***alloc_voronoi(void)
{
	printf("alloc_voronoi->");
	printf("MLADIST = %lf", MLADIST);
	  VORONOI ***voronoi, *vp;
	  int xg, yg;
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
	  shift_y = 1639.5 + sensor_offset_y;
	  shift_x = 1639.5 + sensor_offset_x;
	  voronoi = (VORONOI ***)alloc_3d_array(NUM_OF_GENy+RUN_OVER, NUM_OF_GENx+RUN_OVER, NUM_OF_VGEN, sizeof(VORONOI));

	  for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
		for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
		dblxg = (double)xg;
		dblyg = (double)yg;
        if(yg % 2 == 0){
        
          vp = &voronoi[yg][xg][0];
		   vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) + (dblxg * cent_dist_x/cos1) + sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2));
		   vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
		   vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		   vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          
          vp = &voronoi[yg][xg][1];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						+ (r /2.0);
			vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2 * r);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 1;
          
          vp = &voronoi[yg][xg][2];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						+ r;
            vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
            vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
			vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 2;
          
          vp = &voronoi[yg][xg][3];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						+ (r / 2.0);
            vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2 * r);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 3;
          
          vp = &voronoi[yg][xg][4];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						- (r / 2.0);
            vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2 * r);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 4;
          
          vp = &voronoi[yg][xg][5];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						- r;
            vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 5;
          
          vp = &voronoi[yg][xg][6];
            vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) 
						+ (dblxg * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2)) 
						- (r /2.0);
			vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2 * r);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		    vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 6;
        }
        if(yg % 2 == 1){
          vp = &voronoi[yg][xg][0];
		    vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) + sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2));
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          
          vp = &voronoi[yg][xg][1];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2.0 * r);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 1;
          
          vp = &voronoi[yg][xg][2];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ r;
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 2;
          
          vp = &voronoi[yg][xg][3];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2.0 * r);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 3;
          
          vp = &voronoi[yg][xg][4];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2.0 * r);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 4;
          
          vp = &voronoi[yg][xg][5];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- r;
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 5;
          
          vp = &voronoi[yg][xg][6];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2.0 * r);
			
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 6;
        }
      }
  }
	printf("ok\n");
	return(voronoi);
}



LABEL *alloc_label(int width, int height, int init)
{
    LABEL *label;
	printf("alloc_label->");
    label = (LABEL *)alloc_mem(sizeof(LABEL));
    label->width = width;
    label->height = height;
    label->value = (POINT **)alloc_2d_array(height, width, sizeof(POINT));
	printf("ok\n");
    return (label);
}


RANGE *alloc_range(int num)
{
    RANGE *range;

    range = (RANGE *)alloc_mem(sizeof(RANGE));
    range->vp = (VORONOI **)alloc_mem(sizeof(VORONOI *) * num);
    range->num = num;
    return (range);
}


int set_range(VORONOI *voronoi, LABEL *label, RANGE *range, ENCODER *enc, int id)
{
    int k;

    if (id < 0) {
  		range->tl.x = 0;
  		range->tl.y = 0;
  		range->br.x = enc->width  - 1;
  		range->br.y = enc->height - 1;
  		for (k = 0; k < range->num; k++) {
  			range->vp[k] = &voronoi[k];
  
  		}
    }
    return 0;
}


