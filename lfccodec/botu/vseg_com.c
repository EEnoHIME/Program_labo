#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "mrp1.h"
#include <float.h>
#include "vseg.h"

void fill_label_grobal_all(RBS **, VORONOI ***, VORONOI ***, LABEL *);
void fill_recur_all(VORONOI ***, VORONOI ***, LABEL *);
int set_range(VORONOI *, LABEL *, RANGE *, ENCODER *, int );
RANGE *alloc_range(int );
LABEL *alloc_label(int , int , int );
VORONOI ***alloc_voronoi(void);
VORONOI ***alloc_voronoi2(void);
double comp_dist(double , double , double , double , VORONOI *, VORONOI *);
long comp_dist_cp(int , int , long , long , VORONOI *, VORONOI *);
inline double comp_dist_cp2(int, int, double, double, double, double, double, double);
LABEL *read_label(char *, int, int);



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
}
*/
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


/*クラスをpgmファイルで出力する関数*/
void print_class(char **class, int height, int width, char *filename)
{
    int i, j;
    FILE *fp;

    fp = fileopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
            putc(class[i][j] * 2, fp);
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

void print_label(LABEL *label, int height, int width, char *filename)
{
    int i, j;
    FILE *fp;

    fp = fileopen(filename, "wb");
    //fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
         fprintf(fp, "%d\n",(int)(label->value[i][j].x));
         fprintf(fp, "%d\n",(int)(label->value[i][j].y));
         fprintf(fp, "%d\n",(int)(label->value[i][j].area));
            //putc( (char)(label->value[i][j]), fp); //ここ直さないと
	     }
    }
    fclose(fp);
    return;
}

void print_range(VORONOI ***voronoi, VORONOI ***voronoi2, int num_of_gen, char *filename)
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
    for(l = 1; l < NUM_OF_VGEN; l++){
      for(yg = 0; yg < NUM_OF_GEN_ey; yg++){//探索
        for(xg = 0; xg < NUM_OF_GEN_ex; xg++){
            vp = &voronoi2[yg][xg][l];
            fprintf(fp, "%d\n", (int)vp->tl.x);
            fprintf(fp, "%d\n", (int)vp->tl.y);
            fprintf(fp, "%d\n", (int)vp->br.x);
            fprintf(fp, "%d\n", (int)vp->br.y);
        }
      }
    }
    fclose(fp);
    return;
}

RBS **roughtly_belongs_search_com(int height, int width, VORONOI ***voronoi2, LABEL *label)
{
	int xg, yg, y, x, count;
	VORONOI *vp;
	RBS **rbs;
	int start_y, start_x, end_y, end_x, search_mask;
	search_mask = 20;
	rbs = (RBS **)alloc_2d_array(height, width, sizeof(RBS));
	for(y = 0; y < height; y++){
		for(x = 0; x < width; x++){
			rbs[y][x].count = 0;
		}
	}
	for(yg = 0; yg < NUM_OF_GEN_ey; yg++){
		for(xg = 0; xg < NUM_OF_GEN_ex; xg++){
			vp = &voronoi2[yg][xg][0];
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
					if(label->value[y][x].area == 0) continue;
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

void fill_recur_all3(RBS **rbs, VORONOI ***voronoi2, LABEL *label)
{
	int x, y, height, width, rbsidx, rbsidy, i, count;
	double min_dist, min_dist2, dist2;
	VORONOI *vp, *ivp;
	height = label->height;
	width = label->width;
	int xg, yg, ar;
	printf("h,w=%d,%d\n",height, width);
	for(y = 0; y < height; y++){
		//printf("y = %d\n",y);
		for(x = 0; x < width; x++){
			if(label->value[y][x].area == 0) continue;
			min_dist2 = DBL_MAX;
			min_dist = DBL_MAX;
			count = rbs[y][x].count;
			if(rbs[y][x].count == 1){
				label->value[y][x].x = rbs[y][x].id[0].x;
				label->value[y][x].y = rbs[y][x].id[0].y;
			}else{
				vp = &voronoi2[0][0][0];
				ivp = vp;
				for(i = 0; i < count; i++){
					rbsidy = (int)rbs[y][x].id[i].y;
					rbsidx = (int)rbs[y][x].id[i].x;
					vp = &voronoi2[rbsidy][rbsidx][0];
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
			ivp = &voronoi2[(int)label->value[y][x].y][(int)label->value[y][x].x][0];
			ivp->tl.x = MIN(ivp->tl.x, (double)x);
			ivp->tl.y =	MIN(ivp->tl.y, (double)y);
            ivp->br.x = MAX(ivp->br.x, (double)x);
            ivp->br.y = MAX(ivp->br.y, (double)y);
		}
	}
	for(yg = 0; yg < NUM_OF_GEN_ey; yg++){
		for(xg = 0; xg < NUM_OF_GEN_ex; xg++){
			for(ar = 1; ar < NUM_OF_VGEN; ar++){
				ivp = &voronoi2[yg][xg][ar];
				vp = &voronoi2[yg][xg][0];
				ivp->tl.x = vp->tl.x;
				ivp->tl.y = vp->tl.y;
				ivp->br.x = vp->br.x;
				ivp->br.y = vp->br.y;
			}
		}
	}
}

/***************************画像をボロノイ分割する関数**********************/
int main(int argc, char **argv)
{ //引数未定
  //ENCODER *enc;
  char *infile = "outfile_label_all.txt"; //ラベル読み込むファイル
  char *outfile;
  char *testin=NULL, *testout=NULL;
  char **class;
  outfile = (char *)alloc_mem(sizeof(char));
  int x, y, height, width;
  IMAGE *img;
  //enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
  //enc->height = 3280;
  //enc->width = 3280;
  
  //int k;
  VORONOI ***voronoi, ***voronoi2;
  LABEL *label;
  RBS **rbs;
  
  testin = argv[1];
  testout = argv[2];
  img = read_pgm(testin);
  height = img->height;
  width = img->width;
  printf("main h,w=%d,%d\n",height,width);
  class = (char **)alloc_2d_array(height, width,sizeof(char));
  printf("a\n");
  label = read_label(infile, height, width); //ラベルファイル読み込む
  voronoi = alloc_voronoi();
  voronoi2 = alloc_voronoi2();
  rbs = roughtly_belongs_search_com(height, width, voronoi2, label);
  printf("b\n");

  fill_label_grobal_all(rbs, voronoi, voronoi2, label); //すべての母点に対してラベリング


  for(y = 0;y < height; y++){
    for(x =0; x < width; x++){//ためしにclassにlabel代入してみる
      class[y][x] = (char)(((int)label->value[y][x].x * 10 % 60) + ((int)label->value[y][x].y * 10 % 60));
    }
  }

  sprintf(outfile, "outfile_class_all_re.pgm");
	print_class(class, height, width, outfile);

  sprintf(outfile, "outfile_label_all_re.txt");
  print_label(label, height, width, outfile);

  sprintf(outfile, "outfile_range_re.txt");
  print_range(voronoi, voronoi2, NUM_OF_GEN, outfile);
  print_lens(img, label, testout);
  printf("a\n");

  /*for(y = 0;y < 14; y++){
      printf("%f %f\n",voronoi[y].gen.x,voronoi[y].gen.y);

  }*/
  if (label->value != NULL) {
    free(label->value);
    label->value = NULL;
  }
  printf("a\n");
  if(class != NULL) {
    free(class);
    class = NULL;
  }
  printf("a\n");
  if (label != NULL) {
    free(label);
    label = NULL;
  }
  printf("a\n");
  if (voronoi != NULL) {
    free(voronoi);
    voronoi = NULL;
  }
  printf("a\n");
  if (voronoi2 != NULL) {
    free(voronoi2);
    voronoi2 = NULL;
  }
  printf("a\n");
  //if (enc != NULL) {
  //  free(enc);
    //enc = NULL;
  //}
  printf("a\n");
  if (outfile != NULL) {
    free(outfile);
    outfile = NULL;
  }
  free(rbs);
  //free(label->value);
  //free(enc->class);
  //free(label);
  //free(voronoi);
  //free(voronoi2);
  //free(enc);
//  free(outfile);

  return 0;
}




//voronoiのメモリ確保　９つの領域を1つの領域とするときの母点配置
VORONOI ***alloc_voronoi2(void)
{
  VORONOI ***voronoi, *vp;
  int xg, yg, area;
  double tlx = TOP_LEFT_x, trx = TOP_RIGHT_x, blx = BOTTOM_LEFT_x, brx = BOTTOM_RIGHT_x;
  double tly = TOP_LEFT_y, try = TOP_RIGHT_y, bly = BOTTOM_LEFT_y, bry = BOTTOM_RIGHT_y;
  double r = RADIUS;
  double shift_y, shift_x;
  double sensor_offset_y, sensor_offset_x;
  double rd;
  double dblyg, dblxg;
  double cos1, sin1, cos2, cos3, sin3;
  double cent_dist_x, cent_dist_y;
  double rotated_bly, rotated_blx, rotated_bry, rotated_brx;
  double slope_a,slope_b,slope_c,slope_d;
  rd = ROTATION;
  sensor_offset_y = SENSOR_OFFSET_y;
  sensor_offset_x = SENSOR_OFFSET_x;
  shift_y = 1639.5 + sensor_offset_y;
  shift_x = 1639.5 + sensor_offset_x;
  voronoi = (VORONOI ***)alloc_3d_array(NUM_OF_GEN_ey, NUM_OF_GEN_ex, NUM_OF_VGEN, sizeof(VORONOI));//メモリ確保
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
  
  /*母点の座標すべて設定*/
  for(yg = 0; yg < NUM_OF_GEN_ey; yg++){ //中心の母点の数だけ
      for(xg = 0; xg < NUM_OF_GEN_ex; xg++){ //母点内について
	    dblyg = (double)yg;
		dblxg = (double)xg;
        if(yg % 2 == 0){ //偶数行 x座標左にずらす
          for(area = 0; area < 7; area++){
            vp = &voronoi[yg][xg][area];
            vp->gen.x = tlx + (dblxg * cent_dist_x / cos1 * 3.0) 
						+ sin3 * (tly + (dblyg * cent_dist_y * cos1 * cos3 / cos2 * 3.0));
            vp->gen.y = tly + (dblyg * cent_dist_y * cos1 * cos3 / cos2 * 3.0);
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = area;
          }
        }
        if(yg % 2 == 1){ //奇数行
          for(area = 0; area < 7; area++){
            vp = &voronoi[yg][xg][area];
            vp->gen.x = tlx - (cent_dist_x / cos1 * 1.5) + (dblxg * cent_dist_x / cos1 * 3.0)
						+ sin3 * (tly + (dblyg * cent_dist_y *cos1 *cos3 / cos2 * 3.0));
            vp->gen.y = tly + (dblyg * cent_dist_y *cos1 *cos3 / cos2 * 3.0);
            vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
			vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = area;
          }

        }
      }
  }

  return(voronoi);
}

//周辺部の領域上に母点を配置
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
  //rd = ROTATION * M_PI / 180.0;
  rd = ROTATION;
  //rd = 0.0;
  double cos1, sin1, cos2, cos3, sin3;
  double cent_dist_x, cent_dist_y;
  double rotated_bly, rotated_blx, rotated_bry, rotated_brx;
  double slope_a,slope_b,slope_c,slope_d;
  sensor_offset_y = SENSOR_OFFSET_y;
  sensor_offset_x = SENSOR_OFFSET_x;
  shift_y = 1639.5 + sensor_offset_y;
  shift_x = 1639.5 + sensor_offset_x;
  voronoi = (VORONOI ***)alloc_3d_array(NUM_OF_GENy +RUN_OVER, NUM_OF_GENx +RUN_OVER, NUM_OF_VGEN, sizeof(VORONOI));//メモリ確保 xy座標は少し広めにとってる
  
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
  /*母点の座標すべて設定*/
  for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){ //中心の母点の数だけ
      for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){ //母点内について
		dblxg = (double)xg;
		dblyg = (double)yg;
        if(yg % 2 == 0){ //偶数行 x座標左にずらす
          //中心
          vp = &voronoi[yg][xg][0];
		   vp->gen.x = tlx - (cent_dist_x/cos1 / 2.0) + (dblxg * cent_dist_x/cos1) + sin3 * (tly + ((dblyg -1.0) * cent_dist_y*cos1*cos3/cos2));
		   vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
		   //逆回転後↓
		   vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
		   vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          //右上
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
          //右
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
          //右下
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
          //左下
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
          //左
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
          //左上
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
        if(yg % 2 == 1){ //奇数行
          //中心
          vp = &voronoi[yg][xg][0];
		    vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) + sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2));
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          //右上
          vp = &voronoi[yg][xg][1];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2.0 * r);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly; 
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 1;
          //右
          vp = &voronoi[yg][xg][2];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ r;
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 2;
          //右下
          vp = &voronoi[yg][xg][3];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						+ (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2.0 * r);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 3;
          //左下
          vp = &voronoi[yg][xg][4];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						+ (sqrt(3.0) / 2.0 * r);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 4;
          //左
          vp = &voronoi[yg][xg][5];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- r;
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2);
			//逆回転後↓
			vp->gen.x = (vp->gen.x - tlx) * cos1 + (vp->gen.y - tly) * sin1 + tlx;
			vp->gen.y = -(vp->gen.x - tlx) * sin1 + (vp->gen.y - tly) * cos1 + tly;
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 5;
          //左上
          vp = &voronoi[yg][xg][6];
            vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2))
						- (r /2.0);
		    vp->gen.y = tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2)
						- (sqrt(3.0) / 2.0 * r);
			//逆回転後↓
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




//labelのメモリ確保　値代入
LABEL *alloc_label(int width, int height, int init)
{
    LABEL *label;

    label = (LABEL *)alloc_mem(sizeof(LABEL));
    label->width = width;
    label->height = height;
    label->value = (POINT **)alloc_2d_array(height, width, sizeof(POINT));

    return (label);
}
//rangeのメモリ確保　numを全母点数に
RANGE *alloc_range(int num)
{
    RANGE *range;

    range = (RANGE *)alloc_mem(sizeof(RANGE));
    range->vp = (VORONOI **)alloc_mem(sizeof(VORONOI *) * num);
    range->num = num;
    return (range);
}

//rangeをとりあえず画像の大きさに
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
	/*	    printf("%d:%d\n", k, range->vp[k]->id);*/
  		}
    } //elseの部分はよく分からなかったので、省略してる
    return 0;
}

//画素ごとにラベリング
/*void fill_recur_all(VORONOI ***voronoi, VORONOI ***voronoi2, LABEL *label)
{
    int x, y;
    int xg, yg, area;
    double dx,dy;
    double gex, gey;
    //int id;
    double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp;

    printf("c\n");
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
		for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
			for(area = 1; area < NUM_OF_VGEN; area++){
				vp = &voronoi[yg][xg][area];
				gex = vp->gen.x;
				gey = vp->gen.y;
				ivp = &voronoi2[0][0][area];//ivpの初期化
				dx = DIFFPG(gex,ivp->gen.x); //x - gen.xの絶対値
    			dy = DIFFPG(gey,ivp->gen.y); //今参照してる母点とlabel範囲の点との距離
    			min_dist2 = dx*dx + dy*dy;//2乗の和
    			min_dist = (double)(sqrt((double)min_dist2) + 1.0);
				for (y = 0; y < NUM_OF_GEN_ey; y++) { //label範囲のすべての画素について
					for (x = 0; x < NUM_OF_GEN_ex; x++) {
						vp = &voronoi2[y][x][area];
						dist2 = comp_dist(gex, gey, min_dist, min_dist2, vp, ivp);
						if (dist2 >= 0) {
							ivp = vp;
							min_dist2 = dist2;
							min_dist = (double)(sqrt((double)min_dist2) + 1.0);
						}
					}
				}
				vp = &voronoi[yg][xg][area];

				ivp->tl.x = MIN(ivp->tl.x, vp->tl.x); //母点にラベリングされた画素を含む最小の大きさの四角形の範囲
				ivp->tl.y = MIN(ivp->tl.y, vp->tl.y);
				ivp->br.x = MAX(ivp->br.x, vp->br.x);
				ivp->br.y = MAX(ivp->br.y, vp->br.y);
				for (y = vp->tl.y; y < vp->br.y + 1; y++) { //label範囲のすべての画素について
					for (x = vp->tl.x; x < vp->br.x + 1; x++) {
						if(label->value[y][x].x == xg && label->value[y][x].y == yg && label->value[y][x].area == area){
							label->value[y][x].x = ivp->id.x;
							label->value[y][x].y = ivp->id.y;
						}
					}
				}

			}
		}
    }

    return;
}*/

void fill_recur_all2(VORONOI ***voronoi2, LABEL *label)
{
	int x, y, l;
    int xg, yg;
    //int id;
    double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp;
	double tlx = TOP_LEFT_x, trx = TOP_RIGHT_x, blx = BOTTOM_LEFT_x, brx = BOTTOM_RIGHT_x;
	//double tly = TOP_LEFT_y, try = TOP_RIGHT_y, bly = BOTTOM_LEFT_y, bry = BOTTOM_RIGHT_y;
	double cent_dist_x = ((trx - tlx)/NUM_OF_GENx+(brx - blx)/NUM_OF_GENx)/ 2.0;
	double absolute_flag = cent_dist_x * 3.0 * 0.4;
	printf("absolute_flag = %f\n", absolute_flag);
	int height, width;
	height = label->height;
	width = label->width;
	l=0;
	printf("height, width = %d, %d\n", height, width);
    for (y = 0; y < height; y++) { //label範囲のすべての画素について
	printf("%d\n",(int)y);
  		for (x = 0; x < width; x++) {
			if(label->value[y][x].area == 0) continue;
  			//if (label->value[y][x] >= 0) continue; //ラベル決まってたら以下の処理スキップでforループ繰り返し
		//printf("%d ", x);
        min_dist2 = DBL_MAX;
        min_dist = DBL_MAX;
		
		vp = &voronoi2[0][0][l];
		ivp = vp;
		//printf("num_yg, num_xg = %d,%d", NUM_OF_GENy, NUM_OF_GENx);
		
			for(yg = 0; yg < NUM_OF_GEN_ey; yg++){//探索
				if(min_dist <= absolute_flag) break;
				for(xg = 0; xg < NUM_OF_GEN_ex; xg++){
					vp = &voronoi2[yg][xg][(int)label->value[y][x].area];
					if(fabs((double)y - vp->gen.y) > 25.0 || fabs((double)x - vp->gen.x) > 25.0){
						continue;//xy座標が10より大きかったらスキップ
					}
					
					//dist2 = comp_dist_cp(x, y, min_dist, min_dist2, vp, ivp); //
					dist2 = comp_dist((double)x, (double)y, min_dist, min_dist2, vp, ivp); //
					//printf("E ");
					if (dist2 >= 0.0){
						ivp = vp;
						min_dist2 = dist2;
						min_dist = sqrt((double)min_dist2);
						if(min_dist <= absolute_flag) break;
					}
				}
			}
			
  			label->value[y][x].x = ivp->id.x; //labelのvalueをidに（最も近い母点のid）
			label->value[y][x].y = ivp->id.y; //labelのvalueをidに（最も近い母点のid）
			//label->value[y][x].area = ivp->id.area; //labelのvalueをidに（最も近い母点のid）
  			//ivp->num_pixel++; //母点に対して近い画素一つ決まったからsize++　つまりその母点にラベリングされた画素の数
            ivp->tl.x = MIN(ivp->tl.x, (double)x); //母点にラベリングされた画素を含む最小の大きさの四角形の範囲
            ivp->tl.y =	MIN(ivp->tl.y, (double)y);
            ivp->br.x = MAX(ivp->br.x, (double)x);
            ivp->br.y = MAX(ivp->br.y, (double)y);
			//printf("E");
  		}
	
		
    }
	
}

//レンジ、サイズ、ラベル初期化して再設定
void fill_label_grobal_all(RBS **rbs, VORONOI ***voronoi, VORONOI ***voronoi2, LABEL *label)
{
    int xg, yg, area;
    int tlx, tly, brx, bry;
    FILE *fp;
    fp = fileopen("outfile_range.txt","rb");
    VORONOI *vp;
   printf("ly, lx = %d,%d\n",label->height, label->width);
    for(yg = 0; yg < NUM_OF_GEN_ey; yg++){
      for (xg = 0; xg < NUM_OF_GEN_ex; xg++) {
        for(area = 0; area < 7; area++){
          vp = &voronoi2[yg][xg][area];
      		vp->tl.x = label->width-1;
			vp->tl.y = label->height-1;
      		vp->br.x = vp->br.y = 0;
      		vp->num_pixel = 0; //初期値としてボロノイのサイズ全部　0にする
        }
      }
    }
	area=0;
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
      for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
        //for(area = 0; area < 7; area++){
          vp = &voronoi[yg][xg][area];
          fscanf(fp, "%d\n", &tlx);
          fscanf(fp, "%d\n", &tly);
          fscanf(fp, "%d\n", &brx);
          fscanf(fp, "%d\n", &bry);
          vp->tl.x = (double)tlx;
          vp->tl.y = (double)tly;
          vp->br.x = (double)brx;
      		vp->br.y = (double)bry;

          if(vp->tl.x < 0){vp->tl.x = 0;
					}else if(vp->tl.x > label->width - 1){vp->tl.x = label->width - 1;}
					if(vp->tl.y < 0) {vp->tl.y = 0;
					}else if(vp->tl.y > label->height - 1){vp->tl.y = label->height - 1;}
					if(vp->br.x < 0) {vp->br.x = 0;
					}else if(vp->br.x > label->width - 1){vp->br.x = label->width - 1;}
					if(vp->br.y < 0) {vp->br.y = 0;
					}else if(vp->br.y > label->height - 1){vp->br.y = label->height - 1;}
      		vp->num_pixel = 0; //初期値としてボロノイのサイズ全部　0にする
        //}
      }
    }
	printf("e\n");
	//fill_recur_all(voronoi, voronoi2, label);
	fill_recur_all3(rbs, voronoi2, label);
	printf("grobal->ok\n");
    fclose(fp);
    return;
}

double comp_dist(double x, double y, double min_dist, double min_dist2,
	      VORONOI *vpt, VORONOI *vps)
{  //vpsが参照してる母点 vptを動かす
    double dx=0.0, dy=0.0;
    double dist2=0.0;

    if (vpt == vps) return(-1.0);
    dx = DIFFPG(x, vpt->gen.x);
    if (dx > min_dist) return (-1.0);
    dy = DIFFPG(y, vpt->gen.y);
    if (dy > min_dist) return (-1.0);
    dist2 = (dx*dx + dy*dy);
    if (dist2 > min_dist2) return(-1.0);
    if (dist2 < min_dist2) return(dist2);
    if (vpt->gen.y > vps->gen.y) return(-1.0);
    if (vpt->gen.y < vps->gen.y) return(dist2);
    if (vpt->gen.x > vps->gen.x) return(-1.0);
    if (vpt->gen.x < vps->gen.x) return(dist2);
    //if (vpt->id > vps->id) return(-1L);
    return (dist2);
}

/*long comp_dist_cp(int x, int y, long min_dist, long min_dist2,
	      VORONOI *vpt, VORONOI *vps)
{  //vpsが参照してる母点 vptを動かす
    double dx, dy;
    long double dist2;

    dx = DIFFPG(x, vpt->gen.x);
    if (dx > min_dist) return (-1L);
    dy = DIFFPG(y, vpt->gen.y);
    if (dy > min_dist) return (-1L);
    dist2 = (dx*dx + dy*dy);
    if (dist2 > min_dist2) return(-1L);
    if (dist2 < min_dist2) return(dist2);
    if (vpt->gen.y > vps->gen.y) return(-1L);
    if (vpt->gen.y < vps->gen.y) return(dist2);
    if (vpt->gen.x > vps->gen.x) return(-1L);
    if (vpt->gen.x < vps->gen.x) return(dist2);
    //if (vpt->id > vps->id) return(-1L);
    return (dist2);
}*/

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
    if (dist2 < min_dist2) return(dist2);
	if (vptgeny > vpsgeny) return(-1.0);
	if (vptgeny < vpsgeny) return(dist2);
	if (vptgenx > vpsgenx) return(-1.0);
	if (vptgenx < vpsgenx) return(dist2);
    
    return (dist2);
}

LABEL *read_label(char *filename, int height, int width)
{
  LABEL *label;
  FILE *fp;
  int i, j;
  int x, y, area;

  fp = fileopen(filename, "rb");
  label = alloc_label(width, height, 0);

  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      fscanf(fp, "%d\n", &x);
      fscanf(fp, "%d\n", &y);
      fscanf(fp, "%d\n", &area);
      label->value[i][j].x = (double)x;
      label->value[i][j].y = (double)y;
      label->value[i][j].area = (double)area;
	  //printf("area = %d\n",label->value[i][j].area);
    }
  }

  fclose(fp);
  return (label);
}
