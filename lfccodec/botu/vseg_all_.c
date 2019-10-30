#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <float.h>
#include "mrp1.h"

#include "vseg.h"

void fill_label_grobal_all(VORONOI ***, LABEL *);
void fill_recur_all(VORONOI ***, LABEL *);
int set_range(VORONOI *, LABEL *, RANGE *, ENCODER *, int );
RANGE *alloc_range(int );
LABEL *alloc_label(int , int , int );
VORONOI ***alloc_voronoi(void);
double comp_dist(int , int , double , double , VORONOI *, VORONOI *);
long comp_dist_cp(int , int , long , long , VORONOI *, VORONOI *);
LABEL *read_label(char *, ENCODER *);



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

/*クラスをpgmファイルで出力する関数*/
void print_class(char **class, int height, int width, char *filename)
{
    int i, j;
    FILE *fp;

    fp = fileopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
	     for (j = 0; j < width; j++) {
            putc(class[i][j] * 4, fp);
	     }
    }
    fclose(fp);
    return;
}

void print_area(IMAGE *img, LABEL *label, char *outfile)
{
	printf("***print_area***\n");
    int i, j;
    char name[100];
    FILE *fp;
	int **rorg, **gorg, **borg;
	rorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	gorg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	borg = (int **)alloc_2d_array(img->height, img->width, sizeof(int));
	int cent_pix=0, edge_pix=0;
	double ratio;
    sprintf(name, "./infofile/result/%s_area.ppm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", img->width, img->height);
    for (i = 0; i < img->height; i++) {
			for (j = 0; j < img->width; j++) {
				
				rorg[i][j] = gorg[i][j] = borg[i][j] = img->val[i][j];
				//if(label->value[i][j].area == 0) borg_p[0] = 0;
				//else rorg[i][j] = gorg[i][j] = borg[i][j] = 255;
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
					if(label->value[i][j].area == 1){
						edge_pix++;
						//borg[i][j] /= 2;
						//gorg[i][j] /= 2;
					}else if(label->value[i][j].area == 2){
						edge_pix++;
						//rorg[i][j] /= 2;
						//gorg[i][j] /= 2;
					}else if(label->value[i][j].area == 3){
						edge_pix++;
						//borg[i][j] /= 2;
						//rorg[i][j] /= 2;
					}else if(label->value[i][j].area == 4){
						edge_pix++;
						//gorg[i][j] /= 2;
					}else if(label->value[i][j].area == 5){
						edge_pix++;
						//rorg[i][j] /= 2;
					}else if(label->value[i][j].area == 6){
						borg[i][j] = gorg[i][j] = rorg[i][j] = 255;	
						//borg[i][j] /= 2;
						edge_pix++;
					}else{
						cent_pix++;
					}
					
		            putc(rorg[i][j], fp);//R
					putc(gorg[i][j], fp);//G
					putc(borg[i][j], fp);//B
					
			}
    }
	ratio = (double)edge_pix/(double)cent_pix;
	printf("edge_pix/cent_pix = %f\n", ratio);
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

void print_range(VORONOI ***voronoi, int num_of_gen, char *filename)
{
    FILE *fp;
    VORONOI *vp;
    int xg, yg, l;

    fp = fileopen(filename,"wb");
    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){//探索
      for(xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++){
        for(l = 0; l < NUM_OF_VGEN; l++){
          vp = &voronoi[yg][xg][l];
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



/***************************画像をボロノイ分割する関数**********************/
int main(int argc, char **argv)
{ //引数未定
  ENCODER *enc;
  char *infile = "outfile_label.txt"; //ラベル読み込むファイル
  char *outfile;
  char *testin=NULL, *testout=NULL;
  outfile = (char *)alloc_mem(sizeof(char));
  int x, y;
  IMAGE *img;
  enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
  
  
  //int k;
  testin = argv[1];
  testout = argv[2];
  VORONOI ***voronoi;
  LABEL *label;
  printf("a\n");
  img = read_pgm(testin);
  enc->height = img->height;
  enc->width = img->width;
  label = read_label(infile,enc); //ラベルファイル読み込む

  voronoi = alloc_voronoi();
  
  
  enc->class = (char **)alloc_2d_array(enc->height, enc->width,sizeof(char));
  printf("***testinfile -> %s***\n", testin);
  printf("***testoutfile_name -> %s***\n", testout);
  fill_label_grobal_all(voronoi, label); //すべての母点に対してラベリング


  /*for(k =0;k<NUM_OF_GEN;k++){
    x=(int)(voronoi[k].gen.x +0.5);
    y=(int)(voronoi[k].gen.y+0.5);
    if(x <0 || y<0)continue;
    enc->class[y][x] = 25;
  }*/
  for(y = 0;y < enc->height; y++){
    for(x =0; x < enc->width; x++){//ためしにclassにlabel代入してみる
      enc->class[y][x] = (char)(((int)label->value[y][x].x % 32) + ((int)label->value[y][x].y % 32) + (int)label->value[y][x].area);
    }
  }

  sprintf(outfile, "outfile_class_all.pgm");
	print_class(enc->class, enc->height, enc->width, outfile);

  sprintf(outfile, "outfile_label_all.txt");
  print_label(label, enc->height, enc->width, outfile);

  //sprintf(outfile, "outfile_range.txt");
  //print_range(voronoi, NUM_OF_GEN, outfile);
  print_area(img, label, testout);

  /*for(y = 0;y < 14; y++){
      printf("%f %f\n",voronoi[y].gen.x,voronoi[y].gen.y);

  }*/

  free(label);
  free(voronoi);
  free(enc);
  free(outfile);
  return 0;
}




//voronoiのメモリ確保　４隅の点から母点配置 引数未定void?
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
		    vp->gen.x = tlx + ((dblxg - 1.0) * cent_dist_x/cos1) 
						+ sin3 * (tly + ((dblyg - 1.0) * cent_dist_y*cos1*cos3/cos2));
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
void fill_recur_all(VORONOI ***voronoi, LABEL *label)
{
    int x, y;
    int xg, yg, l;
    double dx,dy;
    int height, width;
	height = label->height;
	width = label->width;
    double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp;
    printf("c\n");
    printf("h,w = %d,%d\n",label->height, label->width);
	ivp = &voronoi[0][0][0];
	vp = &voronoi[0][0][0];
	for (y = 0; y < height; y++) { //label範囲のすべての画素について
  		for (x = 0; x < width; x++) {
			if(label->value[y][x].area == 0) continue;
			xg = (int)label->value[y][x].x;
			yg = (int)label->value[y][x].y;
			min_dist2 = DBL_MAX;
			min_dist = DBL_MAX;
			/* searching */
			ivp = &voronoi[0][0][0];
			for (l = 1;  l < NUM_OF_VGEN; l++) { //マイクロレンズ内の７つの母点でボロノイ分割
				vp = &voronoi[yg][xg][l];
				dist2 = comp_dist(x, y, min_dist, min_dist2, vp, ivp); //x,yと最も近い母点との距
				if (dist2 >= 0) {
					ivp = vp;
					min_dist2 = dist2;
					min_dist = (sqrt((double)min_dist2));
				}
			}
			label->value[y][x].area = ivp->id.area; //labelのvalueをidに（最も近い母点のid）
			ivp->num_pixel++; //母点に対して近い画素一つ決まったからsize++　つまりその母点にラベリングされた画素の数
			ivp->tl.x = MIN(ivp->tl.x, x); //母点にラベリングされた画素を含む最小の大きさの四角形の範囲
			ivp->tl.y = MIN(ivp->tl.y, y);
			ivp->br.x = MAX(ivp->br.x, x);
			ivp->br.y = MAX(ivp->br.y, y);
  		}
		
    }

    return;
}

//レンジ、サイズ、ラベル初期化して再設定
void fill_label_grobal_all(VORONOI ***voronoi, LABEL *label)
{
    int xg, yg, l;
    VORONOI *vp;

    for(yg = 0; yg < NUM_OF_GENy +RUN_OVER; yg++){
      for (xg = 0; xg < NUM_OF_GENx +RUN_OVER; xg++) {
        for(l = 0; l < 7; l++){
          vp = &voronoi[yg][xg][l];
      		vp->tl.x = vp->br.x = vp->gen.x; //初期値としてRANGEのtlとbrを母点に合わせる
      		vp->tl.y = vp->br.y = vp->gen.y;
      		vp->num_pixel = 0; //初期値としてボロノイのサイズ全部　0にする
        }
      }
    }

    /*for (y = 0; y < label->height; y++) { //LABELの範囲のvalueを-1に　初期値？
  		for (x = 0; x < label->width; x++) {
  			label->value[y][x] = -1;
  		}
    }*/

		    fill_recur_all(voronoi, label);


    return;
}

double comp_dist(int x, int y, double min_dist, double min_dist2,
	      VORONOI *vpt, VORONOI *vps)
{  //vpsが参照してる母点 vptを動かす
    double dx, dy;
    double dist2;

    if (vpt == vps) return(-1.0);
    dx = DIFFPG((double)x, vpt->gen.x);
    if (dx > min_dist) return (-1.0);
    dy = DIFFPG((double)y, vpt->gen.y);
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

long comp_dist_cp(int x, int y, long min_dist, long min_dist2,
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
}

LABEL *read_label(char *filename, ENCODER *enc)
{
  LABEL *label;
  FILE *fp;
  int i, j;
  int x, y, area;

  fp = fileopen(filename, "rb");
  label = alloc_label(enc->width, enc->height, 0);

  for (i = 0; i < enc->height; i++) {
    for (j = 0; j < enc->width; j++) {
      fscanf(fp, "%d\n", &x);
      fscanf(fp, "%d\n", &y);
      fscanf(fp, "%d\n", &area);
      label->value[i][j].x = (double)x;
      label->value[i][j].y = (double)y;
	  label->value[i][j].area = area;
    }
  }

  fclose(fp);
  return (label);
}
