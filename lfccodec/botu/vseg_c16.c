#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "mrp.h"

#include "vseg.h"

void fill_label_grobal_all(VORONOI ***, VORONOI ***, LABEL *);
void fill_recur_all(VORONOI ***, VORONOI ***, LABEL *);
int set_range(VORONOI *, LABEL *, RANGE *, ENCODER *, int );
RANGE *alloc_range(int );
LABEL *alloc_label(int , int , int );
VORONOI ***alloc_voronoi(void);
VORONOI ***alloc_voronoi2(void);
VORONOI **alloc_voronoi_center_com(int yg, int xg);
long comp_dist(int , int , long , long , VORONOI *, VORONOI *);
long comp_dist_cp(int , int , long , long , VORONOI *, VORONOI *);
void compo_central16(VORONOI **, LABEL *);
LABEL *read_label(char *);
VORONOI ***read_range(char *filename);



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

LABEL *read_label(char *filename)
{
  LABEL *label;
  FILE *fp;
  int i, j;
  int x, y, area;

  fp = fileopen(filename, "rb");
  label = alloc_label(3280, 3280, 0);

  for (i = 0; i < 3280; i++) {
    for (j = 0; j < 3280; j++) {
      fscanf(fp, "%d\n", &x);
      fscanf(fp, "%d\n", &y);
      fscanf(fp, "%d\n", &area);
      label->value[i][j].x = (double)x;
      label->value[i][j].y = (double)y;
      label->value[i][j].area = (double)area;
    }
  }

  fclose(fp);
  return (label);
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

void print_range(VORONOI **voronoi, VORONOI ***voronoi2, int num_of_gen, char *filename)
{
  printf("print_range start.\n");
  FILE *fp;
  VORONOI *vp;
  int xg, yg, l;

  fp = fileopen(filename,"wb");
  printf("file opened.\n");
	for(yg = 0; yg < 96; yg++){
		for(xg = 0; xg < 83; xg++){
      vp = &voronoi[yg][xg];
      fprintf(fp, "%d\n", (int)vp->tl.x);
      fprintf(fp, "%d\n", (int)vp->tl.y);
      fprintf(fp, "%d\n", (int)vp->br.x);
      fprintf(fp, "%d\n", (int)vp->br.y);
    }
  }
  printf("area0 printed\n");
  for(l = 1; l < NUM_OF_VGEN; l++){
    for(yg = 0; yg < 128; yg++){//探索
      for(xg = 0; xg < 111; xg++){
        vp = &voronoi2[yg][xg][l];
        fprintf(fp, "%d\n", (int)vp->tl.x);
        fprintf(fp, "%d\n", (int)vp->tl.y);
        fprintf(fp, "%d\n", (int)vp->br.x);
        fprintf(fp, "%d\n", (int)vp->br.y);
      }
    }
  }
  printf("area1~6 printed.\n");
  fclose(fp);
  return;
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
RANGE *alloc_range(int num)//この関数は使われていません
{
    RANGE *range;

    range = (RANGE *)alloc_mem(sizeof(RANGE));
    range->vp = (VORONOI **)alloc_mem(sizeof(VORONOI *) * num);
    range->num = num;
    return (range);
}

//rangeをとりあえず画像の大きさに
int set_range(VORONOI *voronoi, LABEL *label, RANGE *range, ENCODER *enc, int id)//この関数は使われていません
{
    int k;

    if (id < 0) {
      range->tl.x = 0;
      range->tl.y = 0;
      range->br.x = enc->width  - 1;
      range->br.y = enc->height - 1;
      for (k = 0; k < range->num; k++) {
        range->vp[k] = &voronoi[k];
  /*      printf("%d:%d\n", k, range->vp[k]->id);*/
      }
    } //elseの部分はよく分からなかったので、省略してる
    return 0;
}

VORONOI ***read_range(char *filename){

  FILE *fp2;
  int area, xg, yg;
  int height, width;
  height = width = 3280;
  VORONOI ***vp, **vp_s, *vpa;
  vp = (VORONOI ***)alloc_3d_array(128, 111, NUM_OF_VGEN, sizeof(VORONOI));
  vp_s = (VORONOI **)alloc_2d_array(NUM_OF_GENy + 2, NUM_OF_GENx + 2, sizeof(VORONOI));

  fp2 = fileopen(filename, "rb");

  area = 0;
  for(yg = 0; yg < NUM_OF_GENy + 2; yg++){
    for(xg = 0; xg < NUM_OF_GENx + 2; xg++){
      fscanf(fp2, "%d\n", &vp_s[yg][xg].tl.x);
      fscanf(fp2, "%d\n", &vp_s[yg][xg].tl.y);
      fscanf(fp2, "%d\n", &vp_s[yg][xg].br.x);
      fscanf(fp2, "%d\n", &vp_s[yg][xg].br.y);
    }
  }
  for(area = 1; area < 7; area++){
    for(yg = 0; yg < 128; yg++){
      for(xg = 0; xg < 111; xg++){
        vpa = &vp[yg][xg][area];
        fscanf(fp2, "%lf\n", &vpa->tl.x);
        fscanf(fp2, "%lf\n", &vpa->tl.y);
        fscanf(fp2, "%lf\n", &vpa->br.x);
        fscanf(fp2, "%lf\n", &vpa->br.y);
        if(vpa->tl.x < 0){
          vpa->tl.x = 0;
        }else if(vpa->tl.x > width - 1){
          vpa->tl.x = width - 1;
        }
        if(vpa->tl.y < 0){
          vpa->tl.y = 0;
        }else if(vpa->tl.y > height - 1){
          vpa->tl.y = height - 1;
        }
        if(vpa->br.x < 0){
          vpa->br.x = 0;
        }else if(vpa->br.x > width - 1){
          vpa->br.x = width - 1;
        }
        if(vpa->br.y < 0){
          vpa->br.y = 0;
        }else if(vpa->br.y > height - 1){
          vpa->br.y = height - 1;
        }
      }
    }
  }

  fclose(fp2);
  return(vp);
}

/***************************画像をボロノイ分割する関数**********************/
int main(void)
{ //引数未定
  //ENCODER *enc;
  char *labelfile = "outfile_label_0007_c1s9.txt"; //ラベル読み込むファイル
  char *rangefile = "outfile_range_0007_c1s9.txt";
  char *outfile;
  char **class;
  outfile = (char *)alloc_mem(sizeof(char));
  int x, y;

  //enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
  //enc->height = 3280;
  //enc->width = 3280;
  class = (char **)alloc_2d_array(3280, 3280, sizeof(char));
  //int k;

  VORONOI ***voronoi;
  VORONOI **voronoi_center_com;//中心部まとめる用
  LABEL *label;
  label = read_label(labelfile); //ラベルファイル読み込む
  voronoi = read_range(rangefile);

  printf("label read\n");

  //voronoi = alloc_voronoi();
  //voronoi2 = alloc_voronoi2();
  voronoi_center_com = alloc_voronoi_center_com(96, 83);
  //printf("voronoi seeds plotted.\n");


  //fill_label_grobal_all(voronoi, voronoi2, label); //すべての母点に対してラベリング

  //printf("label filled.\n");

  compo_central16(voronoi_center_com, label);

  printf("central16 compositted.\n");


  for(y = 0;y < 3280; y++){
    for(x =0; x < 3280; x++){//ためしにclassにlabel代入してみる
      class[y][x] = (char)(((int)label->value[y][x].x * 10 % 60) + ((int)label->value[y][x].y * 10 % 60));
    }
  }

  printf("class[][] completed.\n");
  sprintf(outfile, "outfile_class_cen16.pgm");
	print_class(class, 3280, 3280, outfile);
  printf("class ok\n");
  sprintf(outfile, "outfile_label_cen16.txt");
  print_label(label, 3280, 3280, outfile);
  printf("class ok\n");
  sprintf(outfile, "outfile_range_cen16.txt");
  print_range(voronoi_center_com, voronoi, NUM_OF_GEN, outfile);
  printf("range ok\n");

  /* それぞれの領域の画素数求めるプログラム
  int num_area[7] = {0};

  for(y=0;y<3280;y++){
    for(x=0;x<3280;x++){
      if(label->value[y][x].area == 0){num_area[0]++;}
      if(label->value[y][x].area == 1){num_area[1]++;}
      if(label->value[y][x].area == 2){num_area[2]++;}
      if(label->value[y][x].area == 3){num_area[3]++;}
      if(label->value[y][x].area == 4){num_area[4]++;}
      if(label->value[y][x].area == 5){num_area[5]++;}
      if(label->value[y][x].area == 6){num_area[6]++;}
    }
  }
  int i = 0;
  for(i=0;i<7;i++){
    printf("number of area%d = %d\n", i+1, num_area[i]);
  }
  int sum_sur = 0;
  sum_sur = num_area[1] + num_area[2] + num_area[3] + num_area[4] + num_area[5] + num_area[6];
  printf("number of summary of surrounding area = %d\n", sum_sur);  */

  /*for(y = 0;y < 14; y++){
      printf("%f %f\n",voronoi[y].gen.x,voronoi[y].gen.y);

  }*/
  if (label->value != NULL) {
    free(label->value);
    label->value = NULL;
  }
  printf("free label->value\n");
  if(class != NULL) {
    free(class);
    class = NULL;
  }
  printf("free class\n");
  if (label != NULL) {
    free(label);
    label = NULL;
  }
  printf("free label\n");
  
  if (voronoi != NULL) {
    free(voronoi);
    voronoi = NULL;
  }
  printf("free voronoi\n");
  /*
  if (voronoi2 != NULL) {
    free(voronoi2);
    voronoi2 = NULL;
  }
  printf("free voronoi2\n");
  */
  //if (enc != NULL) {
  //  free(enc);
    //enc = NULL;
  //}
  if (outfile != NULL) {
    free(outfile);
    outfile = NULL;
  }
  printf("free outfile\n");
  //free(label->value);
  //free(enc->class);
  //free(label);
  //free(voronoi);
  //free(voronoi2);
  //free(enc);
//  free(outfile);

  return 0;
}

VORONOI **alloc_voronoi_center_com(int max_yg, int max_xg)
{
  VORONOI **voronoi_center_com, *vp;
  int xg, yg;

  voronoi_center_com = (VORONOI **)alloc_2d_array(max_yg, max_xg, sizeof(VORONOI));//メモリ確保
  /*母点の座標すべて設定*/
	for(yg = 0; yg < max_yg; yg++){
		for(xg = 0; xg < max_xg; xg++){
			vp = &voronoi_center_com[yg][xg];
			vp->tl.x = vp->tl.y = 3279;
      		vp->br.x = vp->br.y = 0;
    }
  }
  return(voronoi_center_com);
}
/*
//voronoiのメモリ確保　９つの領域を1つの領域とするときの母点配置
VORONOI ***alloc_voronoi2(void)
{
  VORONOI ***voronoi, *vp;
  int xg, yg, area;
  int tlx = TOP_LEFT_x;
  int tly = TOP_LEFT_y;
  //int r = RADIUS;

  voronoi = (VORONOI ***)alloc_3d_array(128, 111, NUM_OF_VGEN, sizeof(VORONOI));//メモリ確保
  //母点の座標すべて設定
  for(yg = 0; yg < 128; yg++){ //中心の母点の数だけ
      for(xg = 0; xg < 111; xg++){ //母点内について
        if(yg % 2 == 0){ //偶数行 x座標左にずらす
          for(area = 0; area < 7; area++){
            vp = &voronoi[yg][xg][area];
            vp->gen.x = tlx + (xg * CENT_DIST_x * 3);
            vp->gen.y = tly + (yg * CENT_DIST_y * 3);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = area;
          }
        }
        if(yg % 2 == 1){ //奇数行
          for(area = 0; area < 7; area++){
            vp = &voronoi[yg][xg][area];
            vp->gen.x = tlx - (CENT_DIST_x * 1.5) + (xg * CENT_DIST_x * 3);
            vp->gen.y = tly + (yg * CENT_DIST_y * 3);
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
  VORONOI ***voronoi, *vp;
  int xg, yg;
  int tlx = TOP_LEFT_x;
  int tly = TOP_LEFT_y;
  int r = 4;

  voronoi = (VORONOI ***)alloc_3d_array(NUM_OF_GENy + 2, NUM_OF_GENx + 2, NUM_OF_VGEN, sizeof(VORONOI));//メモリ確保
  //母点の座標すべて設定
  for(yg = 0; yg < NUM_OF_GENy + 2; yg++){ //中心の母点の数だけ
      for(xg = 0; xg < NUM_OF_GENx + 2; xg++){ //母点内について
        if(yg % 2 == 0){ //偶数行 x座標左にずらす
          //中心
          vp = &voronoi[yg][xg][0];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x);
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          //右上
          vp = &voronoi[yg][xg][1];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) + (r / 2.0);
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) - (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 1;
          //右
          vp = &voronoi[yg][xg][2];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) + r;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 2;
          //右下
          vp = &voronoi[yg][xg][3];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) + r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) + (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 3;
          //左下
          vp = &voronoi[yg][xg][4];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) - r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) + (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 4;
          //左
          vp = &voronoi[yg][xg][5];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) - r;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 5;
          //左上
          vp = &voronoi[yg][xg][6];
            vp->gen.x = tlx - (CENT_DIST_x / 2.0) + (xg * CENT_DIST_x) - r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) - (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 6;
        }
        if(yg % 2 == 1){ //奇数行
          //中心
          vp = &voronoi[yg][xg][0];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x);
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 0;
          //右上
          vp = &voronoi[yg][xg][1];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) + (r / 2.0);
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) - (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 1;
          //右
          vp = &voronoi[yg][xg][2];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) + r;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 2;
          //右下
          vp = &voronoi[yg][xg][3];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) + r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) + (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 3;
          //左下
          vp = &voronoi[yg][xg][4];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) - r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) + (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 4;
          //左
          vp = &voronoi[yg][xg][5];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) - r;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 5;
          //左上
          vp = &voronoi[yg][xg][6];
            vp->gen.x = tlx + ((xg - 1) * CENT_DIST_x) - r / 2.0;
            vp->gen.y = tly + ((yg - 1) * CENT_DIST_y) - (sqrt(3) / 2.0 * r);
            vp->id.x = xg;
            vp->id.y = yg;
            vp->id.area = 6;
        }
      }
  }

  return(voronoi);
}

//画素ごとにラベリング
void fill_recur_all(VORONOI ***voronoi, VORONOI ***voronoi2, LABEL *label)
{
    int x, y;
    int xg, yg, area;
    double dx,dy;
    double gex, gey;
    //int id;
    long double min_dist, min_dist2, dist2;
    VORONOI *vp, *ivp;

    for(yg = 0; yg < NUM_OF_GENy + 2; yg++){
		for (xg = 0; xg < NUM_OF_GENx + 2; xg++) {
			for(area = 1; area < 7; area++){
				vp = &voronoi[yg][xg][area];
				gex = vp->gen.x;
				gey = vp->gen.y;
				ivp = &voronoi2[0][0][area];
				dx = DIFFPG(gex,ivp->gen.x); //x - gen.xの絶対値
    			dy = DIFFPG(gey,ivp->gen.y); //今参照してる母点とlabel範囲の点との距離
    			min_dist2 = dx*dx + dy*dy;//2乗の和
    			min_dist = (long double)(sqrt((double)min_dist2) + 1.0);
				for (y = 0; y < 128; y++) { //label範囲のすべての画素について
					for (x = 0; x < 111; x++) {
						vp = &voronoi2[y][x][area];
						dist2 = comp_dist(gex, gey, min_dist, min_dist2, vp, ivp);
						if (dist2 >= 0) {
							ivp = vp;
							min_dist2 = dist2;
							min_dist = (long)(sqrt((double)min_dist2) + 1.0);
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
}

//レンジ、サイズ、ラベル初期化して再設定
void fill_label_grobal_all(VORONOI ***voronoi, VORONOI ***voronoi2, LABEL *label)
{
    int xg, yg, area;
    int tlx, tly, brx, bry;
    FILE *fp;
    fp = fileopen("outfile_range.txt","rb");
    VORONOI *vp;

    for(yg = 0; yg < 128; yg++){
      for (xg = 0; xg < 111; xg++) {
        for(area = 0; area < 7; area++){
          vp = &voronoi2[yg][xg][area];
      		vp->tl.x = vp->tl.y = 3279;
      		vp->br.x = vp->br.y = 0;
      		vp->num_pixel = 0; //初期値としてボロノイのサイズ全部　0にする
        }
      }
    }

    for(yg = 0; yg < NUM_OF_GENy + 2; yg++){
      for (xg = 0; xg < NUM_OF_GENx + 2; xg++) {
        for(area = 0; area < 7; area++){
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
        }
      }
    }

	  fill_recur_all(voronoi, voronoi2, label);

    fclose(fp);
    return;
}

long comp_dist(int x, int y, long min_dist, long min_dist2,
	      VORONOI *vpt, VORONOI *vps)
{  //vpsが参照してる母点 vptを動かす
    double dx, dy;
    long double dist2;

    if (vpt == vps) return(-1L);
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
*/
//中央領域を近傍16ブロックでまとめる(labelの(xg,yg)の書き換え、voronoiの最初と最後の画素の書き換えのみ)
void compo_central16(VORONOI **voronoi, LABEL *label){

  printf("start compo_central16\n");

  int xg, yg, area;
  int x,y;
  double max_y, max_x;

  LABEL *label_cp = label;
  POINT **conv = (POINT **)alloc_2d_array(NUM_OF_GENy + 2, NUM_OF_GENx + 2, sizeof(POINT));
  POINT cp;
  VORONOI *vp;

max_y = 0;
max_x = 0;
  // (xg,yg)と(xg',yg')の対応表(conv)を作る(area=0)
  for(yg = 0; yg < NUM_OF_GENy+2; yg++){
    for(xg = 0; xg < NUM_OF_GENx+2; xg++){

      conv[yg][xg].x = xg / 4;
      conv[yg][xg].y = yg / 4;
      conv[yg][xg].area = 0;
	  if(conv[yg][xg].x > max_x) max_x = conv[yg][xg].x;
	  if(conv[yg][xg].y > max_y) max_y = conv[yg][xg].y;
    }
  }
  printf("中心部を6つにまとめた時の max_y = %lf, max_x = %lf\n", max_y, max_x);
  printf("made conv\n");

  //voronoi = alloc_voronoi_center_com((int)max_y + 1, (int)max_x + 1);
  //label->valueのarea=0のxg,ygの値を対応表を参照しながら書き換える
  for(y = 0; y < 3280; y++){
    for(x = 0; x < 3280; x++){
      if(label->value[y][x].area == 0){
        cp = label_cp->value[y][x];
        label->value[y][x].x = conv[(int)cp.y][(int)cp.x].x;
        label->value[y][x].y = conv[(int)cp.y][(int)cp.x].y;
      }
    }
  }
  printf("rewrote (xg,yg) of label->value\n");

  //voronoiのarea=0のtl,brを更新する
  for(yg = 0; yg < max_y + 1; yg++){
    for(xg = 0; xg < max_x + 1; xg++){
      area = 0;
	  vp = &voronoi[yg][xg];

      //area = 0;

      for(y = 0; y < 3280; y++){
        for(x = 0; x < 3280; x++){
		//cp = label_cp->value[y][x];
          if((int)label->value[y][x].x == xg && 
             (int)label->value[y][x].y == yg &&
             (int)label->value[y][x].area == area){

            vp->tl.x = MIN(x, vp->tl.x);
            vp->tl.y = MIN(y, vp->tl.y);
            vp->br.x = MAX(x, vp->br.x);
            vp->br.y = MAX(y, vp->br.y);
          }
        }
      }
    }
  }
  printf("updated tl,br of voronoi_center_com\n");
}