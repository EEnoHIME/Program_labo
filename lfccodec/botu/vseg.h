//４隅のマイクロレンズの中心座標(4テスト画像のmaxを取った画像のmla4隅に合わせた)
#define TOP_LEFT_x 8.5
#define TOP_LEFT_y 3.0
#define TOP_RIGHT_x 3277.0
#define TOP_RIGHT_y 3.0
#define BOTTOM_LEFT_x 6.0
#define BOTTOM_LEFT_y 3275.0
#define BOTTOM_RIGHT_x 3274.5
#define BOTTOM_RIGHT_y 3274.0

#define PIXELPITCH 1.399999976E-6
#define MLAPITCH 1.389861488E-5
#define MLADIST (double)MLAPITCH / (double)PIXELPITCH
#define ROTATION -0.00028155732434242963791//[rad]//0//-0.0003
#define SCALE_FACTOR_x 1.0007
#define SCALE_FACTOR_y 1.0013//1.0004389286041259766
//#define CENT_DIST_x ((TOP_RIGHT_x - TOP_LEFT_x) / ((double)NUM_OF_GENx - 1.0)) //行の母点間の距離
#define CENT_DIST_x MLADIST * SCALE_FACTOR_x
//#define CENT_DIST_y ((BOTTOM_LEFT_y - TOP_LEFT_y) / ((double)NUM_OF_GENy - 1.0)) //列の母点間の距離
#define CENT_DIST_y MLADIST * sqrt(3.0000) * 0.5 * SCALE_FACTOR_y
//#define CENT_DIST_y MLADIST * 1.0004389
#define SENSOR_OFFSET_y -1.6298500299453736302E-6 / PIXELPITCH
#define SENSOR_OFFSET_x -3.5040323734283452459E-6 / PIXELPITCH
#define NUM_OF_GEN (NUM_OF_CGEN * 7) //すべての母点の数
#define NUM_OF_CGEN (NUM_OF_GENx * NUM_OF_GENy) //マイクロレンズの中心の母点の数
#define NUM_OF_VGEN 7 //マイクロレンズ内の母点の数
#define NUM_OF_GENx 329.0 //行のマイクロレンズ数
#define NUM_OF_GENy 380.0 //列のマイクロレンズ数
#define NUM_OF_GEN_ex 111.0
#define NUM_OF_GEN_ey 128.0
#define RUN_OVER 3

#define RADIUS 7.0//中心から周りの母点までの距離 ここ変えるとマイクロレンズ内の母点の位置かわる

#define MIN(x, y) (((x) < (y))? (x) : (y))
#define MAX(x, y) (((x) > (y))? (x) : (y))

#define DIFFPG(p,g) fabs(p - g) //double型に


//座標構造体
typedef struct {
    double x;
    double y;
    int area;
} POINT;

//ボロノイ構造体
struct voronoi_tag {
    POINT gen;//母点座標
    POINT tl;
    POINT br;
    int num_pixel;//領域の画素数
    POINT id;
};
typedef struct voronoi_tag VORONOI;

typedef struct {
	int count;
	POINT id[10];
} RBS;

//ラベル構造体
typedef struct {
    POINT **value;
    int width;
    int height;
} LABEL;

//探索範囲の構造体
typedef struct {
    int num;
    VORONOI **vp;
    POINT tl;
    POINT br;
} RANGE;

typedef struct {
    int height;
    int width;
    int maxval;
    char **class;
} ENCODER;
