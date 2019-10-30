#define HAVE_CLOCK#define HAVE_64BIT_INTEGER#define QUADTREE_DEPTH	4#define BASE_BSIZE      8#define MAX_BSIZE       32#define MIN_BSIZE       (MAX_BSIZE >> QUADTREE_DEPTH)#define MAX_CLASS		100//63//69//127#define NUM_CLASS_CENTER 100//63//100//20//127#define NUM_CLASS_EDGE	100#define NUM_GROUP       16//#define PRD_ORDER       126//61//max = 96//今回は使っていません。#define COEF_PRECISION  6#define MAX_UPARA       512//2048//5000//512#define UPEL_DIST       4//#define NUM_UPELS		39#define NUM_UPELS       (UPEL_DIST * (UPEL_DIST +1))//126//184//389//(UPEL_DIST * (UPEL_DIST + 1))#define MAX_ITERATION   100//デバック用<2>通常<100>#define EXTRA_ITERATION 10#define NUM_PMODEL      16#define PM_ACCURACY     3#define NUM_KIND_PRD    4          //changed lfc#define MIN_FREQ        1#define MAX_SYMBOL		1024//4096//10000//1024		/* must be >> MAX_UPARA */#define PMCLASS_MAX	16#define PMCLASS_LEVEL	32//#define AUTO_DEL_CL//#define AUTO_PRD_ORDER#define MAX_PRD_ORDER_ALL 96//= MAX_PRD_ORDER#define MAX_PRD_ORDER  	96//389//389//184#define BASE_PRD_ORDER  51//389//389#define NUM_ZMODEL      49#define TOT_ZEROFR      (1 << 10)#define OPT_SIDEINFO#define VLC_MAXLEN      32#define MAGIC_NUMBER    ('M' << 8) + 'R'#define BANNER          "encmrp/decmrp version %.1f (Nov. 2004)"#define VERSION         5#define uint            unsigned int#define img_t           unsigned char#define cost_t          double#ifdef HAVE_64BIT_INTEGER#define RANGE_SIZE 64#if defined(_MSC_VER) || defined(__BORLANDC__)#	define range_t unsigned __int64#else#	define range_t unsigned long long#endif#	define MAX_TOTFREQ (1 << 20)	/* must be < RANGE_BOT */#else#	define RANGE_SIZE 32#	define range_t unsigned int#	define MAX_TOTFREQ (1 << 14)	/* must be < RANGE_BOT */#endif#define RANGE_TOP  ((range_t)1 << (RANGE_SIZE - 8))#define RANGE_BOT  ((range_t)1 << (RANGE_SIZE - 16))//vseg#define NUM_BLOCK 16810 //（仮）領域何個か合体させたブロック//#define NUM_OF_GEN (NUM_OF_CGEN * 7) //母点の数//#define NUM_OF_CGEN (NUM_OF_GENx * NUM_OF_GENy) //マイクロレンズの中心の母点の数//#define NUM_OF_VGEN 7 //マイクロレンズ内の母点の数//#define NUM_OF_VGENc 1//#define NUM_OF_GENx 329.0 //行のマイクロレンズ数(中心部)//330//#define NUM_OF_GENy 380.0 //列のマイクロレンズ数(中心部)//381//#define NUM_OF_GEN_ex 111.0//9個まとめた境界部の大きな6角形の数(行)//111//#define NUM_OF_GEN_ey 128.0//9個まとめら境界部の大きな6角形の数(列)//128//#define RUN_OVER 3.0#define MIN(x, y) (((x) < (y))? (x) : (y))#define MAX(x, y) (((x) > (y))? (x) : (y))//myaza#define NUM_REF_LENS 22//当該レンズも含む//NUM_REF_LENSとlens_center[]とnum_ref[]をよく見て設定する↓#define SEARCH_RANGE_Y 40//pel//必ず偶数#define SEARCH_RANGE_X 70//pel//必ず偶数#define LENS_MASK_MHD 4typedef struct {    int height;//3280    int width;//3280    int maxval;//255    img_t **val;//輝度値0~255} IMAGE;typedef struct {    int size;    int id;    uint *freq;    uint *cumfreq;    float *cost;    float *subcost;    double norm;} PMODEL;typedef struct {    int size;    int max_len;    int *len;    int *index;    int *off;    uint *code;} VLC;typedef struct {    range_t low;    range_t code;    range_t range;} RANGECODER;typedef struct {    int y, x, area;} POINT;typedef struct {    int y, x, p;} POINT_REF;typedef struct{	int y, x, lens_n;	double coef;} POINT_CORRELATION;typedef struct {	int *num_of_geny;	int *num_of_genx;	int max_geny;	int max_genx;	int num_of_vgen;	int num_of_vgenc;    int height;//画像のたての画素数    int width;//画像のよこの画素数    int maxval;//輝度値の最大値    int *num_class;//classの数	int max_num_pixel;    int num_group;//閾値で区切ったグループの数（16グループ）    int *base_prd_order;//予測次数    int *max_prd_order;//最大の予測次数	int max_prd_order_all;    int prd_mhd;//市街地距離の大きさ（使ってない）    int coef_precision;//予測係数の値をbitシフト(何ビットシフトさせるか）    int max_coef;//最大の予測係数//2    int num_pmodel;//確率モデルの数（１６）    int pm_accuracy;//確率モデルの桁落ちを防ぐ    int maxprd;//予測値の最大値（ビットシフト済み）    int f_huffman;//ハフマン符号をするかしないかのフラグ    int quadtree_depth;//四分木構造の何回分割したか    int optimize_loop;//loopが1stか2ndか    int num_kind_prd;//色の数（青、緑奇、赤、緑偶）    int ****predictor;//[cl][co][ar][k]が決まったときに係数をかえす(bitshiftされている６４倍)	int **predict_out;//[co][k]    int sub_cl;	double ***correlation_coef;    cost_t min_cost;//    int ***predictor_center;//    int ****predictor_edge;    int ****nzconv;//非ゼロ係数の位置    int ***num_nzcoef;//非ゼロ係数の数    int ****th;//スレッショルド    int **upara;//各座標の特徴量u    int **prd;//各画素の予測値    int max_class;    int **encval;//各画素の輝度値.予測誤差？（ビットシフト済み）    int **err;//各画素の予測誤差encvalと同じ    int **org;//各画素の輝度値    int *ctx_weight;//芝崎さんのδk    //int ***roff;//参照画素配置の位置[y][x][k]	//double ***roff;//lfc    int qtctx[QUADTREE_DEPTH << 3];    char **qtmap[QUADTREE_DEPTH];    int **class;//画素ごとのクラス選択情報    char **group;//画素ごとのグループ選択情報    int ****uquant;//changed量子化器    int **econv;//econv[org][prd]で予測誤差を返す。    img_t *bconv;//prdを与えるとbitを返す    img_t *fconv;//prdを与えるとbitを返す    PMODEL ***pmodels;    PMODEL **pmlist;    PMODEL spm;    VLC **vlcs;    RANGECODER *rc;//各メンバーに各bit数が入ってる？    double *sigma;//sigma[gr]はcomon.cのsigma_aか他を指す    int *mtfbuf;//move to front一時保存    int *ord2mhd;    int ***num_search;    int *cl_hist;    int *zero_m;    int *zero_fr;    int *coef_m;    cost_t ***coef_cost2;    cost_t **coef_cost;    cost_t *th_cost;    cost_t **class_cost;//vseg（ボロノイ分割の）    cost_t qtflag_cost[QUADTREE_DEPTH << 3];    cost_t **err_cost;    //int ***mask; //vseg    POINT ***gen_tl;	POINT ***gen_br;//vseg。各ボロノイ分割された領域の左上の座標と右下の座標[xg][yg][area][]    int ***num_pixel;//vseg    int ***gen_class;//vseg    POINT **label;//vseg    char **area;//area用配列！    int num_center_class;//vseg 2class    int *num_segment;		//****************符号量算出	int **area_co_err;	int *area_pel;	double *area_class_info;	double *area_predictor_info;	double *area_th_info;		} ENCODER;typedef struct {	  //int ***roff;//lfc	int *num_of_geny;	int *num_of_genx;	int max_geny;	int max_genx;	int num_of_vgen;	int num_of_vgenc;	int **org;//lfc    int version;    int height;    int width;    int maxval;    int num_comp;    int *num_class;    int num_group;    int *max_prd_order;	int max_prd_order_all;    int max_class;    int prd_mhd;    int num_pmodel;    int pm_accuracy;    int maxprd;    int max_coef;    int coef_precision;    int f_huffman;    int quadtree_depth;    int num_kind_prd;    int ****predictor;//四次元へ    int ****nzconv;//そのうち四次元へ    int ***num_nzcoef;    int ****th;    int **err;    int *ctx_weight;    char **qtmap[QUADTREE_DEPTH];    int **class;    int *pm_idx;    PMODEL ***pmodels;    PMODEL spm;    VLC **vlcs;    RANGECODER *rc;    double *sigma;    int *mtfbuf;    int *zero_fr;    int *ord2mhd;    POINT ***gen_tl;	POINT ***gen_br; //vseg    int ***num_pixel;//vseg    int ***gen_class;//vseg    POINT **label;//vseg    char **area;//追加    int num_center_class;//vseg2class} DECODER;/* encmrp.c *//* common.c */FILE *fileopen(char *, char *);void *alloc_mem(size_t);/*void **alloc_2d_array(int, int, int);char ***alloc_3d_array(int, int, int, int);char ****alloc_4d_array(int, int, int, int, int);*/void **alloc_2d_array(int, int, size_t);void ***alloc_3d_array(int, int, int, size_t);void ****alloc_4d_array(int, int, int, int, size_t);void init_array(int *, int, int);//初期化void init_2d_array(int **, int, int, int);void init_3d_array(int ***, int, int, int, int);void init_4d_array(int ****, int, int, int, int, int);IMAGE *alloc_image(int, int, int);int *gen_hufflen(uint *, int, int);void gen_huffcode(VLC *);VLC *make_vlc(uint *, int, int);VLC ***init_vlc(int, int, int *, int);void free_vlc(VLC *);VLC **init_vlcs(PMODEL ***, int, int);PMODEL ***init_pmodels(int, int, int, int *, double *, int);void set_spmodel(PMODEL *, int, int);int *init_ctx_weight(void);int e2E(int, int, int, int);int E2e(int, int, int, int);//void mtf_classlabel(int **, int *, int, int, int, int, int);void mtf_classlabel_vseg(int, int ***, int *, int , int, int, int *);//vseg//void mtf_classlabel_vseg_2class(int ***, int *, int , int, int, int, int);//vseg//int **set_color_position(int, int, int);int color(int, int);double cpu_time(void);/* data put out */void print_predictor(int ****, int, int *, int, int, char *);void print_predict_out(int **, int, int, char *);void print_predict_outall(int **, int, int, char *);void print_threshold(int, int ****, int, int *, int, PMODEL **, int *, char *);void print_class(int **, int, int, char *);void set_blksize(int **, char ***, int, int, int, int, int, int);void print_block_size(img_t **, int, int, int, char ***, char *);void print_color(int, int, int **, char *);void print_error(int **, int, int, char *);void print_prd(int **, int, int, char *);/* rc.c */RANGECODER *rc_init(void);void rc_encode(FILE *, RANGECODER *, uint, uint, uint);void rc_finishenc(FILE *, RANGECODER *);int rc_decode(FILE *, RANGECODER *, PMODEL *, int, int);void rc_startdec(FILE *, RANGECODER *);
