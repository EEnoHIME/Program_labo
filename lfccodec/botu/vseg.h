//�S���̃}�C�N�������Y�̒��S���W(4�e�X�g�摜��max��������摜��mla4���ɍ��킹��)
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
//#define CENT_DIST_x ((TOP_RIGHT_x - TOP_LEFT_x) / ((double)NUM_OF_GENx - 1.0)) //�s�̕�_�Ԃ̋���
#define CENT_DIST_x MLADIST * SCALE_FACTOR_x
//#define CENT_DIST_y ((BOTTOM_LEFT_y - TOP_LEFT_y) / ((double)NUM_OF_GENy - 1.0)) //��̕�_�Ԃ̋���
#define CENT_DIST_y MLADIST * sqrt(3.0000) * 0.5 * SCALE_FACTOR_y
//#define CENT_DIST_y MLADIST * 1.0004389
#define SENSOR_OFFSET_y -1.6298500299453736302E-6 / PIXELPITCH
#define SENSOR_OFFSET_x -3.5040323734283452459E-6 / PIXELPITCH
#define NUM_OF_GEN (NUM_OF_CGEN * 7) //���ׂĂ̕�_�̐�
#define NUM_OF_CGEN (NUM_OF_GENx * NUM_OF_GENy) //�}�C�N�������Y�̒��S�̕�_�̐�
#define NUM_OF_VGEN 7 //�}�C�N�������Y���̕�_�̐�
#define NUM_OF_GENx 329.0 //�s�̃}�C�N�������Y��
#define NUM_OF_GENy 380.0 //��̃}�C�N�������Y��
#define NUM_OF_GEN_ex 111.0
#define NUM_OF_GEN_ey 128.0
#define RUN_OVER 3

#define RADIUS 7.0//���S�������̕�_�܂ł̋��� �����ς���ƃ}�C�N�������Y���̕�_�̈ʒu�����

#define MIN(x, y) (((x) < (y))? (x) : (y))
#define MAX(x, y) (((x) > (y))? (x) : (y))

#define DIFFPG(p,g) fabs(p - g) //double�^��


//���W�\����
typedef struct {
    double x;
    double y;
    int area;
} POINT;

//�{���m�C�\����
struct voronoi_tag {
    POINT gen;//��_���W
    POINT tl;
    POINT br;
    int num_pixel;//�̈�̉�f��
    POINT id;
};
typedef struct voronoi_tag VORONOI;

typedef struct {
	int count;
	POINT id[10];
} RBS;

//���x���\����
typedef struct {
    POINT **value;
    int width;
    int height;
} LABEL;

//�T���͈͂̍\����
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
