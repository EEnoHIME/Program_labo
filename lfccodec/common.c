#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mrp.h"

#if NUM_KIND_PRD == 4
inline int color(int y, int x)
{
	switch(y & 1){
		case 1:
			switch(x & 1){
				case 1:
					return(2);
				default:
					return(3);
			}
		default:
			switch(x & 1){
				case 1:
					return(0);
				default:
					return(1);
			}
	}
}

#else
inline int color(int y, int x)
{	
	int i;
	
	i = y - x;
	if (i < 0) i = -i;
	if (i % 2 == 0) return(0);//color = 0;
	else if (y % 2 == 1 && x % 2 == 0) return(1);//color = 1;
	else return(2);//color = 2;
}	
#endif

/*POINT_REF dyx[NUM_OF_VGEN][126] = {
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	},
	{//area0
		//renzu0
		{0, -1, 3}, {-1, 0, 1},
		{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
		{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
		//k>=10
		{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
		//k>=14
		{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
		//k>=20
		{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
		//k>=25
		{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
		//k>=29
		{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
		//k>=35
		{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
		//k>=40
		{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
		//k>=45
		{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
		//k>=50
		{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
		//k>=55
		{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
		//k>=61
		{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
		//k>=66
		{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
		//k>=71
		{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
		//k>=76
		{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
		//k>=81
		{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
		//k>=86
		{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
		//k>=91
		{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
		//k>=96
		{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
		//k>=101
		{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
		//k>=106
		{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
		//k>=111
		{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
		//k>=116
		{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
		//k>=121
		{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}
	}
};*/

const POINT_REF dyx2[] = 
{
	
	{0, -1, 3}, {-1, 0, 1},
	{0, -2 , 0}, {-1, -1, 2}, {-2, 0, 0}, {-1, 1, 2},
	{0, -9, 3}, {0, -10, 0}, {0, -11, 3}, {-1, -10, 1},
	//k>=10
	{-8, -6 , 0}, {-8, -4, 0}, {-8, 4, 0}, {-8, 6, 0}, 
	//k>=14
	{-9, -5, 2}, {-9, 5, 2}, {-10, -6, 0}, {-10, -4, 0}, {-10, 4, 0}, {-10, 6, 0},
	//k>=20
	{-16, 0, 0}, {-17, 0, 1}, {-17, -1, 2}, {-17, 1, 2}, {-18, 0, 0},
	//k>=25
	{-8, -16, 0}, {-8, -14, 0}, {-8, 14, 0}, {-8, 16, 0}, 
	//k>=29
	{-9, -15, 2}, {-9, 15, 2}, {-10, -16, 0}, {-10, -14, 0}, {-10, 14, 0}, {-10, 16, 0},
	//k>=35
	{-1, 10, 1}, {0, -19, 3}, {0, -20, 0}, {0, -21, 3}, {-1, -20, 1},
	//k>=40
	{-16, -10, 0}, {-17, -9, 2}, {-17, 9, 2}, {-16, 10, 0}, {-17, -10, 1}, 
	//k>=45
	{-17, 10, 1}, {-17, -11, 2}, {-18, -10, 0}, {-18, 10, 0}, {-17, 11, 2},
	//k>=50
	{-1, 20, 1}, {-25, -5, 2}, {-26, -4, 0}, {-26, -5, 3}, {-26, -6, 0}, 
	//k>=55
	{-27, -5, 2}, {-25, 5, 2}, {-26, 4, 0}, {-26, 5, 3}, {-26, 6, 0}, {-27, 5, 2},
	//k>=61
	{0, -29, 3}, {0, -30, 0}, {0, -31, 3}, {-1, -30, 1}, {-1, 30, 1}, 
	//k>=66
	{-8, -26, 0}, {-8, -24, 0}, {-9, -25, 2}, {-10, -26, 0}, {-10, -24, 0},
	//k>=71
	{-8, 26, 0}, {-8, 24, 0}, {-9, 25, 2}, {-10, 26, 0}, {-10, 24, 0}, 
	//k>=76
	{-16, -20, 0}, {-17, -19, 2}, {-17, -20, 1}, {-17, -21, 2}, {-18, -20, 0},
	//k>=81
	{-16, 20, 0}, {-17, 19, 2}, {-17, 20, 1}, {-17, 21, 2}, {-18, 20, 0}, 
	//k>=86
	{-25, -15, 2}, {-26, -14, 0}, {-26, -15, 3}, {-26, -16, 0}, {-27, -15, 2},
	//k>=91
	{-25, 15, 2}, {-26, 14, 0}, {-26, 15, 3}, {-26, 16, 0}, {-27, 15, 2}, 
	//k>=96
	{-34, 0, 0}, {-35, 0, 1}, {-35, -1, 0}, {-35, 1, 0}, {-36, 0, 0},
	//k>=101
	{-34, -10, 0}, {-35, -10, 1}, {-35, -11, 2}, {-35, -9, 2}, {-36, -10, 0}, 
	//k>=106
	{-34, 10, 0}, {-35, 10, 1}, {-35, 9, 2}, {-35, 11, 2}, {-36, 10, 0},
	//k>=111
	{-43, -5, 2}, {-44, -5, 3}, {-44, -6, 0}, {-44, -4, 0}, {-45, -5, 2}, 
	//k>=116
	{-43, 5, 2}, {-44, 5, 3}, {-44, 4, 0}, {-44, 6, 0}, {-45, 5, 2},
	//k>=121
	{-52, 0, 0}, {-53, 0, 1}, {-53, -1, 2}, {-53, 1, 2}, {-54, 0, 0}

};

/*const POINT_REF ref_center[] =
{
	{0, 0, 0, 0, 0, 0}, 
	
	{0, -10, 0, 1, 0, 0}, {-9, -5, 2, 2, 0, 0}, {-9, 5, 2, 3, 0, 0}, {0, 10, 0, 4, 0, 0},
	
	{0, -20, 0, 5, 0, 0}, {-9, -15, 2, 6, 0, 0}, {-17, -10, 1, 7, 0, 0}, {-17, 0, 1, 8, 0, 0}, {-17, 10, 1, 9, 0, 0}, {-9, 15, 2, 10, 0, 0}, {0, 20, 0, 11, 0, 0},
	
	{0, -30, 0, 12, 0, 0}, {-9, -25, 2, 13, 0, 0}, {-17, -20, 1, 14, 0, 0}, {-26, -15, 3, 15, 0, 0}, {-26, -5, 3, 16, 0, 0}, 
	
	{-26, 5, 3, 17, 0, 0}, {-26, 15, 3, 18, 0, 0}, {-17, 20, 1, 19, 0, 0}, {-9, 25, 2, 20, 0, 0}, {0, 30, 0, 21, 0, 0}, 
	
	{-35, 0, 1, 22, 0, 0}, {-35, -10, 1, 23, 0, 0}, {-35, 10, 1, 24, 0, 0}, 
	
	{-44, -5, 3, 25, 0, 0}, {-44, 5, 3, 26, 0, 0},
	
	{-53, 0, 1, 27, 0, 0}
};*/

POINT lens_center[] =
{	//y,x
	//0
	{0, 0},
	//1
	{0, -10},{-9, -5},{-9, 5},{0, 10},
	//√3
	{-9, -15},{-17, 0},{-9, 15},
	//2
	{0, -20},{-17, -10},{-17, 10},{0, 20},
	//√7
	{-26, -5},{-26, 5},{-9, -25},{-17, -20},{-17, 20},{-9, 25},
	//3
	{0, -30},{-26, -15},{-26, 15},{0, 30},
};//必ずnum_refと合わせて設定する

int num_ref[] =
{
	6, 4, 5, 5, 1, 5, 5, 5, 4, 5, 5, 1, 5, 5, 5, 5, 5, 5, 4, 5, 5, 1, 
};
const double sigma_h[] = {0.85, 1.15, 1.50, 1.90, 2.55, 3.30, 4.25, 5.60,
                    7.15, 9.20,12.05,15.35,19.95,25.85,32.95,44.05};
const double sigma_a[] = {0.15, 0.26, 0.38, 0.57, 0.83, 1.18, 1.65, 2.31,
                    3.22, 4.47, 6.19, 8.55,11.80,16.27,22.42,30.89};//WÎ·
const double qtree_prob[] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};
const double zerocoef_prob[NUM_ZMODEL] = {
0.003,0.010,0.020,0.033,0.046,0.062,0.079,0.097,0.116,0.135,0.156,0.177,0.200,
0.222,0.246,0.270,0.294,0.319,0.344,0.370,0.396,0.421,0.447,0.474,0.500,0.526,
0.552,0.578,0.604,0.630,0.656,0.681,0.706,0.730,0.754,0.778,0.800,0.823,0.844,
0.865,0.884,0.903,0.921,0.938,0.954,0.967,0.980,0.990,0.997
};


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

void *alloc_mem(size_t size)//Èãª Á½çt@CÉ«o·@\Â«malloc
{
    void *ptr;
    if ((ptr = (void *)malloc(size)) == NULL) {
        fprintf(stderr, "Can't allocate memory (size = %u)!\n", (int)size);//lfc
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

/*void **alloc_2d_array(int height, int width, int size)
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

char ***alloc_3d_array(int height, int width, int depth, int size)
{
		char ***mat;
		char **ptr2, *ptr;
		int t, k, l;
		t = size * depth;

		mat = (char ***)alloc_mem((sizeof(void *) + sizeof(void **) * width + t * width) * height);
		ptr2 = (char **)(mat + height);
		ptr = (char *)(ptr2 + height * width);

		for(k = 0; k < height; k++){
			mat[k] = ptr2;
			for(l = 0; l < width; l++){
				ptr2[l] = ptr;
				ptr += t;
			}
			ptr2 += width;
		}
		return (mat);
}*/

/*void ***alloc_3d_array(int height, int width, int depth, int size)
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

void ****alloc_4d_array(int height, int width, int depth, int times, size_t size)
{
		void ****mat, ***ptr3, **ptr2;
		void *ptr;
		int t, k, l, m;
		t = size * times;

		mat = (void ****)alloc_mem((sizeof(*mat) + sizeof(**mat) * width
						+ sizeof(***mat) * width * depth + t * width * depth) * height);
		ptr3 = (void ***)(mat + height);
		ptr2 = (void **)(ptr3 + height * width);
		ptr = (void *)(ptr2 + height * width * depth);

		for(k = 0; k < height; k++){
			mat[k] = ptr3;
			for(l = 0; l < width; l++){
				ptr3[l] = ptr2;
				for(m = 0; m < depth; m++){
					ptr2[m] = ptr;
					ptr += t;
				}
				ptr2 += depth;
			}
			ptr3 += width;
		}
		return (mat);
}

/*char ****alloc_4d_array(int height, int width, int depth, int times, int size)
{
		char ****mat, ***ptr3, **ptr2, *ptr;
		int t, k, l, m;
		t = size * times;

		mat = (char ****)alloc_mem((sizeof(void *) + sizeof(void **) * width
						+ sizeof(void ***) * width * depth + t * width * depth) * height);
		ptr3 = (char ***)(mat + height);
		ptr2 = (char **)(ptr3 + height * width);
		ptr = (char *)(ptr2 + height * width * depth);

		for(k = 0; k < height; k++){
			mat[k] = ptr3;
			for(l = 0; l < width; l++){
				ptr3[l] = ptr2;
				for(m = 0; m < depth; m++){
					ptr2[m] = ptr;
					ptr += t;
				}
				ptr2 += depth;
			}
			ptr3 += width;
		}
		return (mat);
}*/

/*void ****alloc_4d_array(int height, int width, int depth, int times, int size)
{
		void ****mat;
		char *ptr;
		int k, l, m;

		mat = (void ****)alloc_3d_array(height, width, depth, sizeof(void *));
		ptr = (char *)alloc_mem(height * width * depth * times * size);
		for (k = 0; k < height; k++) {
			for(l = 0; l < width; l++) {
				for(m = 0; m < depth; m++) {
					mat[k][l][m] = ptr;
					ptr += times * size;
				}
			}
		}
		return(mat);
}*/
//vseg øÌÔÏ¦
/*void ***alloc_3d_array_vseg(int num_gen, int height, int width, int size)
{
    void ***ptr1, **mat;
    char *ptr2;
    int k, l;
		printf("%d\n",sizeof(void *) * height + height * width * num_gen * size);
    mat = (void **)alloc_mem(sizeof(void *) * height + height * width * num_gen * size);
    ptr1 = (void ***)alloc_mem(num_gen * size);
		printf("xyz\n");

		ptr2 = (char *)(mat + height);
    for (k = 0; k < NUM_OF_GEN; k++) {
			ptr1[k] = mat;
			for (l = 0; l < height; l++) {
					ptr1[k][l] = ptr2;
          ptr2 += width * size;
			}
			mat += height * width * size;
    }
    return(ptr1);
}*/

void init_array(int *array, int height, int value){	//1³zñÌú»
	int i=0;

	for(i=0; i<height; i++){
		array[i] = value;
	}
}

void init_2d_array(int **array, int height, int width, int value){
	int i=0, j=0;

	for(i=0; i<height; i++){
		for(j=0; j<width; j++){
			array[i][j] = value;
		}
	}
}

void init_3d_array(int ***array, int height, int width, int depth ,int value){
	int i=0, j=0, k=0;

	for(i=0; i<height; i++){
		for(j=0; j<width; j++){
			for(k=0; k<depth; k++){
				array[i][j][k] = value;
			}
		}
	}
}

void init_4d_array(int ****array, int height, int width, int depth, int time, int value){
	int i=0, j=0, k=0, l=0;

	for(i=0; i<height; i++){
		for(j=0; j<width; j++){
			for(k=0; k<depth; k++){
				for(l=0; l<time; l++){
					array[i][j][k][l] = value;
				}
			}
		}
	}
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

int *gen_hufflen(uint *hist, int size, int max_len)
{
    int i, j, k, l, *len, *index, *bits, *link;

    len = (int *)alloc_mem(size * sizeof(int));
    index = (int *)alloc_mem(size * sizeof(int));
    bits = (int *)alloc_mem(size * sizeof(int));
    link = (int *)alloc_mem(size * sizeof(int));
    for (i = 0; i < size; i++) {
        len[i] = 0;
        index[i] = i;
        link[i] = -1;
    }
    /* sort in decreasing order of frequency */
    for (i = size -1; i > 0; i--) {
	for (j = 0; j < i; j++) {
	    if (hist[index[j]] < hist[index[j + 1]]) {
                k = index[j + 1];
                index[j + 1] = index[j];
                index[j] = k;
	    }
	}
    }
    for (i = 0; i < size; i++) {
        bits[i] = index[i];	/* reserv a sorted index table */
    }
    for (i = size - 1; i > 0; i--) {
        k = index[i - 1];
        l = index[i];
        hist[k] += hist[l];
        len[k]++;
	while (link[k] >= 0) {
            k = link[k];
            len[k]++;
	}
        link[k] = l;
        len[l]++;
        while (link[l] >= 0) {
            l = link[l];
            len[l]++;
	}
	for (j = i - 1; j > 0; j--) {
	    if (hist[index[j - 1]] < hist[index[j]]) {
                k = index[j];
                index[j] = index[j - 1];
                index[j - 1] = k;
	    } else {
                break;
	    }
	}
    }
    /* limit the maximum code length to max_len */
    for (i = 0; i < size; i++) {
	index[i] = bits[i];	/* restore the index table */
        bits[i] = 0;
    }
    for (i = 0; i < size; i++) {
        bits[len[i]]++;
    }
    for (i = size - 1; i > max_len; i--) {
	while (bits[i] > 0) {
            j = i - 2;
            while(bits[j] == 0) j--;
            bits[i] -= 2;
            bits[i - 1]++;
            bits[j + 1] += 2;
            bits[j]--;
	}
    }
    for (i = k = 0; i < size; i++) {
	for (j = 0; j < bits[i]; j++) {
            len[index[k++]] = i;
	}
    }
    free(link);
    free(bits);
    free(index);
    return (len);
}

void gen_huffcode(VLC *vlc)
{
    int i, j, *idx, *len;
    uint k;

    vlc->index = idx = (int *)alloc_mem(vlc->size * sizeof(int));
    vlc->off = (int *)alloc_mem(vlc->max_len * sizeof(int));
    vlc->code = (uint *)alloc_mem(vlc->size * sizeof(int));
    len = vlc->len;
    /* sort in increasing order of code length */
    for (i = 0; i < vlc->size; i++) {
        idx[i] = i;
    }
    for (i = vlc->size -1; i > 0; i--) {
	for (j = 0; j < i; j++) {
	    if (len[idx[j]] > len[idx[j + 1]]) {
                k = idx[j + 1];
                idx[j + 1] = idx[j];
                idx[j] = k;
	    }
	}
    }
    k = 0;
    for (j = 0; j < vlc->max_len; j++) {
	vlc->off[j] = -1;
    }
    j = len[idx[0]];
    for (i = 0; i < vlc->size; i++) {
	if (j < len[idx[i]]) {
	    k <<= (len[idx[i]] - j);
	    j = len[idx[i]];
	}
	vlc->code[idx[i]] = k++;
	vlc->off[j - 1] = i;
    }
    return;
}

VLC *make_vlc(uint *hist, int size, int max_len)
{
    VLC *vlc;

    vlc = (VLC *)alloc_mem(sizeof(VLC));
    vlc->size = size;
    vlc->max_len = max_len;
    vlc->len = gen_hufflen(hist, size, max_len);
    gen_huffcode(vlc);
    return (vlc);
}

void free_vlc(VLC *vlc)
{
    free(vlc->code);
    free(vlc->off);
    free(vlc->index);
    free(vlc->len);
    free(vlc);
    return;
}

VLC **init_vlcs(PMODEL ***pmodels, int num_group, int num_pmodel)
{
    VLC **vlcs, *vlc;
    PMODEL *pm;
    int gr, k;

    vlcs = (VLC **)alloc_2d_array(num_group, num_pmodel, sizeof(VLC));
    for (gr = 0; gr < num_group; gr++) {
			for (k = 0; k < num_pmodel; k++) {
			    vlc = &vlcs[gr][k];
			    pm = pmodels[gr][k];
			    vlc->size = pm->size;
			    vlc->max_len = VLC_MAXLEN;
			    vlc->len = gen_hufflen(pm->freq, pm->size, VLC_MAXLEN);
			    gen_huffcode(vlc);
			}
    }
    return (vlcs);
}

/*
  Natural logarithm of the gamma function
  cf. "Numerical Recipes in C", 6.1
  http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf.html
*/
double lngamma(double xx)
{
    int j;
    double x,y,tmp,ser;
    double cof[6] = {
	76.18009172947146,	-86.50532032941677,
	24.01409824083091,	-1.231739572450155,
	0.1208650973866179e-2,	-0.5395239384953e-5
    };

    y = x = xx;
    tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
    ser = 1.000000000190015;
    for (j=0;j<=5;j++)
	ser += (cof[j] / ++y);
    return (log(2.5066282746310005 * ser / x) - tmp);
}

double calc_ggprob(double beta, double shape, double h, double x)
{
    double p;

    if (x < 0.0) x = -x;
    if (x < 1E-6) {
	p = exp(-pow(beta * h, shape)) + exp(0.0);
    } else {
	p = exp(-pow(beta * (x - h), shape))
	  + exp(-pow(beta * (x + h), shape));
    }
    return (p);
}

void set_freqtable(PMODEL *pm, int size, int ssize,
		   double shape, double sigma, double h, double off)
{
    double beta, norm;
    int i;

    pm->size = size;
    pm->freq = (uint *)alloc_mem((size * 2 + 1) * sizeof(uint));
    pm->cumfreq = &pm->freq[size];
    /* Generalized Gaussian distribution */
    beta = exp(0.5*(lngamma(3.0/shape)-lngamma(1.0/shape))) / sigma;
    norm = 0.0;
    for (i = 0; i < size; i++) {
	norm += calc_ggprob(beta, shape, h, i - off);
    }
    norm = (double)(MAX_TOTFREQ - size * MIN_FREQ) / norm;
    norm += 1E-8;	/* to avoid machine dependent rounding errors */
    pm->norm = norm;
    pm->cumfreq[0] = 0;
    for (i = 0; i < size; i++) {
	pm->freq[i] = norm * calc_ggprob(beta, shape, h, i - off) + MIN_FREQ;
	pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
    }
    if (ssize > 0) {
	pm->cost = (float *)alloc_mem((size + ssize) * sizeof(float));
	pm->subcost = &pm->cost[size];
    }
    return;
}

PMODEL ***init_pmodels(int num_group, int num_pmodel, int pm_accuracy,//KEXªzÌ³ðè
		       int *pm_idx, double *sigma, int size)
{
    PMODEL ***pmodels, *pmbuf, *pm;
    int gr, i, j, num_subpm, ssize;
    double delta_c, c, s, sw, off;

    if (pm_accuracy < 0) {
			num_subpm = 1;
			ssize = 1;
			sw = 0.0;
    } else {
			num_subpm = 1 << pm_accuracy;
			ssize = size;
			size = size + ssize - 1;
			sw = 1.0 / (double)num_subpm;
    }
    delta_c = 3.2 / (double)num_pmodel;
    off = (double)(ssize - 1);
    if (pm_idx != NULL) {
			num_pmodel = 1;
			ssize = 0;
    }
    pmodels = (PMODEL ***)alloc_2d_array(num_group, num_pmodel,
					 sizeof(PMODEL *));
    pmbuf = (PMODEL *)alloc_mem(num_group * num_pmodel * num_subpm
				* sizeof(PMODEL));
    for (gr = 0; gr < num_group; gr++) {
			s = sigma[gr];
			if (pm_accuracy < 0) s *= 2.0;
			for (i = 0; i < num_pmodel; i++) {
			    pmodels[gr][i] = pmbuf;
			    for (j = 0; j < num_subpm; j++) {
						pm = pmbuf++;
						pm->id = i;
						if (num_pmodel > 1) {
						    c = delta_c * (double)(i + 1);
						} else if (pm_idx != NULL) {
						    c = delta_c * (double)(pm_idx[gr] + 1);
						} else {
						    c = 2.0;
						}
						if (c < 0.1) c = 0.1;
						set_freqtable(pm, size, ssize, c, s, sw/2.0, off - sw * j);
			    }
			}
    }
    return (pmodels);
}

/* probaility model for coefficients and thresholds */
void set_spmodel(PMODEL *pm, int size, int m)
{
    int i, sum;
    double p;

    pm->size = size;
    if (m >= 0) {
			p = 1.0 / (double)(1 << (m % 8));
			sum = 0;
			for (i = 0; i < pm->size; i++) {
			    pm->freq[i] = exp(-p * i) * (1 << 10);
			    if (pm->freq[i] == 0) pm->freq[i]++;
			    sum += pm->freq[i];
			}
			if (m & 8) pm->freq[0] = (sum - pm->freq[0]);	/* weight for zero */
    } else {
			for (i = 0; i < pm->size; i++) {
			    pm->freq[i] = 1;
			}
    }
    pm->cumfreq[0] = 0;
    for (i = 0; i < pm->size; i++) {
			pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
    }
    return;
}

int *init_ctx_weight(void)
{
    int *ctx_weight, k;
    //double dy, dx;

    ctx_weight = (int *)alloc_mem(NUM_UPELS * sizeof(int));
    for (k = 0; k < NUM_UPELS; k++) {
	//dy = dyx[k].y;
	//dx = dyx[k].x;
	//ctx_weight[k] = (64.0 / (fabs(dy) + fabs(dx)) + 0.5);
	//ctx_weight[k] = 64.0 / sqrt(dy * dy + dx * dx) + 0.5;//lfc
	ctx_weight[k] = 64.0;//lfc
    }
    return (ctx_weight);
}

/*int **set_color_position(int height, int width, int num_kind_prd)
{
    int x, y, i;
    int **color;

    color = (int **)alloc_2d_array(height, width, sizeof(int));
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
			if (num_kind_prd == 4) {
				if (y % 2 == 0 && x % 2 == 0) color[y][x] = 1;//Â
				else if (y % 2 == 1 && x % 2 == 0) color[y][x] = 0;//Î
				else if (y % 2 == 0 && x % 2 == 1) color[y][x] = 3;//Î
				else color[y][x] = 2;//lfcÔ
			} else {
				i = y - x;
				if (i < 0) i = -i;
				if (i % 2 == 0) color[y][x] = 0;
				else if (y % 2 == 1 && x % 2 == 0) color[y][x] = 1;
				else color[y][x] = 2;
			}
		}
    }
    return(color);
}*/

int e2E(int e, int prd, int flag, int maxval)
{
    int E, th;

    E = (e > 0)? e : -e;
    th = (prd < ((maxval + 1) >> 1))? prd : maxval - prd;
    if (E > th) {
        E += th;
    } else if (flag) {
	E = (e < 0)? (E << 1) - 1 : (E << 1);
    } else {
	E = (e > 0)? (E << 1) - 1 : (E << 1);
    }
    return (E);
}

int E2e(int E, int prd, int flag, int maxval)
{
    int e, th;

    th = (prd < ((maxval + 1) >> 1))? prd : maxval - prd;
    if (E > (th << 1)) {
	e = (prd < ((maxval + 1) >> 1))? E - th : th - E;
    } else if (flag) {
	e = (E & 1)? -((E >> 1) + 1) : (E >> 1);
    } else {
	e = (E & 1)? (E >> 1) + 1 : -(E >> 1);
    }
    return (e);
}

//vsegpÉÂ­é
void mtf_classlabel(int **class, int *mtfbuf, int y, int x,
		    int bsize, int width, int num_class)
{
    int i, j, k, ref[3];
			//ref[0]:ã ref[1]:¶ ref[2]:Eã
    if (y == 0) {
			if (x == 0) {
			    ref[0] = ref[1] = ref[2] = 0;
			} else {
			    ref[0] = ref[1] = ref[2] = class[y][x-1]; //êÔãÌs
			}
    } else {
			ref[0] = class[y-1][x];
			ref[1] = (x == 0)? class[y-1][x] : class[y][x-1]; //êÔ¶Ìñ©Ç¤©
			ref[2] = (x + bsize >= width)?
	    class[y-1][x] : class[y-1][x+bsize];
			if (ref[1] == ref[2]) {
			    ref[2] = ref[0];
			    ref[0] = ref[1];
			}
    }
    /* move to front */
    for (k = 2; k >= 0; k--) {
			if ((j = mtfbuf[ref[k]]) == 0) continue;
			for (i = 0; i < num_class; i++) {
			    if (mtfbuf[i] < j) {
						mtfbuf[i]++;
			    }
			}
			mtfbuf[ref[k]] = 0;
    }
    return;
}

//vseg
void mtf_classlabel_vseg(int max_genx, int ***class, int *mtfbuf, int xg, int yg, int area, int *num_class)
{ 
		int i, j, h, ref[3];
		//ref[0]:ã ref[1]:¶ ref[2]:Eã
		//int upper_class, left_class, upright_class;
		//upper_class = class[area][yg - 1][xg];
		//left_class = class[area][yg][xg - 1];
		//upright_class = class[area][yg - 1][xg + 1];

		if(yg == 0){
			if(xg == 0){
				ref[0] = ref[1] = ref[2] = 0;
			}else{
				if(class[area][yg][xg - 1] < 0){
					ref[0] = ref[1] = ref[2] = 0;
				}else{
					ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];
				}
			}
		}else{
			if(xg == 0){
				if(class[area][yg - 1][xg] < 0){
					if(class[area][yg - 1][xg + 1] < 0){
						ref[0] = ref[1] = ref[2] = 0;
					}else{
						ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg + 1];
					}
				}else if(class[area][yg - 1][xg + 1] < 0){
					ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];
				}else{
					ref[0] = class[area][yg - 1][xg];
					ref[1] = class[area][yg - 1][xg];
					ref[2] = class[area][yg - 1][xg + 1];
				}
			}else if(xg == max_genx-1){
				if(class[area][yg - 1][xg] < 0){
					if(class[area][yg][xg - 1] < 0){
						ref[0] = ref[1] = ref[2] = 0;
					}else{
						ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];
					}
				}else if(class[area][yg][xg - 1] < 0){
					ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];
				}else{
					ref[0] = class[area][yg - 1][xg];
					ref[1] = class[area][yg][xg - 1];
					ref[2] = class[area][yg - 1][xg];
				}
			}else{
				if(class[area][yg - 1][xg] < 0){
					if(class[area][yg][xg - 1] < 0){
						if(class[area][yg - 1][xg + 1] < 0){
							ref[0] = ref[1] = ref[2] = 0;
						}else{
							ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg + 1];
						}
					}else if(class[area][yg - 1][xg + 1] < 0){
						ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];
					}else{
						ref[0] = ref[1] = class[area][yg][xg - 1];
						ref[2] = class[area][yg - 1][xg + 1];
					}
				}else if(class[area][yg][xg - 1] < 0){
					if(class[area][yg - 1][xg + 1] < 0){
						ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];
					}else{
						ref[0] = class[area][yg - 1][xg];
						ref[1] = class[area][yg - 1][xg];
						ref[2] = class[area][yg - 1][xg + 1];
					}
				}else if(class[area][yg - 1][xg + 1] < 0){
					ref[0] = class[area][yg - 1][xg];
					ref[1] = class[area][yg][xg - 1];
					ref[2] = class[area][yg - 1][xg];
				}else{
					ref[0] = class[area][yg - 1][xg];
					ref[1] = class[area][yg][xg - 1];
					ref[2] = class[area][yg - 1][xg + 1];
				}
			}

			if (ref[1] == ref[2]) {
			    ref[2] = ref[0];
			    ref[0] = ref[1];
			}
		}

		/* move to front */
		for (h = 2; h >= 0; h--) {
			if ((j = mtfbuf[ref[h]]) == 0) continue;
			for (i = 0; i < num_class[area]; i++) {
					if (mtfbuf[i] < j) {
						mtfbuf[i]++;
					}
			}
			mtfbuf[ref[h]] = 0;
		}
		return;
}

//vseg
/*void mtf_classlabel_vseg_2class(int ***class, int *mtfbuf, int xg, int yg, int area, int num_class, int num_center_class)
{ 
	int i, j, h, ref[3];
	//ref[0]:ã ref[1]:¶ ref[2]:Eã
	//int upper_class, left_class, upright_class;
	//upper_class = class[area][yg - 1][xg];
	//left_class = class[area][yg][xg - 1];
	//upright_class = class[area][yg - 1][xg + 1];

	if(yg == 0){
		if(xg == 0){//¶ã
			ref[0] = ref[1] = ref[2] = num_center_class;
		}else{//êÔãÌs
			if(class[area][yg][xg - 1] < 0){
				ref[0] = ref[1] = ref[2] = num_center_class;
			}else{
				ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];
			}
		}
	}else{
		if(xg == 0){//êÔ¶Ìñ
			if(class[area][yg - 1][xg] < 0){
				if(class[area][yg - 1][xg + 1] < 0){
					ref[0] = ref[1] = ref[2] = num_center_class;
				}else{
					ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg + 1];//Eãok
				}
			}else if(class[area][yg - 1][xg + 1] < 0){
				ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];//ãok
			}else{
				ref[0] = class[area][yg - 1][xg];//ãAEãok
				ref[1] = ref[2] = class[area][yg - 1][xg + 1];
			}
		}else if(xg == NUM_OF_GENx + 1){//êÔEÌñ
			if(class[area][yg - 1][xg] < 0){
				if(class[area][yg][xg - 1] < 0){
					ref[0] = ref[1] = ref[2] = num_center_class;
				}else{
					ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];//¶ok
				}
			}else if(class[area][yg][xg - 1] < 0){
				ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];//ãok
			}else{
				ref[0] = class[area][yg - 1][xg];//ãA¶ok
				ref[1] = ref[2] = class[area][yg][xg - 1];
			}
		}else{
			if(class[area][yg - 1][xg] < 0){
				if(class[area][yg][xg - 1] < 0){
					if(class[area][yg - 1][xg + 1] < 0){
						ref[0] = ref[1] = ref[2] = num_center_class;
					}else{
						ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg + 1];//Eãok
					}
				}else if(class[area][yg - 1][xg + 1] < 0){
					ref[0] = ref[1] = ref[2] = class[area][yg][xg - 1];//¶ok
				}else{
					ref[0] = ref[1] = class[area][yg][xg - 1];//¶AEãok
					ref[2] = class[area][yg - 1][xg + 1];
				}
			}else if(class[area][yg][xg - 1] < 0){
				if(class[area][yg - 1][xg + 1] < 0){
					ref[0] = ref[1] = ref[2] = class[area][yg - 1][xg];//ãok
				}else{
					ref[0] = class[area][yg - 1][xg];//ãAEãok
					ref[1] = ref[2] = class[area][yg - 1][xg + 1];
				}
			}else if(class[area][yg - 1][xg + 1] < 0){
				ref[0] = class[area][yg - 1][xg];//ãA¶ok
				ref[1] = ref[2] = class[area][yg][xg - 1];
			}else{
				ref[0] = class[area][yg - 1][xg];//Sok
				ref[1] = class[area][yg][xg - 1];
				ref[2] = class[area][yg - 1][xg + 1];
			}
		}

		if (ref[1] == ref[2]) {
				ref[2] = ref[0];
				ref[0] = ref[1];
		}
	}

	// move to front 
	for (h = 2; h >= 0; h--) {
		if ((j = mtfbuf[ref[h]]) == 0) continue;
		for (i = 0; i < num_class; i++) {
				if (mtfbuf[i] < j) {
					mtfbuf[i]++;
				}
		}
		mtfbuf[ref[h]] = 0;
	}
	return;
}*/

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

void swap(int *a, int *b)//a,bðüêÖ¦é
{
	int t = *a;
	*a = *b;
	*b = t;
}

/*void print_predict_out(int **predictor, int prd_order,
                     int num_kind_prd, char *outfile)
{
    int i, j, co, k;
		
	int max_mhd = 100;//10//lfc}nb^fBX^XisXn£j
    for (co = 0; co < num_kind_prd; co++) {
		FILE *fp;
		int **p_prd;
		char name[100];
		char color[10];

		if (co == 0) {//color¶ñpzñÉeFîñðüêéD
				if (num_kind_prd == 4) { sprintf(color, "G_odd");
			} else { sprintf(color, "G");}
			} else if (co == 1) { sprintf(color,"B");
			} else if (co == 2) { sprintf(color,"R");
			} else { sprintf(color, "G_even");
			}
		sprintf(name, "./infofile/predictor/%s_prd_all_%s.csv", outfile, color);//nameÉGNZfile¼ðüêéD
		fp = fileopen(name, "wb");
		p_prd = (int **)alloc_2d_array(max_mhd + 1, max_mhd << 1, sizeof(int));
		for (i = 0; i < max_mhd + 1; i++) {//p_prd[yÀW][xÀW]IÈH
			for (j = 0; j < max_mhd << 1; j++) {//maxmhd<<1Í120
				p_prd[i][j] = INT_MAX;//129
			}	
	    }//dyx[k]ÍkðüêéÆ»ÌÀWðo[x,yÅ¦¹éD¨»ç­»ÌÀWÍ»ÎÛæf©çÌsXn£ÉÈÁÄ¢é
	    for (k = 0; k < prd_order; k++) {//lfc
			p_prd[max_mhd + dyx[k].y][max_mhd + dyx[k].x] = predictor[co][k];	
		}
		//±±©çExelt@CÉ«oµ
	    for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				if (p_prd[i][j] == INT_MAX) fprintf(fp, ",");
				else fprintf(fp, "%d,", p_prd[i][j]);
			}
			fprintf(fp, "\n");
	    }
		fclose(fp);
		free(p_prd);
    }
    return;
}*/
/*void print_predict_outall(int **predictor, int prd_order,
                     int num_kind_prd, char *outfile)
{
    int i, j, co, k;
		
	int max_mhd = 100;//10//lfc}nb^fBX^XisXn£j
    
		FILE *fp;
		int **p_prd, *predict_all;
		char name[100];
		char color[10];

		sprintf(color, "all_color");
		sprintf(name, "./infofile/predictor/%s_prd_all_%s.csv", outfile, color);//nameÉGNZfile¼ðüêéD
		fp = fileopen(name, "wb");
		p_prd = (int **)alloc_2d_array(max_mhd + 1, max_mhd << 1, sizeof(int));
		predict_all = (int *)alloc_mem(prd_order * sizeof(int));
		for(k = 0; k < prd_order; k++){
			predict_all[k] = 0;;
		}
		for (i = 0; i < max_mhd + 1; i++) {//p_prd[yÀW][xÀW]IÈH
			for (j = 0; j < max_mhd << 1; j++) {//maxmhd<<1Í120
				p_prd[i][j] = INT_MAX;//129
			}	
	    }//dyx[k]ÍkðüêéÆ»ÌÀWðo[x,yÅ¦¹éD¨»ç­»ÌÀWÍ»ÎÛæf©çÌsXn£ÉÈÁÄ¢é
	    for (k = 0; k < prd_order; k++) {//lfc
			for (co = 0; co < num_kind_prd; co++) {
				predict_all[k] += predictor[co][k];
			}
		}
		for(k = 0; k < prd_order; k++){
			p_prd[max_mhd + dyx[k].y][max_mhd + dyx[k].x] = predict_all[k];
		}
		//±±©çExelt@CÉ«oµ
	    for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				if (p_prd[i][j] == INT_MAX) fprintf(fp, ",");
				else fprintf(fp, "%d,", p_prd[i][j]);
			}
			fprintf(fp, "\n");
	    }
		fclose(fp);
		free(p_prd);
		free(predict_all);
    return;
}*/
/* data put out */
/*void print_predictor(int ****predictor, int prd_order, int *num_class,
                     int num_kind_prd, int max_coef, char *outfile)
{
    int i, j, cl, co, ar, k;
		int renzu[12];//lfc//QÆæfªzu³êÄ¢é}CNY
		int a[12];//lfc//average?
		int max_mhd = 100;//10//lfc}nb^fBX^XisXn£j
		for (i = 0; i < 12; i++) {//ú»
			renzu[i] = 0;
		}
    for (co = 0; co < num_kind_prd; co++) {
      FILE *fp;
      int **p_prd;
      char name[100];
      char color[10];

      if (co == 0) {
				if (num_kind_prd == 4) { sprintf(color, "G_odd");
	      } else { sprintf(color, "G");
				}
			} else if (co == 1) { sprintf(color,"B");
			} else if (co == 2) { sprintf(color,"R");
			} else { sprintf(color, "G_even");
			}
      sprintf(name, "./infofile/predictor/%s_prd_%s.csv", outfile, color);
      fp = fileopen(name, "wb");
      p_prd = (int **)alloc_2d_array(max_mhd + 1, max_mhd << 1, sizeof(int));
		for(ar = 0; ar < NUM_OF_VGEN; ar++){
			fprintf(fp, "area%d\n", ar);
      for (cl = 0; cl < num_class[ar]; cl++) {
        for (i = 0; i < max_mhd + 1; i++) {
	        for (j = 0; j < max_mhd << 1; j++) {
            p_prd[i][j] = max_coef + 1;//129
					}
	    	}//dyx[k]
	    	for (k = 0; k < prd_order; k++) {//lfc
          p_prd[max_mhd + dyx[k].y][max_mhd + dyx[k].x] = predictor[cl][co][ar][k];
					//printf("p_prd[%d][%d]=%d\n", max_mhd + dyx[k].y, max_mhd + dyx[k].x, predictor[cl][co][k]);
					if (k < 6){
						renzu[0] += abs(predictor[cl][co][ar][k]);
					}
					else if (k < 10){
						renzu[1] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 10 || k == 11 || k == 14 || k == 16 || k == 17){
						renzu[2] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 12 || k == 13 || k == 15 || k == 18 || k == 19){
						renzu[3] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 20 || k == 21 || k == 22 || k == 23 || k == 24){
						renzu[4] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 25 || k == 26 || k == 29 || k == 31 || k == 32){
						renzu[5] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 27 || k == 28 || k == 30 || k == 33 || k == 34){
						renzu[6] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 35){
						renzu[7] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 36 || k == 37 || k == 38 || k == 39){
						renzu[8] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 40 || k == 41 || k == 44 || k == 46 || k == 47){
						renzu[9] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 42 || k == 43 || k == 45 || k == 48 || k == 49){
						renzu[10] += abs(predictor[cl][co][ar][k]);
					}
					else if (k == 50){
						renzu[11] += abs(predictor[cl][co][ar][k]);
					}
				}
				
	    	for (i = 0; i < max_mhd + 1; i++) {
				for (j = 0; j < max_mhd << 1; j++) {
					if (p_prd[i][j] == max_coef + 1) fprintf(fp, ",");
		    		else fprintf(fp, "%d,", p_prd[i][j]);
				}
				fprintf(fp, "\n");
	    	}
        fprintf(fp, "\n");
			}
		}
			
			printf("%d\n", co);
			for (i = 0; i < 12; i++) {
				if (i == 0)a[i]=renzu[i]/6;
				else if (i == 1 || i == 8)a[i] = renzu[i] / 4;
				else if (i == 7 || i == 11)a[i] = renzu[i] ;
				else a[i] = renzu[i] / 5;
			}

			int *p = (int *)a;
			int n = sizeof(a) / sizeof(int);
			int i, j;
			int cpy[sizeof(a) / sizeof(int)];
			int idx[sizeof(a) / sizeof(int)];

			for (i = 0; i < n; i++) {
				cpy[i] = p[i];
				idx[i] = i;
			}
			for (i = 0; i < n - 1; i++) {
				for (j = i + 1; j < n; j++) {
					if (cpy[i] < cpy[j]) {
						swap(&cpy[i], &cpy[j]);
						swap(&idx[i], &idx[j]);
					}
				}
			}
			for (i = 0; i < n; i++) {
				int ii = idx[i] / (sizeof(a[0]) / sizeof(int));
				printf("a[%d] = %d\n", ii, a[ii]);
			}


      fclose(fp);
      free(p_prd);
    }
    return;
}*/

void print_threshold(int num_of_vgen,int ****th, int num_group, int *num_class, int num_kind_prd,
                     PMODEL **pmlist, int *pm_idx, char *outfile)
{
    int cl, co, ar, gr;

    for (co = 0; co <num_kind_prd; co++) {
        char name[100];
        char color[10];
        FILE *fp;

        if (co == 0) {
					if (num_kind_prd == 4) {
						sprintf(color, "G_odd");
	        } else{  sprintf(color, "G");
					}
				} else if (co == 1) { sprintf(color, "B");
				} else if (co == 2) { sprintf(color, "R");
				} else { sprintf(color, "G_even");
				}

        	sprintf(name, "./infofile/threshold/%s_th_%s.csv", outfile, color);
        	fp = fileopen(name, "wb");
					for(ar = 0; ar < num_of_vgen; ar++){
        		for (cl = 0; cl < num_class[ar]; cl++) {
	    				for (gr = 0; gr < num_group - 1; gr++) {
	        			fprintf(fp, "%d,", th[cl][co][ar][gr]);
	    				}
	    				fprintf(fp, "\n");
						}
						fprintf(fp, "\n");
					}
        	fprintf(fp, "\n");
        	if (pmlist != NULL) {
	    			for (gr = 0; gr < num_group; gr++) {
            	fprintf(fp, "%d,", pmlist[gr]->id);
	    			}
            fprintf(fp, "\n");
					}
        	if (pm_idx != NULL) {
	    			for (gr = 0; gr < num_group; gr++) {
            	fprintf(fp, "%d,", pm_idx[gr]);
	    			}
            fprintf(fp, "\n");
					}
        	fclose(fp);
    }
    return;
}

void print_class(int **class, int height, int width, char *outfile)
{
    int i, j;
    char name[100];
    FILE *fp;

    sprintf(name, "./infofile/class/%s_class.pgm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
		            putc(class[i][j] * 4, fp);
			}
    }
    fclose(fp);
    return;
}

void print_error(int **error, int height, int width, char *outfile)
{
    int i, j, e;
    char name[100];
    FILE *fp;

    sprintf(name, "./LOG/%s_error.pgm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				e = error[i][j];
				if(e < 0) e = -1 * e;
				if(e > 255) e = 255;
		            putc(e, fp);
			}
    }
    fclose(fp);
    return;
}

void print_prd(int **prd, int height, int width, char *outfile)
{
    int i, j, prd_buf;
    char name[100];
    FILE *fp;

    sprintf(name, "./LOG/%s_prd.pgm", outfile);
    fp = fileopen(name, "wb");
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				prd_buf = prd[i][j] / 64;
				if(prd_buf < 0) prd_buf = -1 * prd_buf;
				if(prd_buf > 255) prd_buf = 255;
		            putc(prd_buf, fp);
			}
    }
    fclose(fp);
    return;
}

void set_blksize(int **size, char ***qtmap, int tly, int tlx, int blksize, int level,
		 int height, int width)
{
    int y, x, bry, brx;

    if (tly >= height || tlx >= width) return;
    if (level > 0) {
	y = (tly / MIN_BSIZE) >> level;
	x = (tlx / MIN_BSIZE) >> level;
	if (qtmap[level - 1][y][x] == 1) {
	    blksize >>= 1;
	    set_blksize(size, qtmap, tly, tlx, blksize, level - 1,
			height, width);
	    set_blksize(size, qtmap, tly, tlx + blksize, blksize, level - 1,
			height, width);
	    set_blksize(size, qtmap, tly + blksize, tlx, blksize, level - 1,
			height, width);
	    set_blksize(size, qtmap, tly + blksize, tlx + blksize, blksize, level - 1,
			height, width);
	    return;
	}
    }
    if (tly + blksize <= height) bry = tly + blksize;
    else bry = height;
    if (tlx + blksize <= width) brx = tlx +  blksize;
    else brx = width;
    for (y = tly; y < bry; y++) {
	for (x = tlx; x < brx; x++) {
	    size[y][x] = level;
	}
    }
    return;
}

void print_block_size(img_t **val, int height, int width, int qtree_depth, char ***qtmap,
		      char *outfile)
{
    int y, x, bs, blksize, level;
    int **size;
    img_t **img;
    char name[100];
    FILE *fp;

    sprintf(name, "./infofile/block_size/%s_bksize.pgm", outfile);
    fp = fileopen(name, "wb");

    img = (img_t **)alloc_2d_array(height, width, sizeof(img_t));
    size = (int **)alloc_2d_array(height, width, sizeof(int));
    if (qtree_depth > 0) blksize = MAX_BSIZE;
    else blksize = BASE_BSIZE;
    level = qtree_depth;
    for (y = 0; y < height; y += blksize) {
			for (x = 0; x < width; x += blksize) {
			    set_blksize(size, qtmap, y, x, blksize, level, height, width);
			}
    }
    for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
			    bs = MIN_BSIZE << size[y][x];
			    if (y % bs == 0 || x % bs == 0) img[y][x] = 255;
			    else img[y][x] = val[y][x];
			}
    }
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
		            putc(img[y][x], fp);
			}
    }
    free(img);
    free(size);
    fclose(fp);
    return;
}

void print_color(int height, int width, int **color, char *outfile)
{
    int x, y;
    char name[100];
    FILE *fp;

    sprintf(name, "./infofile/%s_color.csv", outfile);
    fp = fileopen(name, "wb");

    for (y = 0; y < height; y ++) {
        for (x = 0; x < width; x++) {
            fprintf(fp, "%d.", color[y][x]);
	}
        fprintf(fp, "\n");
    }
    fclose(fp);
    return;
}
