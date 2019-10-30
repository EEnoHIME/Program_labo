/*********************************************
 * 非線形方程式の解法 ( ニュートン法 )
 *********************************************/
#include "equation.h"
#include <iostream>  // for cout
#include <math.h>    // for fabs()
#include <stdio.h>   // for printf()

using namespace std;

extern value_for_eq coeff;

//double ans;
double a;
double b = 0;

// 方程式定義
double F(double x){
    double y;
	y = (coeff.sqq[1]+coeff.sqq[2]-coeff.sqq[3]-coeff.sqq[0])*x*x*x*x+
        (coeff.sqd[1]+coeff.sqd[2]-coeff.sqd[3]-coeff.sqd[0])*x*x+
        (coeff.c[1]+coeff.c[2]-coeff.c[3]-coeff.c[0]);
    return y;
}
// f(x) の x における傾き ( f(x) を１回微分 )
double  G(double x){
    double y; 
    y = 4*(coeff.sqq[1]+coeff.sqq[2]-coeff.sqq[3]-coeff.sqq[0])*x*x*x+
        2*(coeff.sqd[1]+coeff.sqd[2]-coeff.sqd[3]-coeff.sqd[0])*x;
    return y;
}

/*
 * 計算クラス
 */
class Calc
{
    // 各種定数
    static const double eps = 1e-10;  // 打ち切り精度
    static const int  limit = 10000000;     // 打ち切り回数

    // 各種変数
    double x, dx,ans;  // x, dx 値
    int k;         // LOOP インデックス

    public:
        // 非線形方程式を解く（ニュートン法）
        double calcNonlinearEquation();
};

/*
 * 非線形方程式を解く（ニュートン法）
 */
double Calc::calcNonlinearEquation()
{
  // x 初期値設定
  
  x = 4.0;

  // 打ち切り回数 or 打ち切り誤差になるまで LOOP
  for (k = 1; k <= limit; k++) {
    dx = x;
    x = x - F(x) / G(x);
    if (fabs(x - dx) / fabs(dx) < eps) {
      //printf("x=%f\n", x);
	 
	  return x;
      break;
    }
  }

  // 収束しなかった場合
  if (k > limit)
    cout << "収束しない" << endl;
}


double solve_equation(){
	double ans_a;
    try
    {
        // 計算クラスインスタンス化
        Calc objCalc;

        // 非線形方程式を解く（ニュートン法）
        ans_a = objCalc.calcNonlinearEquation();
    }
    catch (...) {
        cout << "例外発生！" << endl;
        return -1;
    }

    // 正常終了
    return ans_a;
}

/*
 * メイン処理
 
 int main()
{
	 double a;
	 coeff[8] = 3;
	 coeff[7] = -8;
	 coeff[6] = 5;
	 coeff[5] = -9;
	 coeff[4] = 2;
	 coeff[3] = -7;
	 coeff[2] = 3;
	 coeff[1] = -1;
	 coeff[0] = 6;

    a = solve_equation();
	cout << a << endl;
	return 0;
}
*/

 
