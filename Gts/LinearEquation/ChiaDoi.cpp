#include <iostream>
#include <iomanip>
#include <math.h>
#define PI 3.1415926
using namespace std;

// hàm cần tìm nghiệm
double f(double x){
	return 3*sin(x) + x*x*x - 8*x*x + 8*x + 1;
}

int main(){
	cout << setprecision (10);
	// khoảng cách ly nghiệm
	double a,b;
	a = 0;
	b = 2;
	// sai sô
	double eps;
	eps = 0.5e-7;
	// nghiệm cần tìm
	double c;
	// số lần lặp
	int l = 0;
	//cin >> a >> b >> eps;
	do{
		cout << l << '|' << a << '-' << b << '-';
		l++;
		c = (a+b)/2;
		//nghiệm chính xác
		if(f(c) == 0){
			cout << "Exactly solution : " << c << '\n';
			break;
		}
		if(f(c)*f(a) > 0){
			a = c;
			cout << 'a' ;
		}
		if(f(c)*f(b) > 0){
			b = c;
			cout << 'b' ;
		}
		cout << '\n';
	}while(abs(a-b) > eps);
	cout << c << '\n' ;
	double dV = (1/6)*c*c*c*0.5e-7 + 0.5*PI*c*c*eps;
	cout << dV ;
	return 0;
}