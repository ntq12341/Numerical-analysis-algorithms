#include <iostream>
#include <math.h>

#define PI 3.1459

using namespace std;

// hàm cần tìm nghiệm
double f(double x){
	return x - sin(3*x);
}

int main(){
	double x,d,M,m,eps;
	//cin >> x >> d >> M >> m >> eps;
	x = PI/6;
	d = PI/4;
	M = 3;
	m = 1;
	eps = 1e-3;
	int i = 0;
	double y,e;
	do{
		i++;
		y = x;
		x = x - (f(x)*(d - x))/(f(d)-f(x));
		e = (M - m)*abs(y - x)/m;
		cout << i << '|' << x << '\n';
	}while(e > eps);
}