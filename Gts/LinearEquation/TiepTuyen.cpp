#include <iostream>
#include <math.h>
using namespace std;

double f(double x){
	return 1;
}

double f1(double x){
	return 1;
}

int main(){
	double x,M,m,eps;
	cin >> x >> M >> m >> eps;
	double e,y;
	int i = 0;
	do{
		cout << i << '|' << x << '\n';
		y = x;
		x = x - f(x)/f1(x);
		e = (M/(2*m))*abs(x-y);
	}while(e > eps);
}