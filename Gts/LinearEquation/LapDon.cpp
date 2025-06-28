#include <iostream>
#include <math.h>
using namespace std;

double p(double x){
	return 1;
}

int main(){
	double x,q,eps;
	cin >> x >> q >> eps;
	double e,y;
	int i = 0;
	do{
		cout << i << '|' << x << '\n';
		y = x;
		x = p(x);
		e = (q/(1-q))*abs(y-x);
	}while(e > eps);
}