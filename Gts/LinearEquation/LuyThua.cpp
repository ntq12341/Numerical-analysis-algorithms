#include "matrix.h"

matrix Norm1(matrix &a){
	int k = a.row;
	//double s = 0;
	matrix n = a;
	double m = n.mat[0][0];
	for(int i = 0;i < k;i++){
		m = (m > n.mat[i][0]) ? m : n.mat[i][0];
	}
	for(int i = 0;i < k;i++){
		n.mat[i][0] /= m;
	}
	return n;
}

int main(){
	int n,k;
	cin >> n;
	matrix a(n,n);
	matrix x(n,1);
	a.inp();
	x.inp();
	cin >> k;
	// a.out();
	// x.out();
	matrix a1 = a;
	for(int i = 1;i < 2*k;i++){
		a1 = a1*a;
	}
	a1 = a1*x;
	a1.out();
	matrix a2 = a*a1;
	matrix a3 = a*a2;
	a2.out();
	a3.out();
	matrix e1 = Norm1(a1);
	matrix e2 = Norm1(a2);
	matrix e3 = Norm1(a3);
	e1.out();
	e2.out();
	e3.out();
}