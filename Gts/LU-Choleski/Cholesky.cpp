#include "matrix.h"

void Cholesky(matrix &a,matrix &u){
	int n = a.row;
	double s;
	for(int i = 0;i < n;i++){
		for(int k = 0;k < i;k++){
			s += u.mat[k][i]*u.mat[k][i];
		}
		u.mat[i][i] = sqrt(a.mat[i][i] - s);
		s = 0;
		for(int k = i+1;k < n;k++){
			for(int h = 0;h < i;h++){
				s += u.mat[h][i]*u.mat[h][k];
			}
			u.mat[i][k] = (a.mat[i][k] - s)/u.mat[i][i];
			s = 0;
		}
	}
}

int main(){
	int n = 3;
	matrix A(n,n);
	matrix U(n,n);
	A.inp();
	Cholesky(A,U);
	U.out();
}