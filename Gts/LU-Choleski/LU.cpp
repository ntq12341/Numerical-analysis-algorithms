#include "matrix.h"

void LU(matrix &a,matrix &l,matrix &u){
	int n = a.row;
	double s;
	for(int i = 0;i < n;i++){
		for(int j = 0;j < n;j++){
			if(i >= j){
				for(int k = 0;k < j;k++){
					s += l.mat[i][k]*u.mat[k][j];
				}
				l.mat[i][j] = a.mat[i][j] - s;
				s = 0;
			}else{
				for(int k = 0;k < i;k++){
					s += l.mat[i][k]*u.mat[k][j];
				}
				u.mat[i][j] = (a.mat[i][j] - s)/l.mat[i][i] ;
				s = 0;
			}
		}
	}
}
int main(){
	int m,n;
	cin >> m >> n;
	matrix A(m,n);
	matrix L(m);
	matrix U(m);
	A.inp();
	LU(A,L,U);
	A.out();
	L.out();
	U.out();
}