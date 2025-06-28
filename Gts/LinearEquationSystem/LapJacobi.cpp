#include "LapDon.cpp"

matrix LapJacobi(matrix &a,matrix &b,double eps){
	int n = a.row;
	matrix d_i(n);
	for(int i = 0;i < n;i++){
		d_i.mat[i][i] = 1/a.mat[i][i];
	}
	matrix i(n);
	if(a.row_dom()){
		matrix x_0(n,1);
		for(int i = 0;i < n;i++){
			x_0.mat[i][0] = 1;
		}
		matrix k = i - d_i*a;
		matrix h = d_i*b;
		matrix x = lap_don(k,h,x_0,eps,"inp"); 
		return x;
	}else if(a.column_dom()){
		matrix x_0(n,1);
		for(int i = 0;i < n;i++){
			x_0.mat[i][0] = 1;
		}
		matrix d(n);
		double max_ai = 0 , min_ai = -1;
		for(int i = 0;i < n;i++){
			d.mat[i][i] = a.mat[i][i];
			max_ai = (max_ai > fabs(a.mat[i][i])) ? max_ai : fabs(a.mat[i][i]);
			min_ai = (min_ai < fabs(a.mat[i][i])) ? min_ai : fabs(a.mat[i][i]);
		}
		matrix k1 = i - a*d_i;
		matrix k = d_i*k*d;
		matrix h = d_i*b;
		eps *= max_ai/min_ai;
		matrix x = lap_don(k,h,x_0,eps,"one");
		return x;
	}else{
		cout << "Ma tran ko cheo troi\n";
		matrix x(1,1);
		return x;
	}
}

int main(){
	matrix a(2,2);
	a.inp();
	a.out();
}