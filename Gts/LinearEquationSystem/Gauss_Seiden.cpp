#include "LapDon.cpp"

matrix lap_seiden(matrix &a,matrix &b,matrix &x0,double e,const char* s){
	int m = a.row;
	vector<vector<double>> l0(m,vector<double> (m,0));
	vector<vector<double>> u0(m,vector<double> (m,0));
	for(int i = 0;i < m;i++){
		for(int j = 0;j < m;j++){
			if(i > j){
				l0[i][j] = a.mat[i][j];
			}else if(i < j){
				u0[i][j] = a.mat[i][j];
			}
		}
	}
	matrix l(l0);
	matrix u(u0);
	matrix i(m);
	matrix M = (i-l).inverse()*u;
	matrix x = lap_don(M,b,x0,e,"inf");
	return x;
}

matrix gauss_seiden(matrix &a,matrix &b,double e){
	int m = a.row;
	matrix x_0(m,1);
		for(int i = 0;i < m;i++){
			x_0.mat[i][0] = 1;
		}
	vector<vector<double>> d0(m,vector<double> (m,0));
	vector<vector<double>> l0(m,vector<double> (m,0));
	vector<vector<double>> u0(m,vector<double> (m,0));
	for(int i = 0;i < m;i++){
		for(int j = 0;j < m;j++){
			if(i > j){
				l0[i][j] = -a.mat[i][j];
			}else if(i < j){
				u0[i][j] = -a.mat[i][j];
			}else{
				d0[i][j] = a.mat[i][j];
			}
		}
	}
	matrix d(d0);
	matrix l(l0);
	matrix u(u0);
	matrix M(m,m);
	M = (d-l).inverse()*u;
	//cout << M.inf() << endl;
	matrix x = lap_don(M,b,x_0,e,"inf");
	return x;
}

int main(){

}