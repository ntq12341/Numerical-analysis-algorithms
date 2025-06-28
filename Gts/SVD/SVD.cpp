#include "XuongThang.cpp"

void SVD(matrix &a,matrix &U,matrix &V,matrix &Sigma);

int main(){
	// nhập input 
	int m ; 
	int n; 
	cin >> m >> n; // cỡ ma trận 
	matrix a(m,n); // ma trận A
	a.inp();
	a.out();

	matrix Sigma;
	matrix V;
	matrix U;

	SVD(a,U,V,Sigma);
	U.out();
	Sigma.out();
	V.out();
}

void SVD(matrix &a,matrix &U0,matrix &V0,matrix &Sigma0){
	int m = a.row;
	int n = a.column; 
	matrix b = a.t()*a;

	vector<double> e_values1;
	vector<matrix> e_vectors1;
	vector<Complex> e_values_cpx1;
	vector<vector<Complex>> e_vectors_cpx1;

	XuongThang(b,e_values1,e_vectors1,e_values_cpx1,e_vectors_cpx1);

	int r = -1;
	for(int i = b.row-1;i >= 0;i--){
		r = i+1;
		if(fabs(e_values1[i]) > 1e-3){
			break;
		}
	}
	// Tính ma trận giá trị kì dị
	matrix Sigma(r,r);
	for(int i = 0;i < r;i++){
		Sigma.mat[i][i] = sqrt(e_values1[i]);
	} 
	Sigma0 = Sigma;

	// Tính ma trận kì dị phải
	matrix V(n,r);
	for(int col = 0;col < r;col++){
		for(int ro = 0;ro < n;ro++){
			V.mat[ro][col] = e_vectors1[col].mat[ro][0];
		}
	}
	V0 = V;

	b = a*a.t();

	vector<double> e_values2;
	vector<matrix> e_vectors2;
	vector<Complex> e_values_cpx2;
	vector<vector<Complex>> e_vectors_cpx2;

	XuongThang(b,e_values2,e_vectors2,e_values_cpx2,e_vectors_cpx2);

	r = -1;
	for(int i = e_values2.size()-1;i > 0;i--){
		r = i+1;
		if(fabs(e_values2[i]) > 1e-3){
			break;
		}
	}
	// //Tính ma trận giá trị kì dị của phân tích đầy đủ
	// matrix Sigma1(m,n);
	// for(int i = 0;i < r;i++){
	// 	Sigma1.mat[i][i] = sqrt(e_values2[i]);
	// } 
	// Sigma1.out();

	// Tính ma trận kì dị trái
	matrix U(m,r);
	for(int col = 0;col < r;col++){
		for(int ro = 0;ro < m;ro++){
			U.mat[ro][col] = e_vectors2[col].mat[ro][0];
		}
	}
	U0 = U;

	// for(int i = 0;i < e_values1.size();i++){
	// 	cout << e_values1[i] << ' ';
	// }
	// cout << '\n';
	
	// for(int i = 0;i < e_vectors1.size();i++){
	// 	e_vectors1[i].out();
	// }
}