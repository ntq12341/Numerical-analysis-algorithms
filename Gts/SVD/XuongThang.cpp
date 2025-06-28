#include "LuyThua.cpp"

void XuongThang(matrix &a,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values_cpx,vector<vector<Complex>> &e_vectors_cpx);

// int main(){
// 	// nhập input 
// 	int n; cin >> n; // cỡ ma trận 
// 	matrix a(n,n); // ma trận A
// 	a.inp();

// 	vector<double> e_values;
// 	vector<matrix> e_vectors;
// 	vector<Complex> e_values_cpx;
// 	vector<vector<Complex>> e_vectors_cpx;

// 	XuongThang(a,e_values,e_vectors,e_values_cpx,e_vectors_cpx);

// 	for(int i = 0;i < e_values.size();i++){
// 		cout << e_values[i] << ' ';
// 	}
// 	cout << '\n';
	
// 	for(int i = 0;i < e_vectors.size();i++){
// 		e_vectors[i].out();
// 	}
	
// }

void XuongThang(matrix &a,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values_cpx,vector<vector<Complex>> &e_vectors_cpx){
	int n = a.row;

	// số lần lặp k 
	int k = 10;

	int cas = 0;

	int j = 0;

	// vector v là tổ hợp tuyến tính bất kì của các vector riêng
	matrix v(n,1);
	for(int i = 0;i < n;i++){
		v.mat[i][0] = 1;
	}

	vector<double> et_values;
	vector<matrix> et_vectors;
	vector<Complex> et_values_cpx;
	vector<vector<Complex>> et_vectors_cpx;

	matrix q = a;

	bool ngh_kep = false;

	while(j < n){
		LuyThua(n,q,v,k,cas,e_values,e_vectors,e_values_cpx,e_vectors_cpx);
		j += cas;
		if(cas == 1 || cas == 2){
			q = q.t();
			LuyThua(n,q,v,k,cas,et_values,et_vectors,et_values_cpx,et_vectors_cpx);
			q = q.t();
			matrix vj = e_vectors[e_vectors.size() - 1];
			matrix wj = et_vectors[et_vectors.size() - 1];
			matrix tmp_mat = wj.t()*vj;
			double tmp = (e_values[e_values.size()-1]/ tmp_mat.mat[0][0]);
			q = q - vj*wj.t()*tmp;
			if(ngh_kep){
				e_vectors.pop_back();
				e_values.pop_back();
				et_vectors.pop_back();
				et_values.pop_back();
				e_vectors_cpx.pop_back();
				e_values_cpx.pop_back();
				et_vectors_cpx.pop_back();
				et_values_cpx.pop_back();
				j--;
				ngh_kep = false;
			}
			if(cas == 2){
				ngh_kep = true;
			}
		}else{
			// vector<Complex> vj = e_vectors_cpx[e_vectors_cpx.size() - 1];
			// vector<Complex> wj = et_vectors_cpx[e_vectors_cpx.size() - 1];
			// Complex tmp_cpx = Complex(0,0);
			// for(int i = 0;i < vj.size();i++){
			// 	tmp_cpx += vj[i]*wj[i];
			// }
			// Complex tmp_ev = e_values_cpx[e_values_cpx.size()-1];

			// Complex tmp = tmp_ev/tmp_cpx;
			// vector<vector<Complex>> q1(vj.size());
			// for (int i = 0; i < vj.size(); i++) {
			// 	q1[i].resize(vj.size());
			// 	for (int y = 0; y < vj.size(); y++) {
			// 		q1[i][y] = vj[i]*wj[y]*tmp;
			// 		Complex z = Complex(q.mat[i][j],0);
			// 		q1[i][y] = z-q1[i][y];
			// 		cout << q1[i][y] << ' ';
			// 	}
			// 	cout << '\n';
			// }
			cout << "\nNghiem phuc !!!!!\n";
			break;
		}
	}
}