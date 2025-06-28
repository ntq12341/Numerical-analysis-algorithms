#include "../MathHelper.cpp"

matrix akv(matrix &a,matrix &v,int k);
matrix div(matrix &a,matrix &b);
matrix norm_mat(matrix &a);
bool check_vr(matrix &vr);
void case1(matrix &a1,matrix &vr);
void case2(matrix &a1,matrix &a2,matrix &vr);
void case3(matrix &a1,matrix &a2,matrix &a3);
void LuyThua(int n,matrix &a,matrix &v,int k,int &cas,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values1, vector<vector<Complex>> &e_vectors1);

// int main(){
// 	// nhập input 
// 	int n; cin >> n; // cỡ ma trận 
// 	matrix a(n,n); // ma trận A
// 	a.inp();

// 	// số lần lặp k 
// 	int k = 20;
// 	vector<double> e_values;
// 	vector<matrix> e_vectors;
// 	vector<Complex> e_values_cpx;
// 	vector<vector<Complex>> e_vectors_cpx;
// 	// vector v là tổ hợp tuyến tính bất kì của các vector riêng
// 	matrix v(n,1);
// 	for(int i = 0;i < n;i++){
// 		v.mat[i][0] = i+1;
// 	}
// 	int cas ;
// 	LuyThua(n,a,v,k,cas,e_values,e_vectors,e_values_cpx,e_vectors_cpx);

// 	for(int i = 0;i < e_values.size();i++){
// 		cout << e_values[i] << ' ';
// 	}
// 	cout << '\n';
	
// 	for(int i = 0;i < e_vectors.size();i++){
// 		e_vectors[i].out();
// 	}

//  }

matrix akv(matrix &a,matrix &v,int k){
	matrix ans = v;
	for(int i = 0;i < k;i++){
		ans = a*ans;
	}
	return ans;
}

matrix div(matrix &a,matrix &b){
	matrix ans = a;
	for(int i = 0;i < a.row;i++){
		ans.mat[i][0] /= b.mat[i][0]; 
	}
	return ans;
}

matrix norm_mat(matrix &a){
	matrix ans = a;
	double s = 0;
	for(int i = 0;i < ans.row; i++){
		s += ans.mat[i][0]*ans.mat[i][0];
	}
	s = sqrt(s);
	for(int i = 0;i < ans.row; i++){
		ans.mat[i][0] /= s;
	}
	return ans ;
}

bool check_vr(matrix &vr){
	for(int i = 1;i < vr.row;i++){
		if(fabs(vr.mat[i][0] - vr.mat[0][0]) > 1e-3){
			return false;
		}
	}
	return true;
}

void case1(matrix &a1,matrix &vr,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values1,vector<vector<Complex>> &e_vectors1){
	// Tìm Giá trị riêng
	// cout << "e-value : " << '\n';
	// cout << vr.mat[0][0] << '\n';
	e_values.push_back(vr.mat[0][0]);
	// Tìm vector riêng tương ứng TH1 
	matrix e_v = norm_mat(a1);
	// cout << "e-vector : " << '\n';
	// e_v.out();
	e_vectors.push_back(e_v);
	Complex z = Complex(0,0);
	vector<Complex> v_z;
	v_z.push_back(z);

	e_values1.push_back(z);
	e_vectors1.push_back(v_z);
}

void case2(matrix &a1,matrix &a2,matrix &vr,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values1,vector<vector<Complex>> &e_vectors1){
	// Tìm giá trị riêng
	double e_value1 = sqrt(vr.mat[0][0]);
	double e_value2 = -e_value1;
	// cout << "e-value : " << '\n';
	// cout << e_value1 << '\n';
	// cout << e_value2 << '\n';
	e_values.push_back(e_value1);
	e_values.push_back(e_value2);
	// Tìm vector riêng
	matrix e_v1 = a1*e_value1 + a2;
	matrix e_v2 = a1*e_value1 - a2;
	e_v1 = norm_mat(e_v1);
	e_v2 = norm_mat(e_v2);	
	// cout << "e-vector : " << '\n';
	// e_v1.out();
	// e_v2.out();
	e_vectors.push_back(e_v1);
	e_vectors.push_back(e_v2);

	Complex z = Complex(0,0);
	vector<Complex> v_z;
	v_z.push_back(z);
	for(int i = 0;i < 2;i++){
		e_values1.push_back(z);
		e_vectors1.push_back(v_z);
	}
}

void case3(matrix &a1,matrix &a2,matrix &a3,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values1,vector<vector<Complex>> &e_vectors1){
	// Tìm hệ số và giải az^2 + bz +c = 0
	matrix tmp(2,2);
	tmp.mat[0][0] = a2.mat[0][0];
	tmp.mat[1][0] = a2.mat[1][0];
	tmp.mat[0][1] = a1.mat[0][0];
	tmp.mat[1][1] = a1.mat[1][0];
	double a = tmp.det();
	tmp.mat[0][0] = a3.mat[0][0];
	tmp.mat[1][0] = a3.mat[1][0];
	double b = -tmp.det();
	tmp.mat[0][1] = a2.mat[0][0];
	tmp.mat[1][1] = a2.mat[1][0];
	double c = tmp.det();
	matrix x = PTBac2(a,b,c);
	// Tìm 2 giá trị riêng
	Complex e_v1 = Complex(x.mat[0][0],x.mat[1][0]);
	Complex e_v2 = Complex(x.mat[0][1],x.mat[1][1]);
	// cout << "e-value : " << '\n';
	// cout << e_v1 << '\n';
	// cout << e_v2 << '\n' << '\n';
	e_values1.push_back(e_v1);
	e_values1.push_back(e_v2);
	// Tìm 2 vector riêng 
	matrix v1_r = a2 - a1*re(e_v2);
	matrix v1_i = a1*im(e_v2);
	double s = 0;
	vector<Complex> v1(v1_r.row) ;
	for(int i = 0;i < v1_r.row;i++){
		v1[i] = Complex(v1_r.mat[i][0],v1_i.mat[i][0]);
		s += norm(v1[i]); 
	}
	s = sqrt(s);
	for(int i = 0;i < v1_r.row;i++){
		v1[i] = v1[i]/s;
		v1_r.mat[i][0] = re(v1[i]);
		v1_i.mat[i][0] = im(v1[i]);
	}
	matrix v2_r = a2 - a1*re(e_v1);
	matrix v2_i = a1*im(e_v1);
	s = 0;
	vector<Complex> v2(v2_r.row) ;
	for(int i = 0;i < v2_r.row;i++){
		v2[i] = Complex(v2_r.mat[i][0],v2_i.mat[i][0]);
		s += norm(v2[i]); 
	}
	s = sqrt(s);
	for(int i = 0;i < v1_r.row;i++){
		v2[i] = v2[i]/s;
		v2_r.mat[i][0] = re(v2[i]);
		v2_i.mat[i][0] = im(v2[i]);
	}

	// cout << "e-vector : " << '\n';
	// for(int i = 0;i < v1.size();i++){
	// 	cout << v1[i] << '\n';
	// }
	// cout << '\n';
	// for(int i = 0;i < v1.size();i++){
	// 	cout << v2[i] << '\n';
	// }
	e_vectors1.push_back(v1);
	e_vectors1.push_back(v2);

	matrix tmp_1(1,1);
	for(int i = 0;i < 2;i++){
		e_values.push_back(0);
		e_vectors.push_back(tmp_1);
	}
}

void LuyThua(int n,matrix &a,matrix &v,int k,int &cas,vector<double> &e_values,vector<matrix> &e_vectors,vector<Complex> &e_values1,vector<vector<Complex>> &e_vectors1){
	// Tính 2 ma trận A^k*v và A^(k+1)*v 
	matrix a1 = akv(a,v,k) ;
	matrix a2 = akv(a,v,k+1);

	// Tính ma trận v_r = A^(k+1)*v / A^k*v
	matrix v_r = div(a2,a1);
	// Kiểm tra v_r để phân trường hợp
	if(check_vr(v_r)){
		// Trường hợp 1 
		case1(a1,v_r,e_values,e_vectors,e_values1,e_vectors1);
		cas = 1;
	}
	else{
		// Tính ma trận A^(k+2)*v
		matrix a3 = akv(a,v,k+2);
		// Tính ma trận v_r1 = A^(k+2)*v / A^k*v
		matrix v_r1 = div(a3,a1);
		// Kiểm tra v_r1 để phân trường hợp
		if(check_vr(v_r1)){
			// Trường hợp 2
			case2(a1,a2,v_r1,e_values,e_vectors,e_values1,e_vectors1);
			cas = 2; 
		}
		else{
			// Trường hợp 3
			case3(a1,a2,a3,e_values,e_vectors,e_values1,e_vectors1);
			cas = 3;
		}
 	}
}
