#include "../matrix.h"

void case1(matrix &a,int k,matrix &t); 
void case2(matrix &a,int k,matrix &t); 
void case3(matrix &a,int k,matrix &t); 
double ChiaDoi(vector<double> &v,double a,double b);
double f(vector<double> &v);
vector<double> GiaiPhuongTrinhDaThuc(vector<double> &v,double min_x,double max_x);
vector<matrix> e_vecs(vector<double> &e_vals,matrix &t,int n);
matrix norm_mat(matrix &a);

int main(){
	cout << setprecision(5);

	int n ;
	cin >> n;
	matrix a(n,n);
	a.inp();
	matrix t;
	for(int i = 1;i < n;i++){
		if(a.mat[n-i][n-i-1] != 0){
			case1(a,i,t);
		}else{
			case2(a,i,t);
		}
	}
	a.out();

	matrix t_i = t.inverse();
	vector<double> p(n+1);
	for(int i = 1;i < n+1;i++){
		p[i] = -a.mat[0][i-1];
	}
	p[0] = 1;
	vector<double> e_values = GiaiPhuongTrinhDaThuc(p,-20,10);
	vector<matrix> e_vectors = e_vecs(e_values,t_i,n);

	for(int i = 0;i < e_values.size();i++){
		cout << e_values[i] << ' ';
	}
	cout << "\n\n";
	for(int i = 0;i < e_vectors.size();i++){
		e_vectors[i].out();
	}
}

void case1(matrix &a,int k,matrix &t){
	int n = a.row;
	matrix m(n);
	m.mat[n-k-1] = a.mat[n-k];
	matrix m_i(n);
	for(int j = 0;j < n;j++){
		m_i.mat[n-k-1][j] = -a.mat[n-k][j]/a.mat[n-k][n-k-1];
	}
	m_i.mat[n-k-1][n-k-1] = 1/a.mat[n-k][n-k-1];
	a = m*a*m_i;
	if(k == 1){
		t = m;
	}else{
		t = m*t;
	}
	// cout << "-\n";
	// t.out();
}

void case2(matrix &a,int k,matrix &t){
	int n = a.row;
	int h = -1;
	matrix c(n);
	for(int j = n-k-2;j >=0;j--){
		if(a.mat[n-k][j] != 0){
			h = j;
			break;
		}
	}
	if(h == -1){
		case3(a,k,t);
		return;
	}
	c.mat[n-k-1][n-k-1] = 0;
	c.mat[h][n-k-1] = 1;
	c.mat[h][h] = 0;
	c.mat[n-k-1][h] = 1;

	a = c*a*c;
	if(k == 1){
		t = c;
	}else{
		t = c*t;
	}
	case1(a,k,t);
}

void case3(matrix &a,int k,matrix &t){
	int n = a.row;
	matrix s1(n);
	matrix s2(n);
	for(int q = 1;q < k;q++){
		for(int j = 0;j < n-k;j++){
			s1.mat[j][n-k+q] = -a.mat[j][n-k+q-1];
			s2.mat[j][n-k+q] = a.mat[j][n-k+q-1];
		}
		a = s1*a*s2;
		for(int j = 0;j < n-k;j++){
			s1.mat[j][n-k+q] = 0;
			s2.mat[j][n-k+q] = 0;
		}
	}
	matrix b(n-k,n-k);
	for(int i = 0;i < n-k;i++){
		for(int j = 0;j < n-k;j++){
			b.mat[i][j] = a.mat[i][j];
		}
	}
	matrix f(k,k);
	for(int i = 0;i < k;i++){
		for(int j = 0;j < k;j++){
			f.mat[i][j] = a.mat[i][j];
		}
	}
	b.out();
	f.out();
}


double f(vector<double> &v,double x){
	double ans = 0;
	int n = v.size();
	for(int j = 0;j < n;j++){
		ans += v[j]*pow(x,n-j-1);
	}
	return ans;
}

double ChiaDoi(vector<double> &v,double a,double b){
	// sai sô
	double eps;
	eps = 0.5e-3;
	// nghiệm cần tìm
	double c;
	
	do{
		c = (a+b)/2;
		//nghiệm chính xác
		if(f(v,c) == 0){
			return c;
		}
		if(f(v,c)*f(v,a) > 0){
			a = c;
		}
		if(f(v,c)*f(v,b) > 0){
			b = c;
		}
	}while(abs(a-b) > eps);
	return c;
}

vector<double> GiaiPhuongTrinhDaThuc(vector<double> &v,double min_x,double max_x){
	vector<double> ans;
	double step = min_x;
	while(step < max_x){
		if(f(v,step)*f(v,step+0.5) < 0){
			ans.push_back(ChiaDoi(v,step,step+0.5));
		}
		step += 0.5;
	}
	return ans;
}

vector<matrix> e_vecs(vector<double> &e_vals,matrix &t,int n){
	vector<matrix> ans;
	for(int i = 0;i < e_vals.size();i++){
		matrix v(n,1);
		for(int j = 0;j < n;j++){
			v.mat[j][0] = pow(e_vals[i],n-j-1);
		}
		v = t*v;
		v = norm_mat(v);
		ans.push_back(v);
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
