#include "../matrix.h"
using namespace std;

//quy trình thuận(khử)
void gauss_set(matrix &a,matrix &b) {
	int r = a.row;
	int c1 = a.column;
	int c2 = b.column;
	int i = 0, j = 0;
	while (i < r && j < c1)
	{
		if (a.mat[i][j] != 0){
			for (int k = i + 1; k < r; k++) {
				double x = a.mat[k][j] / a.mat[i][j];
				for (int h = 0; h < c1; h++) {
					a.mat[k][h] -= x * a.mat[i][h];
				}
				for (int h = 0; h < c2; h++) {
					b.mat[k][h] -= x * b.mat[i][h];
				}
			}
			i++;
			j++;
		}
		else {
			if (i == r - 1) {
				break;
			}
			int t = i + 1;
			while (a.mat[t][j] == 0 && t < r) {
				t++;
			}
			if (t < r) {
				a.swap_row(t,i);
				b.swap_row(t,i);
			}
			else {
				j++;
				if (j == c1) {
					break;
				}
			}
		}
	}
}

// quy trình nghịch (thế)
matrix gauss_get_inverse(matrix &a,matrix &b) {
	int r = a.row;
	int c1 = a.column;
	int c2 = b.column;
	matrix ans(c1,c2);
	for (int y = 0; y < c2; y++) {
		for (int k = r - 1; k >= 0; k--) {
			double x = 0;
			for (int h = k + 1; h < c1; h++) {
				x += a.mat[k][h] * ans.mat[h][y];
			}
			x = b.mat[k][y] - x;
			ans.mat[k][y] = x / a.mat[k][k];
		}
	}
	return ans;
}

void gauss_jordan_inverse(matrix &a){
    int row = a.row;
    matrix b(a.row);
    vector<vector<int>> lis(2);
    int i,j;
    bool check;
	bool check1;
	double max;
	do {
		max = 0;
		check = false;
		check1 = false;
		for (int q = 0; q < row; q++) {
			if (!search_matrixlib(lis[0], q)) {
				for (int w = 0; w < row; w++) {
					if (!search_matrixlib(lis[1], w)) {
						if (abs(a.mat[q][w]) == 1) {
							i = q;
							j = w;
							lis[0].push_back(i);
							lis[1].push_back(j);
							check = true;
							check1 = true;
							break;
						}
						else {
							if (abs(a.mat[q][w]) > max) {
								max = abs(a.mat[q][w]);
								i = q;
								j = w;
								check = true;
							}
						}
					}
				}
				if (check1) {
					break;
				}
			}
		}
		if (max > 0 && !check1) {
			lis[0].push_back(i);
			lis[1].push_back(j);
		}
		if (check) {
			for (int q = 0; q < row; q++) {
				if (q != i) {
					double x = a.mat[q][j] / a.mat[i][j];
					for (int w = 0; w < row; w++) {
						a.mat[q][w] -= x * a.mat[i][w];
					}
					for (int w = 0; w < row; w++) {
						b.mat[q][w] -= x * b.mat[i][w];
					}
				}
				else {
					double x = a.mat[i][j];
					for (int w = 0; w < row; w++) {
						a.mat[q][w] /= x;
					}
					for (int w = 0; w < row; w++) {
						b.mat[q][w] /= x;
					}
				}
			}
		}
	} while (check);
	vector<int> ind(row,-1);
	for (int k = 0; k < row; k++) {
		for (int h = 0; h < 2*row; h++) {
			if (h < row) {
				if (a.mat[k][h] != 0) {
					ind[k] = h;
					break;
				}
			}
			else {
				if (b.mat[k][h - row] != 0) {
					ind[k] = h;
					break;
				}
			}
		}
	}
	for(int p = 0;p < row - 2;p++){
		bool sort = true;
		for(int q = 0;q < row - p - 1;q++){
			if(ind[q] > ind[q+1]){
				a.swap_row(q,q+1);
				b.swap_row(q,q+1);
				swap(ind[q],ind[q+1]);
				sort = false;
			}
		}
		if(sort) break;
	}
	b.out();
}

void LU (matrix &v,matrix &l,matrix &u){
    int n = v.row;
    for (int i = 0; i < n; i++)
    {
        u.mat[i][i] = 1;
        l.mat[i][0] = v.mat[i][0];
        if (i > 0) {
            u.mat[0][i] = v.mat[0][i] / v.mat[0][0];
        }
    }
    for (int i = 1; i < n; i++)
    {
        for (int f = i; f < n; f++)
        {
            int s1 = 0;
            for (int p = 0; p < i; p++)
            {
                s1 += l.mat[f][p] * u.mat[p][i];
            }
            l.mat[f][i] = v.mat[f][i] - s1;
            if (f > i) {
                int s2 = 0;
                for (int p = 0; p < n; p++)
                {
                    s2 += l.mat[i][p] * u.mat[p][f];
                }
            u.mat[i][f] = (v.mat[i][f] - s2) / l.mat[i][i];
            }
        }
    }
}

void LU_inverse(matrix &a){
	int n = a.row;
    matrix b(n);
    matrix l(n,n);
    matrix u(n,n);
    LU(a,l,u);
    matrix y(n,n);
    gauss_set(l,b);
    y = gauss_get_inverse(l,b);
    gauss_get_inverse(u,y).out();
}

void choleski(matrix &v,matrix &u){
	if(v.mat[0][0] < 0){
		cout << "ko the phan tach" << endl;
		return ;
	}
    int n = v.row;
    double sum;
    u.mat[0][0] = sqrt(v.mat[0][0]);
    for(int i = 1;i < n;i++){
        u.mat[0][i] = v.mat[0][i]/u.mat[0][0];
    }
    for(int i = 1;i < n;i++){
        sum = 0;
        for(int p = 0;p < i;p++){
            sum += (u.mat[p][i]*u.mat[p][i]);
        }
        if(v.mat[i][i] < sum){
        	cout << "ko the phan tach" << endl;
        	return;
        }
        u.mat[i][i] = sqrt(v.mat[i][i] - sum);
        for(int k = i+1;k < n;k++){
            sum = 0;
            for(int j = 0;j < i;j++){
                sum += (u.mat[j][i]*u.mat[j][k]);
            }
            u.mat[i][k] = (v.mat[i][k] - sum)/u.mat[i][i];
        }
    }
}

void choleski_inverse(matrix &a){
	int n = a.row;
	matrix b(n);
	if(!(a == a.t())){
		b = a.t()*b;
		a = a.t()*a;
	}
	matrix u(n,n);
	choleski(a,u);
	matrix y(n,n);
	matrix ut(n,n);
	ut = u.t();
	gauss_set(ut,b);
	y = gauss_get_inverse(ut,b);
	gauss_get_inverse(u,y).out();
}

void lap_don(matrix &a,matrix &b,matrix &x0,double e,const char* s){
	int m = a.row;
	int n = a.column;
	int p = b.column;
	if(s == "inf"){
		double q = a.inf();
		if(q >= 1){
			cout << "ko tm dk lap" << endl;
			return ;
		}
		matrix x(n,p);
		double d;
		double si;
		int i = 0;
	    cout << "q : " << q << endl;
	    do{
		    i++;
			x = a*x0 + b;
			d = abs((x - x0).inf())*(q/(1-q));
			si = (d/x.inf())*100;
			x0 = x;
	    }while(d > e);
	    cout << i << '\n' << '-' << d << '-' << endl;
		cout << '-' << si << '%' << '-' << endl;
		x.out();
	}else if(s == "one"){
		double q = a.one();
		if(q >= 1){
			cout << "ko tm dk lap" << endl;
			return ;
		}
		matrix x(n,p);
		double d;
		double si;
		int i = 0;
	    cout << "q : " << q << endl;
	    do{
		    i++;
			x = a*x0 + b;
			d = abs((x - x0).one())*(q/(1-q));
			si = (d/x.one())*100;
			x0 = x;
	    }while(d > e);
	    cout << i << '\n' << '-' << d << '-' << endl;
		cout << '-' << si << '%' << '-' << endl;
		x.out();
	}
}

void gauss_seiden_inverse(matrix &a,double e){
	int m = a.row;
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
	matrix n(m,m);
	M = (d-l).inverse()*u;
	cout << M.inf() << endl;
	n = (d-l).inverse();
	lap_don(M,n,a,e,"inf");
}

void vien_quanh(matrix &a){
	int n = a.row;
	// m = a^t * a 
	matrix m(n,n);
	m = a.t()*a;
	int i = 1;
	matrix m1;
	while(i <= n){
		matrix m2(i,i);
		if(i == 1){
			m2.mat[0][0] = 1/m.mat[0][0];
		}else{
			matrix a1(i-1,1); // ma trận cột ngoài của A
			matrix a2(1,i-1); // ma trận hàng dưới của A 
			for(int j = 0;j < i-1;j++){
				a1.mat[j][0] = m.mat[j][i-1];
			}
			a2 = a1.t();
			double a3 = m.mat[i-1][i-1]; // phần tử dưới cùng bên phải của A

			double b3 = 1/(a3 - ((a2*m1*a1).mat[0][0])); // phần tử dưới cùng bên phải của A^-1
			matrix b0(i-1,i-1); // ma trận trên cùng bên trái của A^-1
			matrix b1(i-1,1); // ma trận cột ngoài của A^-1
			matrix b2(1,i-1); // ma trận hàng dưới của A^-1
			matrix e(i-1); // ma trận đơn vị cấp i-1
			// b1 = -(m1*a1)/b3;
			// b2 = -(a2*m1)/b3;
			// b0 = m*(e + (a1*a2*m1)/b3);

			b0 = m1*(e + (a1*a2*m1)*b3);
			b1 = (m1*a1)*(-b3);
			b2 = (a2*m1)*(-b3);
			for(int j = 0;j < i-1;j++){
				m2.mat[j][i-1] = b1.mat[j][0];
				m2.mat[i-1][j] = b2.mat[0][j];
				for(int k = 0;k < i-1;k++){
					m2.mat[j][k] = b0.mat[j][k];
				}
			}
			m2.mat[i-1][i-1] = b3;
		}
		m1 = m2;
		i++;
	}
	matrix ans(n,n);
	ans = m1*a.t();
	ans.out();
}

int main(){
	matrix a(2,2);
	a.inp();
	matrix b = a;
	matrix e(2);
	gauss_set(b,e);
	gauss_get_inverse(b,e).out();
	a.inverse().out();
	return 0;
}