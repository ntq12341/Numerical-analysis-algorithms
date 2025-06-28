#include "../matrix.h"

matrix lap_don(matrix &a,matrix &b,matrix &x0,double e,const char* s){
	int m = a.row;
	int n = a.column;
	int p = b.column;
	matrix x(n,p);
	if(s == "inf"){
		double q = a.inf();
		if(q >= 1){
			cout << "ko tm dk lap" << endl;
		}
		double d;
		double si;
		int i = 0;
	    // cout << "q : " << q << endl;
	    do{
		    i++;
			x = a*x0 + b;
			d = abs((x - x0).inf())*(q/(1-q));
			si = (d/x.inf())*100;
			x0 = x;
	    }while(d > e);
	    // cout << i << '\n' << '-' << d << '-' << endl;
		// cout << '-' << si << '%' << '-' << endl;
	}else if(s == "one"){
		double q = a.one();
		if(q >= 1){
			cout << "ko tm dk lap" << endl;
		}
		matrix x(n,p);
		double d;
		double si;
		int i = 0;
	    // cout << "q : " << q << endl;
	    do{
		    i++;
			x = a*x0 + b;
			d = abs((x - x0).one())*(q/(1-q));
			si = (d/x.one())*100;
			x0 = x;
	    }while(d > e);
	    // cout << i << '\n' << '-' << d << '-' << endl;
		// cout << '-' << si << '%' << '-' << endl;
	}
	return x;
}

// int main(){
// 	matrix a(2,2);
// 	a.inp();
// 	a.out();
// }