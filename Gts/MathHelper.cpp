#include "matrix.h"
#include "complex.h"
#include "complex.cpp"

matrix PTBac2(double a,double b,double c){
	matrix ans(1,2);
	if(a == 0){
		cout << "Ko phai pt bac 2" << '\n';
		return ans;
	}else{
		double delta = b*b - 4*a*c;
		if(delta >= 0){
			ans.mat[0][0] = (-b + sqrt(delta))/(2*a);
			ans.mat[0][1] = (-b - sqrt(delta))/(2*a);
			return ans;
		}else{
			matrix ans1(2,2);
			Complex delta1 = Complex(delta,0);
			Complex sq_delta = sqrt(delta1);
			Complex b1 = Complex(b,0);
			Complex a1 = Complex(a,0);
			Complex x1 = (-b1 + sq_delta)/(2*a1);
			Complex x2 = (-b1 - sq_delta)/(2*a1);
			ans1.mat[0][0] = re(x1);
			ans1.mat[1][0] = im(x1);
			ans1.mat[0][1] = re(x2);
			ans1.mat[1][1] = im(x2);
			return ans1;
		}
	}
}
