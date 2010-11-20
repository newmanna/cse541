// Author: Nathan Newman
// Date: 
// CSE 541, LAB 2, Malkiman MW 5:30

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;
//------------------
typedef double item;
typedef std::vector<item> Array;
typedef std::vector<Array> Matrix;
//------------------
void Swap (item& a, item& b);
void Trapezoid (item f(item), item a, item b, item h, int N, item& sum);
void Romberg (ostream& outs, item f(item), item a, item b, int N, std::vector< std::vector<item> >& r);
item Simpson (ostream& out, item f(item), item a, item b, item epsilon, int level, int level_max);
//------------------
// function 1 is g(x) = 2^x,
//               g'(x) = 2^x * log(2),
//               g''(x) = 2^x * (log(2))^2, used in 1.a
inline item g(item x) {return (pow(2, x));}
inline item frst_deriv_g(item x) {return (pow(2,x) * pow(log(2), 1));}
inline item sec_deriv_g(item x) {return (pow(2,x) * pow(log(2), 2));}

// function 2.a is a helper function to save typing 
// e(x) = e^-((ln(x))^2/2)
inline item e(item x) {return(exp(-1*pow(log(x), 2)/2));}
// function 2 is h(x) = e(x)/x,
//               h'(x) = ,
//               h''(x) = ,
inline item h(item x) {return (e(x)/x);}
inline item frst_deriv_h(item x) {return (-1*(e(x)/pow(x,2))*(log(x) - 1));}
inline item sec_deriv_h(item x) {return ((e(x)/pow(x,3))*(pow(log(x),2)+log(x)-1));}

// function 3 is l(x) = x^(1/3),
inline item l(item x) {return (pow(x, 1.0/3.0));}

// function 4 is m(x) = t^2/(t^4 + 1)
inline item m(item x) {return (pow(x, 2)/(pow(x,4)+1));}
inline item frst_deriv_m(item x) {return(((2*x)/(pow(x,4)+1)-((4*pow(x,5))/pow(pow(x,4)+1,2))));}
inline item sec_deriv_m(item x) {return((-28*pow(x,4)/pow(pow(x,4)+1,2))+(2/pow(x,4)+1))+(32*pow(x,8)/pow(pow(x,4)+1,3));}

//function 5 is p(x) = (1/x)*ln(x)*(1/(t+1))
inline item p(item x) {return((1/x)*log(x)*(1/(x+1)));}
int main ()
{
	ofstream outfile("test.txt");
	item upper_bound, lower_bound, sum, n;
	int N;
	//---------------
	cout << "Problem #1, Part a..." << endl;
	// First check g''(x) = 0, and then compute max (g''(upper_bound), g''(lower_bound), g''(xi)):
	// Since g'''(x) != 0, xi doesn't exist, therefore just need to look at upper_bound
	//	and lower bound. Since g''(0) < g''(2) then xi = 2 which is the upper bound.

	upper_bound = 2; 
	lower_bound = 0; 
	item epsilon = pow(10, -8);

	item h = sqrt(12 / (fabs((upper_bound - lower_bound) * sec_deriv_g(upper_bound))) * epsilon);

	n = (upper_bound - lower_bound) / h;
	n = ceil (n);
	N = (int) n;
	Trapezoid (g, lower_bound, upper_bound, h, N, sum);
	cout << "Integeral approximation of g(x) = 2^x from 0 to 2 = " << sum << ", by Trapeziod Rule." << endl;
	cout << endl << endl;
    	//--------------------------------------
	cout << "Problem #1, Part b..." << endl;
	// After the appropiate manipulations and subtitutions of x = -ln(t) you come upwith the function
	// h(x) = e(x) / x where e(x) = e^-(ln(x)^2/2).  and the a change a boundaries implies that the
	// new lower bound  = 0 and the new upper bound = 1.  Knowing that the lower bound of 0 when plugged 
	// into the function would yield and undefined answer it is safe to use the limit of h(x) as x->0
	// which is 0.  So you can use a number very close to zero instead of zero for the lower bound. 
	// Something like lower_bound = .0000000001.  Also using the third derivative it is ieasy to see that
	// h''(x) has a maximu at approx x = .0999809, therefore xi = .0999809
	upper_bound = 1; 
	lower_bound = pow(10, -11); 
	epsilon = pow(10, -10);

	h = sqrt(12 / (fabs((upper_bound - lower_bound) * sec_deriv_h(.0999809))) * epsilon);

	n = (upper_bound - lower_bound) / h;
	n = ceil (n);
	N = (int) n;
	Trapezoid (g, lower_bound, upper_bound, h, N, sum);
	cout << "Integeral approximation of h(x) = e^-((ln(x))^2/2)/x from 0 to 1 = " << sum << ", by Trapeziod Rule." << endl;
	cout << endl << endl;
    	//--------------------------------------
	cout << "Problem #1, Part c..." << endl;
	// Pretty explanatory when you look at output.
	epsilon = pow(10, -4);
	upper_bound = 1;
	lower_bound = 0;
  	N = 6;
  	Matrix r(N+1, Array(N+1, 0.0));
  	Romberg(cout, l, lower_bound, upper_bound, N, r);
  	cout << "Using Romberg Algorithm (with " << N*N << " runs) for function l(x) = x^(1/3): " << endl;

  	cout << endl << endl;
    	//--------------------------------------
	cout << "Problem #2, Part a..." << endl;
	// Similar to the other problems using the trapezoidal formula looking at the graph of this function
    	// m(x) = t^2/(t^4+1) and the graph of its second derivtive shows a definite maximum at approx 
    	// x = .5401828 which means xi = .5401828
   	upper_bound = 1; 
	lower_bound = 0; 
	epsilon = pow(10, -12);

	h = sqrt(12 / (fabs((upper_bound - lower_bound) * sec_deriv_m(.5401828))) * epsilon);

	n = (upper_bound - lower_bound) / h;
	n = ceil (n);
	N = (int) n;
	Trapezoid (m, lower_bound, upper_bound, h, N, sum);
	cout << "Integeral approximation of m(x) = t^2/(t^4+1) from 0 to 1 = " << sum << ", by Trapeziod Rule." << endl;
	cout << endl << endl;
    	//---------------------------------------
	cout << "Problem #2, Part b..." << endl;
	// Using adaptive simpson's rule it is shown that the integral is similar
	// in approximation to the #2 part a
  	epsilon = 0.0000005; 
  	int level = 0;
  	int maxlevel = 20;
  	item result;
  
  	result = Simpson(cout, m, lower_bound, upper_bound, epsilon, level, maxlevel);
  	cout.precision(10);
  	cout << "Using Adaptive Simpson's (with " << maxlevel << " runs)" << endl; 
  	cout << "Approximate integral = " << setw(15) << result << endl;	
  	cout << endl << endl;
    	//-------------------------------
	cout << "Problem #3..." << endl;
	//  After a proper substitution of x = -ln(t) the new function 
	// p(x) = (1/x)*ln(x)/(1+x) then using adaptive simpson's rule it is easy
	// to approx the integral.  Also p(x) is undefined at x = 0 which means which
	// must consider the limit like in #1.b so we consider a lower bound 
	// very close to zero
	upper_bound = 1;
	lower_bound = .1;
  	epsilon = 0.0000005; 
  	level = 0;
  	maxlevel = 20;
  
  	result = Simpson(cout, p, lower_bound, upper_bound, epsilon, level, maxlevel);
  	cout.precision(10);
  	cout << "Using Adaptive Simpson's (with " << maxlevel << " runs)" << endl; 
  	cout << "Approximate integral = " << setw(15) << result << endl;
	cout << endl << endl;
	//--------------------------------
	system("PAUSE");
	return EXIT_SUCCESS;
}
//------------------
void Swap (item& a, item& b)
{
	item temp = a;
	a = b;
	b = temp;
}
void Trapezoid (item f(item), item a, item b, item h, int N, item& sum)
{
	item x;
	sum = 0.5 * (f(a) + f(b));
	for (int i = 1; i < N; i++)
	{
		x = a + i * h;
		sum = sum + f(x);
	}
	sum = sum * h;
}
void Romberg (ostream& outs, item f(item), item a, item b, int N, std::vector< std::vector<item> >& r)
{
	int i, j, k;
	item h, sum;
	outs.precision(10);
	h = b - a;
	r[0][0] = 0.5 * h * (f(a) + f(b));
	outs << "r[0][0] = 0.5 = " << setw(15) << r[0][0] << endl;
	for (i = 1; i <= N; i++)
	{
		h = h * 0.5;
		sum = 0.0;
		for (k = 1; k <= pow(2, i) - 1; k = k + 2)
		{
			sum = sum + f(a + k * h);
		}
		r[i][0] = 0.5 * r[i-1][0] + sum * h;
		outs << "r[" << i << "][0] = " << setw(15) << r[i][0] << endl;
		for (j = 1; j <= i; j++)
		{
			r[i][j] = r[i][j-1] + (r[i][j-1] - r[i-1][j-1]) / (pow(4,j)-1);
			outs << "r["<< i << "][" << j << "] = "<< setw(15) << r[i][j] << endl;
		}
	}
}
item Simpson (ostream& out, item f(item), item a, item b, item epsilon, int level, int level_max)
{
  	item c, d, e, h;
  	item one_simpson, 
      	two_simpson, 
      	simpson_result,
      	left_simpson,
      	right_simpson;
 
  	level++;
  	h = b - a;
  	c = 0.5 * (a + b);

  	one_simpson = h * (f(a) + 4 * f(c) + f(b)) / 6.0;
  	d = 0.5 * (a + c);
  	e = 0.5 * (c + b);
  	two_simpson = h * (f(a) + 4 * f(d) + 2 * f(c) + 4 * f(e) + f(b)) / 12.0;
  	if (level >= level_max)
  	{
  	  	simpson_result = two_simpson;
		out << "maximum level reached" << endl;
  	}
  	else
  	{
    		if ((fabs(two_simpson - one_simpson))  < 15.0 * epsilon)
      	{	
			simpson_result = two_simpson + (two_simpson - one_simpson) / 15.0;
		}
    		else
      	{
			left_simpson = Simpson(out, f, a, c, epsilon/2, level, level_max);
			right_simpson = Simpson(out, f, c, b, epsilon/2, level, level_max);
			simpson_result = left_simpson + right_simpson;
      	}
   	}
	return simpson_result;
}
