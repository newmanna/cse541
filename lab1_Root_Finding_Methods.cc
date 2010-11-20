//.cc file

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <fstream>
#include "lab1_Root_Finding_Methods.h"

using namespace std;

typedef double item;

/*
         Declarations and Definitions of functions used in Problems 1-4, can be interchanged and altered
         easily in order to calculate different roots for different functions
*/
// function 1 is x^3 - 2x + 1 - x^2, used in problem 1.a
inline item f(item x) { return (pow(x,3) - 2*x + 1 - pow(x, 2)); }

//inline item deriv_f(item x) { return (3*x*x - 2 - 2*x); }
// function 2 is x, used in problem 1.b
inline item g(item x) { return (x); }

//function 3 is 2 - x^2*e^(-.385*x), used in 1.c
inline item h(item x) { return (2 - pow(x, 2) * exp(-.385 * x)); }

//function 4 is sin(x) - 0.6*x*cos(x)-1, used in 1.d
inline item k(item x) { return (sin(x) - .6 * x * cos(x) -1); }

//function 5 is tan(x) - x and its derivative sec^2(x) - 1, used in 2.a
inline item l(item x) { return (tan(x) - x);}
inline item deriv_l(item x) { return (pow(cos(x), -2) - 1); }

//function 6 is (x-1)^2 - 0.5*e^-x and its derivative is 2(x-1) + .5*e^-x, used in 2.b
inline item m(item x) { return (pow(x-1, 2) - .5*exp(-1*x));}
inline item deriv_m(item x) { return (2*(x-1) + .5 *exp(-1*x)); }

//function 7 is x^3-x^2-.00501*x-.0000077, and its derivative is 3*x^2-2*x-.00501, used in 2.c
inline item n(item x) { return (pow(x, 3) - pow(x, 2) - .00501*x - .000077);}
inline item deriv_n(item x) { return (3*pow(x, 2) - 2*x - .00501); }

//function 8 is 5x - 8ln(x) - 8, used in 3
inline item p(item x) { return (5*x - 8 * log(x) - 8);}
inline item deriv_p(item x) { return (5 - 8 * pow(x, -1));}

//function 9 is sin(x+ (1.2 - x^2)^(1/2)) - 1.4x + .2, which comes from solving the equation
//x^2+y^2 - 1.2 = 0 for y and substituting it into the equation sin (x+y) - 1.4x +2 = 0, used in problem 4
inline item q(item x) { return (sin(x + pow(1.2 - pow(x, 2), .5)) - 1.4*x +.2);}
inline item deriv_q(item x) { return ((1 - x * pow(1.2 - pow(x, 2), -.5)) * cos(x + pow(1.2 - pow(x, 2), .5)) - 1.4);}
inline item y(item x) { return (pow(1.2 - pow(x, 2), .5));}

//function 10 is -(8/9)*x^2+(25/9)*x - 1 which is used from problem 5 of homework 3 and its derivative is
//-(16/9)*x+(25/9)
inline item z(item x) { return (-(8.0/9)*x*x + (25.0/9)*x - 1);}
inline item deriv_z(item x) { return (-(16.0/9)*x + (25.0/9));}
/*
        Operation definitions follwing the information found in the .h declarations
        
        Also, for the sake of brevity, all output statements within the operations are commented out but
        can be un commented in order to allow them to print to the outfile.
*/
item Bisection(ofstream& outs, item f(item), item a, item b, int MAX, item delta, item epsilon)
{
        item root, w;
        item fa = f(a), 
        fb = f(b), 
        e = b - a;
        
        //outs << "\tLower Bound (a): " << lower << endl << "\tUpper Bound (b): " << upper 
             //<< endl << "\tf(a) = " << fa << endl << "\tf(b) = " << fb << endl;
        if (fa*fb < 0)
        {
                  for (int i = 0; i < MAX; i++)
                  {
                      e = e / 2;
                      root = a + e;
                      w = f(root);
                      //outs << "\ti = " << i << ", (c) = " << root << ", f(c) = " 
                           //<< w << ", (b - a) = " << e << endl;
                      if (fabs(e) < delta || fabs(w) < epsilon)
                      {
                                  outs << "\tNumber of steps: " << i << endl;
                                  return root;
                      }
                      if (w * fa < 0)
                      {
                               b = root;
                               fb = w;
                      }
                      else 
                      {
                               a = root;
                               fa = w;
                      }
                  }
        }
        outs << "\tNumber of steps > MAX" << endl;
        return root;
}
item Newton (ofstream& outs, item f(item), item deriv_f(item), item initial_pt, int MAX, item delta, item epsilon)
{
     item v, root = initial_pt;
     v = f(root);
     //outs << "\t" << 0 << ", " << initial_pt << ", " << v << endl;
     if (fabs(v) >= epsilon)
     {
                 for (int i = 0; i < MAX; i++)
                 {
                     item x_1;
                     x_1 = root - (v / deriv_f(root));
                     v = f(x_1);
                     //outs << "\t" << i << ", " << root << ", " << v << endl;
                     if ( fabs(x_1 - root) < delta || fabs(v) < epsilon)
                     {
                          outs << "\tNumber of steps: " << i << endl;
                          return root;
                     }
                     root = x_1;
                 }
     } 
     outs << "\tNumber of steps > MAX" << endl;
     return root;
}
void Swap (item& a, item& b)
{
     item temp;
     temp = a;
     a = b;
     b = temp;
}
item Secant (ofstream& outs, item f(item), item a, item b, int MAX, item delta, item epsilon)
{
     item s, root = a;
     item fa = f(root), fb = f(b);
     //outs << '\t' << 0 << ", " << root << ", " << fa << endl;   
     //outs << '\t' << 1 << ", " << b << ", " << fb << endl;
     for (int i = 2; i < MAX; i++)
     {
         if (fabs(fa) > fabs(fb))
         {
                      Swap(root, b);
                      Swap(fa, fb);
         }
         s = (b - root) / (fb - fa);
         b = root;
         fb = fa;
         root = root - (fa * s);
         fa = f(root);
         //outs << '\t' << i << ", " << root << ", " << fa << endl;
         if (fabs(fa) < epsilon || fabs(b- root) < delta)
         {
                      outs << "\tNumber of steps: " << i << endl;
                      return root;
         }
     }
     outs << "\tNumber of steps > MAX" << endl;
     return root;
}
