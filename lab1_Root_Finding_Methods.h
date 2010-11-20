//.h file

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

#ifndef _lab1_Root_Finding_Methods_
#define _lab1_Root_Finding_Methods_

typedef double item;

/*
        Declarations of defined operations Bisection, Newton, Swap, Secant needed
        in order to calculate the root estimations for the given functions/equations 
        
        Bisection:
                  given the upper and lower (upper > lower) bounds of some function f, will estimate the 
                  roots according to the bisection method 
        Newton:
               given some initial point of some function f, will estimate the roots according 
               to the Newton method
        Swap:
             given two variables a and b, swaps the two values
        Secant:
               Similar to Bisection method
               
        For all functions, the variable MAX is used in case a function does not converge.
*/
item Bisection (ofstream& outs, item f(item), item a, item b, int MAX, item delta, item epsilon);
item Newton (ofstream& outs, item f(item), item deriv_f(item), item initial_pt, int MAX, item delta, item epsilon);
void Swap (item& a, item& b);
item Secant (ofstream& outs, item f(item), item a, item b, int MAX, item delta, item epsilon);

#include "lab1_Root_Finding_Methods.cc"
#endif
