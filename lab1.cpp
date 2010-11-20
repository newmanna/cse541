//.cpp file

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <fstream>
#include "lab1_Root_Finding_Methods.h"

using namespace std;

typedef double item;

/*
       Main Body is ultimately a test driver given function def/decl and operation def/decl and
       is formatted to produce output to file below in the form of a solution set for Lab 1
*/
int main()
{
    //comments and output are self explanatory and should be easy to read
    ofstream outfile("cse 541 Lab #1 Newman.txt");
    item upper, lower, delta, epsilon;
    int MAX = 30;
    outfile << "Author: Nathan Newman\nDate: " << endl << "Lab 1 Problem 1:" << endl << "------------------------" << endl
            << "Part a) " << endl;
    //Looking at the graph of the function x^3 - 2x + 1 - x^2 shows roots between 
    //(-1.5, -.5), (.25, .75), (1.5, 2)
            
    upper = -.5;
    lower = -1.5;
    delta = pow(10, -4);
    epsilon = pow(10, -5);
    item answer1 = Bisection(outfile, f, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of f(x) = x^3 - 2x + 1 - x^2 at x = " << answer1 << endl;
    
    upper = .75;
    lower = .25;
    delta = pow(10, -4);
    epsilon = pow(10, -5);
    item answer2 = Bisection(outfile, f, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of f(x) = x^3 - 2x + 1 - x^2 at x = " << answer2 << endl;
    
    upper = 2;
    lower = 1.5;
    delta = pow(10, -4);
    epsilon = pow(10, -5);
    item answer3 = Bisection(outfile, f, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of f(x) = x^3 - 2x + 1 - x^2 at x = " << answer3 << endl << endl;
    
    outfile << "These three roots imply that there exists three points of intersection\n\t"
         << "for the functions x^3-2x+1 and x^2 and they are the roots shown \n\tabove." << endl;
    //---------------------------------------------     
    outfile << "Part b) " << endl;
    
    int months_in_a_year = 12;
    item interest = 0.0675, interest_per_month = interest / months_in_a_year;
    
    
    //---------------------------------------------
    outfile << "Part c) " << endl;
    //Looking at the graph of the given function shows there exists roots between 
    //(-2, -1), (2, 3), (10, 11)
    upper = -1;
    lower = -2;
    delta = pow(10, -4);
    epsilon = pow(10, -5);
    answer1 = Bisection(outfile, h, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of h(x) = 2 - x^2*e^(-.385*x) at x = " << answer1 << endl;
    
    upper = 3;
    lower = 2;
    delta = pow(10, -4);
    epsilon = pow(10, -5);
    answer2 = Bisection(outfile, h, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of h(x) = 2 - x^2*e^(-.385*x) at x = " << answer2 << endl;
    
    upper = 11;
    lower = 10;
    delta = pow(10, -5);
    epsilon = pow(10, -6);
    answer3 = Bisection(outfile, h, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of h(x) = 2 - x^2*e^(-.385*x) at x = " << answer3 << endl;
    //-------------------------------------------------
    outfile << "Part d) " << endl;
    //looking at the graph of the function sin x - .6 x cos x - 1 implies that the first two positive roots lie between
    // (1.5, 2), (3.5, 4)
    upper = 2;
    lower = 1.5;
    delta = pow(10, -5);
    epsilon = pow(10, -6);
    answer1 = Bisection(outfile, k, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of k(x) = sin(x) - .6*x*cos(x) - 1 at x = " << answer1 << endl;
    
    upper = 4;
    lower = 3.5;
    delta = pow(10, -5);
    epsilon = pow(10, -6);
    answer2 = Bisection(outfile, k, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of k(x) = sin(x) - .6*x*cos(x) - 1 at x = " << answer2 << endl;
    //----------------------------------------------------
    
    
    outfile << "\nLab 1 Problem 2" << endl <<"-------------------------" << endl;
    
    //for problem 5 from homework 3
    outfile << "-------------------------" << endl;
    outfile << "PROBLEM 5 FROM HOMEWORKK 3" << endl;
    outfile << "-------------------------" << endl;
    item initial_x_value = .5;
    delta = .00001;
    epsilon = .000001;
    answer1 = Newton(outfile, z, deriv_z, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of z(x) at x = " << answer1 << endl;
    outfile << "-------------------------" << endl;
    outfile << "END OF PROBLEM 5 FROM HOMEWORKK 3" << endl;
    outfile << "-------------------------" << endl;    
    outfile << "Part a)" << endl;      
    //looking at the graph of the function tan(x) - x and given the boundary [4, 7.5], it appears there might be 
    //roots at x = 4.5
    initial_x_value = 4.5;
    delta = .0001;
    epsilon = .00001;
    answer1 = Newton(outfile, l, deriv_l, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of l(x) = tan(x) - x at x = " << answer1 << endl;
    //----------------------------------------------------    
    outfile << "Part b)" << endl;
    //looking at the graph of the function (x-1)^2 - 0.5*e^-x and given the boundary x > 0, it appears there might be 
    //roots at x = .5 and x = 1.25
    initial_x_value = .5;
    delta = .0001;
    epsilon = .00001;
    answer1 = Newton(outfile, m, deriv_m, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of m(x) = (x-1)^2 - 0.5*e^-x at x = " << answer1 << endl;
    //----------------------------------------------------   
    initial_x_value = 4.5;
    delta = .0001;
    epsilon = .00001;
    answer2 = Newton(outfile, m, deriv_m, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of m(x) = (x-1)^2 - 0.5*e^-x at x = " << answer2 << endl;
    //----------------------------------------------------  
    outfile << "Part c)" << endl;    
    //looking at the graph of the function x^3-x^2-.00501*x-.0000077 it appears there might be roots at x~0 
    //and x=1
    initial_x_value = 0;
    delta = .0001;
    epsilon = .00001;
    answer1 = Newton(outfile, n, deriv_n, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of n(x) = x^3-x^2-.00501*x-.0000077 at x = " << answer1 << endl;
    //----------------------------------------------------    
    initial_x_value = 1;
    delta = .5*pow(10, -5);
    epsilon = pow(10, -6);
    answer2 = Newton(outfile, n, deriv_n, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of n(x) = x^3-x^2-.00501*x-.0000077 at x = " << answer2 << endl;
    //----------------------------------------------------
    outfile << "\nLab 1 Problem 3" << endl <<"-------------------------" << endl;
    outfile << "\t\t\tUsing Secant Method:" << endl << endl;
    //looking at the graph of the function 5x - 8ln(x) - 8 it appears there might be 
    //roots at x = .5 and x = 3.5 
    upper = .25;
    lower = .75;
    delta = pow(10, -10);
    epsilon = pow(10, -17);
    answer1 = Secant (outfile, p, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of p(x) = 5x - 8ln(x) - 8 at x = " << answer1 << endl;
    //----------------------------------------------------  
    upper = 3.25;
    lower = 3.75;
    delta = pow(10, -10);
    epsilon = pow(10, -17);
    answer2 = Secant (outfile, p, lower, upper, MAX, delta, epsilon);
    outfile << "There exists a root of p(x) = 5x - 8ln(x) - 8 at x = " << answer2 << endl;
    //----------------------------------------------------  
    outfile << "\n\t\t\tUsing Newton Method:" << endl << endl;
        initial_x_value = .25;
    delta = pow(10, -10);
    epsilon = pow(10, -17);
    answer1 = Newton(outfile, p, deriv_p, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of p(x) = 5x - 8ln(x) - 8 at x = " << answer1 << endl;
    //----------------------------------------------------    
    initial_x_value = 3.25;
    delta = pow(10, -10);
    epsilon = pow(10, -17);
    answer2 = Newton(outfile, p, deriv_p, initial_x_value, MAX, delta, epsilon);
    outfile << "There exists a root of p(x) = 5x - 8ln(x) - 8 at x = " << answer2 << endl;
    
    outfile << "\nTherefore, its easy to see that for this particular function and initial bounds/point for the\n"
            << "Secant/Newton methods (respectively) that the amount of steps for the Newton Method are roughly\n"
            << "2/3's the number of steps for the Secant Method.  Since the number of steps for the secant method\n"
            << "is greater than that for newton's method then the convergence is not quadratic." << endl;
     //----------------------------------------------------
    outfile << "\nLab 1 Problem 4" << endl <<"-------------------------" << endl;
    outfile << "Part a)" << endl;
    //looking at the graph of the function sin(x+ (1.2 - x^2)^(1/2)) - 1.4x + .2 implies
    //that the function has one root around x = .75         
    initial_x_value = .75;
    delta = .5 * pow(10, -4);;
    epsilon = pow(10, -4);;
    answer1 = Newton(outfile, q, deriv_q, initial_x_value, MAX, delta, epsilon);
    //Now you need to plug answer1 into function y in order to find the y-value of the 
    //solution to the sytem of the equations stated in problem
    answer2 = y(answer1);
    outfile << "There exists a solution to the set of equations: " 
            << "\n\tsin(x+y)-1.4x+.2=0\n\tx^2+y^2-1.2=0 "
            << "\nat the points x = " << answer1 << ", and y = " << answer2 << endl;
        
    //system("PAUSE");
    return EXIT_SUCCESS;
}
