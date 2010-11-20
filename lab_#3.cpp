// Author: Nathan Newman
// Date:
// CSE 541, Lab #3, Malkiman, MW 5:30

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;
//-------------------------------------------------------------------------
typedef double item;
typedef std::vector<item> Array;
typedef std::vector<Array> Matrix;

typedef std::vector<int> IntArray;
//-------------------------------------------------------------------------
// function declarations
//Helper operations
void outputArr (ofstream& outs, int n, std::vector<item> a);
void outputIntArr(ofstream& outs, int n, std::vector<int> l);
void outputMatrix (ofstream& outs, int n, std::vector< std::vector<item> > a);
void ClearArr (int n, std::vector<item>& a);
void ClearMatrix (int n, std::vector< std::vector<item> >& a);
//Matix Operations
void Forw_Subs (int n, std::vector< std::vector<item> > l, std::vector<item>& z, std::vector<item> b);
void Back_Subs (int n, std::vector< std::vector<item> > u, std::vector<item>& x, std::vector<item> z);
void Inverse3x3 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A);
void Inverse4x4 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A);
void Inverse5x5 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A);
item Determinant (int n, std::vector< std::vector<item> > a);
//Lower-Upper Decomposition Operation
void LUDecomp(int n, std::vector< std::vector<item> > a, std::vector< std::vector<item> >& l, std::vector< std::vector<item> >& u);
//Tridiagonal
void Tridiagonal(int n, std::vector<item> a, std::vector<item> d, std::vector<item> c, std::vector<item> b, std::vector<item>& x);
//Gauss/Scaled Pivoting
void Gauss (int n, std::vector< std::vector<item> >& a, std::vector<int>& l);
void Solve (int n, std::vector< std::vector<item> > a, std::vector<int> l, std::vector<item>& b, std::vector<item>& x);
//Gauss-Seidel Method
void GaussSeidel(ofstream& outs, int n, std::vector< std::vector<item> > a, std::vector<item> b, std::vector<item>& x);
bool epsilon_comp(int size, std::vector<item> x, std::vector<item> y);
//-------------------------------------------------------------------------

int main ()
{
    ofstream outfile("CSE_541_Lab_#3_Newman.txt");
    outfile << " Author: Nathan Newman \n Date:8/20/2010 \n Lab:#3" << endl << endl;
    //--------------------------------------------------------------------------
    // Problem 1, Part a...
    outfile << "Problem #1, part a..." << endl;
    //It is clear to see that we have 3x3 matrices so therefore our n = 3.
    int n = 3;
    //intialize Matrix A
    Matrix Aa(n+1, Array(n+1, 0.0));
    Aa[1][1] = 3.1; Aa[1][2] = -0.6; Aa[1][3] = 0.3;
    Aa[2][1] = -0.3; Aa[2][2] = 6.1; Aa[2][3] = -0.3;
    Aa[3][1] = -0.7; Aa[3][2] = 0.6; Aa[3][3] = -2.9;
    //initialize Matrix B
    Array Ba(n+1, 0.0);  
    Ba[1] = .4; Ba[2] = .8; Ba[3] = .6;
    Matrix La(n+1, Array(n+1, 0.0));
    Matrix Ua(n+1, Array(n+1, 0.0));
    //LU decmposition on Matrix A
    LUDecomp(n, Aa, La, Ua);
    //output matrix A
    outfile << "\nMatrix A:" << endl;
    outputMatrix (outfile, n, Aa);
    
    //output Matrix L
    outfile << "\nMatrix L:" << endl;
    outputMatrix (outfile, n, La);
    
    //output Matrix U
    outfile << "\nMatrix U:" << endl;
    outputMatrix (outfile, n, Ua);
    
    //Forward Substitution implies...Lz=b
    Array za(n+1, 0.0);
    Forw_Subs(n, La, za, Ba);
    
    //output Matrix z
    outfile << "\nMatrix z:" << endl;
    outputArr (outfile, n, za);
    
    //Backward Substitution implies...Ux = z
    Array xa(n+1, 0.0);
    Back_Subs (n, Ua, xa, za);
    
    //output Matrix x
    outfile << "\nMatrix x: (Corresponding solutions: x1, x2, x3 to the given system of equations)" << endl;
    outputArr (outfile, n, xa);   
    
    //Find Determinant of Matrix A
    item DetA = Determinant(n, La);
    outfile << "\nDeterminant of A = " << DetA << endl;
      
    //Find Inverse of A
    Matrix inv_Aa(n+1, Array(n+1, 0.0));
    Inverse3x3 (n, La, Ua, inv_Aa);
    outfile << "\nMatrix Inverse A:" << endl;
    outputMatrix (outfile, n, inv_Aa);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;
    
    // Problem 1, Part b...
    outfile << "Problem #1, part b..." << endl;
    //It is clear to see that we have 3x3 matrices so therefore our n = 3.
    n = 4;
    //intialize Matrix A
    Matrix Ab(n+1, Array(n+1, 0.0));
    Ab[1][1] = 1.0671; Ab[1][2] = .0216; Ab[1][3] = -.02995; Ab[1][4] = .0101;
    Ab[2][1] = .0101; Ab[2][2] = .5671; Ab[2][3] = .0216; Ab[2][4] = -.02995;
    Ab[3][1] = -.02995; Ab[3][2] = .0101; Ab[3][3] = 1.1671; Ab[3][4] = .0216;
    Ab[4][1] = 1.0216; Ab[4][2] = -.02995; Ab[4][3] = .0101; Ab[4][4] = .6671;
    //initialize Matrix B
    Array Bb(n+1, 0.0);  
    Bb[1] = 1.941; Bb[2] = -.430; Bb[3] = 1.941; Bb[4] = -.430;
    Matrix Lb(n+1, Array(n+1, 0.0));
    Matrix Ub(n+1, Array(n+1, 0.0));
    //LU decmposition on Matrix A
    LUDecomp(n, Ab, Lb, Ub);
    //output matrix A
    outfile << "\nMatrix A:" << endl;
    outputMatrix (outfile, n, Ab);
    
    //output Matrix L
    outfile << "\nMatrix L:" << endl;
    outputMatrix (outfile, n, Lb);
    
    //output Matrix U
    outfile << "\nMatrix U:" << endl;
    outputMatrix (outfile, n, Ub);
    
    //Forward Substitution implies...Lz=b
    Array zb(n+1, 0.0);
    Forw_Subs(n, Lb, zb, Bb);
    
    //output Matrix z
    outfile << "\nMatrix z:" << endl;
    outputArr (outfile, n, zb);
    
    //Backward Substitution implies...Ux = z
    Array xb(n+1, 0.0);
    Back_Subs (n, Ub, xb, zb);
    
    //output Matrix x
    outfile << "\nMatrix x: (Corresponding solutions: x1, x2, x3, x4 to the given system of equations)" << endl;
    outputArr (outfile, n, xb);   
    
    //Find Determinant of Matrix A
    DetA = Determinant(n, Lb);
    outfile << "\nDeterminant of A = " << DetA << endl;
      
    //Find Inverse of A
    Matrix inv_Ab(n+1, Array(n+1, 0.0));
    Inverse4x4 (n, Lb, Ub, inv_Ab);
    outfile << "\nMatrix Inverse A:" << endl;
    outputMatrix (outfile, n, inv_Ab);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;  
    
    // Problem 1, Part c...
    outfile << "Problem #1, part c..." << endl;
    //It is clear to see that we have 3x3 matrices so therefore our n = 3.
    n = 5;
    //intialize Matrix A
    Matrix Ac(n+1, Array(n+1, 0.0));
    Ac[1][1] = 10.25; Ac[1][2] = -16.94; Ac[1][3] = 11.06; Ac[1][4] = -16.52; Ac[1][5] = 10.01;
    Ac[2][1] = -17.43; Ac[2][2] = 11.51; Ac[2][3] = -20.58; Ac[2][4] = 13.93; Ac[2][5] = 13.16;
    Ac[3][1] = 0; Ac[3][2] = 12.53; Ac[3][3] = -18.03; Ac[3][4] = 8.82; Ac[3][5] = -17.85;
    Ac[4][1] = 0; Ac[4][2] = 0; Ac[4][3] = 8.54; Ac[4][4] = -14.11; Ac[4][5] = 9.17;
    Ac[5][1] = 0; Ac[5][2] = 0; Ac[5][3] = 0; Ac[5][4] = 7.77; Ac[5][5] = -14.53;
    //initialize Matrix B
    Array Bc(n+1, 0.0);  
    Bc[1] = -1.20; Bc[2] = -2.33; Bc[3] = -.14; Bc[4] = -.09; Bc[5] = 2.19;
    Matrix Lc(n+1, Array(n+1, 0.0));
    Matrix Uc(n+1, Array(n+1, 0.0));
    //LU decmposition on Matrix A
    LUDecomp(n, Ac, Lc, Uc);
    //output matrix A
    outfile << "\nMatrix A:" << endl;
    outputMatrix (outfile, n, Ac);
    
    //output Matrix L
    outfile << "\nMatrix L:" << endl;
    outputMatrix (outfile, n, Lc);
    
    //output Matrix U
    outfile << "\nMatrix U:" << endl;
    outputMatrix (outfile, n, Uc);
    
    //Forward Substitution implies...Lz=b
    Array zc(n+1, 0.0);
    Forw_Subs(n, Lc, zc, Bc);    
    //output Matrix z
    outfile << "\nMatrix z:" << endl;
    outputArr (outfile, n, zc);
    
    //Backward Substitution implies...Ux = z
    Array xc(n+1, 0.0);
    Back_Subs (n, Uc, xc, zc);
    
    //output Matrix x
    outfile << "\nMatrix x: (Corresponding solutions: x1, x2, x3, x4, x5 to the given system of equations)" << endl;
    outputArr (outfile, n, xc);   
    
    //Find Determinant of Matrix A
    DetA = Determinant(n, Lc);
    outfile << "\nDeterminant of A = " << DetA << endl;
      
    //Find Inverse of A
    Matrix inv_Ac(n+1, Array(n+1, 0.0));
    Inverse5x5 (n, Lc, Uc, inv_Ac);
    outfile << "\nMatrix Inverse A:" << endl;
    outputMatrix (outfile, n, inv_Ac);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;   
    
    // Problem 2...
    outfile << "Problem #2..." << endl << endl;
    n = 3;
    //Clear matrices needed for calculations
    ClearMatrix(n, Aa);
    ClearArr(n, Ba);
    ClearArr(n, xa);
    
    //Initialize A
    Aa[1][1] = 15.035; Aa[1][2] = 5.45; Aa[1][3] = 7.46;
    Aa[2][1] = -3.3775; Aa[2][2] = 23.06; Aa[2][3] = 11.06;
    Aa[3][1] = 7.2225; Aa[3][2] = -1.1475; Aa[3][3] = 10.2475;
    
    //Initialize B
    Ba[1] = 50.91; Ba[2] = 50.17; Ba[3] = 42.68;

    //Use Gauss-Seidel to solve system
    GaussSeidel(outfile, n, Aa, Ba, xa);
    
    //Output Matrix x, which is the solution set for the system
    outfile <<"\nMatrix x:" << endl;
    outputArr(outfile, n, xa);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;  
    
    //Problem 3, part a...
    outfile << "Problem #3, part a..." << endl;   
    n = 3;
    //Clear Matrices
    ClearMatrix(n, Aa);
    ClearArr(n, Ba);
    ClearArr(n, xa);
    
    Aa[1][1] = 3.3; Aa[1][2] = 3.5; Aa[1][3] = -3;
    Aa[2][1] = 4; Aa[2][2] = 1.7; Aa[2][3] = 5.3;
    Aa[3][1] = -5; Aa[3][2] = 2.5; Aa[3][3] = 3.7;
    Ba[1] = 21.07; Ba[2] = 1.21; Ba[3] = -16.25;
    IntArray l1(n+1, 0);
    
    //output matrix A and Matrix B
    outfile << "\nMatrix A:" << endl;
    outputMatrix(outfile, n, Aa);
    outfile << "\nMatrix B:" << endl;
    outputArr(outfile, n, Ba);
    
    Gauss(n, Aa, l1);
    
    //output matrix A and Matrix l
    outfile << "\nMatrix A:" << endl;
    outputMatrix(outfile, n, Aa);
    outfile << "\nMatrix l:" << endl;
    outputIntArr(outfile, n, l1);
    
    Solve (n, Aa, l1, Ba, xa);
    
    //output matrix x
    outfile << "\nMatrix x:" << endl;
    outputArr(outfile, n, xa);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;  
    
    //Problem 3, part b..
    outfile << "Problem #3, part b..." << endl;   
    n = 4;
    //Clear Matrices
    ClearMatrix(n, Ab);
    ClearArr(n, Bb);
    ClearArr(n, xb);
    
    Ab[1][1] = 5.1; Ab[1][2] = 2.7; Ab[1][3] = 4.8; Ab[1][4] = 7.5;
    Ab[2][1] = -9.4; Ab[2][2] = 3.5; Ab[2][3] = 12.4; Ab[2][4] = 14.3;
    Ab[3][1] = 0; Ab[3][2] = 6; Ab[3][3] = 8.7; Ab[3][4] = 2.6;
    Ab[4][1] = 0; Ab[4][2] = 0; Ab[4][3] = -10; Ab[4][4] = -12;
    Bb[1] = 4; Bb[2] = -13; Bb[3] = 1; Bb[4] = -5;
    IntArray l2(n+1, 0);
    
    //output matrix A and Matrix B
    outfile << "\nMatrix A:" << endl;
    outputMatrix(outfile, n, Ab);
    outfile << "\nMatrix B:" << endl;
    outputArr(outfile, n, Bb);
    
    Gauss(n, Ab, l2);
    
    //output matrix A and Matrix l
    outfile << "\nMatrix A:" << endl;
    outputMatrix(outfile, n, Ab);
    outfile << "\nMatrix l:" << endl;
    outputIntArr(outfile, n, l2);
    
    Solve (n, Ab, l2, Bb, xb);
    
    //output matrix x
    outfile << "\nMatrix x:" << endl;
    outputArr(outfile, n, xb);
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    outfile << "*--------------------------------------------------------------*" << endl;
    outfile << endl << endl;  
    
    //Problem 4..
    outfile << "Problem #4..." << endl; 
    
    
	//system("PAUSE");
	return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------

void outputArr (ofstream& outs, int n, std::vector<item> a)
{
     for (int i = 1; i <= n; i++)
     {
        outs << "|" << a[i] << "|" << endl;
     }
}
void outputIntArr(ofstream& outs, int n, std::vector<int> l)
{
     for (int i = 1; i <= n; i++)
     {
        outs << "|" << l[i] << "|" << endl;
     }
}
void outputMatrix (ofstream& outs, int n, std::vector< std::vector<item> > a)
{
    int i, j;
    for (i = 1; i <= n; i++)
    {
        
        outs << "|";
        for (j = 1; j <= n; j++)
        {
              outs << a[i][j] << " ";
        }               
        outs << "|\n";
    }     
}
void ClearArr (int n, std::vector<item>& a)
{
     for (int i = 1; i <= n; i++)
     {
         a[i] = 0.0;
     }     
}
void ClearMatrix (int n, std::vector< std::vector<item> >& a)
{
     int i, j;
     for (i = 1; i <= n; i++)
     {
         for (j = 1; j <= n; j++)
         {
             a[i][j] = 0;
         }
     }
}
void Forw_Subs (int n, std::vector< std::vector<item> > l, std::vector<item>& z, std::vector<item> b)
{
     z[1] = b[1]/l[1][1];
     int i, s;
     for (i = 2; i <= n; i++)
     {
         item temp_sum = 0;
         for (s = 1; s <= i-1; s++)
         {
             temp_sum = temp_sum + l[i][s]*z[s];
         }
         z[i] = (b[i] - temp_sum)/l[i][i];
     }     
}
void Back_Subs (int n, std::vector< std::vector<item> > u, std::vector<item>& x, std::vector<item> z)
{
     x[n] = z[n];
     int i, s;
     for (i = n-1; i >= 1; i--)
     {
         item temp_sum = 0;
         for (s = i+1; s <= n; s++)
         {
             temp_sum = temp_sum + u[i][s]*x[s];
         }
         x[i] = z[i] - temp_sum;
     }
}
void Inverse3x3 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A)
{
     Array b1(n+1, 0.0), b2(n+1, 0.0), b3(n+1, 0.0), 
           z1(n+1, 0.0), z2(n+1, 0.0), z3(n+1, 0.0), 
           x1(n+1, 0.0), x2(n+1, 0.0), x3(n+1, 0.0);  
     b1[1] = 1; b1[2] = 0; b1[3] = 0;
     b2[1] = 0; b2[2] = 1; b2[3] = 0;
     b3[1] = 0; b3[2] = 0; b3[3] = 1;
     //solve for x1
     Forw_Subs(n, l, z1, b1);
     Back_Subs(n, u, x1, z1);
     //solve for x2
     Forw_Subs(n, l, z2, b2);
     Back_Subs(n, u, x2, z2); 
     //solve for x3
     Forw_Subs(n, l, z3, b3);
     Back_Subs(n, u, x3, z3); 
     //produce inv_A
     for (int i = 1; i <= n; i++)
     {
         for (int j = 1; j <= n; j++)
         {
             if (i == 1)
             {
                   inv_A[j][i] = x1[j];
             }
             else if (i == 2)
             {
                  inv_A[j][i] = x2[j];
             }
             else if (i == 3)
             {
                  inv_A[j][i] = x3[j];
             }
         }
     }
}
void Inverse4x4 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A)
{
     Array b1(n+1, 0.0), b2(n+1, 0.0), b3(n+1, 0.0), b4(n+1, 0.0), 
           z1(n+1, 0.0), z2(n+1, 0.0), z3(n+1, 0.0), z4(n+1, 0.0),
           x1(n+1, 0.0), x2(n+1, 0.0), x3(n+1, 0.0), x4(n+1, 0.0);  
     b1[1] = 1; b1[2] = 0; b1[3] = 0, b1[4] = 0;
     b2[1] = 0; b2[2] = 1; b2[3] = 0, b2[4] = 0;
     b3[1] = 0; b3[2] = 0; b3[3] = 1, b3[4] = 0;
     b4[1] = 0; b4[2] = 0; b4[3] = 0, b4[4] = 1;
     //solve for x1
     Forw_Subs(n, l, z1, b1);
     Back_Subs(n, u, x1, z1);
     //solve for x2
     Forw_Subs(n, l, z2, b2);
     Back_Subs(n, u, x2, z2); 
     //solve for x3
     Forw_Subs(n, l, z3, b3);
     Back_Subs(n, u, x3, z3); 
     //solve for x4
     Forw_Subs(n, l, z4, b4);
     Back_Subs(n, u, x4, z4);
     //produce inv_A
     for (int i = 1; i <= n; i++)
     {
         for (int j = 1; j <= n; j++)
         {
             if (i == 1)
             {
                   inv_A[j][i] = x1[j];
             }
             else if (i == 2)
             {
                  inv_A[j][i] = x2[j];
             }
             else if (i == 3)
             {
                  inv_A[j][i] = x3[j];
             }
             else if (i == 4)
             {
                  inv_A[j][i] = x4[j];
             }
         }
     }       
}
void Inverse5x5 (int n, std::vector< std::vector<item> > l, std::vector< std::vector<item> > u, std::vector< std::vector<item> >& inv_A)
{
     Array b1(n+1, 0.0), b2(n+1, 0.0), b3(n+1, 0.0), b4(n+1, 0.0), b5(n+1, 0.0),
           z1(n+1, 0.0), z2(n+1, 0.0), z3(n+1, 0.0), z4(n+1, 0.0), z5(n+1, 0.0),
           x1(n+1, 0.0), x2(n+1, 0.0), x3(n+1, 0.0), x4(n+1, 0.0), x5(n+1, 0.0);  
     b1[1] = 1; b1[2] = 0; b1[3] = 0, b1[4] = 0, b1[5] = 0;
     b2[1] = 0; b2[2] = 1; b2[3] = 0, b2[4] = 0, b2[5] = 0;
     b3[1] = 0; b3[2] = 0; b3[3] = 1, b3[4] = 0, b3[5] = 0;
     b4[1] = 0; b4[2] = 0; b4[3] = 0, b4[4] = 1, b4[5] = 0;
     b5[1] = 0; b5[2] = 0; b5[3] = 0, b5[4] = 0, b5[5] = 1;
     //solve for x1
     Forw_Subs(n, l, z1, b1);
     Back_Subs(n, u, x1, z1);
     //solve for x2
     Forw_Subs(n, l, z2, b2);
     Back_Subs(n, u, x2, z2); 
     //solve for x3
     Forw_Subs(n, l, z3, b3);
     Back_Subs(n, u, x3, z3); 
     //solve for x4
     Forw_Subs(n, l, z4, b4);
     Back_Subs(n, u, x4, z4);
     //solve for x5
     Forw_Subs(n, l, z5, b5);
     Back_Subs(n, u, x5, z5);
     //produce inv_A
     for (int i = 1; i <= n; i++)
     {
         for (int j = 1; j <= n; j++)
         {
             if (i == 1)
             {
                   inv_A[j][i] = x1[j];
             }
             else if (i == 2)
             {
                  inv_A[j][i] = x2[j];
             }
             else if (i == 3)
             {
                  inv_A[j][i] = x3[j];
             }
             else if (i == 4)
             {
                  inv_A[j][i] = x4[j];
             }
             else if (i == 5)
             {
                  inv_A[j][i] = x5[j];
             }
         }
     }
}
item Determinant (int n, std::vector< std::vector<item> > a)
{
    item det = 1;
    int i;
    for (i = 1; i <= n; i++)
    {
        det = det * a[i][i];
    }
    return det;      
}
void LUDecomp (int n, std::vector< std::vector<item> > a, std::vector< std::vector<item> >& l, std::vector< std::vector<item> >& u)
{
     int i, j, s;
     for (i = 1; i <= n; i++)
     {
         l[i][1] = a[i][1];
     }
     for (i = 1; i <= n-1; i++)
     {
   		u[i][i] = 1;
		for (j = i+1; j <= n; j++)
		{
            item temp_sum = 0;
            for (s = 1; s <= i-1; s++)
            {
                temp_sum = temp_sum + l[i][s]*u[s][j];
            }
            u[i][j] = (a[i][j] - temp_sum)/l[i][i];
		}
		for (j = i+1; j <= n; j++)
		{
            item temp_sum = 0;
            for (s = 1; s <= i; s++)
            {
                temp_sum = temp_sum + l[j][s]*u[s][i+1];
            }
			l[j][i+1] = a[j][i+1] - temp_sum;
		}
	} 
	u[n][n] = 1;
}
//---------------------------------------------------------------------
void Tridiagonal (int n, std::vector<item> a, std::vector<item> d, std::vector<item> c, std::vector<item> b, std::vector<item>& x)
{
 	int i;
 	item xmult;

  	for (i = 2; i <= n; i++)
  	{
    		xmult = a[i-1]/d[i-1];
    		d[i] = d[i] - xmult * c[i-1];
    		b[i] = b[i] - xmult * b[i-1];
  	}
  	x[n] = b[n]/d[n];
  	for (i = n - 1; i >= 1; i--)
   	{
		x[i] = (b[i]- c[i]* x[i+1]) / d[i];
	}
}
//--------------------------------------------------------------------------
// Gaussian-Elimination with scaled row pivoting
#define max(a,b) (a > b) ? a : b
void Gauss (int n, std::vector< std::vector<item> >& a, std::vector<int>& l)
{
	Array s(n+1, 0.0);
	int i, j, k, temp_k;
	item r, rmax, smax, xmult;
  
	for (i = 1; i <= n; i++)
	{
		l[i] = i;
		smax = 0.0;
		for (j = 1; j <= n; j++)
			smax = max(smax, fabs(a[i][j]));
		s[i] = smax;
	}

	for (k = 1; k <= n - 1; k++)
	{
		rmax = 0.0;
		for (i = k; i <= n; i++)
		{
			r = fabs(a[l[i]][k]/ s[l[i]]);
			if (r > rmax)
			{
				rmax = r;
				j = i;
			}
		}
		temp_k = l[j];
		l[j] = l[k];
		l[k] = temp_k;

		for (i = k + 1; i <= n; i++)
		{
			xmult = a[l[i]][k] / a[l[k]][k];
			a[l[i]][k] = xmult;
			for (j = k + 1; j <= n; j++)
				a[l[i]][j] -= xmult * a[l[k]][j];
		}
	}
}
void Solve (int n, std::vector< std::vector<item> > a, std::vector<int> l, std::vector<item>& b, std::vector<item>& x)
{
	int i, j, k;
	item sum;
  
	for (k = 1; k <= n - 1; k++)
    {
		for (i = k + 1; i <= n; i++)
		{
			b[l[i]] = b[l[i]] - a[l[i]][k] * b[l[k]];
        }
    }
	x[n] = b[l[n]] / a[l[n]][n];
	for (i = n - 1; i >= 1; i--)
	{
		sum = b[l[i]];
		for (j = i + 1; j <= n; j++)
        {
			sum = sum - a[l[i]][j] * x[j];
        }
		x[i] = sum / a[l[i]][i];
	}   
}
//--------------------------------------------------------------------------
// Gauss-Seidel
bool epsilon_comp(int size, std::vector<item> x, std::vector<item> y)
{
     const item epsilon = 0.00005;
     int i, count;
  
     count = 0;
     for (i = 1; i <= size; i++)
     {
         if (fabs(x[i]-y[i]) < epsilon)
         {
            count++;
         }
     }
     if (count == size)
     {
         return true;
     }
     else
     {
         return false;
     }
}
void GaussSeidel(ofstream& out, int n, std::vector< std::vector<item> > A, std::vector<item> b, std::vector<item>& x)
{
	const int k_max = 100;

  	Array y(n+1, 0.0);
  	item temp_sum;

  	for (int k = 1; k <= k_max; k++)
	{
		y = x;
      	for (int i = 1; i <= n; i++)
		{
	  		temp_sum = b[i];
	  		for (int j = 1; j <= n; j++)
   		    {
				if (j != i)
	      		{
					temp_sum = temp_sum - A[i][j]*y[j];
				}
			}
	  		x[i] = temp_sum/A[i][i];     		 
		}
		out << "Iteration " << k << ": ";
		for (int i = 1; i <= n; i++)
		{
            out << "\tx[" << i << "]"<< "=" << x[i] << "\t";
        }
        out << endl;
      	if (epsilon_comp(n, x, y))
      	{
	  		return;
		}
 	}
  	out << "maximum iterations reached" << endl;
}


