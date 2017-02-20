// g++ -Wall agF.cpp -o a -lgsl -lgslcblas
/*
Calculate arithmetic-geometric mean of two input numbers ag(x,y) and evaluate 
the identity F=1/ag(1,sqrt(1-r^2)) proven by Gauss in 1799, where F is 
the hypergeometric function 2_F_1(a,b;c,d). F used is from the GSL.
 */

#include <iostream>
#include <math.h> 
#include <stdlib.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;

// Function for calculaing arithmetic-geometric mean
double ag(double a, double b, int n)
{
  double an, bn;
  int i=0;
  do
    {
      an = (a+b)/2;
      bn = sqrt(a*b);

      a = an;
      b = bn;
      
      i++;
    }
  while(i<n);
  
  return an;
}

int main()
{
  double N, RHS, LHS;
  
  cout<<"Difference between the right hand side and left hand side of the";
  cout<<"identity for k=1..19:"<<endl;
  for(int k=1;k<20;k++)
    {
      // Calculating the right hand side of the identity
      N=ag(1, sqrt(1-pow(0.05*k,2)), 100000);
      RHS = 1/N;

      // Calculating the left hand side of the identity
      LHS = gsl_sf_hyperg_2F1(0.5, 0.5, 1, pow(0.05*k, 2));

      cout<<LHS-RHS<<endl;
    } 
    
}
