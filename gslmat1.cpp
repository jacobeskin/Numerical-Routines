// g++ -Wall gslmat1.cpp -o gslmat1.exe -lgsl -lgslcblas
/*
Least squares solution for Ax=b, A is mxn matrix where m>n using  
two solution methods, GSL QR decomposition and direct GSL LSQ function.
y(t)=x1t+x2 with data {(2,1),(5,2),(7,3),(8,3)}.
*/

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>

using namespace std;

int main()
{

  
  // First solve with QR decomposition
  
  gsl_matrix *A;
  gsl_vector *b, *x, *tau, *residual;
  int i;
  double sum;

  // Allocate space for the matricies
  A = gsl_matrix_alloc(4,2);
  b = gsl_vector_alloc(4);
  x = gsl_vector_alloc(2);
  residual = gsl_vector_alloc(4);
  tau = gsl_vector_alloc(2);

  // Put the data into the matricies
  gsl_matrix_set(A, 0, 0, 2.0);
  gsl_matrix_set(A, 0, 1, 1.0);
  gsl_matrix_set(A, 1, 0, 5.0);
  gsl_matrix_set(A, 1, 1, 1.0);
  gsl_matrix_set(A, 2, 0, 7.0);
  gsl_matrix_set(A, 2, 1, 1.0);
  gsl_matrix_set(A, 3, 0, 8.0);
  gsl_matrix_set(A, 3, 1, 1.0);

  gsl_vector_set(b, 0, 1.0);
  gsl_vector_set(b, 1, 2.0);
  gsl_vector_set(b, 2, 3.0);
  gsl_vector_set(b, 3, 3.0);

  // Solution and results
  gsl_linalg_QR_decomp(A, tau);
  gsl_linalg_QR_lssolve(A, tau, b, x, residual);
  cout<<" "<<endl;
  cout<<"x1 and x2 from QR decomposition:"<<endl;
  cout<<gsl_vector_get(x,0)<<"  "<<gsl_vector_get(x,1)<<endl;
  cout<<" "<<endl;
  cout<<"And the sum of the squared residuals is:"<<endl;
  sum = 0.0;
  for(i=0;i<4;i++)
    {
      sum = gsl_vector_get(residual,i)* gsl_vector_get(residual,i)+sum;
    }
  cout<<sum<<endl;
  cout<<" "<<endl;

  // Free the space
  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(tau);
  gsl_vector_free(residual); 

  
  // Solution with gsl_fit_linear(), without gsl_matricies for funzies

  
  int n = 4;
  double t[4] = {2.0, 5.0, 7.0, 8.0};
  double y[4] = {1.0, 2.0, 3.0, 3.0};
  double x1, x2, cov00, cov01, cov11, sumsq;

  // Solution and results
  gsl_fit_linear(t, 1, y, 1, n, &x1, &x2, &cov00, &cov01, &cov11, &sumsq);
  cout<<"x1 and x2 from linear fit funcion:"<<endl;
  cout<<x2<<"  "<<x1<<endl;
  cout<<" "<<endl;
  cout<<"And the sum of the residuals squared is:"<<endl;
  cout<<sumsq<<endl;
  cout<<" "<<endl;

}


  

  

  

  
