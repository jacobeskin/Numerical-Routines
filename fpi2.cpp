// g++ -Wall fpi2.cpp -o b
/*
Fixed point iteration solution to exp(-x)=x
 */

#include <iostream>
#include <math.h> 
#include <stdlib.h>

using namespace std;

int main()
{
  
  double x, y, e ,tol=1.0*pow(10,-10);
  int count = 0;
  
  // Enter staring point for iteration
  cout<<"Enter starting point for iteration:"<<endl;
  cin>>x;
  cout<<" "<<endl;
  
  // Begin the iteration
  do
    {
      y = exp(-x);
      e = fabs(x-y);
      x = y;
      count++;
      if (count>100)
	{
	  cout<<"Max iterations reached!"<<endl;
	  exit(1);
	}
    }
  while(e>tol);
  
  // Printing the root
  cout<<"Root:"<<" "<<y<<endl;
  cout<<"Number of iterations required:"<<" "<<count<<endl;
  cout<<"|f(x)-x|="<<e<<endl;
}
  
