// g++ -Wall mt19937-64.c MCint.cpp -o MCint.exe 
/* 
Calculate area of the curves y = sin(x) and y = cos(x) in the rectangle
0<x<2pi, -1<y<1 and then the area between these curve using Monte Carlo 
integration. Uses Mersenne Twister mt19937-64.c.
*/

//#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "mt64.h"

using namespace std;

int main()
{
  float x, y, area;
  int seed, i;

  // I/O

  cout<<"Give integer for RNG seed:"<<endl;
  cin>>seed;
  cout<<" "<<endl;

  init_genrand64(seed);
  
  // Calculate first y = sin(x)
  area = 0.0;
  for(i=0; i<=1000000; i++)
    {
      x = genrand64_real3()*2.0*M_PI;
      y = sin(x);
      area = area+y;
    }

  cout<<"Area for y=sin(x) is: "<<area/1000000<<endl;
  cout<<" "<<endl;

  // Calculate then y = cos(x)
  area = 0.0;
  for(i=0; i<=100000; i++)
    {
      x = genrand64_real3()*2.0*M_PI;
      y = cos(x);
      area = area+y;
    }

  cout<<"Area for y=cos(x) is: "<<area/1000000<<endl;
  cout<<" "<<endl;

  // Calculate then area between the curves
  area = 0.0;
  for(i=0; i<=1000000; i++)
    {
      x = genrand64_real3()*2.0*M_PI;
      y = 2.0*genrand64_real3()-1.0;
      if( ((y>sin(x))&&(y<cos(x))) || ((y>cos(x))&&(y<sin(x))) )
	area++;
    }

  cout<<"Area between the curves is: "<<(4.0*M_PI*area/1000000)<<endl;
  cout<<" "<<endl;

}
      
      
  
