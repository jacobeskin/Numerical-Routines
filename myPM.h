// Park-Miller RNG with Schrage trick

#ifndef myPM_included
#define myPM_included

#include <stdlib.h>
#include <vector>

using namespace std;

vector<float> PM(int n, int seed)
{
  unsigned int m=2147483647;
  unsigned int a=16807;
  unsigned int q, r; 
  unsigned int I0, I1, d;
  int i;
  float I;
  vector<float> N(n);
  
  // Defining variables...
  q = m/a;
  r = m%a;

  // Generating n random numbers into array N
  I0 = seed;
  for(i=0;i<n;i++)
    {
      d = I0/q;
      I1 = a*I0-a*d*q-r*d;
      I = I1;
      N[i] = I/m;
      I0 = I1;
    }

  return N;
}

#endif
