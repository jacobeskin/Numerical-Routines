// Linear congruential RNG

#ifndef myLCG_included
#define myLCG_included

#include <stdlib.h>
#include <vector>

using namespace std;

// LCG. 
vector<float> LCG(int n, int seed) // n is the # of random numbers desired
{
  unsigned int m=12386880; // Prime factors are 2, 3 , 5, 11, 17 and 23
  unsigned int c=982451653;       // Satisfies rule #1 
  int a=258061;   // Satisfies rules #2 and #3
  int i; 
  unsigned int I0, I1;
  vector<float> N(n);
  float I;

  // Generating n random numbers into array N
  I0 = seed;
  for(i=0;i<n;i++)
    {
      I1 = (a*I0+c)%m;
      // if(I1<0) {I1+=m;} // Found this solution online for negative results of %
      I = I1;           // Type conversion, otherwise doesn't seem to work...
      N[i] = I/m;
      I0 = I1;
    }

  return N;
}

#endif
