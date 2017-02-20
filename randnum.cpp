/* 
Calls 3 different types of RNG's; LCG, Park-Miller and   Mersenne twister. 

Compilation: g++ -Wall -std=c++11 randnum.cpp mt19937-64.c -o a

Example how to run the code:

$ ./a
Give number of RNG's and a seed:
3
856
 
Random numbers by LCG:
0.147294
0.193513
0.597149
 
Random numbers by PM:
0.00669937
0.596347
0.798156
 
Random numbers by MT:
0.453326
0.00888841
0.975014

*/

 

#include <iostream>
#include <stdlib.h>
#include <vector>
// Including RNG header files
#include "myLCG.h"
#include "myPM.h"
#include "mt64.h"

using namespace std;
   
int main()
{
  int n, seed, i;
  vector<float> lcg(n);
  vector<float> pm(n);

  // I/O
  cout<<"Give number of RNG's and a seed:"<<endl;
  cin>>n;
  cin>>seed;
  cout<<" "<<endl;

  // Calling my linear RNG's
  lcg = LCG(n, seed);
  pm = PM(n, seed);

  // Initializing Mersenne Twister
  init_genrand64(seed);

  // Printing out the random numbers
  cout<<"Random numbers by LCG:"<<endl;
  for(i=0;i<n;i++){cout<<lcg[i]<<endl;}

  cout<<" "<<endl;
  
  cout<<"Random numbers by PM:"<<endl;
  for(i=0;i<n;i++){cout<<pm[i]<<endl;}

  cout<<" "<<endl;

  cout<<"Random numbers by MT:"<<endl;
  for(i=0;i<n;i++){cout<<genrand64_real3()<<endl;}
  
  cout<<" "<<endl;
}
    
  

  
