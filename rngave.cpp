/* Calculating how the average value of random numbers from my RNG's. Create vectors of random numbers of different lengths and print these on a file for plotting. 

Compilation:g++ rngave.cpp mt19937-64.c -o a

 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
// Including RNG header files
#include "myLCG.h"
#include "myPM.h"
#include "mt64.h"

using namespace std;

int main()
{
  unsigned int nlcg=1, npm=1, nmt=1, seed=8, i, b=10000;
  vector<float> lcg(1);
  vector<float> pm(1);
  vector<double> mt(1);
  float ave, sum;

  // Vectors for averages
  vector<float> avelcg(b);
  vector<float> avepm(b);
  vector<float> avemt(b);
   

  // LCG
  sum = 0;
  ofstream aveLCG;
  aveLCG.open("aveLCG.txt");
  do
    {
      lcg = LCG(nlcg, seed);
      for(i=0;i<nlcg;i++){sum += lcg[i];}
      ave = sum/nlcg;
      aveLCG<<nlcg<<" "<<ave<<endl;
      lcg.clear();
      nlcg++;
      lcg.resize(nlcg);
      sum = 0;
    }
  while(nlcg<b);
  aveLCG.close();

  // PM
  sum = 0;
  ofstream avePM;
  avePM.open("avePM.txt");
  do
    {
      pm = PM(npm, seed);
      for(i=0;i<npm;i++){sum += pm[i];}
      ave = sum/npm;
      avePM<<npm<<" "<<ave<<endl;
      pm.clear();
      npm++;
      pm.resize(npm);
      sum = 0;
    }
  while(npm<b);
  avePM.close();

  // MT
  init_genrand64(seed);
  sum = 0;
  ofstream aveMT;
  aveMT.open("aveMT.txt");
  do
    {
      for(i=0;i<nmt;i++){sum += genrand64_real3();} 
      ave = sum/nmt;
      aveMT<<nmt<<" "<<ave<<endl;
      nmt++;
      sum = 0;
    }
  while(nmt<b);
  aveMT.close();

}
  

