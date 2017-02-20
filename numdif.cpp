// g++ -Wall numdif.cpp -o a
/* 
Numerically approximate dy/dx from the data in h013.dat 
*/
 
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
  double x[21], y[21], y1;
  int i;
  ifstream read("h013.dat");
  
  // Throwing away the first line by calling getline() once.
  string line;
  getline(read,line);
  
  // Reading the data into arrays x and y.
  for (i=0; i<21; i++)
    {
      read>>x[i]>>y[i];
    }

  // Doing the approximation y'(x0) = (y(x)-y(x0))/(x-x0) and writing
  // (x,y')- pairs on the screen.
  for (i=0;i<20;i++)
    {
      y1=(y[i+1]-y[i])/(x[i+1]-x[i]);
      cout<<x[i]<<" "<<y1<<endl;
    }

}
