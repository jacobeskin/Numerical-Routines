// g++ -Wall -std=c++11 ls.cpp -o a 
/*
Least squares fit line ax+b=y to the data in vectors x and y below.
 */
#include <iostream>
#include <vector>
#include <valarray>
#include <math.h>

using namespace std;

int main()
{

  vector<double> x {1.0,2.0,3.0,4.0,5.0,6.0} ;
  vector<double> y {0.8969, 0.9005, 0.8961, 0.8919, 0.8793, 0.8818};
  vector<double> xy(6), z(6);
  double a, b, avex=21/6, sumx=21;
  int i;

  // Calculating a and b by "hand"
  for (i=0;i<6;i++)
    {
      xy[i]=x[i]*y[i];
      z[i]=pow(x[i]-avex,2);
    }

  // Calculating the sums needed for the formulas
  valarray<double> sumxy(&xy[0],6);
  valarray<double> sumy(&y[0],6);
  valarray<double> sumz(&z[0],6);

  // Calculating a and b
  a = (sumxy.sum()-((1/6)*sumx*sumy.sum()))/sumz.sum();
  b = (sumy.sum()-(a*21))/6;
  
  cout<<"Coefficients a and b:"<<endl;
  cout<<"a="<<a<<" "<<"b="<<b<<endl;
 
}
    

  

  
