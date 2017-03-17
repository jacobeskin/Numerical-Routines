# Numerical-Routines
Programs and modules of numerical methods 

Modules and source code files of numerical solutions to different problems. I will try and have a detailed description of the routines as comments inside the files, but will add short descritions here as well.

tridiag.f90: Three routines for solving second order differential equation that reduces to solving tridiagonal matrix equation; SOR, Thomas algorithm and Numerov method.

SE_H_bisectio.f90: Bisection method solution for the Schrödinger equation for Hydrogen atom. 

The C++ scripts are much shorter and easier to read. I did NOT write the Mersenne twister code and header, but I did write the other RNG's.

agF.cpp:
Calculate arithmetic-geometric mean of two input numbers ag(x,y) and evaluate 
the identity F=1/ag(1,sqrt(1-r^2)) proven by Gauss in 1799, where F is 
the hypergeometric function 2_F_1(a,b;c,d). F used is from the GSL.

fpi1.cpp:
Fixed point iteration solution to equation cos(x)=x

fpi2.cpp:
Fixed point iteration solution to exp(-x)=x

fpi3.cpp:
Fixed point iteration solution to 1-cosh(x)=x

ls.cpp:
Least squares fit line ax+b=y to the data in vectors x and y below

myLGC.h:
Linear congruential random number generator

myPM.h:
Park-Miller random number generator with Schrage trick

numdif.cpp:
Numerically approximate dy/dx from the data in h013.dat 

randnum.cpp:
Calls 3 different types of RNG's; LCG, Park-Miller and   Mersenne twister. 

rngave.cpp:
Calculating how the average value of random numbers from my RNG's. Create vectors of random numbers of different lengths 
and print these on a file for plotting.

MCint.cpp:
Calculate area of the curves y = sin(x) and y = cos(x) in the rectangle
0<x<2pi, -1<y<1 and then the area between these curve using Monte Carlo 
integration. Uses Mersenne Twister mt19937-64.c.

gslmat1.cpp:
Least squares solution for Ax=b, A is mxn matrix where m>n using  
two solution methods, GSL QR decomposition and direct GSL LSQ function.
y(t)=x1t+x2 with data {(2,1),(5,2),(7,3),(8,3)}.

