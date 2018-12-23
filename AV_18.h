//header file for the V_18 potential

#ifndef AV_18_H
#define AV_18_H
#endif

#define INVHBARC	0.197328 //this is inverse (hbar*c) in GeV

#include <iostream>
#include <math.h>
 using namespace std;
//prototypes of the functions that calculates the wf and partial waves
float fd_AV18(float , int);
float U_V18float); //"s" partial wave

float W_V18(float); //"d" partial wave
 
float C1[12], D1[12], BM1[12]; //arrays declared

//definitions of the functions start here

float fd_AV18(float X, int i) {
  float ans;
 if(i==1) {
	 C1[12] = { 0.706699e+00, -0.169743e+00, 0.112368e+01, -0.852995e+01, 0.195033e+02,-0.757831e+02, 0.283739e+03,
	   -0.694734e+03, 0.885257e+03, -0.720739e+03, 0.412969e+03, -0.103336e+03 };

	 D1[12] = {0.176655e-01, -0.124551e+00, -0.108815e+01, 0.384848e+01, -0.852442e+01, 0.209435e+02, -0.490728e+02,
	  0.577382e+02, -0.127114e+01, -0.628361e+02, 0.581016e+02, -0.177062e+02 };

	BM1[12] = {0.2316, 1.0, 1.5, 2.0, 2.5, 3.5, 4.5, 5.5, 6.5, 8.0, 9.5, 11.0 };
        ans = 0.0;
	return ans;
	}
 else{
   ans = (powf(U_18(X/INVHBARC,2.0) + powf(W_18(X/INVHBARC, 2.0))/powf(INVHBARC, 3.0);
  return ans;

  } 
} //end of the function

//"s" partial waves
 
 float U_V18(float X) {
  float A, F, ans;
  int j;	
  A = 0.0; F = 1.0;	
   for(j = 0;j<=11;j++) {
	A += D1[j]/(powf(x,2.0) + powf(BM1[j],2.0));
 } //end of for
 ans = A*F/ sqrt(4.0*M_PI);
 return ans;
} //end of function

//"d" partial waves

 float W_V18(float X) {
  float A, F, ans;
  int j;
  A = 0.0; F = 1.0;
  for(j = 0;j<=11;j++) {
    A += D1[j]/(powf(x,2.0) + powf(BM1[j],2.0));
  }
  ans = A*F/sqrt(4.0*M_PI);
 return ans;
} //end of funciton
