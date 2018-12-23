//header file for the Bonn potential

#ifndef BONN_H
#define BONN_H
#endif


#include <iostream>
#include <math.h>


//function prototypes

float fd_bonn (float, int);
float UUB(float);
float WWB(float);

float C2[11], D2[11], BM2[11];

using namespace std;
//function definitions begin here

float fd_bonn (float X, int i) {
 
 float tm0, tm1, tm2;	
 float q, w, e, q1, w1, e1, q2, w2, e2;
  float a, ans;
   if( i == 1) {
     C2[11] = {0.88472985, -0.26408759, -0.44114404e-01, -0.14397512e+02, 0.85591256e+02, -0.31876761e+03,
	       0.70336701e+03, -0.90049586+03, 0.66145441e+03, -0.25958894e+03 };
  	a = 0.0;
    for(int j = 0; j<=9; j++) {
	 a+= C2[j];
       }      
	C2[10] = -a;
    
	D2[11] = {0.22623762e-01, -0.50471056e+00, 0.56278897e+00, -0.16079764e+02, 0.11126803e+03,
		  -0.44667490e+03, 0.10985907e04, -0.16114995e+04   };
    for(int j = 0;j<=10; j++) {
	BM2[j] = 0.2315380 + (float) j * 0.9;
      }
     for(int k = 0; k<=7;++) {
	tm0+= D2[k];
	tm1+= D2[k]/powf(BM2[k], 2.0);
	tm2+= D2[k]* powf(BM2[k], 2.0);
      }
	q = (powf(BM2[10], 2.0) - powf(BM2[8],2.0)) * (powf(BM2[9],2.0) - powf(BM2[8],2.0));
	w = -powf(BM2[9],2.0)*powf(BM2[10],2.0)*tm1; e = (powf(BM2[9],2.0) + powf(BM2[10],2.0))*tm0 - tm2;
	D2[8] = powf(BM2[8],2.0)/q *(w + e);  

	q1 = (powf(BM2[8], 2.0) - powf(BM2[9],2.0)) * (powf(BM2[10],2.0) - powf(BM2[9],2.0));
	w1 = -powf(BM2[10],2.0)*powf(BM2[8],2.0)*tm1; e1 = (powf(BM2[10],2.0) + powf(BM2[8],2.0))*tm0 - tm2;
	D2[9] = powf(BM2[9]/q1 * (w1 + e1);

	q2 = (powf(BM2[9], 2.0) - powf(BM2[10],2.0)) * (powf(BM2[8],2.0) - powf(BM2[10],2.0));
	w2 = -powf(BM2[8],2.0)*powf(BM2[9],2.0)*tm1; e2 = (powf(BM2[8],2.0) + powf(BM2[9],2.0))*tm0 - tm2;
	D2[10] =  powf(BM2[10], 2.0) /q2 *(w2 + e2);						
   }//end of if
  else {
  ans = 0.0;
   }
 return ans;
} //end of fd_bonn function


//s waave definition
float UUB(float q) {
  float u, ans;
    u = 0.0;	 
   for(int k = 0; k<=10; k++) {
	u+= C2[k]/(powf(q,2.0) + powf(BM2[k],2.0));
      }
   ans = 1.0/(sqrt(2.0) * acos(-1.0)) * u;
   return ans;
 } //end of UUB function


//d wave definition
 float WWB(float q) {
  float u, ans;
   u = 0.0;
  for(int k = 0; k<=10;k++) {
   u+= D2[k]/(powf(q, 2.0) + powf(BM2[k],2.0));
  }
  ans = 1.0/(sqrt(2.0)*acos(-1.0) * u;
  return ans;	
}// end of function WWB

