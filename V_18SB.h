//sabina's version of V18 potential

#ifndef V_18SB_H
#define V_18SB_H
#endif

#include <iostream>
#include <math.h>

//function prototypes

float fd_V18SB(float, int);
float U_V18SB(float);
float W_V18SB(float);

float C3[12], D3[12], AM3[12];
  
//definitions

float fd_V18SB (float X, int i) {
 float sp2, sm, sm2;
  float a, b, cc, a1, b1, cc1, a2, b2, cc2, a3, b3, cc3;
   float tmp, tmp1, tmp2; 
  float ans;
if( i ==1) {
 AM3[12] = {0.232500e+00, 0.500000e+00, 0.800000e+00, 0.120000e+01, 0.160000e+01, 0.200000e+01,
	   0.4000e+01, 0.6000e+01, 0.1000e+02, 0.140000e+02, 0.180000e+02, 0.22000e+02	};
  
C3[12] = {0.105252223e+02, 0.124352529e+02, -0.687541641e+02, 0.239111042e+03, -0.441014422e+03, 0.300140328e+03,
	  -0.230639939e+03, 0.409671540e+03, -0.73345611e+03, 0.123506081e+04, -0.120520606e+04, 0.0  };

 for(int k = 0; k<=10; k++) {
  C3[12]-= C3[k];
  }
 
D3[12] = {0.280995496e+00, 0.334117629e-01, -0.727192237e+00, -0.302809607e+01, -0.903824982e+01, 0.496045967e+01,
	  -0.271985613e+02, 0.125334598e+03, -0.346742235e+03, 0.0, 0.0, 0.0 };	


  sp2 = sm = sm2 = 0.0;
 for(int lp = 0; lp<=0; lp++){
    sp2+= D3[lp]/powf(AM3[lp],2.0);
    sm+= D3[lp];	
    sm2+= D3[lp]*powf(AM3[lp],2.0);
  }
    a = powf(AM3[11], 2.0); b = powf(AM3[10], 2.0); cc = powf(AM3[9],2.0);
   tmp = (a - cc) * (b - cc); tmp1 = -b*a*sp2; tmp2 = (b + a)*sm - sm2;
  D3[9] =cc/tmp * (tmp1 + tmp2);

   a1 = powf(AM3[9],2.0); b1 = powf(AM3[11],2.0); c1 = powf(AM3[10],2.0);
    tmp = (a1 - cc1) * (b1 - cc1); tmp1 = -b1*a1*sp2; tmp2 = (b1 + a1)*sm - sm2;	
  D3[10] = cc/tmp *(tmp1 + tmp2);

  a3 = powf(AM3[10],2.0); b3 = powf(AM3[9],2.0); c3 = powf(AM3[11],2.0);
  tmp = (a3 - cc3) * (b3 - cc3); tmp1 = -b3*a3*sp2; tmp2 = (b3 + a3)*sm - sm2;	
  D3[11] = cc/tmp * (tmp1 + tmp2);

 for(int j = 0; j<=11;j++) {
 C3[j]/= 4.0*M_PI * sqrt(M_PI/2.0);
 D3[j]/=  4.0*M_PI * sqrt(M_PI/2.0);
 }
ans = 0.0;
return ans;
} //end of if
else {
float fg = 0.197328;
 ans = (powf(U_V18SB(X/fg),2.0) + powf(W_V18SB(X/fg),2.0))/powf(fg,3.0);
return ans;
}
}// end of function fd_V18SB

//s wave definition

float U_V18SB(float X) {
 float A = 0.0;
 for(int j = 0; j<=11; j++) {
  A+= C3[j]/(powf(X,2.0) + powf(AM3[j],2.0));
 float F = 1.0;
 ans = A*F/sqrt(4.0 * 3.14159265);
 return ans;
} //end of S wave function


float W_V18SB(float X) {
 float A = 0.0;
 for(int j = 0; j<=11; j++) {
  A+= D3[j]/(powf(X,2.0) + powf(AM3[j],2.0));
 float F = 1.0;
 ans = A*F/sqrt(4.0 * 3.14159265);
 return ans;
} //end of D wave function



