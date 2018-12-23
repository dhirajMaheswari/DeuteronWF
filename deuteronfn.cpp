#include "Paris.h"
//#include "AV_18.h"
#include <stdlib.h>

#define INVHBARC	0.197328 //this is inverse (hbar*c) in GeV

 struct COMPLEX {
   float real;
   float imag; };

//functions prototypes

float uu(float); //s waves
float wd(float); //d waves
COMPLEX DeuteronWf(int, int, int, int, float, float, float, int,int); //calculates deuteron wavefunction
/* in this function, the first parameter takes values of 1, 0 ,-1 : deuteron spin projection
  second parameter takes values of 1,-1: struck out nucleon isospin [1:proton;-1:neutron]	
  third parameter {1,-1}: spin projection of proton
  fourth parameter {1,-1} : spin projection of neutron
  fifth param: momentum in GeV/c sixth param: polar angle in radians
 seventh param: azimuthal angle in radians
 eighth param: {1,2,3,4}::{Paris, V18, cdbonn, V14} used to select the type of potential 
ninth param: {1,0}::{1->initializing, 0->calculating}
 the function returns a complex number
*/

float MomentumDistribution(float); //calculates the momentum distribution
//functions definitions

float uu(float p) {
float ans;
  ans = U(p/0.197328)/sqrt(powf(0.197328,3.0));
 return ans;
} //end of uu

float wd(float p) {
 float ans; 
  //try switching off the d waves	
 ans  = W(p/0.197328)/sqrt(powf(0.197328,3.0));
 return ans;
}//end of wd

COMPLEX DeuteronWf(int kj, int kt, int ksp, int ksn, float p, float theta, float phi, int iw, int ini) {
   
  float wf_re, wf_im;
  float xxx;  
  float fis; //fis==1 implies proton was struck, -1 implies neutron was struck 
 COMPLEX result;
   wf_re = wf_im = 0.0;
   if(ini == 1) {
	xxx = fd_Paris(0.0,1);
	}
     fis = 1.0; //assume proton was struck
   if(kt == -1) fis = -1.0; //neutron was struck
	
   if(kj == 1) { //A:   when deuteron spin projection is +1
 
	 //consider the case when both proton and neutron had their spins up	
	if((ksp == 1) && (ksn == 1)) { //B
	 wf_re = uu(p) + wd(p)/sqrt(8.0)*(3.0*powf(cos(theta),2.0)-1.0);
      wf_im = 0.0;
	} //end of B
  //consider the case when proton is up and neutron is down
 else if((ksp == 1) && (ksn == -1)) {//C
	wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi);   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi);   
 	}//end of C
 //consider the case when proton is down and neutron is up
 else if((ksp ==-1) && (ksn == 1)) { //D
	wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi);   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi);   
      }//end of D
 //consider the case when proton and neutron both are down
 else { //if(ksp ==-1 && ksn ==-1) { //E
	wf_re = wd(p)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*cos(2.0*phi);
      wf_im = wd(p)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*sin(2.0*phi);
      }//end of E
    } //A	 
else if(kj ==0) {//A1: deuteron spin projection is 0
	//consider the case when both protn and neutron are up
	if((ksp ==1) && (ksn ==1)) { //B1
	wf_re = wd(p) * 3.0 /2.0 * cos(theta)*sin(theta)*cos(phi);
	wf_im = - wd(p) * 3.0 /2.0 * cos(theta)*sin(theta)*sin(phi);
	} //end of B1
	//consider the case when proton is up neutron is down
	else if((ksp == 1) && (ksn ==-1)) { //C1
	wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0 * (3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	} //end of C
	//consider the case when proton is down neutron is up
	else if((ksp == -1) && (ksn == 1) ) { //D1
	wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0*(3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	} //end of D1
	//when proton and neutron both are down
	else {//if(ksp == -1 && ksn == -1) { //E1
	wf_re = - wd(p) * 3.0 /2.0 * cos(theta)*sin(theta)*cos(phi);
	wf_im = - wd(p) * 3.0 /2.0 * cos(theta)*sin(theta)*sin(phi);
	} //end of E1
     } //end of A1

   else  { //A2: deuteron spin projection is -1
	//consider when both spins are up
	if((ksp == 1) && (ksn == 1)) { //B2
	wf_re = wd(p)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*cos(2.0*phi);
	wf_im = - wd(p)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*sin(2.0*phi);
	} //end of B2
 	//proton up neutron down
    else if( (ksp == 1) && (ksn ==-1)) { //C2
	wf_re = - wd(p)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*cos(phi);
	wf_im = wd(p)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*sin(phi);
	}//end of C2
	//proton down neutron up
    else if((ksp == -1) && (ksn == 1)) { //D2
	wf_re = - wd(p)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*cos(phi);
	wf_im = wd(p)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*sin(phi);
	}// end of D2
   	//proton down neutron down
   else {//if( ksp == -1 && ksn == -1) { //E2
	wf_re = uu(p) + wd(p)/sqrt(8.0) * (3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	}//end of E2 
 	
} //end of A2



result.real = fis*wf_re;
result.imag = fis*wf_im;
 return result;
} //end of DeuteronWf function

float MomentumDistribution(float p) {
 COMPLEX WF;
 float res;
  float pi = acos(-1.0);
   int ini = 0; int kt = 1;
   int iw;	
   float theta = pi/10.0; float phi = pi/7.0;
   float sum = 0.0;
   int kj, ksp, ksn;
     for(kj = -1; kj<=1; kj++) { //1
	for(ksp = -1; ksp<=1; ksp+=2) { //2
	   for(ksn = -1; ksn <=1; ksn+=2) { //3
		WF = DeuteronWf(kj, kt, ksp, ksn, p, theta, phi, iw, ini);
		sum += powf(WF.real, 2.0) + powf(WF.imag, 2.0);       
	}//3		
     }//2
   } //1

 res = sum / 3.0; //averaged
 return res;
} //end of function


