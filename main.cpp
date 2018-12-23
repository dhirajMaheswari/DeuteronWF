
 #include "deuteronfn.cpp"
#include <fstream>
#include <iomanip>
#include "gadap.cpp"

/*
extern "C" {
 //gaussian adaptive integration routine 
 //integrates funp from a to b with tolerance of eps and result is stored in sun_md!!
//the third parameter is the pointer to a function with one argument, that returns a floating point pointer
void  gadap_(float* a,float* b,float (*funp)(float*),float* eps,float* sun_md);
 }
*/

double p2d3p(double);

int main() {

 int kj, kt,ksp, ksn, IW, INI;
 float Momenta,  pp;
 float THETA, PHI;
 int j;
 float pi;
  double a, b, eps, sumd, normalized_integral;
  a = 0.250; b = 2.0; eps = 0.001;
 COMPLEX WF, ans;
 fstream datafile;
 pi = acos(-1.0);
 IW = 1; INI = 1;
 datafile.open("S_D_wf.dat", ios::out);
 WF = DeuteronWf(kj, kt,ksp, ksn, Momenta, THETA, PHI, IW, INI);
  for(j = 0;j<=1000;j+=50) { 
    Momenta = float (j)/1000.0;
    THETA = 90.0*pi/180.0;
     PHI = 45.0 * pi/180.0;	
    kj = 1; kt = 1; ksp = 1; ksn = 1; INI = 0;
 	ans = DeuteronWf(kj, kt, ksp, ksn, Momenta, THETA,PHI, IW,INI);
	cout<<Momenta<<"\t"<<ans.real<<"\t"<<ans.imag<<endl;	
	datafile<<Momenta<<"\t"<<uu(Momenta)<<"\t"<<wd(Momenta)<<endl;	 //"s" and "d" waves are written in the dat file
 } //end of for
 cout<<"***Momentum Distributions***"<<endl;
 //calculating the momentum distribution function for a range of momenta
	for(j = 360;j<=600; j+=20) { //A
	sumd = 0.0;
   	   pp = (float)j/1000.0; //momentum changed into GeV {also could be started in GeV units from beginning}
		//gadap_(&a, &b, p2d3p, &eps, &sumd);
		sumd = gadap(a, b, &p2d3p, eps);
		normalized_integral = 4.0 * 3.141952*sumd;
	     cout<<pp<<"\t"<<MomentumDistribution(pp)<<"\t"<<"Normalization check: "<<normalized_integral<<endl;
	} //A		  
 datafile.close();
 return 0;
} //end of main


double p2d3p(double x) {
  //This function is to be integrated to check the normalization ! 
 double ans;
  ans = MomentumDistribution(x)*powf(x,2.0);
  return ans;
 } 

