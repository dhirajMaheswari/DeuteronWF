//This is the header file that calculates the Deuteron wf with Paris potential 
//also have functions that calculates the "s" and "d" partial waves


#ifndef PARIS_H
#define PARIS_H
#endif

#define INVHBARC	0.197328 //this is inverse (hbar*c) in GeV

#include <iostream>
#include <cmath>
 using namespace std;
//prototypes of the functions that calculates the wf and partial waves
float fd_Paris(float , int);
float U(float); //"s" partial wave

float W(float); //"d" partial wave
 
float C[13], D[13], BM[13]; //arrays declared

//definitions of the functions start here

 float fd_Paris(float X, int i) {
  float A = 0.0;
  float AA, B, CC;
   AA = B = CC = 0.0;
  float z, z1, y, y1, s, s1; //dummy variables
 int j; 
 float ans;
if(i==1) {
   C[0] =0.88688076; C[1] = -0.34717093; C[2]= -3.050238; C[3]= 56.207766; C[4] = -749.57334;
	C[5] = 5336.5279; C[6] = -22706.863;
	C[7]= 60434.4690; C[8] = -102920.58; C[9] = 112233.57; C[10] =-75925.226; C[11] = 29059.715; //note that the element at C[12] is zero here
  								 //will be calculated later below	
  for(j=0;j<=11;j++) { //label this for as F1
    A += C[j];
    } //end of for F1 
  C[12] = -A;
 D[0] = 0.023135193; D[1] = -0.85604572; D[2] = 5.6068193; D[3] = -69.462922; D[4] = 416.31118;
D[5] = -1254.6621; D[6] = 1238.783; D[7] = 3373.9172; D[8] = -13041.151; D[9] = 19512.524; //remaining 3 elements are zero here..will be calculated below 

  for(j = 0;j<=12;j++) { //this is F2
   BM[j] = 0.23162461 + j; //fill the BM array {in FORTRAN code, we had (j-1). We use j here to account for the fact
 				//that array index starts with 0 instead of 1}
  } //end of F2

  for(j = 0;j<=9;j++) { //this is F3
   AA = AA + D[j]/(powf(BM[j],2.0));
   B = B + D[j];
   CC = CC + D[j]*(powf(BM[j],2.0));
     } //end of F3
	//calculation of D[10]
 	z = (powf(BM[12],2.0) - powf(BM[10],2.0))* (powf(BM[11],2.0) - powf(BM[10],2.0));
 	z1 = -powf(BM[11],2.0) * powf(BM[12],2.0) * AA + (powf(BM[11],2.0) + powf(BM[12],2.0)) * B -CC;	
	D[10] = powf(BM[10],2.0)/z * z1; 
   	//calculation of D[11]
	y = (powf(BM[10],2.0) - powf(BM[11],2.0))* (powf(BM[12],2.0) - powf(BM[11],2.0));
	y1 = -powf(BM[12],2.0) * powf(BM[10],2.0) * AA + (powf(BM[12],2.0) + powf(BM[10],2.0)) * B -CC;
	D[11] = powf(BM[11], 2.0)/y * y1;
	//calculation of D[12]
	s = (powf(BM[11],2.0) - powf(BM[12],2.0))* (powf(BM[10],2.0) - powf(BM[12],2.0));
	s1 = -powf(BM[10],2.0) * powf(BM[11],2.0) * AA + (powf(BM[10],2.0) + powf(BM[11],2.0)) * B -CC;
	D[12] = powf(BM[12],2.0)/s * s1;
 	//printing of the matrices suppressed	
	cout<<"C"<<"\t";
	for(j=0;j<=12;j++) {
	cout<<C[j]<<"\t";
	}
	cout<<endl<<"D"<<"\t";
	for(j = 0;j<=12;j++) {
	cout<<D[j]<<"\t";
	}
	cout<<endl<<"BM"<<"\t";
	for(j = 0; j<=12;j++) {
	cout<<BM[j]<<"\t";
	}
	cout<<endl;
} //end of if
 ans = (powf(U(X/INVHBARC),2.0) + powf(W(X/INVHBARC),2.0))/powf(INVHBARC,3.0);
 return ans;	
} //end of function fd_Paris(X,i)


// "s" partial waves
 float U(float X) {
  float A,F, ans;
  int j;
 F = 0.79788456; A = 0.0;
	for(j=0;j<=12;j++) {
	  A = C[j] /(powf(X,2.0) + powf(BM[j],2.0)) + A;	
	}	
   ans = A * F /(sqrt(4.0*3.14159265));
	return ans;		
} //end of function U(X)

//"d" partial waves
 float W(float X) {
  float A,F, ans;
  int j;
  A = 0.0; F = 0.79788456;
  for(j=0;j<=12;j++) {
	A = D[j]/(powf(X,2.0) + powf(BM[j],2.0)) + A;
   }
  ans = A * F /(sqrt(4.0*3.14159265));
 return ans;
 }// end of funtion W(X)

