//This code here is the translation of the gadap routines written in FORTRAN into C++. 


using namespace std;
double  gadap (double , double , double (*pfun)(double), double );
double gadap2(double ,double , double (*pfl)(double), double (*pfu)(double), double (*pf)(double,double),double ); 
double dsum(double ,double , double , double , double ); 
double  fgadap (double x,double a0,double b0,double (*pfun)(double,double),double eps);
//**********************//

double gadap(double a0, double b0, double (*pfun)(double), double eps) {
// **********************************************************************    
//   
//   THE FOLLOWING INTEGRATION ROUTINES WHERE OBTAINED FROM THE  
//   LUND UNIVERSITY COMPUTER CENTER.    
//   It  was transferred to C++ on November 15
//   FIU, Miami, FL
// **********************************************************************    
//.......................................................................    
//   
//   PURPOSE           - INTEGRATE A FUNCTION pfun(X)   
//   METHOD            - ADAPTIVE GAUSSIAN   
//   USAGE             - sum = GADAP(A0,B0,F,EPS) 
//   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)    
//               B0    - UPPER LIMIT (INPUT,REAL)    
//               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE 
//                       SUPPLIED BY THE USER. (INPUT, double REAL FUNCTION) 
//               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS  
//                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT, double REAL) 
//               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT, double REAL) 
//   PRECISION         - double  
//   REQ'D PROG'S      - C++   
//   AUTHOR            - THOMAS JOHANSSON, LDC,1973  
//   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71    
//   Adaptaion to C++  - Misak Sargsian, FIU, Miami, FL  Nov, 2004
//.......................................................................    
    
      double a[300], b[300], f1[300], f2[300], f3[300], s[300];
      int    n[300], l, i, ifu;
      double red, sum, c;
      double w1,u2,ss,sold;
      
      red = 1.3;
      l=1;
      i=1;
      sum=0;
      c = sqrt(15.0)/5.;
      a[1]  = a0;
      b[1]  = b0;
      f1[1] = pfun(0.5*(1+c)*a0+0.5*(1-c)*b0); 
      f2[1] = pfun(0.5*(a0+b0));  
      f3[1] = pfun(0.5*(1-c)*a0+0.5*(1+c)*b0);    
      ifu   = 3; 
      s[1]  = dsum(f1[1],f2[1],f3[1],a0,b0);
      //      100 continue
    label100:
      l    =  l+1; 
      n[l] =  3;    
      eps=eps*red;   
      a[i+1] = a[i]+c*(b[i]-a[i]); 
      b[i+1] = b[i];   
      a[i+2] = a[i]+b[i]-a[i+1];   
      b[i+2] = a[i+1]; 
      a[i+3] = a[i];   
      b[i+3] = a[i+2]; 

      w1=a[i]+(b[i]-a[i])/5.;
      u2=2.*w1-(a[i]+a[i+2])/2.;

        f1[i+1] = pfun(a[i]+b[i]-w1);
	f2[i+1] = f3[i];
        f3[i+1] = pfun(b[i]-a[i+2]+w1);
	f1[i+2] = pfun(u2) ;
	f2[i+2] = f2[i] ;
	f3[i+2] = pfun(b[i+2]+a[i+2]-u2)  ; 
	f1[i+3] = pfun(a[i]+a[i+2]-w1) ;
	f2[i+3] = f1[i] ;
	f3[i+3] = pfun(w1); 
	ifu = ifu + 6; 

	//      IF(IFU.GT.5000) GOTO 130
	if(ifu>5000) goto label130;
	s[i+1]=  dsum(f1[i+1],f2[i+1],f3[i+1],a[i+1],b[i+1]);  
	s[i+2]=  dsum(f1[i+2],f2[i+2],f3[i+2],a[i+2],b[i+2]);  
	s[i+3]=  dsum(f1[i+3],f2[i+3],f3[i+3],a[i+3],b[i+3]);  
	ss=s[i+1]+s[i+2]+s[i+3];   
	i=i+3; 

	//      IF(I.GT.300)GOTO 120
	if(i>300) goto label120;
	sold=s[i-3];   
	//      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
	if(fabs(sold-ss)>eps*(1.+fabs(ss))/2.) goto label100;  
	sum=sum+ss;    
	i=i-4; 
	n[l]=0;    
	l=l-1; 	
	//  110 CONTINUE 
    label110:
	//      IF(L.EQ.1) GOTO 130   
	if(l==1) goto label130;
	n[l]=n[l]-1;   
	eps=eps/red;   
	//      IF(N(L).NE.0) GOTO 100    
	if(n[l]!=0)goto label100;
	i=i-1; 
	l=l-1; 
	//      GOTO 110
	goto label110;
	//  120 CONTINUE
     label120:
	//C      WRITE(6,1)    
	// 130  RETURN    
	//      END   

	//double result =1.0;
      //     result = pfun(a0);
    label130:
	//     return result;
     return sum;
     }

 //************************************//
double  dsum(double f1f,double f2f,double f3f,double aa,double bb) {

return 5./18.*(bb-aa)*(f1f+1.6*f2f+f3f);
} 

//++++++++++++++++++++++++++//

double gadap2(double a0,double b0, double (*pfl)(double), double (*pfu)(double), double (*pf)(double,double),double eps) {  
     //       int iifu;

//   PURPOSE           - INTEGRATE A FUNCTION F(X,Y) OF TWO VARIABLES    
//   METHOD            - ADAPTIVE GAUSSIAN IN BOTH DIRECTIONS    
//   USAGE             - CALL GADAP2(A0,B0,FL,FU,F,EPS,SUM)  
//   PARAMETERS  A0    - LOWER X-LIMIT (INPUT,REAL)  
//               B0    - UPPER X-LIMIT (INPUT,REAL)  
//               FL    - USER SUPPLIED FUNCTION FL(X) GIVING THE LOWER   
//                      Y-LIMIT FOR A GIVEN X-VALUE 
//                       (INPUT,REAL FUNCTION)   
//               FU    - USER SUPPLIED FUNCTION FU(X) GIVING THE UPPER   
//                       Y-LIMIT FOR A GIVEN X-VALUE 
//                       (INPUT,REAL FUNCTION)   
//               F     - USER SUPPLIED FUNCTION F(X,Y) TO BE INTEGRATED  
//                       (INPUT,REAL FUNCTION)   
//               EPS   - DESIRED ACCURACY (INPUT,REAL)   
//               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL) 
//   PRECISION         - SINGLE  
//   REQ'D PROG'S      - FL,FU,F,FGADAP  
//   AUTHOR            - THOMAS JOHANSSON, LDC,1973  
//   
//.......................................................................    

      double a[300], b[300], f1[300], f2[300], f3[300], s[300];
      int   iifu,  n[300],l,i;
      double red,sum,c,x,ay,by,w1,u2,ss,sold;      
      
      red=1.4;   
      l   = 1;   
      i   = 1;   
      sum = 0.0;    
      c   = sqrt(15.0)/5.;
      a[1] = a0;   
      b[1] = b0;   
      x    = 0.5*(1+x)*a0+0.5*(1-c)*b0;   
      ay   = pfl(x);  
      by   = pfu(x);  
      f1[1]= fgadap(x,ay,by,pf,eps);   
      x    = 0.5*(a0+b0); 
      ay   = pfl(x);  
      by   = pfu(x);  
      f2[1]= fgadap(x,ay,by,pf,eps);   
      x    = 0.5*(1-c)*a0+0.5*(1+c)*b0;   
      ay   = pfl(x);  
      by   = pfu(x);  
      f3[1]= fgadap(x,ay,by,pf,eps);   
      iifu = 3; 
      s[1] = dsum(f1[1],f2[1],f3[1],a0,b0);  
 label100:
      l    = l+1; 
      n[l] = 3;    
      eps  = eps*red;   
      a[i+1] = a[i]+c*(b[i]-a[i]); 
      b[i+1] = b[i];   
      a[i+2] = a[i]+b[i]-a[i+1];   
      b[i+2] = a[i+1]; 
      a[i+3] = a[i];   
      b[i+3] = a[i+2]; 
      w1     = a[i]+(b[i]-a[i])/5.;    
      u2     = 2.*w1-(a[i]+a[i+2])/2.; 
      x      = a[i]+b[i]-w1;    

      ay     = pfl(x);  
      by     = pfu(x);  
      f1[i+1]= fgadap(x,ay,by,pf,eps); 
      f2[i+1]= f3[i]; 
      x      = b[i] - a[i+2] + w1;  
      ay     = pfl(x);  
      by     = pfu(x);  
      f3[i+1]= fgadap(x,ay,by,pf,eps); 
      x      = u2;  
      ay     = pfl(x);  
      by     = pfu(x);  
      f1[i+2]= fgadap(x,ay,by,pf,eps); 
      f2[i+2]= f2[i]; 

      x      = b[i+2]+a[i+2]-u2;    
      ay     = pfl(x);  
      by     = pfu(x);  
      f3[i+2]= fgadap(x,ay,by,pf,eps); 
      x      = a[i]+a[i+2]-w1;  
      ay     = pfl(x);  
      by     = pfu(x);  
      f1[i+3]= fgadap(x,ay,by,pf,eps); 
      f2[i+3]= f1[i]; 
      x      = w1;  
      ay     = pfl(x);  
      by     = pfu(x);  
      f3[i+3]= fgadap(x,ay,by,pf,eps); 
      iifu   = iifu+6;
      if(iifu>5000) goto label130;  

      s[i+1] = dsum(f1[i+1],f2[i+1],f3[i+1],a[i+1],b[i+1]);  
      s[i+2] = dsum(f1[i+2],f2[i+2],f3[i+2],a[i+2],b[i+2]);  
      s[i+3] = dsum(f1[i+3],f2[i+3],f3[i+3],a[i+3],b[i+3]);  
      ss     = s[i+1]+s[i+2]+s[i+3];   
      i      = i+3; 
      if(i>300)goto label120;  
      sold=s[i-3];   
      if(fabs(sold-ss)>eps*(1.+fabs(ss))/2.) goto label100;  

      sum  = sum+ss;    
      i    = i-4; 
      n[l] = 0;  
      l    = l-1; 
  label110:
      if(l==1) goto label130;   
      n[l] = n[l]-1;   
      eps  = eps/red;   
      if(n[l]!=0) goto label100;    
      i=i-1; 
      l=l-1; 
	goto label110;  
    label120:
    label130:
      return sum;
}

      double  fgadap (double x,double a0,double b0,double (*pfun)(double,double),double eps) {
//      COMMON/GADAP_2/ NUM,IFU    
//      EXTERNAL F    

      double a[300], b[300], f1[300], f2[300], f3[300], s[300];
      int    iifu,n[300],l,i;
      double red,sum,c,ay,by,w1,u2,ss,sold,fgdp;      
      //      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      //      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      red=1.4;   
      l     = 1;   
      i     = 1;   
      sum   = 0.;    
      c     = sqrt(15.)/5.;    
      a[1]  = a0;   
      b[1]  = b0;   
      f1[1] = pfun(x,0.5*(1+c)*a0+0.5*(1-c)*b0);  
      f2[1] = pfun(x,0.5*(a0+b0));    
      f3[1] = pfun(x,0.5*(1-c)*a0+0.5*(1+c)*b0);  
      iifu  = 3; 
      s[1]=  dsum(f1[1],f2[1],f3[1],a0,b0);  
	label100:
      l=l+1; 
      n[l]=3;    
      eps=eps*red;   
      a[i+1]=a[i]+c*(b[i]-a[i]); 
      b[i+1]=b[i];   
      a[i+2]=a[i]+b[i]-a[i+1];   
      b[i+2]=a[i+1]; 
      a[i+3]=a[i];   
      b[i+3]=a[i+2]; 

      w1=a[i]+(b[i]-a[i])/5.;    
      u2=2.*w1-(a[i]+a[i+2])/2.; 

      f1[i+1]=pfun(x,a[i]+b[i]-w1); 
      f2[i+1]=f3[i]; 
      f3[i+1]=pfun(x,b[i]-a[i+2]+w1);   
      f1[i+2]=pfun(x,u2);   
      f2[i+2]=f2[i]; 
      f3[i+2]=pfun(x,b[i+2]+a[i+2]-u2); 
      f1[i+3]=pfun(x,a[i]+a[i+2]-w1);   
      f2[i+3]=f1[i]; 
      f3[i+3]=pfun(x,w1);   
      iifu=iifu+6; 
	if(iifu>5000) goto label130;  
	s[i+1]=  dsum(f1[i+1],f2[i+1],f3[i+1],a[i+1],b[i+1]);  
	s[i+2]=  dsum(f1[i+2],f2[i+2],f3[i+2],a[i+2],b[i+2]);  
	s[i+3]=  dsum(f1[i+3],f2[i+3],f3[i+3],a[i+3],b[i+3]);  
	ss=s[i+1]+s[i+2]+s[i+3];   
	i=i+3; 
	if(i>300)goto label120;  
	sold=s[i-3];   
	if(fabs(sold-ss)>eps*(1.+fabs(ss))/2.) goto label100;  
	sum=sum+ss;    
	i=i-4; 
	n[l]=0;    
	l=l-1; 
     label110:
	if(l==1) goto label130;   
	n[l]=n[l]-1;   
	eps=eps/red;   
	if(n[l]!=0) goto label100;    
	i=i-1; 
	l=l-1; 
	goto label110;  
	label120:
	label130:
      fgdp = sum;
      eps=eps/red;   
      return fgdp;    
	}




