#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"main.h"
//#include "Fifthlopp.h"
//#include "xmalloc.h"
//#include "xmalloc.c"
//#include "linspace.h"
//#include "linspace.c"


double *Ones(int n);// Ones function




void fifthLoop(int TNum, int DSNum, double OmegaR, double Gamma, double C_0, double *MaX, double *MaY, double *MaZ, double DX, double DY, double DZ, double DR, double OX, double OY,\
	 double OZ, double Theta, double *DataXR, double *DataYR, double *DataZR, double *DOrX, double *DOrY, double *DOrZ, double *DOr, double *DORStar, double *DOR,\
		 double *DORStarX, double *DORStarY, double *DORStarZ, double *DORX, double *DORY, double *DORZ, double *RGama)
{
	//double Tint = 1.0/30;
	

	

  for(int i=0;i<TNum;i++)
    {

		DataXR[i] = DR*cos(OmegaR*Tint*linSpace(TNum)[i]+atan2(DY,DX)+Theta);
        DataYR[i] = DR*sin(OmegaR*Tint*linSpace(TNum)[i]+atan2(DY,DX)+Theta);
		DataZR[i] = DZ*Ones(TNum)[i];
		DOrX[i] = OX-DataXR[i];
		DOrY[i] = OY-DataYR[i];
		DOrZ[i] = OZ-DataZR[i];
		
		DOr[i] = sqrt(pow(DOrX[i],2)+pow(DOrY[i],2)+pow(DOrZ[i],2));
		DORStar[i] = sqrt(pow(DOr[i],2)+pow(Gamma,2)*pow((MaX[0]*DOrX[i]+MaY[0]*DOrY[i]+MaZ[0]*DOrZ[i]),2))/Gamma;
		DOR[i] = pow(Gamma,2)*(DORStar[i]-(MaX[0]*DOrX[i]+MaY[0]*DOrY[i]+MaZ[0]*DOrZ[i]));
		
		DORStarX[i] = (DOrX[i]+pow(Gamma,2)*(MaX[0]*DOrX[i]+MaY[0]*DOrY[i]+MaZ[0]*DOrZ[i])*MaX[0])/(pow(Gamma,2)*DORStar[i]);
		DORStarY[i] = (DOrY[i]+pow(Gamma,2)*(MaX[0]*DOrX[i]+MaY[0]*DOrY[i]+MaZ[0]*DOrZ[i])*MaY[0])/(pow(Gamma,2)*DORStar[i]);
		DORStarZ[i] = (DOrZ[i]+pow(Gamma,2)*(MaX[0]*DOrX[i]+MaY[0]*DOrY[i]+MaZ[0]*DOrZ[i])*MaZ[0])/(pow(Gamma,2)*DORStar[i]);
		
		DORX[i] = pow(Gamma,2)*(DORStarX[i]-MaX[0]);
		DORY[i] = pow(Gamma,2)*(DORStarY[i]-MaY[0]);
		DORZ[i] = pow(Gamma,2)*(DORStarZ[i]-MaZ[0]);
		
		RGamma[i] = Tint*linSpace(TNum)[i]+DOR[i]/C_0;

		
	  
    }

}



/*function definition ofr ones 
double *Ones(int n)
{
	double *Bir;
	make_vector(Bir,n);
	return Bir;
}
*/
double *Ones(int Tnum)
{
  double *One;
  One = make_vector(One,Tnum);
  for(int i=0;i<Tnum;i++)
  {
	  One[i] = 1.0;
  }
  return One;
}

  
