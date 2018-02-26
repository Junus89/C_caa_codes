#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include "xmalloc.h"
//#include "xmalloc.c"
//#include "linspace.h"
//#include "linspace.c"


double *Ones(int n);// Ones function
//void *fifthLoop(int TNum, int DSNum, double OmegaR, double DX, double DY, double DZ, double DR, double OX, double OY,\
	 double OZ, double Theta, double *DataXR, double *DataYR, double *DataZR, double *DOrX, double *DOrY, double *DOrZ, double *DOr);



void fifthLoop(int TNum, int DSNum, double OmegaR, double Gamma, double *MaX, double *MaY, double *MaZ, double DX, double DY, double DZ, double DR, double OX, double OY,\
	 double OZ, double Theta, double *DataXR, double *DataYR, double *DataZR, double *DOrX, double *DOrY, double *DOrZ, double *DOr, double *DORStar, double *DOR,\
		 double *DORStarX, double *DORStarY, double *DORStarZ)
{
	double Tint = 1.0/30;
	

	

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

		
	  
    }
	//printf("DataXR = %g, DataYR = %g, DataZR = %g, DOrX = %g, DOrY = %g, DOrZ = %g, DOr = %g\n", DataXR[1], DataYR[1], DataZR[1], DOrX[1],DOrY[1],DOrZ[1],DOr[1]);
  



	

	//free_vector(DataXR);free_vector(DataYR);free_vector(DataZR);
	//free_vector(DOrX);free_vector(DOrY);free_vector(DOrZ);
    //printf("DataXR[3433] = %g, DataYR = %g, DataZR = %g\n",DataXR[3444],DataYR[3444],DataZR[3444]);
}



//function definition ofr ones 
double *Ones(int n)
{
	double *Bir;
	make_vector(Bir,n);
	return Bir;
}


  
