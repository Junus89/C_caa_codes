#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */

#include"FDSolver.h"

/* ---------- function definitions --------- */


/*double **make_dmatrix(size_t m, size_t n) //2D array 
{
  double **a;
  make_vector(a, m+1);
  for(size_t i=0;i<m;i++)
    make_vector(a[i],n);
  a[m] = NULL;
  return a;
}

void free_dmatrix(double **a, size_t n)
{
  if(a!=NULL){
    for(size_t i=0;i<n;i++)
      free_vector(a[i]);
    free(a);
  }
}

/* 3D matrix definition *
double ***make_3Ddmatrix(size_t p, size_t q, size_t r)
{
  double ***a;
  make_vector(a,p+1); //make the pointer vector
  for(size_t i=0;i<p;i++)
    make_vector(a[i],q);
  for(size_t i=0;i<p;i++) // make the row vecotrs 
    for(size_t j=0;j<q;j++)
      make_vector(a[i][j],r);
  a[p]=NULL;
  return a;
}

/* free_3Ddmatrix *
void free_3Ddmatrix(double ***a, size_t q, size_t r)
{
  if(a!=NULL){
    for(size_t i=0;i<r;i++)
      for(size_t j=0;j<r;j++)
		  free_vector(a[i][j]);
    for(size_t i=0;i<q;i++)
      free_vector(a[i]);
    free(a);
  }
}

/* 4D array function definition *
double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r)
{
  double ****a;
  make_vector(a,o+1);//make the pointer vector
  for(size_t i=0;i<o;i++)
    make_vector(a[i],p);
  for(size_t i=0;i<o;i++)
    for(size_t j=0;j<p;j++)
      make_vector(a[i][j],q);
  for(size_t i=0;i<o;i++)
    for(size_t j=0;j<p;j++)
      for(size_t k=0;k<q;k++)
		  make_vector(a[i][j][k],r);
  a[o]=NULL;
  return a;
}

void free_4Ddmatrix(double ****a, size_t p, size_t q, size_t r)
{
	if(a!=NULL){
		for(int i=0;i<r;i++)
		    for(int j=0;j<r;j++)
		      for(int k=0;k<r;k++)
		        free(a[i][j][k]);
		for(int i=0;i<q;i++)
		    for(int j=0;j<q;j++)
		      free(a[i][j]);
		for(int i=0;i<p;i++)
		    free(a[i]);
		free(a);
	}
}*/

double *Ones(int Tnum)
{
  double *One;
  One = make_vector(One,Tnum);
  for(int i=0;i<Tnum;i++)
  {
	  One[i] = 1.0;
	  
  }
  free(One);
  return One;
}

//fifthloop function
/*double complex SP1T=Inicomp; double complex SP1T_sum=Inicomp; double complex SP1B=Inicomp; double complex SP1E=Inicomp;
double complex SP2T=Inicomp; double complex SP2T_sum=Inicomp; double complex SP2B=Inicomp; double complex SP2E=Inicomp;
double complex SP3T=Inicomp; double complex SP3T_sum=Inicomp; double complex SP3B=Inicomp; double complex SP3E=Inicomp; */




void fifthLoop(double OmegaR,double Omega, double MaX, double MaY, double MaZ, double ka, double pF[FNum][1])
{
	
	double DORStar[TNum],RGamma[TNum], Q[TNum], FRStarM[TNum], FRStarP[TNum], FRM[TNum], FRP[TNum];
	
	double complex SP1T=0.0+0.0*I; double complex SP1B=0.0+0.0*I; double complex SP1E=0.0+0.0*I;
	double complex SP2T=0.0+0.0*I; double complex SP2B=0.0+0.0*I; double complex SP2E=0.0+0.0*I;
	double complex SP3T=0.0+0.0*I; double complex SP3B=0.0+0.0*I; double complex SP3E=0.0+0.0*I;
	
    //FILE *fwrite2,*fwrite3;
    //fwrite2 = fopen("cexpVSexp.txt","w");
	//fwrite3 = fopen("SP1T.txt","w");
	
	

	
	
  for(int i=0;i<TNum;i++)
    {
		
		
		DataXR = DR*cos(OmegaR*Time[i]+atan2(DY,DX)+Theta);
        DataYR = DR*sin(OmegaR*Time[i]+atan2(DY,DX)+Theta);
		DataZR = DZ*Ones(TNum)[i];
		DOrX = OX-DataXR;
		DOrY = OY-DataYR;
		DOrZ = OZ-DataZR;
		DOr = sqrt(pow(DOrX,2)+pow(DOrY,2)+pow(DOrZ,2));

		//printf("DOR[%i] = %4.9f\n",i,DOR[i]);
		
		DORStarX = (DOrX+pow(Gamma,2)*(MaX*DOrX+MaY*DOrY+MaZ*DOrZ)*MaX)/(pow(Gamma,2)*DORStar[i]);
		DORStarY = (DOrY+pow(Gamma,2)*(MaX*DOrX+MaY*DOrY+MaZ*DOrZ)*MaY)/(pow(Gamma,2)*DORStar[i]);
		DORStarZ = (DOrZ+pow(Gamma,2)*(MaX*DOrX+MaY*DOrY+MaZ*DOrZ)*MaZ)/(pow(Gamma,2)*DORStar[i]);
		DORStar[i] = sqrt(pow(DOr,2)+pow(Gamma,2)*pow((MaX*DOrX+MaY*DOrY+MaZ*DOrZ),2))/Gamma;
		
		DORX = pow(Gamma,2)*(DORStarX-MaX);
		DORY = pow(Gamma,2)*(DORStarY-MaY);
		DORZ = pow(Gamma,2)*(DORStarZ-MaZ);
		DOR = pow(Gamma,2)*(DORStar[i]-(MaX*DOrX+MaY*DOrY+MaZ*DOrZ));
		
		RGamma[i] = Time[i]+DOR/C_0;
		//printf("RGamma[%i] = %4.9f\n",i,RGamma[i]);
		
		Vx = -DR*OmegaR*sin(OmegaR*Time[i]+atan2(DY,DX)+Theta);
		Vy = DR*OmegaR*cos(OmegaR*Time[i]+atan2(DY,DX)+Theta);
		Vz = 0*Ones(TNum)[i];
		
		Q[i] = cos(OmegaM*Time[i])*A;
		//printf("Q[%i] = %4.9f\n",i,Q[i]);
		
		Lx = 0*A*Ones(TNum)[i];
		Ly = 0*A*Ones(TNum)[i];
		Lz = 0*A*Ones(TNum)[i];
		
		FxM = Lx;
		FyM = Ly;
		FzM = Lz;
		
		FxP = -C_0*MaX*Q[i]; // here MaX, MaY, and MaZ are zero, is given from the file so the remaining result is zero
		FyP = -C_0*MaY*Q[i];
		FzP = -C_0*MaZ*Q[i];
		
		FRStarM[i] = FxM*DORStarX+FyM*DORStarY+FzM*DORStarZ;
		FRStarP[i] = FxP*DORStarX+FyP*DORStarY+FzP*DORStarZ;
		
		FRM[i] = FxM*DORX+FyM*DORY+FzM*DORZ;
		FRP[i] = FxP*DORX+FyP*DORY+FzP*DORZ;
		
		//--------------following the time and spacial integration --------------//
		

		SP1T += DT*(z1*Omega*Q[i]/DORStar[i])*cexp(-z1*Omega*RGamma[i]);
		SP1B = (DT/2)*(z1*Omega*Q[0]*cexp(-z1*Omega*RGamma[0]))/DORStar[0];
		SP1E = (DT/2)*(z1*Omega*Q[TNum-1]*cexp(-z1*Omega*RGamma[TNum-1]))/DORStar[TNum-1];
		SP1 = SP1T-SP1B-SP1E;
		
		SP2T += DT*(z1*ka*FRM[i]/DORStar[i]+FRStarM[i]/pow(DORStar[i],2)*cexp(-z1*Omega*RGamma[i]));
		SP2B = (DT/2)*(z1*ka*FRM[0]/DORStar[0]+FRStarM[0]/pow(DORStar[0],2))*cexp(-z1*Omega*RGamma[0]);
		SP2E = (DT/2)*(z1*ka*FRM[TNum-1]/DORStar[TNum-1]+FRStarM[TNum-1]/pow(DORStar[TNum-1],2))*cexp(-z1*Omega*RGamma[TNum-1]);
		SP2 = SP2T-SP2B-SP2E;
		
		
		SP3T += DT*(z1*ka*FRP[i]/DORStar[i]+FRStarP[i]/pow(DORStar[i],2)*cexp(-z1*Omega*RGamma[i]));
		SP3B = (DT/2)*(z1*ka*FRP[0]/DORStar[0]+FRStarP[0]/pow(DORStar[0],2))*cexp(-z1*Omega*RGamma[0]);
		SP3E = (DT/2)*(z1*ka*FRP[TNum-1]/DORStar[TNum-1]+FRStarP[TNum-1]/pow(DORStar[TNum-1],2))*cexp(-z1*Omega*RGamma[TNum-1]);
		SP3 = SP3T-SP3B-SP3E;
		
	  
    }
	//fclose(fwrite2);
	//fprintf(fwrite3,"%4.8f %4.8f\n",creal(SP1T),cimag(SP1T)); 
	
	//fclose(fwrite3);

}


