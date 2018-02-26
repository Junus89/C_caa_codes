#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */


/* ----- including -----header defined functions ------*/
#include "xmalloc.h"
#include "xmalloc.c"
#include "nextpower.h"
#include "nextpower_c.c"
#include "FDSolver.h"
#include "FDSolver.c"


int main()
{
	
	/*reading ConstantVariable file including file names */
    FILE *fp;
    char fdName[128], sdName[128], odName[128];
	//int BNum, TNum, FNum;
	//double MaX, MaY, MaZ;
    fp = fopen("ConstantV.txt","r");
    while(1 == fscanf(fp,"%s%*[^\n] %s%*[^\n] %s%*[^\n] %d%*[^\n] %lf%*[^\n] %d%*[^\n] %d%*[^\n] %lf%*[^\n] %lf%*[^\n] %lf%*[^\n]\n",fdName, sdName, odName,\
    	 &BNum,&R,&TNum,&FNum,&MaX,&MaY,&MaZ)){
    }
	
	fclose(fp);
	
	
  // opening and reading flowdata files
  FILE *fp_f;
  fp_f = fopen("flowdata.txt","r");
  //storing flowdata 
  flowData FlowD[]={0};
  for(int i =0;i<sizeof(FlowD);i++)
    {
      fscanf(fp_f,"%lf %d %lf",&FlowD[i].rho_0,&FlowD[i].c_0,&FlowD[i].RMaTip);
      //    printf(" Rho_0 = %4.4f, c_0 = %d, RMaTip = %4.4f\n",(double)FlowD[i].rho_0,FlowD[i].c_0,(double)FlowD[i].RMaTip);
    }
  printf("--------------FlowData------------\n\n");
  printf(" Rho_0 = %4.4f, c_0 = %d, RMaTip = %4.4f\n\n",(double)FlowD[0].rho_0,FlowD[0].c_0,(double)FlowD[0].RMaTip);
  double Rho_0 = FlowD[0].rho_0;
  C_0 = FlowD[0].c_0;
  double RMaTip = FlowD[0].RMaTip;
  fclose(fp_f);


  double OmegaR = RMaTip*C_0/R;
  double fR = OmegaR/(2*PI);         // rotation frequency
  double TR = 1/fR;
  OmegaM = 150*2*PI;         
  double fM = 150;                  //pulsation frequency
  printf("---------Matematical relations ------------\n\n");
  printf("OmegaR = %4.4f [rad/s], fR = %4.4f [hz], TR = %4.4f [s], OmegaM = %4.4f [rad/s]\n\n",OmegaR, fR, TR, OmegaM);
  

  /* Mathematical relations of TNum data */
  double Tint = TR;
  DT = Tint/(TNum-1);          //...discrete time or sampling time? or sampling time step?


  //double Time[TNum];
  for(int i=0;i<TNum;i++)
  {
	  Time[i]=i/(TNum*1.0-1.0)*Tint;
	  printf("Time[%d] = %4.6f\n",i,(double)Time[i]);
  }
    
  
  int NFFT = nextPower2(FNum);
  printf("the next power of FNum->NFFT = %d\n",NFFT);
  double ODT = Tint/NFFT;
  printf("ODT = %lf\n",ODT);
  /* for construction of OTime */
  //double *OTime;
  
  OTime = make_vector(OTime,NFFT-1); /* calling forOTime function and initializing it */
  for(int i=0;i<NFFT-1;i++)
    {
		OTime[i] = ODT*i;
      printf("OTime[%d] = %4.4f\n",i,OTime[i]);
    }
  printf("----checking ----> OTime[5] = %4.4f\n",OTime[5]); /*checking for OTime */
  double DF = ((1/ODT)/2)*(1/(1.0*FNum/2)); /* here as FNum is int type, should multiply with 1.0 to get the double type DF result */
  printf("----checking ----> DF = %4.9f\n",DF);
  
  /* reading WriteDataSMeshBlades.dat data */
  FILE *fp_meshB;
  fp_meshB = fopen("WriteDataSMeshBlade.txt","r");
  //storing the data into the structure
  WDataSMB BladeMData[]={0};
  for(int i=0;i<sizeof(BladeMData);i++)
    {
      fscanf(fp_meshB,"%lf %lf %lf %lf %lf %lf",&BladeMData[i].DataX, &BladeMData[i].DataY, &BladeMData[i].DataZ, &BladeMData[i].nxData,&BladeMData[i].nyData,&BladeMData[i].nzData);

    }
  
  printf("----checking ---> DataX = %lf, DataY= %lf, DataZ = %lf\n, nxData = %lf\n, nyData = %lf\n, nzData = %lf\n",BladeMData[0].DataX, BladeMData[0].DataY, BladeMData[0].DataZ,BladeMData[0].nxData,BladeMData[0].nyData,BladeMData[0].nzData);
  int DSNum = 1;
  fclose(fp_meshB);


  /* reading WriteObserverSMesh.dat data  */
  double XO_temp;
  double YO_temp;
  double ZO_temp;
  double *XO = NULL;
  double *YO = NULL;
  double *ZO = NULL;
  int index, OSNum = 0, j = 0, i = 1;

  FILE *fp_ObMesh;
  fp_ObMesh = fopen("WriteObserverSMesh.txt","r");
  //storing the data into the structure
  while(fscanf(fp_ObMesh,"%lf %lf %lf",&XO_temp, &YO_temp, &ZO_temp)!=EOF)
    {
      if(XO == NULL && YO == NULL && ZO == NULL)
	{

	  XO = malloc(sizeof(XO_temp));
	  YO = malloc(sizeof(YO_temp));
	  ZO = malloc(sizeof(ZO_temp));
	  *XO = XO_temp;
	  *YO = YO_temp;
	  *ZO = ZO_temp;
	}/* in case the memory is not enough, it hase to be reallocated */
      else
	{
	  i++;
	  XO = realloc(XO,sizeof(XO)*i);
	  YO = realloc(YO,sizeof(YO)*i);
	  ZO = realloc(ZO,sizeof(ZO)*i);
	  index = i-1;
	  *(XO+index)=XO_temp;
	  *(YO+index)=YO_temp;
	  *(ZO+index)=ZO_temp;
	}
      /* showing and checking the result */
    }
  XO_temp = 0.0;
  YO_temp = 0.0;
  ZO_temp = 0.0;
  while(index>=0)
    {
      printf("[%d]: XO = %lf, YO = %lf, ZO = %lf\n",j, XO[j], YO[j], ZO[j]);
      index--;
      OSNum = 1+j;
      j++;
      XO_temp++;
      YO_temp++;
      ZO_temp++;
    }
  printf("----- checking ----> OSNum = %d\n",OSNum);
  printf("----- checking ----> [53]: XO[53] = %lf, YO[53] = %lf, ZO[53] = %lf\n", XO[53], YO[53], ZO[53]);
  printf("----- checking ----> [99]: XO[99] = %lf, YO[99] = %lf, ZO[99] = %lf\n", XO[99], YO[99], ZO[99]);
  printf("----- checking ----> [149]: XO[149] = %lf, YO[1499] = %lf, ZO[149] = %lf\n", XO[149], YO[149], ZO[149]);
  fclose(fp_ObMesh);
  /* freeing the allocated memory for 'WritesObserverSMeshData.dat */
  


  Gamma = sqrt(1/(1-(pow(MaX,2)+pow(MaY,2)+pow(MaZ,2)))); 
  printf("Gamma = %lf\n",Gamma);
  /* reading data for ObserverSmeshparameter.dat data */
  
  int ObsrSMeshD, *ObsrSMeshD_ptr, Counter_for_ObsrSMeshD=1, idx_for_ObsrSMeshD, k_for_ObsrSMeshD;
  FILE *fp_ObsrSMeshD;
  fp_ObsrSMeshD = fopen("ObserverMeshParameter.txt","r");
  
  while(fscanf(fp_ObsrSMeshD,"%d\n",&ObsrSMeshD)!=EOF)
    {
      if(ObsrSMeshD_ptr==NULL)
	{
	  ObsrSMeshD_ptr = malloc(sizeof(ObsrSMeshD));
	  *ObsrSMeshD_ptr = ObsrSMeshD;
	}
      else 
	{
	  Counter_for_ObsrSMeshD++;
	  ObsrSMeshD_ptr = realloc(ObsrSMeshD_ptr,sizeof(ObsrSMeshD)*Counter_for_ObsrSMeshD);
	  idx_for_ObsrSMeshD = Counter_for_ObsrSMeshD-1;
	  *(ObsrSMeshD_ptr+idx_for_ObsrSMeshD)=ObsrSMeshD;
	}
    }
  ObsrSMeshD = 0;
  int ObsrSthetaNum = ObsrSMeshD_ptr[1];
  int ObsrSFaiNum = ObsrSMeshD_ptr[2];
  double HalfNum = (ObsrSthetaNum-1)*(ObsrSFaiNum+1)/2;
  printf("\n-------->Checking:\n ObsrSthetaNum = %d\n ObsrSFeinum = %d\n HalfNum = %lf\n",ObsrSthetaNum, ObsrSFaiNum, HalfNum);
  
  
    
  /*----------------------------------
  double **Area,**UnitNormal;
  Area = make_dmatrix(DSNum,1);
  UnitNormal = make_dmatrix(DSNum,3);
  for(int m=1;m<=DSNum;m++)
    {
      Area[m][0] = sqrt(pow(BladeMData[m].nxData,2)+pow(BladeMData[m].nyData,2)+pow(BladeMData[m].nzData,2));
      UnitNormal[m][0] = (double)BladeMData[m].nxData/Area[m][0];
      UnitNormal[m][1] = (double)BladeMData[m].nyData/Area[m][0];
      UnitNormal[m][2] = (double)BladeMData[m].nzData/Area[m][0];
      printf("UnitNormal[%d][0] = %g\n",UnitNormal[m][0]);
    }
  double **DataSArea, **DataSVector;
  DataSArea = Area;
  DataSVector = UnitNormal;
  
  ---------------------------------*/

  PBD1 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  
  PBD2 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD3 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  

  PB1 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB2 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB3 = make_3Ddmatrix(FNum,OSNum,BNum);
  

  P1 = make_dmatrix(FNum,OSNum);
  P2 = make_dmatrix(FNum,OSNum);
  P3 = make_dmatrix(FNum,OSNum);
  
  


  make_vector(DataXR,TNum);
  make_vector(DataYR,TNum);
  make_vector(DataZR,TNum);

  
  make_vector(DOrX,TNum);
  make_vector(DOrY,TNum);
  make_vector(DOrZ,TNum);
  make_vector(DOr,TNum);
  
  make_vector(DORStar,TNum);
  make_vector(DOR,TNum);
  
  make_vector(DORStarX,TNum);
  make_vector(DORStarY,TNum);
  make_vector(DORStarZ,TNum);
  
  make_vector(DORX,TNum);
  make_vector(DORY,TNum);
  make_vector(DORZ,TNum);
  make_vector(RGamma,TNum);
  
  make_vector(Vx,TNum);
  make_vector(Vy,TNum);
  make_vector(Vz,TNum);
  make_vector(Q,TNum);
  
  make_vector(Lx,TNum);
  make_vector(Ly,TNum);
  make_vector(Lz,TNum);
  make_vector(FxM,TNum);
  make_vector(FyM,TNum); 
  make_vector(FzM,TNum); 
  make_vector(FxP,TNum); 
  make_vector(FyP,TNum); 
  make_vector(FzP,TNum);
  make_vector(FRStarM,TNum); 
  make_vector(FRStarP,TNum);
  make_vector(FRM,TNum); 
  make_vector(FRP,TNum);
  
  pF = make_dmatrix(FNum,1);
  
  //double **pXOYM,**Op,*OF;
  Op = make_dmatrix(FNum,1);
  pXOYM = make_dmatrix(FNum,1);

 
  for(int n =0; n< FNum;n++)
  {
	  double Omega= OmegaM+BNum*(n-5)*OmegaR;
	  double ka = Omega/C_0;
	  //printf("Omega[%d]=%4.4g   ka[%d] = %4.4g\n",n,Omega[n],n,ka[n]);
  
	  double complex Transa1 = 0.0 + 0.0*I;
	  double complex Transb1 = 0.0 + 0.0*I;
	  double complex Transc1 = 0.0 + 0.0*I;

	  double complex Transa2 = 0.0 + 0.0*I;
	  double complex Transb2 = 0.0 + 0.0*I;
	  double complex Transc2 = 0.0 + 0.0*I;
  
	  
	  for(int m=20;m<21;m++)//OSNum=
	  {
		  OX = XO[m];//the first value of XO vector
		  OY = YO[m];// the first value of YO vector
		  OZ = ZO[m];// the first value of ZO vector

		  
		  for(int k=0;k<BNum;k++)
		  {
			  Theta = 2*((k+1)-1)*PI/BNum; // Theta = 2*(k-1)*PI/BNum;
			  
		  
		  		for(int j=0;j<DSNum;j++)
		  	  	{
		  			printf("Harmonic number %d, Observer number %d, Blade number %d, Source number %d\n",n,m+1,k+1,j+1); // since m, k, j starts from zero, add it with 1
					
					DX=BladeMData[j].DataX;
					DY=BladeMData[j].DataY;
					DZ=BladeMData[j].DataZ;
					DR=sqrt(pow(DX,2)+pow(DY,2));
					

					
					/*--------DVX, DVY, DVZ; about unitnormal *----------*/
					
					//double A = 1.0;//Area[j][1];
				
					// calling fifthLoop function
					fifthLoop(OmegaR, Omega, MaX, MaY, MaZ, ka, pF);
					
					
					PBD1[n][m][k][j] = SP1;
					PBD2[n][m][k][j] = SP2;
					PBD3[n][m][k][j] = SP3;
					
					Transa1 += PBD1[n][m][k][j];
					Transb1 += PBD2[n][m][k][j];
					Transc1 += PBD3[n][m][k][j];
				
					printf("\n------checking-----: SP1 = %g + %gi\n",creal(SP1),cimag(SP1));
					printf("the output is PBD1[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD1[n][m][k][j]),cimag(PBD1[n][m][k][j]));
					
		  	  	}
				printf("the output is PBD1[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD1[n][m][k][j]),cimag(PBD1[n][m][k][j]));
				PB1[n][m][k] = Transa1;
				PB2[n][m][k] = Transb1;
				PB3[n][m][k] = Transc1;
				
				Transa2 += PB1[n][m][k];
				Transb2 += PB2[n][m][k];
				Transc2 += PB3[n][m][k];
				printf("PB1[%d][%d][%d] = %g + %g\n",n,m,k,creal(PB1[n][m][k]),cimag(PB1[n][m][k]));
				
			}
			
			P1[n][m] = Transa2;
			P2[n][m] = Transb2;
			P3[n][m] = Transc2;
			printf("P1[%d][%d] = %g + %g\n",n,m,creal(P1[n][m]),cimag(P1[n][m])); // m+1 just for show its index not zero but 1
			printf("P2[%d][%d] = %g + %g\n",n,m,creal(P2[n][m]),cimag(P2[n][m]));
			printf("P3[%d][%d] = %g + %g\n",n,m,creal(P3[n][m]),cimag(P3[n][m]));
			
			// pressure final
		     // OSNum = 1;
		    pF[n][m] = (P1[n][m]+P2[n][m]+P3[n][m])/(4*PI*Tint);
			printf("pF[%d][%d]= %g + %g\n",n,m,creal(pF[n][m]),cimag(pF[n][m]));// m+1 just for show its index not zero but 1
			
			
	  }
	  

  }


  
		 
for(int i=0;i<FNum;i++)
	{
		for(int j=0;j<1;j++)
			 {
			 	pXOYM[i][j] = pF[i][j];
				printf("pXOYM[%d][%d]= %g + %g\n",i,j,creal(pXOYM[i][j]),cimag(pXOYM[i][j]));
			 }
	}

  for(int i=0;i<FNum;i++)
  {
  	//for(int j=(HalfNum-(ObsrSthetaNum-1))+1;j<HalfNum;j++) // try to use with while loop
	for(int j=20;j<21;j++)
	{
	  //pXOYM[i][j] = pF[i][j];
		Op[i][j] = 2*sqrt(pow(creal(pXOYM[i][j]),2)+pow(cimag(pXOYM[i][j]),2));
	  //Op[i][j] = 2*fabs(pF[i][j]);
	  printf("Op[%d][%d] = %g\n",i,j,Op[i][j]);
	  
	}

  }
  
  int rr;
  OF = make_vector(OF,FNum);
  for(rr=0;rr<FNum;rr++)
  {
	  OF[rr] = (OmegaM+BNum*(rr-5)*OmegaR)/(2*PI);
	  printf("OF[%d] = %g\n",rr,OF[rr]);
  }
  
		  
		 
  //OF = (OmegaM+BNum*(-5:FNum-6)*OmegaR)/(2*PI);
  //

  
  rmdir('FD_Spectrum','s');
  mkdir('FD_Spectrum');
  
  FILE *fwrite1;
  fwrite1 = fopen("FDPressureSpectrum.txt","w");
  for(int k=20;k<21;k++)
  {
	  
	  for(int j=0;j<FNum;j++)
	  {
		  fprintf(fwrite1,"%12.6f\t %12.6f\n",OF[j],Op[j][k]);
	  }
	  
  }
  fclose(fwrite1);
  
  
  
  
  
  
  
  
  
  
  
  
  //------checking variable values-------
  printf("----- checking ----: fR: %g\n",fR);
  printf("PBD1[%d][%d][%d][%d] = %g + %gi\n",0,0,0,0,creal(PBD1[0][0][0][0]),cimag(PBD1[0][0][0][0]));
  printf("----- checking ----: XO = %g, YO = %g, ZO = %g\n", OX, OY, OZ);
  printf("----- checking ----: Theta: %g\n",Theta);
  printf("------checking-----: DataXR[3] = %g\n",DataXR[3]);
  printf("------checking-----: DataYR[32] = %g\n",DataYR[32]);
  printf("------checking-----: DataZR[3] = %g\n",DataZR[34]);
  printf("----- checking ----: DX = %g, DY = %g, DZ = %g, DR = %g\n",DX,DY,DZ,DR);
  printf("------checking-----: DORStar[3] = %g\n",DORStar[343]);
  printf("------checking-----: DOrX[43] = %g\n",DOrX[143]);
  printf("------checking-----: DOrY[43] = %g\n",DOrY[143]);
  printf("------checking-----: DORZ[43] = %g\n",DOrZ[143]);
  printf("------checking-----: DOr[43] = %g\n",DOr[143]);
  printf("------checking-----: DOR[43] = %g\n",DOR[143]);
  printf("------checking-----: DORStarX[143] = %g\n",DORStarX[143]);
  printf("------checking-----: DORStarY[243] = %g\n",DORStarY[243]);
  printf("------checking-----: DORStarZ[343] = %g\n",DORStarZ[343]);
  printf("------checking-----: DORX[143] = %g\n",DORX[41]);
  printf("------checking-----: DORY[243] = %g\n",DORY[43]);
  printf("------checking-----: DORZ[343] = %g\n",DORZ[243]);
  printf("------checking-----: RGamma[0] = %g\n",RGamma[0]);
  printf("------checking-----: RGamma[TNum-1] = %g\n",RGamma[TNum-1]);
  printf("------checking-----: Vx[65] = %g\n",Vx[65]);
  printf("------checking-----: Vy[243] = %g\n",Vy[243]);
  printf("------checking-----: Vz[65] = %g\n",Vz[35]);
  printf("------checking-----: Q[0] = %g\n",Q[0]);
  printf("------checking-----: Q[TNum-1] = %g\n",Q[TNum-1]);
  printf("------checking-----: OmegaM = %g\n",OmegaM);
  printf("------checking-----: Lx[25] = %g\n",Lx[25]);
  printf("------checking-----: Ly[63] = %g\n",Ly[63]);
  printf("------checking-----: Lz[44] = %g\n",Lz[44]);
  printf("------checking-----: FxM[255] = %g\n",FxM[255]);
  printf("------checking-----: FyM[165] = %g\n",FyM[165]);
  printf("------checking-----: FzM[165] = %g\n",FzM[165]);
  printf("------checking-----: FxP[5] = %g\n",FxP[5]); // here MaX, MaY, and MaZ are zero, is given from the file so the remaining result is zero
  printf("------checking-----: FyP[311] = %g\n",FyP[311]);
  printf("------checking-----: FzP[265] = %g\n",FzP[265]);
  printf("------checking-----: FRStarM[185] = %g\n",FRStarM[185]);
  printf("------checking-----: FRStarP[355] = %g\n",FRStarP[355]);
  printf("------checking-----: FRM[171] = %g\n",FRM[191]);
  printf("------checking-----: FRP[335] = %g\n",FRP[335]);

  printf("------checking-----: Gamma = %g\n",Gamma);
  printf("------checking-----: DT = %g\n",DT);
  printf("------checking-----: Tint = %g\n",Tint);
  printf("------checking-----: ODT = %g\n",ODT);
  printf("------checking-----: Q[0] = %g\n",Q[0]);
  printf("------checking-----: Q[TNum-1] = %g\n",Q[TNum-1]);
  
  printf("------checking-----: PBD1[34][1][1][1] = %g\n",PBD1[34][0][0][0]);
  printf("------checking-----: PBD1[155][1][1][1] = %g\n",PBD2[155][0][0][0]);
  printf("------checking-----: PBD1[238][1][1][1] = %g\n",PBD1[238][0][0][0]);
  
  printf("------checking-----: PB1[34][1][1] = %g\n",PB1[34][0][0]);
  printf("------checking-----: PB1[155][1][1] = %g\n",PB2[155][0][0]);
  printf("------checking-----: PB1[238][1][1] = %g\n",PB1[238][0][0]);
  
  printf("------checking-----: P1[34][1] = %g\n",P1[34][0]);
  printf("------checking-----: P1[155][1] = %g\n",P2[155][0]);
  printf("------checking-----: P1[238][1] = %g\n",P3[238][0]);
  printf("------checking-----: pF[78][1] = %g\n",pF[78][0]);

  
  printf("\n------------ Pressure Prediciton is DONE!----------------\n");
  
  
  

 
  
  /*--------------------freeing all used memories--------------------*/

  free(XO);
  free(YO);
  free(ZO);
  
  /* freeing the pointer */
  free(ObsrSMeshD_ptr);
  free_vector(OTime);

  
  free_vector(DataXR);
  free_vector(DataYR);
  free_vector(DataZR);
  free_vector(DOrX);
  free_vector(DOrY);
  free_vector(DOrZ);
  free_vector(DOr);
  free_vector(DORStar);
  free_vector(DOR);
  free_vector(DORStarX);
  free_vector(DORStarY);
  free_vector(DORStarZ);
  free_vector(DORX);
  free_vector(DORY);
  free_vector(DORZ);
  free_vector(RGamma);
  free_vector(Vx); 
  free_vector(Vy); 
  free_vector(Vz); 
  free_vector(Q); 
  free_vector(Lx);
  free_vector(Ly);
  free_vector(Lz);
  free_vector(FxM);
  free_vector(FyM); 
  free_vector(FzM); 
  free_vector(FxP); 
  free_vector(FyP); 
  free_vector(FzP);
  free_vector(FRStarM); 
  free_vector(FRStarP);
  free_vector(FRM); 
  free_vector(FRP);
  free_vector(OF);
  
  free_4Ddmatrix(PBD1,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD2,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD3,OSNum,BNum,DSNum);
  free_3Ddmatrix(PB1,OSNum,BNum);
  free_3Ddmatrix(PB2,OSNum,BNum);
  free_3Ddmatrix(PB3,OSNum,BNum); 
  free_dmatrix(P1,OSNum);
  free_dmatrix(P2,OSNum);
  free_dmatrix(P3,OSNum);
  free_dmatrix(pF,1);
  free_dmatrix(Op,ObsrSthetaNum-1);
  free_dmatrix(pXOYM,HalfNum);
  

  


  return 0;
}


