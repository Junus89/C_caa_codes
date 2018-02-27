#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>


/* ----- including -----header defined functions ------*/
#include "xmalloc.h"
#include "xmalloc.c"
#include "linspace.h"
#include "linspace.c"
#include "nextpower.h"
#include "nextpower_c.c"
#include "OTime.h"
#include "OTime.c"
#include "main.h"
#include "main.c"
//#include "FifthLoop.h"
//#include "FifthLoop.c"






/* structure for storing flow data */
typedef struct flowData{
  double rho_0;
  int c_0;
  double RMaTip;
}flowData;

/* structure for storing geometry data */
typedef struct geomData{
  int BNum;         // Blade Number
  double R;         // Radius of propeller blade [m]
}geomData;

/*structure for storing intTnum data */
typedef struct intTnumData{
  int TNum;         // ...
}intTnumData;

/*structure for storing FNum data */
typedef struct FNumData{
  int FNum;
}FNumData;

/*structure for storing WriteDataSMeshBlades.dat data  */
typedef struct WDataSMB{
  double DataX;
  double DataY;
  double DataZ;
  double nxData;
  double nyData;
  double nzData;
}WDataSMB;





int main(){
  // opening and reading flowdata file
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
  
  // opening and reading geometry data
  FILE *fp_g;
  fp_g = fopen("GeometryData.txt","r");
  //storing geometry data
  geomData GeomD[]={0};
  for(int i=0; i<sizeof(GeomD);i++)
    {
      fscanf(fp_g,"%d %lf",&GeomD[i].BNum, &GeomD[i].R);
    }
  printf("----------GeometryData------------\n\n");
  printf("Bnum = %d, R = %f\n\n",GeomD[0].BNum, GeomD[0].R);
  double BNum = GeomD[0].BNum;
  double R = GeomD[0].R;
  fclose(fp_g);


  double OmegaR = RMaTip*C_0/R;
  double fR = OmegaR/(2*PI);         // rotation frequency
  double TR = 1/fR;
  double OmegaM = 150*2*PI;         
  double fM = 150;                  //pulsation frequency
  printf("---------Matematical relations ------------\n\n");
  printf("OmegaR = %4.4f [rad/s], fR = %4.4f [hz], TR = %4.4f [s], OmegaM = %4.4f [rad/s]\n\n",OmegaR, fR, TR, OmegaM);
  
  /* import the intTnum.dat data */
  FILE *fp_Tnum;
  fp_Tnum = fopen("intTNum.txt","r");
  //storing intTnum data
  intTnumData TNumD[]={0};
  for(int i = 0; i<sizeof(TNumD);i++)
    {
      fscanf(fp_Tnum,"%d",&TNumD[i].TNum);
    }
  int TNum = TNumD[0].TNum;
  printf("------intTNumData-------\n\n");
  printf("TNum = %d\n",TNum);
  fclose(fp_Tnum);
  /* Mathematical relations of TNum data */
  double Tint = TR;
  double DT = Tint/(TNum-1);          //...discrete time or sampling time? or sampling time step?
  double *Time = NULL;
  Time = linSpace(TNum); 
  for(int i=0;i<TNum;i++)
    {
      Time[i]=Tint*linSpace(TNum)[i]; 
      printf("Time[%d] = %4.9f\n",i,(double)Time[i]);
    } 
  printf("-----checking------>Time[3409] = %4.4f\n",(double)Time[3409]);   
    


  FILE *fp_Fnum;
  fp_Fnum = fopen("FNum.txt","r");
  //storing FNum data
  FNumData FNumD[]={0};
  for(int i=0;i<sizeof(FNumD);i++)
    {
      fscanf(fp_Fnum,"%d",&FNumD[i].FNum);
    }
  int FNum = FNumD[0].FNum;
  printf("---------FNumData------\n\n");
  printf("FNum = %d\n",FNum); /* checking for FNUM */

  
  int NFFT = nextPower2(FNum);
  printf("the next power of FNum->NFFT = %d\n",NFFT);
  double ODT = Tint/NFFT;
  printf("ODT = %lf\n",ODT);
  /* for construction of OTime */
  double *OTime=NULL;
  OTime = forOTime(NFFT-1); /* calling forOTime function and initializing it */
  for(int i=0;i<NFFT-1;i++)
    {
		OTime[i] = ODT*OTime[i];
      //OTime[i] = ODT*forOTime(NFFT-1)[i];
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
  


  /* reading 'MachParameter.dat' data */
  double MaX_temp;
  double MaY_temp;
  double MaZ_temp;
  //double Gamma;
  double *MaX = NULL;
  double *MaY = NULL;
  double *MaZ = NULL;
  int idx, r=1, k=0;/*r=i,k=j */
  
  FILE *fp_MMa;
  fp_MMa = fopen("MachParameter.txt","r");
  
  while(fscanf(fp_MMa,"%lf %lf %lf", &MaX_temp, &MaY_temp,&MaZ_temp)!=EOF)
    {
      if(MaX==NULL&&MaY==NULL&&MaZ==NULL)
	{
	  MaX = malloc(sizeof(MaX_temp));
	  MaY = malloc(sizeof(MaY_temp));
	  MaZ = malloc(sizeof(MaZ_temp));
	  *MaX = MaX_temp;
	  *MaY = MaY_temp;
	  *MaZ = MaZ_temp;
	  //  printf(" MaX = %lf\n MaY = %lf\n MaZ = %lf\n",MaX,MaY,MaZ); 
	}
      else
	{
	  r++;
	  MaX = realloc(MaX,sizeof(MaX)*r);
	  MaY = realloc(MaY,sizeof(MaY)*r);
	  MaZ = realloc(MaZ,sizeof(MaZ)*r);
	  idx = r-1;
	  *(MaX+idx) = MaX_temp;
	  *(MaX+idx) = MaY_temp;
	  *(MaX+idx) = MaZ_temp;
	  //  printf("MaX = %lf\n MaY = %lf\n MaZ = %lf\n",MMa[i],MMa[i],MMa[i]); 
	}
    }
  MaX_temp = 0.0;
  MaY_temp = 0.0;
  MaZ_temp = 0.0;

  printf("----checking for the result-----\n");
  while(idx>=0)
    {
      printf(" MaX = %lf\n MaY = %lf\n MaZ = %lf\n",MaX[k],MaY[k],MaZ[k]);
      idx--;
      k++;
    }
	Gamma = sqrt(1/(1-(pow(*MaX,2)+pow(*MaY,2)+pow(*MaZ,2)))); // check this part, it seems strange!
	printf("Gamma = %lf\n",Gamma);
  fclose(fp_MMa);
  
  /* freeing the allocated memory for MachParameters.dat data */
  
  
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

  /* freeing the pointer */
  
  
    
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
  double ****PBD1,****PBD2,****PBD3;
  PBD1 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD2 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD3 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  
  double ***PB1,***PB2,***PB3;
  PB1 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB2 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB3 = make_3Ddmatrix(FNum,OSNum,BNum);
  
  double **P1,**P2,**P3;
  P1 = make_dmatrix(FNum,OSNum);
  P2 = make_dmatrix(FNum,OSNum);
  P3 = make_dmatrix(FNum,OSNum);
  
  
  double *Omega,*ka;
  make_vector(Omega,FNum);
  make_vector(ka,FNum);
  double OX,OY,OZ;
  double Theta;
  double DX,DY,DZ,DR;
  /*double *DataXR,*DataYR,*DataZR; */
  make_vector(DataXR,TNum);
  make_vector(DataYR,TNum);
  make_vector(DataZR,TNum);
  //double *DOrX, *DOrY, *DOrZ, *DOr;
  
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
  
  

  
  for(int n =0; n< FNum;n++)
  {
	  Omega[n]= OmegaM+BNum*(n-6)*OmegaR;
	  ka[n] = Omega[n]/C_0;
	  //printf("Omega[%d]=%4.4g   ka[%d] = %4.4g\n",n,Omega[n],n,ka[n]);
	  
	  for(int m=0;m<1;m++)//OSNum=
	  {
		  OX = XO[0];//the first value of XO vector
		  OY = YO[0];// the first value of YO vector
		  OZ = ZO[0];// the first value of ZO vector
		  
		  for(int k=0;k<BNum;k++)
		  {
			  Theta = 2*((k+1)-1)*PI/BNum; // Theta = 2*(k-1)*PI/BNum;
			  
		  
		  		for(int j=0;j<DSNum;j++)
		  	  	{
		  			printf("Harmonic number %d, Observer number %d, Blade number %d, Source number %d\n",n,m+1,k+1,j+1); // since m, k, j starts from zero, add it with 1
					
					DX=BladeMData[0].DataX;
					DY=BladeMData[0].DataY;
					DZ=BladeMData[0].DataZ;
					DR=sqrt(pow(DX,2)+pow(DY,2));
					//printf("DX=%g, DY=%g, DZ= %g, DR = %g\n",DX,DY,DZ,DR);
					
					/*--------DVX, DVY, DVZ; about unitnormal *----------*/
					
					double A = 1.0;//Area[j][1];

					
					
		  	  	}
			}
	  }
	  
    
	 
  
 

  
  
  }
  
  fifthLoop(3600, 1, OmegaR, Gamma, C_0, MaX, MaY, MaZ, DX, DY, DZ, DR, OX, OY ,OZ, Theta, DataXR, DataYR, DataZR, DOrX, DOrY, DOrZ, DOr, DORStar, DOR, DORStarX,\
		 DORStarY, DORStarZ, DORX, DORY, DORZ, RGamma);
  
  
  //checking for OX, OY, OZ
  printf("----- checking ----: XO = %g, YO = %g, ZO = %g\n", OX, OY, OZ);
  printf("----- checking ----: Theta: %g\n",Theta);
  printf("------checking-----: DataXR[3] = %g\n",DataXR[3]);
  printf("------checking-----: DataYR[32] = %g\n",DataYR[32]);
  printf("------checking-----: DataZR[3334] = %g\n",DataZR[3334]);
  printf("----- checking ----: DX = %g, DY = %g, DZ = %g, DR = %g\n",DX,DY,DZ,DR);
  printf("------checking-----: DORSar[3] = %g\n",DORStar[343]);
  printf("------checking-----: DOR[543] = %g\n",DOR[543]);
  printf("------checking-----: DORStarX[1543] = %g\n",DORStarX[1543]);
  printf("------checking-----: DORStarY[2543] = %g\n",DORStarY[2543]);
  printf("------checking-----: DORStarZ[3543] = %g\n",DORStarZ[3543]);
  printf("------checking-----: DORX[1543] = %g\n",DORX[41]);
  printf("------checking-----: DORY[2543] = %g\n",DORY[543]);
  printf("------checking-----: DORZ[3543] = %g\n",DORZ[2343]);
  printf("------checking-----: RGamma[333] = %g\n",RGamma[333]);
  printf("------checking-----: RGamma[333] = %g\n",RGamma[3333]);
  
  //printf("DataXR[3433] = %g, DataYR = %g, DataZR = %g\n",DataXR[3444],DataYR[3444],DataZR[3444]);
  
  /*--------------------freeing all used memories--------------------*/
  /* freeing the allocated memory for MachParameters.dat data */
  //free(Time);
  free(OTime);
  free(XO);
  free(YO);
  free(ZO);
  free(MaX);
  free(MaY);
  free(MaZ);
  
  /* freeing the pointer */
  free(ObsrSMeshD_ptr);
  
  free_vector(Omega);
  free_vector(ka); 
  
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
  
  free_4Ddmatrix(PBD1,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD2,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD3,OSNum,BNum,DSNum);
  free_3Ddmatrix(PB1,OSNum,BNum);
  free_3Ddmatrix(PB2,OSNum,BNum);
  free_3Ddmatrix(PB3,OSNum,BNum); 
  free_dmatrix(P1,OSNum);
  free_dmatrix(P2,OSNum);
  free_dmatrix(P3,OSNum);

  


  return 0;
}


