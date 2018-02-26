#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/* ----- including -----header defined functions ------*/
#include "xmalloc.h"
#include "xmalloc.c"
#include "linspace.h"
#include "linspace.c"
#include "nextpower.h"
#include "nextpower_c.c"
#include "OTime.h"
#include "OTime.c"



/*------ constants definition------*/
#define PI 3.1416 

/* --------Function prototype----------*/
double **make_dmatrix(size_t m, size_t n);// construct 2D array type of double
void free_dmatrix(double **a, size_t n);// freeing it

double ***make_3Ddmatrix(size_t p, size_t q, size_t r);// construct 3D array type of double
void free_3Ddmatrix(double ***a, size_t q, size_t r); // freeing it

double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r); // construct 4D array type of double
void free_4Ddmatrix(double ****a, size_t p, size_t q, size_t r);




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
  double C_0 = FlowD[0].c_0;
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
  double *Time;
  Time = linSpace(TNum); /* initialize Time */
  for(int i=0;i<TNum;i++)
    {
      Time[i]=Tint*linSpace(TNum)[i]; /* update Time */
      printf("Time[%d] = %4.9f\n",i,(double)Time[i]);
    } 
  printf("-----checking------>Time[3409] = %4.4f\n",(double)Time[3409]);    /*checking Time */
     

  /*import the FNum.dat data */
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
  double *OTime;
  OTime = forOTime(NFFT-1); /* calling forOTime function and initializing it */
  for(int i=0;i<NFFT-1;i++)
    {
      OTime[i] = ODT*forOTime(NFFT-1)[i];
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
  double Gamma;
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
	Gamma = sqrt(1/(1-(pow(*MaX,2)+pow(*MaY,2)+pow(*MaZ,2))));
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
  double *OX,*OY,*OZ;
  make_vector(OX,OSNum);
  make_vector(OY,OSNum);
  make_vector(OZ,OSNum);
  double Theta;
  double *DX,*DY,*DZ,*DR;
  make_vector(DX,DSNum);
  make_vector(DY,DSNum);
  make_vector(DZ,DSNum);
  make_vector(DR,DSNum);
  double *DataXR,*DataYR,*DataZR;
  make_vector(DataXR,DSNum);
  make_vector(DataYR,DSNum);
  make_vector(DataZR,DSNum);
  

  
  for(int n =1; n<=FNum;n++)
  {
	  Omega[n]= OmegaM+BNum*(n-6)*OmegaR;
	  ka[n] = Omega[n]/C_0;
	  printf("Omega[%d]=%4.4g   ka[%d] = %4.4g\n",n,Omega[n],n,ka[n]);
	  
	  for(int m=1;m<=OSNum;m++)
	  {
		  OX[m] = XO[m];
		  OY[m] = YO[m];
		  OZ[m] = ZO[m];
		  
		  /*for(int k=1;k<=BNum;)
		  {
			  Theta = 2*(k-1)*PI/BNum;
			  
		  
		  		for(int j=1;j<=DSNum;j++)
		  	  	{
		  			//printf("Harmonic number %d, Observer number %d, Blade number %d, Source number %d\n",n,m,k,j);
					
					DX[j]=BladeMData[j].DataX;
					DY[j]=BladeMData[j].DataY;
					DZ[j]=BladeMData[j].DataZ;
					DR[j]=sqrt(pow(DX[j],2)+pow(DY[j],2));
					//printf("%g %g %g %g\n",DX[j],DY[j],DZ[j],DR[j]);
					
					//DVX, DVY, DVZ; about unitnormal
					
					double A = 1.0;//Area[j][1];
					DataXR[j] = DR[j]*cos(OmegaR*Time[j]+atan2(DY[j],DX[j])+Theta);
					DataYR[j] = DR[j]*sin(OmegaR*Time[j]+atan2(DY[j],DX[j])+Theta);
					
					//printf("%g\n",DataXR[j]);
					
		  	  	}
			}*/
	  }
	  
    
	 
  
 

  
  
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //checking for OX, OY, OZ
  printf("----- checking ----> [53]: XO[53] = %g, YO[53] = %g, ZO[53] = %g\n", OX[53], OY[53], OZ[53]);
  printf("----- checking ----> [99]: XO[99] = %g, YO[99] = %g, ZO[99] = %g\n", OX[99], OY[99], OZ[99]);
  printf("----- checking ----> [149]: XO[149] = %g, YO[149] = %g, ZO[149] = %g\n", OX[149], OY[149], OZ[149]);
  printf("----- checking ----: Theta: %g\n",Theta);
  
  
  /*--------------------freeing all used memories--------------------*/
  /* freeing the allocated memory for MachParameters.dat data */
  free(OTime);
  free(XO);free(YO);free(ZO);
  free(MaX);free(MaY);free(MaZ);
  
  /* freeing the pointer */
  free(ObsrSMeshD_ptr);
  
  free_vector(Omega);free_vector(ka); 
  free_vector(OX);free_vector(OY);free_vector(OZ);
  free_vector(DX);free_vector(DY);free_vector(DZ);free_vector(DR);
  free_vector(DataXR);free_vector(DataYR);free_vector(DataZR);
  free_4Ddmatrix(PBD1,OSNum,BNum,DSNum);free_4Ddmatrix(PBD2,OSNum,BNum,DSNum);free_4Ddmatrix(PBD3,OSNum,BNum,DSNum);
  free_3Ddmatrix(PB1,OSNum,BNum);free_3Ddmatrix(PB2,OSNum,BNum);free_3Ddmatrix(PB3,OSNum,BNum); 
  free_dmatrix(P1,OSNum);free_dmatrix(P2,OSNum);free_dmatrix(P3,OSNum);

  


  return 0;
}


/* ---------- function definitions --------- */
/* ---------- function definitions --------- */
double **make_dmatrix(size_t m, size_t n) /*2D array */
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

/* 3D matrix definition */
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

/* free_3Ddmatrix */
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

/* 4D array function definition */
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
}


