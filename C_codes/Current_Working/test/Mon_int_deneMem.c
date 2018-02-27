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
double Theta, OX,OY,OZ, DX,DY,DZ,DR;
double *Omega,*ka,*DataXR,*DataYR,*DataZR;
double ****PBD1,****PBD2,****PBD3;



/* ----- function prototyping----- */

double **make_dmatrix(size_t m, size_t n);// construct 2D array type of double
void free_dmatrix(double **a);// freeing it

double ***make_3Ddmatrix(size_t p, size_t q, size_t r);// construct 3D array type of double
void free_3Ddmatrix(double ****a, size_t p, size_t q); // freeing it

//double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r); // construct 4D array type of double
//void free_4Ddmatrix(double ****a);


double ****make_4ddmatrix(size_t o, size_t p, size_t q, size_t r);
void free_4Ddmatrix(double ****a, size_t p, size_t q, size_t r);


double *Ones(int Tnum);

/* structure for storing flow data */
typedef struct flowData{
  float rho_0;
  int c_0;
  float RMaTip;
}flowData;

/* structure for storing geometry data */
typedef struct geomData{
  int BNum;         // Blade Number
  float R;         // Radius of propeller blade [m]
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
  float DataX;
  float DataY;
  float DataZ;
  float nxData;
  float nyData;
  float nzData;
}WDataSMB;





int main(){
  // opening and reading flowdata file
  FILE *fp_f;
  fp_f = fopen("flowdata.txt","r");
  //storing flowdata 
  flowData FlowD[]={0};
  for(int i =0;i<sizeof(FlowD);i++)
    {
      fscanf(fp_f,"%f %d %f",&FlowD[i].rho_0,&FlowD[i].c_0,&FlowD[i].RMaTip);
      //    printf(" Rho_0 = %4.4f, c_0 = %d, RMaTip = %4.4f\n",(double)FlowD[i].rho_0,FlowD[i].c_0,(double)FlowD[i].RMaTip);
    }
  printf("--------------FlowData------------\n\n");
  printf(" Rho_0 = %4.4f, c_0 = %d, RMaTip = %4.4f\n\n",(double)FlowD[0].rho_0,FlowD[0].c_0,(double)FlowD[0].RMaTip);
  float Rho_0 = FlowD[0].rho_0;
  float C_0 = FlowD[0].c_0;
  float RMaTip = FlowD[0].RMaTip;
  fclose(fp_f);
  
  // opening and reading geometry data
  FILE *fp_g;
  fp_g = fopen("GeometryData.txt","r");
  //storing geometry data
  geomData GeomD[]={0};
  for(int i=0; i<sizeof(GeomD);i++)
    {
      fscanf(fp_g,"%d %f",&GeomD[i].BNum, &GeomD[i].R);
    }
  printf("----------GeometryData------------\n\n");
  printf("Bnum = %d, R = 4.4%f\n\n",GeomD[0].BNum, GeomD[0].R);
  float BNum = GeomD[0].BNum;
  float R = GeomD[0].R;
  fclose(fp_g);


  float OmegaR = RMaTip*C_0/R;
  float fR = OmegaR/(2*PI);         // rotation frequency
  float TR = 1/fR;
  float OmegaM = 150*2*PI;         
  float fM = 150;                  //pulsation frequency
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
  printf("-----checking------>Time[3409] = %g\n",(double)Time[3409]);    /*checking Time */
     

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
  float ODT = Tint/NFFT;
  printf("ODT = %f\n",ODT);
  /* for construction of OTime */
  double *OTime;
  OTime = forOTime(NFFT-1); /* calling forOTime function and initializing it */
  for(int i=0;i<NFFT-1;i++)
    {
      OTime[i] = ODT*forOTime(NFFT-1)[i];
      printf("OTime[%d] = %g\n",i,OTime[i]);
    }
  printf("----checking ----> OTime[5] = %4.4f\n",OTime[5]); /*checking for OTime */
  double DF = ((1/ODT)/2)*(1/(1.0*FNum/2)); /* here as FNum is int type, should multiply with 1.0 to get the float type DF result */
  printf("----checking ----> DF = %g\n",DF);
  
  /* reading WriteDataSMeshBlades.dat data */
  FILE *fp_meshB;
  fp_meshB = fopen("WriteDataSMeshBlade.txt","r");
  //storing the data into the structure
  WDataSMB BladeMData[]={0};
  for(int i=0;i<sizeof(BladeMData);i++)
    {
      fscanf(fp_meshB,"%f %f %f %f %f %f",&BladeMData[i].DataX, &BladeMData[i].DataY, &BladeMData[i].DataZ, &BladeMData[i].nxData,&BladeMData[i].nyData,&BladeMData[i].nzData);

    }
  
  printf("----checking ---> DataX = %f, DataY= %f, DataZ = %f\n, nxData = %f\n, nyData = %f\n, nzData = %f\n",BladeMData[0].DataX, BladeMData[0].DataY, BladeMData[0].DataZ,BladeMData[0].nxData,BladeMData[0].nyData,BladeMData[0].nzData);
  int DSNum = 1;
  fclose(fp_meshB);


  /* reading WriteObserverSMesh.dat data  */
  float XO_temp;
  float YO_temp;
  float ZO_temp;
  float *XO = NULL;
  float *YO = NULL;
  float *ZO = NULL;
  int index, OSNum = 0, j = 0, i = 1;

  FILE *fp_ObMesh;
  fp_ObMesh = fopen("WriteObserverSMesh.txt","r");
  //storing the data into the structure
  while(fscanf(fp_ObMesh,"%f %f %f",&XO_temp, &YO_temp, &ZO_temp)!=EOF)
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
      printf("[%d]: XO = %f, YO = %f, ZO = %f\n",j, XO[j], YO[j], ZO[j]);
      index--;
      OSNum = 1+j;
      j++;
      XO_temp++;
      YO_temp++;
      ZO_temp++;
    }
  printf("----- checking ----> OSNum = %d\n",OSNum);
  printf("----- checking ----> [53]: XO[53] = %f, YO[53] = %f, ZO[53] = %f\n", XO[53], YO[53], ZO[53]);
  printf("----- checking ----> [99]: XO[99] = %f, YO[99] = %f, ZO[99] = %f\n", XO[99], YO[99], ZO[99]);
  printf("----- checking ----> [149]: XO[149] = %f, YO[1499] = %f, ZO[149] = %f\n", XO[149], YO[149], ZO[149]);
  fclose(fp_ObMesh);
  /* freeing the allocated memory for 'WritesObserverSMeshData.dat */
  free(XO);free(YO);free(ZO);


  /* reading 'MachParameter.dat' data */
  float MaX_temp;
  float MaY_temp;
  float MaZ_temp;
  float Gamma;
  float *MaX = NULL;
  float *MaY = NULL;
  float *MaZ = NULL;
  int idx, r=1, k=0;/*r=i,k=j */
  
  FILE *fp_MMa;
  fp_MMa = fopen("MachParameter.txt","r");
  
  while(fscanf(fp_MMa,"%f %f %f", &MaX_temp, &MaY_temp,&MaZ_temp)!=EOF)
    {
      if(MaX==NULL&&MaY==NULL&&MaZ==NULL)
	{
	  MaX = malloc(sizeof(MaX_temp));
	  MaY = malloc(sizeof(MaY_temp));
	  MaZ = malloc(sizeof(MaZ_temp));
	  *MaX = MaX_temp;
	  *MaY = MaY_temp;
	  *MaZ = MaZ_temp;
	  //  printf(" MaX = %f\n MaY = %f\n MaZ = %f\n",MaX,MaY,MaZ); 
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
	  //  printf("MaX = %f\n MaY = %f\n MaZ = %f\n",MMa[i],MMa[i],MMa[i]); 
	}
    }
  MaX_temp = 0.0;
  MaY_temp = 0.0;
  MaZ_temp = 0.0;

  printf("----checking for the result-----\n");
  while(idx>=0)
    {
      printf(" MaX = %f\n MaY = %f\n MaZ = %f\n",MaX[k],MaY[k],MaZ[k]);
      idx--;
      k++;
    }
	Gamma = sqrt(1/(1-(pow(*MaX,2)+pow(*MaY,2)+pow(*MaZ,2))));
	printf("Gamma = %f\n",Gamma);
  fclose(fp_MMa);
  
  /* freeing the allocated memory for MachParameters.dat data */
  free(MaX);free(MaY);free(MaZ);
  
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
  printf("\n-------->Checking:\n ObsrSthetaNum = %d\n ObsrSFeinum = %d\n HalfNum = %f\n",ObsrSthetaNum, ObsrSFaiNum, HalfNum);

  /* freeing the pointer */
  free(ObsrSMeshD_ptr);
  
    
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
  DataSVector = UnitNormal; */
  
 // --------------------------------- 
  double ****PBD1,****PBD2,****PBD3;
  
  /*PBD1= make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD2= make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD3= make_4Ddmatrix(FNum,OSNum,BNum,DSNum);*/
  
  double ***PB1,***PB2,***PB3;
  PB1 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB2 =make_3Ddmatrix(FNum,OSNum,BNum);
  PB3 =make_3Ddmatrix(FNum,OSNum,BNum);
  
  double **P1,**P2,**P3;
  P1 = make_dmatrix(FNum,OSNum);
  P2 = make_dmatrix(FNum,OSNum);
  P3 = make_dmatrix(FNum,OSNum);
  
  
  PBD1 = make_4ddmatrix(FNum,1,BNum,DSNum);

  for(int n=0;n<FNum;n++)
    {
		Omega = make_vector(Omega,FNum);
		ka = make_vector(ka,FNum);
  	  	Omega[n]= OmegaM+BNum*(n-6)*OmegaR;
  	  	ka[n] = Omega[n]/C_0;
		//printf("Omega[%d]=%4.4g   ka[%d] = %4.4g\n",n,Omega[n],n,ka[n]);
      
	  for(int m=0;m<1;m++)
        {
  		  OX = XO[m];//the first value of XO vector
  		  OY = YO[m];// the first value of YO vector
  		  OZ = ZO[m];// the first value of ZO vector
		  
		  for(int k=0;k<BNum;k++)
		  {
		  	Theta = 2*((k+1)-1)*PI/BNum; // Theta = 2*(k-1)*PI/BNum;
			
			for(int l=0;l<DSNum;l++)
			{
	  			printf("Harmonic number %d, Observer number %d, Blade number %d, Source number %d\n",n,j+1,k+1,l+1); // since m, k, j starts from zero, add it with 1
				
				DX=BladeMData[l].DataX;
				DY=BladeMData[l].DataY;
				DZ=BladeMData[l].DataZ;
				DR=sqrt(pow(DX,2)+pow(DY,2));
				
				
				/*--------DVX, DVY, DVZ; about unitnormal *----------*/
				
				double A = 1.0;//Area[j][1];
				
				DataXR = make_vector(DataXR,TNum);
				DataYR = make_vector(DataYR,TNum);
				DataZR = make_vector(DataZR,TNum);
				
				
			    for(int i=0;i<100;i++)
			      {

					  DataXR[i] = DR*cos(OmegaR*Tint*linSpace(TNum)[i]+atan2(DY,DX)+Theta);
			          DataYR[i] = DR*sin(OmegaR*Tint*linSpace(TNum)[i]+atan2(DY,DX)+Theta);
					  DataZR[i] = DZ*Ones(TNum)[i];
				  }
				  
				PBD1[n][m][k][l]=1.0/(n+m+k+l+1.0);
				  
				PB1[n][m][k] += PBD1[n][m][k][l];
				printf("PB1[%d][%d][%d] = %g\n",n,m,k,PB1[n][m][k]);
				
			}
		  }
          /*for(int k=0;k<1;i++)
            {
              PBD1 = make_4ddmatrix(FNum,1,BNum,DSNum);
              for(int i=0;i<FNum;i++)
                for(int j=0;j<1;j++)
                  for(int k=0;k<BNum;k++)
                    for(int l=0;l<DSNum;l++)
                      {
                        PBD1[i][j][k][l]=1.0/(1+i+j+k+l);
                        printf("PBD1[%d][%d][%d][%d]=%g\n",i,j+1,k+1,l+1,PBD1[i][j][k][l]);
                      }
            } */
        }

    }
  //printf("-------> checking:H[2][2][2][4]=%g\n",PBD1[24][232][32][41]);
  //printf("-------> checking:H[3][3][3][4]=%g\n",PBD1[33][53][63][214]);
  
  
  
  
  
  
  	printf("----- checking ----: XO = %g, YO = %g, ZO = %g\n", OX, OY, OZ);
	printf("----- checking ----: Theta: %g\n",Theta);
	printf("----- checking ----: DX = %g, DY = %g, DZ = %g, DR = %g\n",DX,DY,DZ,DR);
    printf("------checking-----: DataXR[3] = %g\n",DataXR[3]);
    printf("------checking-----: DataYR[32] = %g\n",DataYR[32]);
    printf("------checking-----: DataZR[3334] = %g\n",DataZR[3334]);
	//free(Omega);
	//free(ka);
	free_vector(Omega);
	free_vector(ka); 
	free_vector(DataXR);
	free_vector(DataYR);
	free_vector(DataZR);
	
  	free_4Ddmatrix(PBD1,OSNum,BNum,DSNum);
  	free_3Ddmatrix(PB1,OSNum,BNum);
  
  
  

























  return 0;
}

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

void free_dmatrix(double **a)
{
  if(a!=NULL){
    for(size_t i=0;a[i] !=NULL;i++)
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
    for(size_t j=0;j<p;j++)
      make_vector(a[i][j],r);
  a[p]=NULL;
  return a;
}

/* free_3Ddmatrix 
void free_3Ddmatrix(double ***a)
{
  if(a!=NULL){
    for(size_t i=0;a[i]!=NULL;i++)
      for(size_t j=0;a[i][j]!=NULL;j++)
	free_vector(a[i][j]);
    for(size_t i=0;a[i]!=NULL;i++)
      free_vector(a[i]);
    free(a);
  }
}
*/


void free_3Ddmatrix(double ****a, size_t p, size_t q)
{
  for(int i=0;i<q;i++)
    for(int j=0;j<q;j++)
        free(a[i][j]);
  for(int i=0;i<p;i++)
    free(a[i]);
  free(a);
}




// 4D deneme
double ****make_4ddmatrix(size_t o, size_t p, size_t q, size_t r)
{
  double ****a;
  make_vector(a,o+1); // make the pointer vector                                                                        
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



//function definition ofr ones 
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



/* 4D array function definition 
double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r)
{
  double ****a;
  make_vector(a,o+1);//make the pointer vector
  for(size_t i=0;i<o;i++)
    make_vector(a[i],p);
  for(size_t i=0;i<o;i++)
    for(size_t j=0;j<o;j++)
      make_vector(a[i][j],q);
  for(size_t i=0;i<o;i++)
    for(size_t j=0;j<o;j++)
      for(size_t k=0;k<o;k++)
	make_vector(a[i][j][k],r);
  a[o]=NULL;
  return a;
}
//freeing 4Ddmatrix 
void free_4Ddmatrix(double ****a)
{
  if(a!=NULL){
    for(size_t i=0;a[i]!=NULL;i++)
      for(size_t j=0;a[i][j]!=NULL;j++)
	for(size_t k=0;a[i][j][k]!=NULL;k++)
	  free_vector(a[i][j][k]);
    for(size_t i=0;a[i]!=NULL;i++)
      for(size_t j=0;a[i][j]!=NULL;j++)
	free_vector(a[i][j]);
    for(size_t i=0;a[i]!=NULL;i++)
      free_vector(a[i]);
    free(a);
  }
}
*/
   
