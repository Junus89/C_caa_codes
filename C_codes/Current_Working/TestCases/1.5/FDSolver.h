#ifndef FDSolver_H_
#include<complex.h> /* Standard Library of Complex Numbers */
#define FDSolver_H_
#define PI 3.1416 


//double Tint,TR,fR,OmegaR;
double DT;
int DSNum, BNum, FNum, TNum;
double R;
double OX,OY,OZ, MaX, MaY, MaZ;
double Theta;
double DX,DY,DZ,DR;
double Time[194401];

double complex ****PBD1,****PBD2,****PBD3;
double complex ***PB1,***PB2,***PB3;
double complex **P1,**P2,**P3;

double complex **pF, **pXOYM;
double  **Op, *OF;

double A = 1.0;
double *OTime;
double OmegaM, OmegaR, Gamma, C_0,  DX, DY, DZ,  DR, OX, OY,\
	 OZ, Theta, *DataXR, *DataYR, *DataZR, *DOrX, *DOrY,  *DOrZ, *DOr, *DORStar, *DOR,\
		 *DORStarX, *DORStarY, *DORStarZ, *DORX, *DORY,*DORZ, *RGamma, *Vx, *Vy, *Vz, *Q,\
			 *Lx, *Ly, *Lz, *FxM, *FyM, *FzM, *FxP, *FyP, *FzP;
double *FRStarM, *FRStarP, *FRM, *FRP;
double complex SP1;
double complex SP2;
double complex SP3;
double complex z1 = 0.0 + 1.0*I;
double complex Inicomp = 0.0 + 0.0*I;



double complex try;





// structure definition

/* structure for storing flow data */
typedef struct flowData{
  double rho_0;
  int c_0;
  double RMaTip;
}flowData;


/*structure for storing WriteDataSMeshBlades.dat data  */
typedef struct WDataSMB{
  double DataX;
  double DataY;
  double DataZ;
  double nxData;
  double nyData;
  double nzData;
}WDataSMB;



double **make_dmatrix(size_t m, size_t n);// construct 2D array type of double
void free_dmatrix(double **a, size_t n);// freeing it

double ***make_3Ddmatrix(size_t p, size_t q, size_t r);// construct 3D array type of double
void free_3Ddmatrix(double ***a, size_t q, size_t r); // freeing it

double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r); // construct 4D array type of double
void free_4Ddmatrix(double ****a, size_t p, size_t q, size_t r);


// ----fifthloop prepreation-----



//void fifthLoop();
void fifthLoop(double OmegaR, double Omega, double MaX, double MaY, double MaZ, double ka, double **pF);
//void fifthLoop(double *DataXR, double *DataYR, double *DataZR, double *DOrX, double *DOrY, double *DOrZ, double *DOr, double *DORStar, double *DOR, double *DORStarX,\
		 double *DORStarY, double *DORStarZ, double *DORX, double *DORY, double *DORZ, double *RGamma);


double *Ones(int n);// Ones function


#endif


