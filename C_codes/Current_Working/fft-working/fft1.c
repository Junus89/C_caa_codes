#include <stdio.h>
#include <math.h>
#include <complex.h>

typedef double complex cplx;
 
#define PI atan2(1, 1) * 4

 
void fftRB(cplx buf[], cplx out[], int n, int step, int isign)
{
  if (step < n) {
    fftRB(out, buf, n, step * 2,isign);
   	fftRB(out + step, buf + step, n, step * 2,isign);
 
    for (int i = 0; i < n; i += 2 * step) {
      cplx t = cexp(-I * PI * i / (n*isign)) * out[i + step];
      buf[i / 2]     = out[i] + t;
      buf[(i + n)/2] = out[i] - t;
    }
  }
}

 
void fft(cplx buf[], int n)
{
  cplx out[n];
  for (int i = 0; i < n; i++) out[i] = buf[i];
 
  fftRB(buf, out, n, 1,1);
}


 
/*void ifft(cplx buf[], int n)
{
  cplx out[n];
  for (int i = 0; i < n; i++) out[i] = buf[i];
  fftRB(buf, out, n, 1,-1);
  //fft(buf,n);
  
}*/

void ifft(cplx buf[],int n)
{
	cplx out[n];
	for(int i=0;i<n;i++) out[i] = buf[i];
	out[n]=conj(buf[n]);
	fft(buf,n);
	out[n]=conj(buf[n]);
	out[n]/=n;
}
 
void Print(const char * s, cplx buf[],int N) {
  printf("%s", s);
  for (int i = 0; i < N; i++)
    if (!cimag(buf[i]))
      printf("%g, 0.000i ", creal(buf[i]));
    else
      printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}
 
int main()
{

  cplx test[] = {1, 1, 1, 1, 0, 0, 0, 0}; 
  //cplx test[] = {1, 1, 1, 1,2,2,2,2, 0, 0, 0, 0,1,2,3,4};
  int N = 8;
 
  Print("Data : ", test,N);
  fft(test, N);
  Print("\nFFT : ", test,N);
  ifft(test,N);
  Print("\niFFT :", test,N);
  
 
  return 0;
}
