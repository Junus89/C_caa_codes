/* This file calculating the FFT of a given set of complex
numbers assuming array size is power of 2 */

#include<stdio.h>
#include<math.h>
#include<complex.h>

double PI;
typedef double complex cplx;

//#define PI atan2(1,1)*4;



/* fftRB function definition: fft recursion and bit reverse definition */
void fftRB(cplx buf[], cplx out[], int n, int step)
{
	if(step<n) {
		fftRB(out, buf, n, step*2);
		fftRB(out+step, buf+step,n,step*2);
		
		for(int i=0;i<n;i += 2*step)
		{
			cplx t = cexp(-I*PI*i/n)*out[i+step];
			buf[i/2] = out[i]+t; /*fe the even part */
			buf[(i+n)/2] = out[i]-t; /*fo the odd part */
		}
	}
}

/* fft function definition */
void fft(cplx buf[], int n)
{
	cplx out[n];
	for(int i=0;i<n;i++) out[i] = buf[i];
	/*calling fftBR */
	fftBR(buf,out,n,1);
}

/* print function defition: print out the result */
void Print(const char *s, int N, cplx buf[])
{
	printf("%s",s);
	for (int i=0;i<N;i++)
	{
		if(!cimag(buf[i]))
			printf("%g\n",creal(buf[i]));
		else
			printf("(%g,%g)",creal(buf[i]),cimag(buf[i]));
	}
}

/*testing function */


int main()
{
	 PI = atan2(1,1)*4;
	int N=8;
	cplx buf[] = {1, 1, 1, 1, 0, 0, 0, 0};
	Print("Data: ",buf,N);
	fft(buf,N);
	Print("\nFFT: ",buf,N);
	
	return 0;
}

