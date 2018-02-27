#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "fft2.c"

#define DEBUGFFT 0
#define PI acos(-1.0)

void realfft(double *data, unsigned long n, int isign)
{
	unsigned long i,j;
	double theta, wpr, wpi, wr, wi, wtemp, c1, c2, h1i,h1r,h2i,h2r;
	
	#if DEBUGFFT
	for (i=0;i<n;i+=2){
		printf("\n realfft: (%f, %f)",data[i],data[i+1]);
	}
	#endif
	theta = PI/((double)n/2.);
	c1 = 0.5;
	if (isign == 1)
	{
		c2 = -0.5;
		fft(data,n,1);
		
	}else {
		c2=0.5;
		theta = -theta;
	}
	#if DEBUGFFT
	printf("\n\n");
	for(i=0;i<n;i+=2)
	{
		printf("\n realfft after fft: (%f, %f)",data[i],data[i+1]);
	}
	#endif
	wpr = -2*pow(sin(0.5*theta),2.);
	wpi=sin(theta);
	wr = 1.+wpr;
	wi = wpi;
	#if DEBUGFFT
	printf("\n wr=%f wi = %f\n",wr,wi);
	#endif
	for(i=2;i<n/2-1;i+=2)
	{
		/*devide----*/
		j=n-1;
		#if DEBUGFFT
		printf("\n i1=%d 02 = %d i3=%d i4 = %d",i,i+1,j,j+1);
		#endif
		
		h1r = c1*(data[i]+data[j]);
		h1i = c1*(data[i+1]-data[j+1]);
		h2r = -c2*(data[i+1]+data[j+1]);
		h2i = c2*(data[i]-data[j]);
		
		/* and conquer ---*/
		
		data[i] = h1r+wr*h2r-wi*h2i;
		data[i+1] = h1i+wr*h2i+wi*h2r;
		data[j] = h1r-wr*h2r+wi*h2i;
		data[j+1]=  -h1i+wr*h2i+wi*h2r;
		#if DEBUGFFT
			printf("\n h1r = %f h1i = %f h2r = %f h2i = %f",h1r,h1i,h2r,h2i);
			printf("\n data[%d] = %f data[%d] = %f data[%d] = %f data[%d] = %f",i,data[i],i+1,data[i+1],j,data[j],j+1,data[j+1]);
		#endif
			
			/*trigonometric recursion */
			
			wtemp = wr;
			wr = wr*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1)
	{
		h1r = data[0];
		data[0] = h1r+data[1];
		data[1] = h1r-data[1];
		
	}
	else {
		h1r = data[0];
		data[0] = c1*(h1r+data[1]);
		data[1] = c1*(h1r-data[1]);
		fft(data,n,-1);
	}
	#if DEBUGFFT
	printf("\n\n");
	for(i=0;i<n;i+=2)
	{
		printf("\n realfft: (%f %f)",data[i],data[i+1]);
	}
	#endif
	return;
	
}

int main()
{
	int N=8;
	double test[]={1,1,1,1,0,0,0,0};
	double *tes;
	tes = &test[0];
	printf("\n\n----Original Data---\n\n");
	
    for (int i = 0; i < N; i++)

        printf("[%g] ", tes[i]);
	
	printf("\n\n----fft below----\n\n");
	realfft(tes,N,1);
	
    for (int i = 0; i < N; i++)

        printf("[%g] ", tes[i]);
	
	printf("\n\n----ifft below----\n\n");
	realfft(tes,N,-1);
	
    for (int i = 0; i < N; i++)

        printf("[%g] ", tes[i]);
	
	
	
 	/*
    fft(test,8,1);
    for (int i = 0; i < N; i++)

        printf("[%g ]", tes[i]);
	printf("\n\n----ifft below----\n\n");
    fft(tes,8,-1);
	
    for (int i = 0; i < N; i++)

        printf("[%g ]", tes[i]); */
    
	
	return 0;
}


