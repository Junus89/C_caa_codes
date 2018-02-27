#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define DEBUGFFT 0
#define PI acos(-1.0)

/* function prototypes */
void fft(double *x, unsigned long length, int isign);
int bit_reverse(int val, int nbits);
void swap(double *x, double *y);

/*
int main()
{
	int N=8;
	double test[]={1,1,1,1,0,0,0,0};
	double *tes;
	tes = &test[0];
    for (int i = 0; i < N; i++)

        printf("[%g] ", tes[i]);
	printf("\n\n----fft below----\n\n");
	
	
 
    fft(test,8,1);
    for (int i = 0; i < N; i++)

        printf("[%g ]", tes[i]);
	printf("\n\n----ifft below----\n\n");
    fft(tes,8,-1);
	
    for (int i = 0; i < N; i++)

        printf("[%g ]", tes[i]);
    
	
	return 0;
}

*/


void fft(double *x, unsigned long length, int isign)
{
	unsigned long i, j, nbits, mmax, istep, m,n;
	double theta, wpr, wpi, wr, wi, wtemp, tempi, tempr;
	
    #if DEBUGFFT
	for (i=0;i<length;i+=2){
		printf("(%f, %f)\n",x[i],x[i+1]);
	}
    #endif
	/* CALCULATES THE BITLENGTH OF THE DATA */
	n = length/2+1;
	nbits = (int)(log((double)(n))/log(2.0)); 
	
	/* Check if the length of the data is power of 2 else exit */
	
	if (2*pow(2.,nbits)!=length){
		printf("FFT ERROR: input data length must be a power of 2.\n");
		exit(1);
	}
	
	/* callling the bit reverse function to get the order of bit reversal permutation */
	for (i=1;i<n-1;i++)
	{
		if ((j=bit_reverse(i,nbits))>i){ /* calling the bit reverse function */
			swap(&x[2*i],&x[2*j]); //even
			swap(&x[2*i+1],&x[2*j+1]); // odd
			
		}
		#if DEBUGFFT
		printf("i=%d j=%d\n",i,j);
		#endif
	}
	
	/* Danielson-Lanczos section */
	
	mmax = 2;
	while(length>mmax){
		istep = 2*mmax;
		theta = 2*PI/(isign*mmax);
		wpr = -2.*pow(sin(0.5*theta),2.);
		wpi = sin(theta);
		wr = 1.;
		wi = 0.;
		
		for (m=1;m<=mmax;m+=2)
		{
			for (i=m-1;i<(length-1);i+=istep){
				
				j = i+mmax;
				tempr = wr*x[j]-wi*x[j+1];
				tempi = wr*x[j+1]+wi*x[j];
				x[j] = x[i]-tempr;
				x[j+1] = x[i+1]-tempi;
				x[i] = x[i]+tempr;
				x[i+1] = x[i+1]+tempi;
			}
			/* Trigonometric Recurrence */
			wtemp = wr;
			wr = wr*wpr-wi*wpi*wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
		
	}
	
	#if DEBUGFFT
	for(i=0;i<length;i+=2)
	{
		printf("\n (%f, %f)",x[i],x[i+1]);
	}
	#endif
	return;
}

int bit_reverse(int val, int nbits)
{
	int mask_low, mask_high;	/* masks for bits to interchange */
	int bit_low, bit_high;		/*values of the bits in question */
	int i;
	mask_low =1;
	mask_high = 1 << (nbits-1);
	for(i=0;i<nbits/2;i++)
	{
		bit_low = (val & mask_low)>>i;
		bit_high = (val & mask_high) >> (nbits-i-1);
		mask_low <<=1;
		mask_high >>=1;
		if(bit_low != bit_high){
			val ^= (1<< i)^(1<<(nbits-i-1));
		}
	}
	return (val);
	
}

void swap(double *x, double *y)
{
	double temp;
	temp =*x;
	*x = *y;
	*y=temp;
}

/*void Print(const char * s, buf[],int N) {
  printf("%s", s);
  for (int i = 0; i < N; i++)

      printf("(%g, %g) ", buf[i]);
}*/

