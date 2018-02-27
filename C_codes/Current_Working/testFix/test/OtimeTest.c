#include <stdio.h>
#include <stdlib.h>


#include "OTime.h"
#include "xmalloc.h"
#include "xmalloc.c"

int main(){
	
	
	double *Otime;
	int n = 12;
	Otime = make_vector(Otime,n);

	for(int k=0;k<n;k++)
	{
		Otime[k] = k;
		printf("Otime[%d]  = %g\n",k,Otime[k]);
	}
	printf("size of Otime = %lu\n",sizeof(Otime));
	
	FILE *fp;
	fp = fopen("OtimeTest.txt","w");
	for(int k=0;k<sizeof(Otime);k++)
	{
		fprintf(fp,"%12.4f\n",Otime[k]);
		
	}
	fclose(fp);
	return 0;
}

