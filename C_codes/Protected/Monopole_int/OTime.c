
#include "OTime.h"

double *forOTime(int s)
{
  double *OTime=NULL,x=0.0;
  int i;
  //OTime = (double*)calloc(s,sizeof(double));
  for(i=0;i<s;i++)
    {
		OTime = (double*)calloc(s,sizeof(double));
      	x = x+1;
      	OTime[i]=x-1;
    }
	
	/*double *OTime, x=0.0;
	int i, j=1, idx;
	while(j<s)
	{
		if(OTime==NULL)
		{
			OTime = (double*)calloc(s,sizeof(double));
	      	x = x+1;
	      	*OTime=x-1;
		}
		else
		{
			j++;
			OTime = realloc(OTime, sizeof(OTime)*j);
			idx = j-1;
			*(OTime+idx) = x-1;
			
		}
	}*/
  return OTime;
}
