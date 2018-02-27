
#include "OTime.h"

double *forOTime(int s)
{
  double *OTime,x=0.0;
  int i;
  //OTime = (double*)calloc(s,sizeof(double));
  for(i=0;i<s;i++)
    {
		OTime = (double*)calloc(s,sizeof(double));
      	x = x+1;
      	OTime[i]=x-1;
    }
	

  return OTime;
}
