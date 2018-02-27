#include <stdio.h>
#include <stdlib.h>
#include "linspace.h"


double *linSpace(int TNum )
{
	
  double Delta = 1.0/TNum;
  double *newArr=NULL, x;
  
  int i;
  /*
  newArr = (double*)malloc(TNum*sizeof(double));
  for(i = 0; i <TNum; i++)
    {
      x += Delta;
      newArr[i]= x;
	  //Time[i] = Tint*newArr[i];
      //printf("linspace[]=%4.4f\n",(double) linspace[i]);
    }*/
	newArr = (double*)realloc(newArr,TNum*sizeof(double));
	for(i = 0;i<TNum;i++)
	{
		x += Delta;
		newArr[i] = x;
	}

  return newArr;
}

