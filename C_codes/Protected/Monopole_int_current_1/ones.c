#include<stdio.h>
#include<stdlib.h>
#include "xmalloc.h"
#include "xmalloc.c"
#include "math.h"

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

int main()
{
  
  //  int TNum = 12;
  double *Bir;
  //make_vector(One,TNum);
  Bir = Ones(12);
  for(int i=0;i<12;i++)
    {
      //Bir[i] = 1.0*i;
      printf("Bir[1][%d] = %g\n",i,Bir[i]);
    }
  
  return 0;
}
