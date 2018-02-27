#include <stdio.h>
#include "xmalloc.h"

void *malloc_or_exit(size_t(nbytes), const char *file, int line)
{
  void *x;
  if((x = malloc(nbytes))==NULL)
    {
      fprintf(stderr,"%s:line %d: malloc() of %zu bytes failed\n",file,line,nbytes);
      exit(EXIT_FAILURE);
    }else
        return x;
}

/** int main(void){
  int n=13;
  double *x;
  x =xmalloc(n*sizeof(*x));
  //x = malloc_or_exit(n*sizeof(*x),__FILE__, __LINE__);
  for(int i=0;i<n;i++)
    {
      x[i] = 1.0/(i+1);
      printf("x[%d] = %g\n",i,x[i]);
    }
  printf("it uses memory of %lu bytes\n",n*sizeof(*x));
  free(x);
  return 0;
  } **/
