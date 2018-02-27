//#include <stdio.h>
#include <math.h>
#include "nextpower.h"

float nextPower2(int x){
  int P;
  P = log(x)/log(2);
  return P+1;
}
