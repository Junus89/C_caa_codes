#include <stdio.h>
#include <string.h>


int main(){
  FILE *fp;
  char fdName[128], sdName[128], odName[128];
  int BNum, TNum, FNum;
  double R, MaX, MaY, MaZ;

  fp = fopen("ConstantV.txt","r");
  while(1 == fscanf(fp,"%s%*[^\n] %s%*[^\n] %s%*[^\n] %d%*[^\n] %lf%*[^\n] %d%*[^\n] %d%*[^\n] %lf%*[^\n] %lf%*[^\n] %lf%*[^\n]\n",fdName, sdName, odName,\
  	 &BNum,&R,&TNum,&FNum,&MaX,&MaY,&MaZ)){
  }
  //printf(" %s\n %s\n %s\n %d\n %4.4f\n %d\n %d\n %4.4f\n %4.4f\n %4.4f\n",fdName,sdName,odName,BNum,R,TNum,FNum,MaX,MaY,MaZ);
  //int cc = BNum+TNum;
  //printf("cc = BNum + TNum =%d\n",cc);
  //printf("the flow data name is: %s\n",fdName);
  fclose(fp);
  return 0;
}
