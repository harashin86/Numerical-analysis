#include <stdio.h>
#include <math.h>
#define N 3
#define EPSILON 1.0E-5 /* 行列の次数 */
float a[N][N]={{ 5.0,2.0,1.0},{1.0,4.0,2.0},{2.0,1.0, 4.0}};
float b[N]={17.0,16.0,11.0};
int k;/* 繰り返し回数*/
float x[N];
void jacobi(void);
int main(void){
  printf("\n\n  ====係数行列と定数=====\n\n");
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      printf("%10.4f",a[i][j] );
      printf("  %10.4f\n",b[i] );
    }
  }
  jacobi();
  printf("\n\n  ===== ヤコビ法による連立一次程式の解　=====\n\n");
  printf("繰り返し回数＝ %dn\n\n",k );
  for(int i=0; i<N; i++){
    printf("x[%2d] = %10.4f \n",i+1,x[i] );
  }

  return 0;
}
void jacobi(void){
  float Nk,xx,xdif;
  float x0[N];
  FILE *outfile;
  if((outfile = fopen("result.csv", "w")) == NULL){
    fprintf(stderr, "Can't open file \n");
    exit(2);
  }
  k=0;
  for(int i=0; i<N; i++)x[i]=0.0;
  do {
    k+=1;
    Nk=0.0;
    for(int i = 0; i < N; i++)x0[N]=x[i];
    for(int j=0;j<N;j++){
      xx=b[j];
      for(int i=0;i<N;i++){
        if(i!=j){
          xx=xx-a[j][i]*x0[i];
        }
      }
      x[j]=xx/a[j][j];
      xdif=x[j]-x0[j];
      if(xdif<0.0)xdif=-xdif;
      Nk+=xdif;
    }
    fprintf(outfile, "%d, %e\n", k, Nk);
  } while(Nk>EPSILON);
  fclose(outfile);
}
