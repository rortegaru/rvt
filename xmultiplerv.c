#include <stdio.h>
#include <stdlib.h>

extern void  multiplerv_(float *, float *);
extern float  ran0(long *);

int main(){
float dist,magnitude,magn;
long seed;
int i,imagnitude;
seed=986782;
dist=40.0;
for(i=0;i<1000;i++){
magnitude=ran0(&seed)*700; 
imagnitude=(int) magnitude;
magn=(float) imagnitude / 100;
magnitude=5.21;
printf("magnitude = %f \n",magn);
multiplerv_(&magn,&dist); 
}

} 
