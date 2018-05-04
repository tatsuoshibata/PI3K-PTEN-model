#include "reactions.h"

X1=X0;
Y1=Y0;
Z1=Z0;
do{
  X0=X1;
  Y0=Y1;
  Z0=Z1;
  react[1]=PTENEnzymaticReaction(Dt,X0,Y0,Z0);
  react[2]=PI3KEnzymaticReaction(Dt,X0,Y0,Z0);
  react[4]=Dt*X0;
  react[5]=Dt*1.0;
  react[6]=PTENBinding(Dt,X0,Y0,Z0,Z0);
  react[7]=Dt*lambda*Z0;
  react[11]=Dt*PIP2degradation*Y0;
  
  X1=X0
    -react[1]
    +react[2]
    -react[4];
  
  Y1=Y0
    +react[1]
    -react[2]
    +react[5]
    -react[11];
  
  Z1=Z0
    +react[6]-react[7];
  //printf("%lf %lf\n",Y0,Y1);
 }while(fabs(Y1-Y0)>1e-16);

X0=X1;
Y0=Y1;
Z0=Z1;

printf("%lf %lf %lf\n",X0,Y0,Z0);
