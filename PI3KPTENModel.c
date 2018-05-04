#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include "cell.h"

#define PTENEnzymaticReaction(delta,x,y,z) ((delta)*v1*(z)*(x)/(K1+(x)))
#define PI3KEnzymaticReaction(delta,x,y,z) ((delta)*(v2)*(y)/(K2+(y)))
#define PLCEnzymaticReaction(delta,x,y,z)  ((delta)*(kdash-1.0)*(y)/(K3+(y)))
#define PTENBinding(delta,x,y,z,ztot) ((delta)*v0*(1-(ztot))*((y)/(K5+(y)))*(K4+alpha*(x))/(K4+(x)))
#define PIP2DegradationReaction(delta,x,y,z) ((delta)*(PIP2degradation)*(y))

double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
void init_by_array(unsigned long init_key[], int key_length);
void init_genrand(unsigned long s);

double kappa=0.0;
double v1=0.0;
double K1=0.0;
double v2=0.0;
double alpha=0.0;
double v0=0.0;
double K4=0.0;
double D=0.0;
double D_PTEN=0.0;
double K2=0.0;
double K3=0.0;
double printTime=0.0;
double K5=0.0;
double lambda=0.0;
double PIP2degradation=0.001;
double perturbTime=-1.0;

double sizeOfSystem=2.0*M_PI*5.;
double Dx;
double Dt=.005;

char outputDir[256]="out";
char outputFile[256];
char pip2File[256];
char pip3File[256];
char ptenFile[256];

#define DELIMITER   "= *"
#define numOfReact 12

int paramcmp(char *str0, double *x, char *paramname)
{
  char *token;
  char str[80];
  double y;

  strcpy(str,str0);
  token = strtok(str,DELIMITER );
  if(strcmp(str,paramname)==0){
    //    printf("%s %s\n",str,paramname);
    token = strtok(NULL, DELIMITER);
    //    printf("%s\n",token);
    if(token!=NULL){
      y=atof(token);
      if(*x!=y)
      {
        *x=y;
        printf("%s=%lf\n",paramname,*x);
        return 1;
      }
    }
  }
  return 0;
}

void write_param(double t, double *x, char *paramname)
{
  FILE *FP_HIST;
  static double time=-10.0;

  if((FP_HIST=fopen(outputFile,"a"))==NULL) {
    printf("cannot open param_hist file \n");
    exit(1);
  }

  if(t!=time){
    time=t;
    fprintf(FP_HIST,"#time=%lf\n",t);
  }
  printf("%s=%lf\n",paramname,*x);

  fprintf(FP_HIST,"%s=%lf\n",paramname,*x);

  if(fclose(FP_HIST)!=0){
    printf("cannot close param_hist file \n");
  }
}

int read_param(char *parameterFile)
{
  int i;
  char str[80];
  char changed=0;
  FILE *FP_PARAM;

  if((FP_PARAM=fopen(parameterFile,"r")) == NULL){
    printf("cannot read parameter file \"%s\"\n",parameterFile);
    exit(1);
  }

  while((fscanf(FP_PARAM,"%s\n",str))!=EOF) {
    changed+=paramcmp(str,&kappa,"kappa");
    changed+=paramcmp(str,&v1,"v1");
    changed+=paramcmp(str,&K1,"K1");
    changed+=paramcmp(str,&v2,"v2");
    changed+=paramcmp(str,&alpha,"alpha");
    changed+=paramcmp(str,&v0,"v0");
    //    changed+=paramcmp(str,&V_dash,"V_dash");
    changed+=paramcmp(str,&K4,"K4");
    changed+=paramcmp(str,&D,"D");
    changed+=paramcmp(str,&D_PTEN,"D_PTEN");
    changed+=paramcmp(str,&K2,"K2");
    //    changed+=paramcmp(str,&kdash,"kdash");
    changed+=paramcmp(str,&PIP2degradation,"PIP2degradation");
    changed+=paramcmp(str,&K3,"K3");
    changed+=paramcmp(str,&printTime,"printTime");
    changed+=paramcmp(str,&K5,"K5");
    changed+=paramcmp(str,&lambda,"lambda");
    changed+=paramcmp(str,&tau,"tau");
    changed+=paramcmp(str,&perturbTime,"perturbTime");
  }

  fclose(FP_PARAM);

  return changed;
}


void outputToFile(double x[numGrids+2], double y[numGrids+2], double z[numGrids+2])
{
  int i;

  double datax[numGrids];
  double datay[numGrids];
  double dataz[numGrids];

  FILE *FP_PIP3;
  FILE *FP_PIP2;
  FILE *FP_PTEN;

  if((FP_PIP3=fopen(pip3File,"a")) == NULL
  || (FP_PIP2=fopen(pip2File,"a")) == NULL
  || (FP_PTEN=fopen(ptenFile,"a")) == NULL) {
    printf("open error\n");
    exit(1);
  }

  for(i=1;i<=numGrids;i++){
    datax[i-1]=x[i];
    datay[i-1]=y[i];
    dataz[i-1]=z[i];
  }
  fwrite(datax,sizeof(double),numGrids,FP_PIP3);
  fwrite(datay,sizeof(double),numGrids,FP_PIP2);
  fwrite(dataz,sizeof(double),numGrids,FP_PTEN);

  fclose(FP_PIP3);
  fclose(FP_PIP2);
  fclose(FP_PTEN);
}

void
gc(char *parameterFile)
{
  int i,j;
  double x0[numGrids+2];
  double x1[numGrids+2];

  double y0[numGrids+2];
  double y1[numGrids+2];

  double z0[numGrids+2];
  double z1[numGrids+2];

  //  double X0=0.4,Y0=120.0,Z0=0.5;
  double X0=0.4,Y0=120.0,Z0=0.5;

  double diff_x[numGrids+2];
  double diff_y[numGrids+2];
  double diff_z[numGrids+2];

  double Ztot;
  double time;

  double temp;

  double react[numOfReact];

  /* init */
  tau=1.0;

  read_param(parameterFile);

  time=0;

  write_param(time,&kappa,"kappa");
  write_param(time,&v1,"v1");
  write_param(time,&K1,"K1");
  write_param(time,&v2,"v2");
  write_param(time,&alpha,"alpha");
  write_param(time,&v0,"v0");
  write_param(time,&K4,"K4");
  write_param(time,&D,"D");
  write_param(time,&D_PTEN,"D_PTEN");
  write_param(time,&K2,"K2");
  write_param(time,&PIP2degradation,"PIP2degradation");
  write_param(time,&K3,"K3");
  write_param(time,&printTime,"printTime");
  write_param(time,&K5,"K5");
  write_param(time,&lambda,"lambda");
  write_param(time,&tau,"tau");
  write_param(time,&perturbTime,"perturbTime");

  sizeOfSystem=2.0*M_PI*5.0;
  Dx=(sizeOfSystem/(double)numGrids);
  Dt/=tau;
  printf("sizeOfSystem=%lf\n",sizeOfSystem);

  for(i=0;i<=numGrids+1;i++){
    x0[i]=X0+.1*genrand_real1();
    y0[i]=Y0+.1*genrand_real1();
    z0[i]=Z0+.1*genrand_real1();
  }

  j=0;


  while(time<1000.0){
    if(j%(int)(1./Dt)==0){//every 1 second
      outputToFile(x0,y0,z0);
      read_param(parameterFile);
    }

    j++;
    time=(double)j*Dt;


    Ztot=0.0;
    for(i=1;i<=numGrids;i++){
      Ztot+=z0[i];
    }
    Ztot=Ztot/(double)numGrids;

    for(i=1;i<=numGrids;i++){
      diff_x[i]=(x0[i+1]+x0[i-1]-2.0*x0[i])/(Dx*Dx);
      diff_y[i]=(y0[i+1]+y0[i-1]-2.0*y0[i])/(Dx*Dx);
      diff_z[i]=(z0[i+1]+z0[i-1]-2.0*z0[i])/(Dx*Dx);
    }

    for(i=1;i<=numGrids;i++){
      react[1]=PTENEnzymaticReaction(Dt,x0[i],y0[i],z0[i]);
      react[2]=PI3KEnzymaticReaction(Dt,x0[i],y0[i],z0[i]);
      //  react[3]=PLCEnzymaticReaction(Dt,x0[i],y0[i],z0[i]);
      react[4]=Dt*x0[i];
      react[5]=Dt*1.0;
      react[6]=PTENBinding(Dt,x0[i],y0[i],z0[i],Ztot);
      //react[6]=Dt*v0*(1-Z0-Ztot)*y0[i]/(K5+y0[i])*(K4+alpha*x0[i])/(K4+x0[i]);
      react[7]=Dt*lambda*z0[i];
      react[11]=Dt*PIP2degradation*y0[i];

      x1[i]=x0[i]+Dt*diff_x[i] -react[1] +react[2] -react[4];
      y1[i]=y0[i]+Dt*diff_y[i] +react[1] -react[2] +react[5] -react[11];
      z1[i]=z0[i]+Dt*D_PTEN*diff_z[i] +react[6]-react[7];
    }

    x1[0]=x1[numGrids];
    x1[numGrids+1]=x1[1];

    y1[0]=y1[numGrids];
    y1[numGrids+1]=y1[1];

    z1[0]=z1[numGrids];
    z1[numGrids+1]=z1[1];

    for(i=0;i<=numGrids+1;i++){
      x0[i]=x1[i];
      y0[i]=y1[i];
      z0[i]=z1[i];
    }

    // if(time==perturbTime){
    //   for(i=numGrids/2-50;i<=numGrids/2+50;i++){
    //     x0[i+50]+=5.0;
    //     z0[i]+=3.0+2.0*genrand_real1();
    //   }
    //
    // }
    //
    if(time==perturbTime){
      for(i=numGrids/2-50;i<=numGrids/2+50;i++){
        x0[i]+=5.0;
            z0[i+100]+=3.0+2.0*genrand_real1();
        // z0[i]+=0.0;
      }
    }
  }




}

int main(int argc, char **argv)
{
  int c;
  int digit_optind = 0;
  time_t         *tp = (time_t *) 0;
  time_t         timer;
  char parameterFile[256];


  sprintf(parameterFile, "param.dat");

  while (1)
  {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] =
    {
      {"file", 1, 0, 'f'},
      {"out", 1, 0, 'o'},
      {0, 0, 0, 0}
    };

    c = getopt_long_only (argc, argv, "",
    long_options, &option_index);

    if (c == -1)
    break;
    switch (c)
    {
      case 'f':
      //printf("%s\n",optarg);
      strcpy(parameterFile,optarg);
      break;
      case 'o':
      //printf("%s\n",optarg);
      strcpy(outputDir,optarg);

      strcpy(outputFile,outputDir);
      strcat(outputFile,"/param_hist.txt");

      strcpy(pip2File,outputDir);
      strcat(pip2File,"/pip2.dat");

      strcpy(pip3File,outputDir);
      strcat(pip3File,"/pip3.dat");

      strcpy(ptenFile,outputDir);
      strcat(ptenFile,"/pten.dat");

      if (mkdir(outputDir,
        S_IRUSR | S_IWUSR | S_IXUSR |         /* rwx */
        S_IRGRP | S_IWGRP | S_IXGRP |         /* rwx */
        S_IROTH | S_IXOTH | S_IXOTH) == 0){
          /* rwx */
          printf("make %s\n", outputDir);
        }
        else {
          printf("cannnot make %s\n",outputDir);
          perror("");
        }
        break;
        default:
        printf ("?? getopt returned character code 0%o ??\n", c);
      }
    }

    //printf("%s\n",parameterFile);
    init_genrand((unsigned long)time(&timer));

    gc(parameterFile);
  }
