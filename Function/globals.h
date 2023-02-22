/**************Globals**************/
CD OneI=0.0+1I;
#define PI 3.141592653589793238
 double h=0.01;
#define ohm0 1.0
int NDelay;
double **Inputs,***Outputs;
double **Inputs1,**Outputs1;
double *Preprocessing_Mask,***eDelay;
double PHIG;
double input_scaling;
double s,s1,s2,FE;
int step=34;

double T_clock,temp,tau,tau_Max;

double kap_min;
double kap;
double delta,*taubar;
int *delayNum, *delayIndex;
