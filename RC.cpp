#include<stdio.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include <stdlib.h>
// #include <fftw3.h>
#include<complex>
#include <random>
#include <armadillo>

using namespace std;
using namespace arma;

typedef complex<double> CD;
#define Discard 490
#define NTrain 5000
#define NTest 2500
#define NThetas 200
#define Nsamples 4 //always is 1
#define Ninputs 2
#define MaxDelayPerLayer 500
int NElements=Discard+NTrain+NTest;

#include "Function/definition.h"
#include "Function/Distributions.h"
#include "Function/globals.h"
#include "Function/preprocessing.h"
#include "Function/Helpers.h"
#include "Function/Tasks.h"
#include "Function/FMG.h"
#include "Function/Weights.h"
/***********************************Delay initializing*************************************/

void initializeDelay(double a[]){
    Definitions obj1;
    int dn=tau_Max/h;
    obj1.def(eDelay,Nsamples,MaxDelayPerLayer,dn);
    for(int i=0;i<Nsamples;i++){
        for(int j=0;j<NDelay;j++){
            delayIndex[j] = 0;
            for(int k = 0; k < dn; k++){
                eDelay[i][j][k] = a[i];

            }

        }

    }

}

void Ikeda(double x[],double y[], double b[],double by[],double V[]){
    int DI;
    double kf=kap/NDelay;

    /**************First Layer************/
    for(int i=0;i<Nsamples;i++){

        double outloop=input_scaling*V[i];

        double var;
        FE=0.0;

        for(int j=0;j<NDelay;j++){
            DI=delayIndex[j];
            var=sin(eDelay[i][j][DI]+PHIG+outloop);
            FE+=var*var;
        }
        FE*=kf;
        b[i] = -x[i]-delta*y[i]+FE;
        by[i] = x[i];

    }

}


void Euler(double a[],double ay[], double h, double t,double V[]){
    double x[Nsamples], b[Nsamples],y[Nsamples], by[Nsamples];
    for(int j = 0; j < Nsamples; j++){
        x[j] = a[j];y[j] = ay[j];
    }
    Ikeda(x,y, b,by,V);

    for(int i = 0; i < Nsamples; i++){
        a[i] += h * b[i];
        ay[i] += h * by[i];
    }
    for(int i = 0; i < Nsamples; i++){
        for(int j=0;j<NDelay;j++){
            int DI=delayIndex[j];
            eDelay[i][j][DI] = a[i];


        }

    }
    for(int j=0;j<NDelay;j++){
        int DI=delayIndex[j];
        delayIndex[j] = (DI + 1) % delayNum[j];
    }
}
void rungeKutta(double a[],double ay[], double h, double t,double V[]){
double x[Nsamples], b[4][Nsamples],y[Nsamples], by[4][Nsamples];
for(int i = 0; i < 4; i++){
    for(int j = 0; j < Nsamples; j++){
        if(i == 0){ x[j] = a[j];y[j] = ay[j];}
        if(i == 1){ x[j] = a[j] + h * b[0][j] / 2.0;y[j] = ay[j] + h * by[0][j] / 2.0;}
        if(i == 2){ x[j] = a[j] + h * b[1][j] / 2.0;y[j] = ay[j] + h * by[1][j] / 2.0;}
        if(i == 3){ x[j] = a[j] + h * b[2][j];y[j] = ay[j] + h * by[2][j];}
    }
    Ikeda(x,y, b[i],by[i],V);
}

for(int i = 0; i < Nsamples; i++){
   a[i] += h * (b[0][i] + 2.0 * b[1][i]+ 2.0 * b[2][i] + b[3][i]) / 6.0;
   ay[i] += h * (by[0][i] + 2.0 * by[1][i]+ 2.0 * by[2][i] + by[3][i]) / 6.0;
}
for(int j=0;j<NDelay;j++){
        int DI=delayIndex[j];
        delayIndex[j] = (DI + 1) % delayNum[j];
    }
    
}
int main(){
    ofstream outputz1("Data/Samples/sample1T4.txt");
    ofstream outputz2("Data/Samples/sample2T4.txt");
    ofstream outputz3("Data/Samples/sample3T4.txt");
    ofstream outputz4("Data/Samples/sample4T4.txt");
    ofstream outputz5("Data/Samples/sample5T4.txt");
    ofstream outputz6("Data/Samples/sample6T4.txt");
    ofstream outtrain("Data/train4.txt");
    ofstream outtest("Data/test4.txt");
    Definitions obj1;
    double a[Nsamples],ay[Nsamples],t;

    obj1.def(Inputs,NTrain+Discard+NTest,Nsamples);
    obj1.def(Outputs,NTrain+NTest,Ninputs+1,Nsamples);
    obj1.def(Preprocessing_Mask,NThetas);
    obj1.def(Inputs1,NTrain+Discard+NTest,Nsamples);
    obj1.def(Outputs1,NTrain+NTest,Ninputs+1);
    double theta;
    ofstream out("test");
    int NDS=210,NDS_gamma=1;
//       string Case="Spatial";
       string Case="temporal";
    
    if(Case=="Spatial"){NDS=30;NDS_gamma=50;}
    s=1.0;
    s1=sqrt(s),s2=sqrt(1.0-s);
    string datastyle="Narma10";
    choose_task(datastyle);
    /************************First loop************************/
    for(int nth_gamma=0;nth_gamma<NDS_gamma;nth_gamma++){

        PHIG=2.0*PI/5.0;
        NDelay=1;
        cout<<"NDE="<<NDelay<<endl;
        
        
        CreateSamples(Preprocessing_Mask);
        /************************Tasks************************/
        

        T_clock=40.0+nth_gamma*0.0;
        theta =T_clock*1.0/NThetas;
        
//      tau=(nth_gamma*0.0*theta)+40.4;
        if(Case=="Spatial"){input_scaling=0.02+(nth_gamma)*0.02;temp=79.8;}
        if(Case=="temporal"){input_scaling=0.22;temp=(nth_gamma*1.0*theta)+36.0;}
        
        cout<<"input_scaling="<<input_scaling<<endl;
        for(int nth_kap=0;nth_kap<NDS;nth_kap++){
              
              

            if(Case=="Spatial"){kap_min=0.2+(nth_kap*1.0)*0.05;tau=40.4;}
            if(Case=="temporal"){kap_min=1.45;tau=(nth_kap*1.0*theta)+76.0;}
            //            kap_min=2.5+(nth_kap*1.0)*0.1;
            cout<<"/**********************"<<nth_gamma+1<<"/"<<NDS_gamma<<"**********************/"<<endl;
            cout<<"/**********************"<<nth_kap+1<<"/"<<NDS<<"**********************/"<<endl;
            cout<<"s="<<s<<endl;
            cout<<"T_clock="<<T_clock<<endl;

            cout<<"theta="<<theta<<endl;
               h=theta/4;
            int NTrainh=getMaxindex(NTrain*T_clock);
            int NTesth=getMaxindex(NTest*T_clock);
            int Discardh=getMaxindex(Discard*T_clock);
            int NElementsh=getMaxindex(NElements*T_clock);
            int determin=NElementsh/NElements;
            int kmask=0;
            int indexmax=theta/h;
            int ntheta=NThetas+1;

            delaytimes(NDelay);
            calcParameter(h);
            for(int j=0;j<Nsamples;j++){
                a[j] = 0.01+j*0.01;
                ay[j] = kap_min/delta*sin(PHIG)*sin(PHIG)+0.01+j*0.01;
                if(nth_gamma==0 && nth_kap==0 && j==Nsamples-1)cout<<"a["<<j+1<<"]="<<a[j]<<'\t'<<"ay["<<j+1<<"]="<<ay[j]<<endl;
            }
            initializeDelay(a);

            int cycle=NTrainh/NTrain;
             int thetah=cycle/NThetas;
             cout<<"beta="<<kap_min<<'\t'<<"temp="<<temp<<'\t'<<"ohm0="<<ohm0<<endl;
             cout<<"thetah="<<thetah<<'\t'<<theta/h<<'\t'<<"Ntheta="<<NThetas<<endl;


            for(int i = 0; i < 200*cycle; i++){
                double external[Nsamples];
                for(int j=0;j<Nsamples;j++){
                    external[j]=0.0;
                }
                Euler(a,ay, h, t,external);

            }
            /*********************Discard***************************/
            // calculation for transient
            int add=1;
            int inpindex=0;
            kmask=0;
            for(int i = 0; i < Discardh; i++){
                double external[Nsamples];
                for(int j=0;j<Nsamples;j++){
                    external[j]=Inputs[kmask][j]*Preprocessing_Mask[inpindex];

                }
                Euler(a,ay, h, t,external);
                if((i+add)%thetah==0){inpindex=(inpindex+1)%NThetas;}
                if((i+add)%cycle==0) {kmask++;}

            }
            if(nth_gamma==0 && nth_kap==0)cout<<"Discard is finished"<<endl;
            /***********************Training Data*****************************/

             double ***Xtrain_aug,**new_Xtrain_aug,***Xtest_aug,**new_Xtest_aug;
             obj1.def(Xtrain_aug,Nsamples,NTrain,NThetas);
             obj1.def(new_Xtrain_aug,Nsamples*NTrain,ntheta);
             int train_index=0;
             cout<<"inpindex="<<inpindex<<'\t'<<cycle<<endl;
             for(int i = 0; i < NTrainh; i++){
                 double external[Nsamples];
                 for(int j=0;j<Nsamples;j++){
                     external[j]=Inputs[kmask][j]*Preprocessing_Mask[inpindex];
                 }
                 Euler(a,ay, h, t,external);
                 if((i+add)%thetah==0){for(int j=0;j<Nsamples;j++){Xtrain_aug[j][train_index][inpindex]=a[j];}inpindex=(inpindex+1)%NThetas;}
                 if((i+add)%cycle==0) {kmask++;train_index++;}
             }
             /*******************weights***************************/
             double CCC=1e-3;
             double **weights,**Ytrain_aug,**_Ybar_train,**Ytest_aug,**_Ybar_test;
             obj1.def(weights,Ninputs+1,ntheta);
             obj1.def(Ytrain_aug,Ninputs+1,Nsamples*NTrain);
             obj1.def(_Ybar_train,Ninputs+1,Nsamples*NTrain);
             obj1.def(_Ybar_test,Ninputs+1,Nsamples*NTest);
             obj1.def(Ytest_aug,Ninputs+1,Nsamples*NTest);
             for(int i=0;i<NTrain;i++){
                 for(int j=0;j<Nsamples;j++){
                     for(int k=0;k<Ninputs+1;k++){
                         Ytrain_aug[k][j*NTrain+i]=Outputs[i][k][j];

                     }

                 }

             }
            new_x(Xtrain_aug,new_Xtrain_aug,NTrain,CCC,ntheta,Nsamples);
            Normal_weights_Brunner(new_Xtrain_aug,Ytrain_aug,weights,ntheta,Nsamples);

             for(int i=0;i<NTrain*Nsamples;i++){
                for(int k=0;k<Ninputs+1;k++){
                    _Ybar_train[k][i]=0;
                    for(int j=0;j<ntheta;j++){
                        _Ybar_train[k][i]+=new_Xtrain_aug[i][j]*weights[k][j];
                    }
                    if(nth_gamma==0 && nth_kap==0){
                        if(k==0)outputz1<<i<<'\t'<<Ytrain_aug[k][i]<<'\t'<<_Ybar_train[k][i]<<endl;
                        if(k==1)outputz2<<i<<'\t'<<Ytrain_aug[k][i]<<'\t'<<_Ybar_train[k][i]<<endl;
                    }
                }
            }

            double mean_nrmse[Ninputs+1];
            cout<<"mean_nrmse=";

            for(int k=0;k<1;k++){
                mean_nrmse[0]=0.0;
                for(int i=0;i<Nsamples;i++){
                    double tip=sqrt(NRMSE(_Ybar_train[k],Ytrain_aug[k],i*NTrain,(i+1)*NTrain));
                    if(tip>1.0)tip=1.0;
                    mean_nrmse[k]+=tip;

                }
            }
            printf("%f",mean_nrmse[0]/Nsamples);                       
            outtrain<<NDelay<<'\t'<<kap_min<<'\t'<<input_scaling<<'\t'<<mean_nrmse[0]/Nsamples<<endl;

            obj1.def(Xtest_aug,Nsamples,NTest,NThetas);
            obj1.def(new_Xtest_aug,Nsamples*NTest,ntheta);

            inpindex=0;
            train_index=0;
            for(int i = 0; i < NTesth; i++){
                 double external[Nsamples];
                 for(int j=0;j<Nsamples;j++){
                     external[j]=Inputs[kmask][j]*Preprocessing_Mask[inpindex];
                 }
                 Euler(a,ay, h, t,external);
                 if((i+add)%thetah==0){for(int j=0;j<Nsamples;j++){Xtest_aug[j][train_index][inpindex]=a[j];}inpindex=(inpindex+1)%NThetas;}
                 if((i+add)%cycle==0) {kmask++;train_index++;}
             }

            new_x(Xtest_aug,new_Xtest_aug,NTest,CCC,ntheta,Nsamples);

            for(int i=0;i<NTest;i++){
                 for(int j=0;j<Nsamples;j++){
                     for(int k=0;k<Ninputs+1;k++){
                         Ytest_aug[k][j*NTest+i]=Outputs[i+NTrain][k][j];

                     }

                 }

             }
             for(int i=0;i<NTest*Nsamples;i++){
                for(int k=0;k<Ninputs+1;k++){
                    _Ybar_test[k][i]=0;
                    for(int j=0;j<ntheta;j++){
                        _Ybar_test[k][i]+=new_Xtest_aug[i][j]*weights[k][j];
                    }
                    if(nth_gamma==0 && nth_kap==0){
                        if(k==0)outputz4<<i<<'\t'<<Ytest_aug[k][i]<<'\t'<<_Ybar_test[k][i]<<endl;
                        if(k==1)outputz5<<i<<'\t'<<Ytest_aug[k][i]<<'\t'<<_Ybar_test[k][i]<<endl;
                    }
                }
            }


            for(int k=0;k<2;k++){
                mean_nrmse[k]=0.0;
                for(int i=0;i<Nsamples;i++){
                    double tip=sqrt(NRMSE(_Ybar_test[k],Ytest_aug[k],i*NTest,(i+1)*NTest));
                    if(tip>1.0)tip=1.0;
                    mean_nrmse[k]+=tip;

                }
            }
            cout<<endl<<"mean_nrmse1=";
            printf("%f\t%f",mean_nrmse[0]/Nsamples,mean_nrmse[1]/Nsamples);
            cout<<endl;
            if(Case=="Spatial"){
                outtest<<NDelay<<'\t'<<input_scaling<<'\t'<<kap_min<<'\t'<<mean_nrmse[0]/Nsamples<<'\t'<<mean_nrmse[1]/Nsamples<<endl;
            }
            
            if(Case=="temporal"){
                outtest<<NDelay<<'\t'<<temp<<'\t'<<tau<<'\t'<<mean_nrmse[0]/Nsamples<<'\t'<<mean_nrmse[1]/Nsamples<<endl;
            }
            obj1.dd(taubar);
            obj1.dd(delayIndex);
            obj1.dd(delayNum);
            obj1.dd(eDelay,Nsamples,MaxDelayPerLayer);
            obj1.dd(weights,Ninputs+1);
            obj1.dd(Ytrain_aug,Ninputs+1);
            obj1.dd(_Ybar_train,Ninputs+1);
            obj1.dd(Xtest_aug,Nsamples,NTest);
            obj1.dd(new_Xtest_aug,Nsamples*NTest);
            obj1.dd(_Ybar_test,Ninputs+1);
            obj1.dd(Ytest_aug,Ninputs+1);
            obj1.dd(Xtrain_aug,Nsamples,NTrain);
            obj1.dd(new_Xtrain_aug,Nsamples*NTrain);


        }


    }
    obj1.dd(Inputs,NTrain+Discard+NTest);
    obj1.dd(Outputs,NTrain+NTest,Ninputs+1);
    obj1.dd(Preprocessing_Mask);
    obj1.dd(Inputs1,NTrain+Discard+NTest);
    obj1.dd(Outputs1,NTrain+NTest);


}
