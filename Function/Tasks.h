
void Narma10(int N,double **Inputs,double ***Outputs,string o,string cc,int kth_sample){

    Definitions obj1;
    if(kth_sample==0)cout<<"s1="<<s1<<'\t'<<"s2="<<s2<<endl;
    uniform_real_distribution<double> distribution(0.0,0.5);
    random_device device;
    mt19937 generator(kth_sample+5);


    int warmup=2;
    /**********Loading second Input**********/
    double deleting11[N+warmup+10][Nsamples];
    ifstream input1("mg17");
    for (int i = 0; i < N+warmup+10; i++) {
        for(int j=0;j<Nsamples;j++){
        input1 >>deleting11[i][j];
        }
    }
    /**********Initializing Inputs**********/
    double **y_base,**u;
    obj1.def(y_base,Ninputs,N+warmup+10);
    obj1.def(u,Ninputs,N+warmup+10);
    for(int i=0;i<N+warmup+10;i++){
        for(int k=0;k<Ninputs;k++){
            u[k][i]=0.0;
            y_base[k][i]=0.0;
        }
    }

    for(int i=0;i<N+warmup+10;i++){
        u[0][i]=distribution(generator);
        u[1][i]=deleting11[i][kth_sample];
    }
    /**********Initializing Outputs**********/
    double SUM;
    for(int i=10;i<N+warmup+10;i++){
        SUM=0.0;
        for(int j=1;j<11;j++){
            SUM+= y_base[0][i-j];

        }
        double t=i*5e-3;
        y_base[0][i]=(0.3 * y_base[0][i-1] + 0.05 * y_base[0][i-1] * SUM + 1.5 * u[0][i-1] * u[0][i-10] + 0.1);

    }
//        normalizeData1(u[0],N+warmup+10,0.0,1.0);
//        normalizeData1(y_base[0],N+warmup+10,0.0,1.0);
      normalizeData1(u[1],N+warmup+10,0.0,1.0);
int fin=10-step;
    normalizeData1(u[1],N+warmup+10,0.0,1.0);

    for(int i=0;i<N+warmup+fin;i++){
        y_base[1][i]=u[1][i+step-1];
    }
    for(int i=0;i<N;i++){
      Inputs1[i][0]=s1*u[0][i+warmup+fin]+s2*u[1][i+warmup+fin];
    }
    for(int k=0;k<Ninputs+1;k++){
    for(int i=Discard;i<N;i++){
        if(k<Ninputs)Outputs1[i-Discard][k]=y_base[k][i+warmup+fin+step];
        if(k==Ninputs)Outputs1[i-Discard][k]=s1*Outputs1[i-Discard][0]+s2*Outputs1[i-Discard][1];
    }




    }
    obj1.dd(y_base,Ninputs);
    obj1.dd(u,Ninputs);
        for(int i=0;i<N;i++){
            Inputs[i][kth_sample]=Inputs1[i][0];
        }
        for(int k=0;k<Ninputs+1;k++){
            for (int i = 0; i < N-Discard; i++) {
                Outputs[i][k][kth_sample]=Outputs1[i][k];
//                if(k==0 && kth_sample==3) printf("%e\t%d\n",Outputs[i][k][kth_sample],kth_sample);
            }
    }
//  if(kth_sample==3) exit(1);
}


void SANTA(int N,double **Inputs,double ***Outputs,int kth_sample,string o,string cc){

    double deleting1[4][NTrain+Discard+NTest+step],deleting2[4][NTrain+Discard+NTest+step],deleting3[4][NTrain+Discard+NTest+step];

    ofstream outg1("Data/extra/general1T1.txt");
    Definitions obj1;

//   ifstream input1("Lorenz/Data/X");
    ifstream input2("Lorenz/Y");
    ifstream input3("mg14");
    ifstream input4("Lorenz/Data/Z");
     ifstream input1("mg17");
//  ifstream input3("Santa1.txt");
    double omega;
    for (int i = 0; i < (N+step); i++) {
        for(int j=0;j<4;j++){
        input1 >>deleting1[j][i];
        input2 >>deleting2[j][i];
        if(cc=="MG")input3 >>deleting3[j][i];
        if(cc=="unique") input4 >>deleting3[j][i];

        }
    }

    double *u;
    obj1.def(u,N+step);
    /*****************************************************//*

    double desired_mean=MEANx(deleting1,N+step,0,N+step);
    double desired_sd=calculateSDx(deleting1,N+step,0,N+step);

    normalizeData1(deleting3,N+step,desired_mean,desired_sd);*/
    /*****************************************************/

        normalizeData1(deleting1[kth_sample],N+step,0.0,1.0);
        normalizeData1(deleting3[kth_sample],N+step,0.0,1.0);

    for(int i=0;i<N+step;i++){
        u[i]=s1*deleting1[kth_sample][i]+s2*deleting3[kth_sample][i];

    }
    if(kth_sample==0)cout<<N+step<<endl;

//     normalizeData1(u,N,0.0,1.0);
    for(int i=0;i<N;i++){
        Inputs1[i][0]=u[i];

    }
    if(kth_sample==0)cout<<MEANx(u,N,0,N)<<'\t'<<calculateSDx(u,N,0,N)<<'\t'<<s1<<'\t'<<s2<<endl;

    for(int i=0;i<N;i++){
        Inputs[i][kth_sample]=Inputs1[i][0];

    }
    for(int k=0;k<3;k++){
        for(int i=Discard;i<N;i++){
            if(k==0)Outputs1[i-Discard][k]=deleting1[kth_sample][i+step];
            if(k==1)Outputs1[i-Discard][k]=deleting3[kth_sample][i+step];
            if(k==2)Outputs1[i-Discard][k]=u[i+step];


        }

    }
    obj1.dd(u);
    for(int k=0;k<3;k++){
        for (int i = 0; i < N-Discard; i++) {
            Outputs[i][k][kth_sample]=Outputs1[i][k];

        }

    }


}

void choose_task(string datastyle){
    string datastyle1="osc1";
    string datastyle2="MG";
    string datastyle3="unique";
    for(int kth_sample=0;kth_sample<Nsamples;kth_sample++){
        if(datastyle=="Narma10"){
            Narma10(NTrain+Discard+NTest,Inputs,Outputs,datastyle,datastyle1,kth_sample);
            if(kth_sample==0)cout<<"input data is Narma10"<<endl;
        }
        if(datastyle=="SANTA"){
            SANTA(NTrain+Discard+NTest,Inputs,Outputs,kth_sample,datastyle,datastyle3);
            if(kth_sample==0)cout<<"input data is SANTA"<<endl;
        }

    }

}
