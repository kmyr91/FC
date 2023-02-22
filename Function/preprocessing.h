void CreateSamples(double X[]){
    double offset=0.0;
    random_device device;
    mt19937 generator(1);
    string o="4levels";
    if(o=="2levels"){
    uniform_int_distribution<int> distribution(0,1);
    double x[2]={0.0,1.0};
    for (int i = 0; i < NThetas; i++){
        int y=distribution(generator);
        X[i]=((2.0*distribution(generator))-1.0) ;
//         X[i]*=input_scaling;

    }
    }
    if(o=="4levels"){
    uniform_int_distribution<int> distribution(0,3);
    for (int i = 0; i < NThetas; i++){
        int random_index = distribution(generator);
        switch(random_index) {
            case 0:
                X[i] = -1.0;
                break;
            case 1:
                X[i] = -1.0/3.0;
                break;
            case 2:
                X[i] = 1.0/3.0;
                break;
            case 3:
                X[i] = 1.0;
                break;
            
                break;

           
            default:
                break;
        }
    }
}
if(o=="6levels"){
    uniform_int_distribution<int> distribution(0,5);
    for (int i = 0; i < NThetas; i++){
        int random_index = distribution(generator);
        switch(random_index) {
            case 0:
                X[i] = -1.0;
                break;
            case 1:
                X[i] = -3.0/5.0;
                break;
            case 2:
                X[i] = -1.0/5.0;
                break;
            case 3:
                X[i] = 1.0/5.0;
                break;
            case 4:
                X[i] = 3.0/5.0;
                break;
            case 5:
                X[i] = 1.0;
                break;

           
            default:
                break;
        }
    }
}
}

void new_x(double ***X,double **new_X,int n,double CCC,int ntheta,int nsamples){
        if(ntheta==NThetas+1){
            for(int k=0;k<nsamples;k++){
                for(int i=0;i<n;i++){
                    new_X[i+k*n][0]=CCC;
                    for(int j=1;j<ntheta;j++){
                        new_X[i+k*n][j]=X[k][i][j-1];
                    }
                }
            }
        }
        else {
            for(int k=0;k<nsamples;k++){
                for(int i=0;i<n;i++){
                    for(int j=0;j<ntheta;j++){
                        new_X[i+k*n][j]=X[k][i][j];
                    }
                }
            }
        }
}
