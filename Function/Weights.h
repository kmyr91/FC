void Normal_weights_Brunner(double **X,double **Ytrain_aug,double **weights,int ntheta,int nsamples){

    mat A(NTrain*nsamples,ntheta,fill::zeros);
    mat target(NTrain*nsamples,Ninputs+1,fill::zeros);
    mat I=eye(ntheta,ntheta);
    mat WEIGHTS(ntheta,Ninputs+1,fill::zeros);
    mat var(ntheta,ntheta,fill::zeros);

    double lambda=1e-8;
    cout<<"lambda="<<lambda<<endl;

    for(int i=0;i<NTrain*nsamples;i++){

        for(int k=0;k<Ninputs+1;k++){
            target(i,k)=Ytrain_aug[k][i];
        }
        for(int j=0;j<ntheta;j++){
            double D=0.0*randn(0.0,1.0)*lambda;
            A(i,j)=X[i][j]+D;

        }

    }
    for(int i=0;i<ntheta;i++){
        for(int j=0;j<ntheta;j++){
            I(i,j)*=(lambda);

        }
    }
    var=((A.t()*A)+I);
    WEIGHTS=pinv(var)*(A.t()*target);
    for(int i=0;i<ntheta;i++){
        for(int k=0;k<Ninputs+1;k++){
            weights[k][i]=WEIGHTS(i,k);

        }
    }
}
