
 void delaytimes(int NDE,double theta){
     srand(0);
     Definitions obj1;
     obj1.def(taubar,MaxDelayPerLayer);

     int ND=NDelay;
     
     
     
     for(int j=1; j<NDelay+1; j++){
         double min_dist=0.0,max_dist=temp;
         double rng=unifRand(min_dist,max_dist);
         rng=rng-remainder(rng,theta);
         if(j==1)rng=0.0;
         taubar[j-1]=tau+(j-1)*temp;
//            cout<<j<<'\t'<<Tau[i][j-1]<<endl;
	}
	sort(taubar,NDelay);
        cout<<"average of the delays="<<MeaN(taubar,NDelay)<<'\t'<<taubar[0]<<'\t'<<taubar[NDelay-1]<<endl;
        tau_Max=taubar[NDelay-1]+1.0;

}
void calcParameter(double h){
    Definitions obj1;
    obj1.def(delayNum,MaxDelayPerLayer);
    obj1.def(delayIndex,MaxDelayPerLayer);
    /**************************************************************************************************/
    delta=ohm0;
    kap=kap_min;
    cout<<"delaynum=";
    for(int j=0; j<NDelay; j++){
        delayNum[j] = getMaxindex(taubar[j]);
        if(j==0)cout<<delayNum[j]<<'\t';
        if(j!=0)cout<<delayNum[j]<<"("<<delayNum[j]-delayNum[j-1]<<")"<<'\t';
    }
    cout<<endl;
    cout<<"input_scaling="<<input_scaling<<endl;
    s1=sqrt(s),s2=sqrt(1.0-s);
}
