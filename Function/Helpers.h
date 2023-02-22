double MEAN(double data[],int n,int begin,int end){
    double sum = 0.0, mean, standardDeviation = 0.0;
    int i;
    for(i = begin; i < end; ++i){
        sum += data[i];
    }
    mean = sum/n;

    return mean;
}
 int getMaxindex(double inputD){
         double MAXxx=inputD*(1.0/h);
         int outputI;
         MAXxx = MAXxx + 0.5 - (MAXxx<0);
         outputI = (int)(MAXxx);
     return outputI;
}
double MeaN(double *data,int n){
    double sum=0.0;
    for(int i=0;i<n;i++){
        sum+=data[i];

    }
    return sum/n;
}
 void sort(double a[],int n){
    int i,j;
    double TEMP;
   for(i=0;i<n;i++){
for(j=0;j<n-1;j++){
if(a[j]>a[j+1]){
TEMP=a[j];
a[j]=a[j+1];
a[j+1]=TEMP;
}
}
}
}
double find_max(double list[],int MAX_delay){
    double highNum=list[0];
    double index=0;
    for (int m = 0 ; m < MAX_delay ; m++)
    {
        if (list[m] > highNum){
            highNum = list[m];
        index=highNum;
        }
//             cout << list[m];
    }
    return index;
}
double find_min(double list[],int MAX_delay){
    double highNum=list[0];
    double index=0;
    for (int m = 0 ; m < MAX_delay ; m++)
    {
        if (list[m] < highNum){
            highNum = list[m];
        index=highNum;
        }
//             cout << list[m];
    }
    return index;
}
double calculateSD(double data[],int n,int begin,int end){
    double sum = 0.0, mean, standardDeviation = 0.0;
    int i;
    for(i = begin; i < end; ++i){
        sum += data[i];
    }
    mean = sum/n;
    for(i = begin; i < end; ++i){
        standardDeviation += pow(data[i] - mean, 2);
    }
    return sqrt(standardDeviation / n);
}
double calculateSDx(double *data,int n,int begin,int end){
    double sum = 0.0, mean, standardDeviation = 0.0;
    int i;
    for(i = begin; i < end; ++i){
        sum += data[i];
    }
    mean = sum/n;
    for(i = begin; i < end; ++i){
        standardDeviation += pow(data[i] - mean, 2);
    }
    return sqrt(standardDeviation / n);
}
double MEANx(double *data,int n,int begin,int end){
    double sum = 0.0, mean, standardDeviation = 0.0;
    int i;
    for(i = begin; i < end; ++i){
        sum += data[i];
    }
    mean = sum/n;
    return mean;
}
void normalizeData(double *data,int ndata){
    double meanx=MEANx(data,ndata,0,ndata),sdx=calculateSDx(data,ndata,0,ndata);
    for(int i=0;i<ndata;i++){
        data[i]=(data[i]-meanx)/(sdx);
    }
}
void normalizeData1(double *data,int ndata,double desired_mean,double desired_sd){
    double meanx=MEANx(data,ndata,0,ndata),sdx=calculateSDx(data,ndata,0,ndata);
    for(int i=0;i<ndata;i++){
        data[i]=desired_mean+(data[i]-meanx)/(sdx/desired_sd);
    }
}
void normalizeData2(double *data,int ndata){
    double max=find_max(data,ndata),min=find_min(data,ndata);
    cout<<min<<endl;
    for(int i=0;i<ndata;i++){
        data[i]=2.0*((data[i]-min)/(max-min))-1.0;
    }
}

double NRMSE(double x[],double y[],int begin,int end){
    double mean=0.0;
    double var=0.0;
    int length=end-begin;
    var=calculateSD(y,length,begin,end);
    var*=var;
    mean=0.0;
    for(int i=begin;i<end;i++){
        mean+=(((x[i])-(y[i]))*((x[i])-(y[i])));

    }
    return /*sqrt()*/mean/((end-begin)*var);
}

double setDigits(double _number, int _digits){
    double tenth = pow((double)10,_digits);
    _number *= tenth;
    _number = floor(_number);
    _number /= tenth;

    return _number;
}
