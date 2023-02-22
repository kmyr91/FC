class Definitions { 
    public: 
      
    void def(double *&x,int n){
        x=new double[n];
        for(int i=0;i<n;i++){
        x[i]=0.0;
    }
}
void def(double **&x,int n,int m){
x=new double*[n];
for(int i=0;i<n;i++){
x[i]=new double[m];
}
 for(int i=0;i<n;i++){
     for(int j=0;j<m;j++){
        x[i][j]=0.0;
     }
    }
}
void def(double ***&x,int n,int m,int l){
x=new double**[n];
for(int i=0;i<n;i++){
x[i]=new double*[m];
for(int j=0;j<m;j++){
x[i][j]=new double[l];
}
}
}

void def(int *&x,int n){
x=new int[n];
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }
}
void def(int **&x,int n,int m){
x=new int*[n];
for(int i=0;i<n;i++){
x[i]=new int[m];
}
 for(int i=0;i<n;i++){
     for(int j=0;j<m;j++){
        x[i][j]=0.0;
     }
    }
}
void def(int ***&x,int n,int m,int l){
x=new int**[n];
for(int i=0;i<n;i++){
x[i]=new int*[m];
for(int j=0;j<m;j++){
x[i][j]=new int[l];
}
}
}
void dd(double *x){

    delete[] x;
}
void dd(double **x,int n){
for(int i=0;i<n;i++){
    delete[] x[i];
}
    delete[] x;
}
void dd(double ***x,int n,int m){
    
for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
    delete[] x[i][j];
    }
    delete[] x[i];
    
}
    delete[] x;
}

void dd(int *x){

    delete[] x;
}
void dd(int **x,int n){
for(int i=0;i<n;i++){
    delete[] x[i];
}
    delete[] x;
}
void dd(int ***x,int n,int m){
for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
    delete[] x[i][j];
    }
    delete[] x[i];
    
}
    delete[] x;
}
}; 

