#include<iostream>
#include <fstream>
#include<math.h>
#include<string>

using namespace std;
typedef double (*Fun) (double*);

//用w更新a
void awithw(double a[6],double w[2]){
    double p1[6] = {-0.173254752956534  , 0.495427762543929 , -0.249163520610373 , 
                    -0.492528483867111, 0.617110244172297 , -0.197591249282208};
    double p2[6] = {-0.098312411528301  , 0.418552568664483 , -0.691086159374925 , 
                    0.545067181420883, -0.199524101733420 , 0.025302922551279};
    double p3[6] = {0.2000,-7./12.,         0.,         0.,    1./3.,    0.0500,};
    for(int i=0;i<6;i++) a[i] = p1[i]*w[1] +p2[i]*w[2] +p3[i];
}
double kr(double alp,double a[6]){
    double ret = 0.;
    for(int i=0;i<6;i++){
        ret += a[i]*cos((i-3)*alp);
    }
    return ret;
}
double ki(double alp,double a[6]){
    double ret = 0.;
    for(int i=0;i<6;i++){
        ret += a[i]*sin((i-3)*alp);
    }
    return ret;
}
double Ja(double a[6]){
    double ret = (2*a[0])/3 - (a[1] - a[5]) + (a[0]*a[0])/2 + (2*a[2] - 2*a[4]) 
    + ((a[1] - a[5])*(a[1] - a[5]))/2 + ((a[2] - a[4])*(a[2] - a[4]))/2+
    2*((a[0]*a[0])/2 + ((a[1] + a[5])*(a[1] + a[5]))/2 + ((a[2] + a[4])*(a[2] + a[4]))/2); 
    return ret;
}
double Jw(double w[2]){
    double a[6];
    awithw(a,w);
    return Ja(a);
}


void outputcurve(string filename, double a[]){
    filename += "alpha,kr,ki.csv";
    ofstream csvfile;
    csvfile.open(filename, ios::out | ios::trunc);
    csvfile << "alpha" << "," << "kr" << "," << "ki" <<","<<"alpha"<< endl;
    for(int i=0;i<20;i++){
        csvfile << M_PI*i/20 << "," << kr(M_PI*i/20,a) <<","<<ki(M_PI*i/20,a)<<","<< M_PI*i/20<< endl;
    }
    csvfile.close();
}
double Newscheme(double u[],int usize,double a[],int asize,int idx){
    double ux=0.;
    for(int i=0;i<asize;i++){
        ux += u[((idx+i-3)+(usize-1))%(usize-1)]*a[i];
    }
    return ux;
}
void updateu(double u[],int usize,double a[],int asize,double dx,double dt,int Nt){
    double* unext = new double[usize];
    double* u1 = new double[usize];
    double* u2 = new double[usize];
    int i;
    double cfl = dt/dx;
    double *tmp;
    for(int j=1;j<=Nt;j++){
        for(i=0;i<usize;i++){
            u1[i] = u[i] - cfl*Newscheme(u,usize,a,asize,i);    
        }
        u2[0] =0.75*u[0] + 0.25*(u1[0] - cfl*(u1[0]-u1[usize-1])) ;
        for(i=0;i<usize;i++){
            u2[i] = 0.75*u[i] + 0.25*(u1[i] - cfl*Newscheme(u1,usize,a,asize,i));
        }
        unext[0] =u[0]/3.0 + 2.0/3.0* (u2[0] - cfl*(u2[0]-u2[usize-1]));
        for(i=0;i<usize;i++){
            unext[i] = u[i]/3.0 + 2.0/3.0*(u2[i] - cfl*Newscheme(u2,usize,a,asize,i));
        }
        tmp = u;
        u = unext;
        unext = tmp;   
    }
    delete[] u1;
    delete[] u2;
    delete[] unext;
}
int main(){
    double a[6]={-0.03333333,0.25,-1,0.333333333,0.5,-0.05};
    outputcurve("5order",a);
    double w[2] = { -4.,-4.} ;
    

    //开始向时间步推进
    double dt = 0.01;
    int Nx = 11;
    double dx = 2*M_PI/(Nx-1);
    
    double *u = new double[Nx];
    int Nt = 1000;
    for(int i=0;i<Nx;i++){
        u[i] = sin(i*dx);
    }
    updateu(u,Nx,a,6,dx,dt,Nt);

    awithw(a,w);
    outputcurve("",a);
    for(int i=0;i<6;i++){
        cout << a[i] << ",";
    }
    double *uf = new double[Nx];
    for(int i=0;i<Nx;i++){
        uf[i] = sin(i*dx);
    }
    updateu(uf,Nx,a,6,dx,dt,Nt);
    
    ofstream csvfile;
    csvfile.open("calc_exact.csv", ios::out | ios::trunc);
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx << "," <<  u[i] <<","<< uf[i]<<","<< sin(i*dx-dt*Nt) << endl;
    }
    csvfile.close();
    delete[] u;
    delete[] uf;
    return 0;
}