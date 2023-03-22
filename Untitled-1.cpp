#include<iostream>
#include <fstream>
#include<math.h>
#include<string>

using namespace std;
typedef double (*Fun) (double);

//优化函数：黄金分割法
double golds(Fun f,double a,double b,double eps=1e-7,double delta=1e-7){
    double p = a + 0.381966*(b-a);
    double q = a + 0.618034*(b-a);
    double fp = f(p),fq = f(q);
    while((q-p)>delta || abs(fp-fq)>eps){
        if(fp<fq){
            b = q;
            q = p;
            fq = fp;
            p = a + 0.381966*(b-a);
            fp = f(p);
        }
        else{
            a = p;
            p = q;
            fp = fq;
            q = a + 0.618034*(b-a);
            fq = f(q);
        }
    }
    return (fp<fq ? p : q);
}

//用w更新a
void awithw(double a[6],double w){
    double p1[6] = {0.1,-0.5,1.,-1,0.5,-0.1};
    double p2[6] = {0.,0.08333333333,-0.6666666667,0.,0.6666666667,-0.08333333333};
    for(int i=0;i<6;i++) a[i] = p1[i]*w+p2[i];
}
double Ja(double a[6]){
    double ret = (2*a[0])/3 - (a[1] - a[5]) + (a[0]*a[0])/2 + (2*a[2] - 2*a[4]) 
    + ((a[1] - a[5])*(a[1] - a[5]))/2 + ((a[2] - a[4])*(a[2] - a[4]))/2+
    2*((a[0]*a[0])/2 + ((a[1] + a[5])*(a[1] + a[5]))/2 + ((a[2] + a[4])*(a[2] + a[4]))/2); 
    return ret;
}
double Jw(double w){
    double a[6];
    awithw(a,w);
    return Ja(a);
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
    if(idx<=2 || idx>=usize-3){
        ux = u[idx] - u[idx-1];
    }
    else{
        for(int i=0;i<asize;i++){
            ux += u[idx+i-3]*a[i];
        }
    }
    return ux;
}
void updateu(double u[],int usize,double a[],int asize,double dx,double dt,int Nt){
    double* unext = new double[usize];
    double* u1 = new double[usize];
    double* u2 = new double[usize];
    int i;
    double cfl = dt/dx;
    double uleft;
    double *tmp;
    ofstream csvfile;
    csvfile.open("u.csv", ios::out | ios::trunc);
    for(int j=1;j<=Nt;j++){
        
        u1[0] = u[0] - cfl*(u[0]-u[usize-1]);
        for(i=1;i<usize;i++){
            u1[i] = u[i] - cfl*Newscheme(u,usize,a,asize,i);    
        }
        u2[0] =0.75*u[0] + 0.25*(u1[0] - cfl*(u1[0]-u1[usize-1])) ;
        for(i=1;i<usize;i++){
            u2[i] = 0.75*u[i] + 0.25*(u1[i] - cfl*Newscheme(u1,usize,a,asize,i));
        }
        unext[0] =u[0]/3.0 + 2.0/3.0* (u2[0] - cfl*(u2[0]-u2[usize-1]));
        for(i=1;i<usize;i++){
            unext[i] = u[i]/3.0 + 2.0/3.0*(u2[i] - cfl*Newscheme(u2,usize,a,asize,i));
        }
        tmp = u;
        u = unext;
        unext = tmp;
        
        
        for(int i=0;i<usize;i++){
            csvfile << u[i] << "," ;
        }
        csvfile << endl;
        
    }
    csvfile.close();
    delete[] u1;
    delete[] u2;
    delete[] unext;
}
int main(){
    double a[6]={-0.03333333,0.25,-1,0.333333333,0.5,-0.05};
    outputcurve("5order",a);
    double w ;
    //while(w<2.){
    //    awithw(a,w);
    //    cout << w <<"," << Ja(a) << endl;
    //    w += 0.01;
    //}
    w = golds(Jw,-2.,2.);
    //cout << w << endl;    
    awithw(a,w);
    //cout << a;
    //outputcurve("",a);
    //cout << Ja(a);
    
    //开始向时间步推进
    double dt = 0.01;
    int Nx = 21;
    double dx = 2*M_PI/(Nx-1);
    
    double *u = new double[Nx];
    int Nt = 2000;
    for(int i=0;i<Nx;i++){
        u[i] = sin(i*dx);
    }
    for(int i=0;i<6;i++){
        cout << a[i] << ",";
    }
    updateu(u,Nx,a,6,dx,dt,Nt);
    
    delete[] u;
    
    //double s = 0;
    //for(int i=0;i<6;i++){
    //    s += a[i]*(i-3)*(i-3)*(i-3)*(i-3)*(i-3);
    //    cout << a[i] << " , ";
    //}
    //cout << endl;
    //cout << s <<endl;
    return 0;
}