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
    double p1[6] = {  3./44,-2./33,-3./11,0.,-37./132,6./11};
    double p2[6] = {-0.033349806441991,-0.011116602147330,0.377964473009228,
                    -0.733695741723794,0.544713505219180,-0.144515827915293};
    for(int i=0;i<6;i++) a[i] = p1[i]*w+p2[i];
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
    double alp = M_PI/3.;
    double ki2s = 0;
    double kr2s = 0;
    for(int i=1;i<=10;i++){
        ki2s += (ki((alp*i)/10,a)-(alp*i)/10)*(ki((alp*i)/10,a)-(alp*i)/10);
        kr2s += kr((alp*i)/10,a)*kr((alp*i)/10,a);
    }
    return 2.*ki2s+5.*kr2s;  
}
double Jw(double w){
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
double Newscheme(double u[],int usize,double a[6],int idx){
    double ux=0.;
    for(int i=0;i<6;i++){
        ux += u[((idx+i-4)+(usize-1))%(usize-1)]*a[i];
    }
    return ux;
}
void updateu(double u[],int usize,double a[6],double dx,double dt,int Nt){
    double* unext = new double[usize];
    double* u1 = new double[usize];
    double* u2 = new double[usize];
    int i;
    double cfl = dt/dx;
    double *tmp;
    for(int j=1;j<=Nt;j++){
        for(i=0;i<usize;i++){
            u1[i] = u[i] - cfl*Newscheme(u,usize,a,i);    
        }
        for(i=0;i<usize;i++){
            u2[i] = 0.75*u[i] + 0.25*(u1[i] - cfl*Newscheme(u1,usize,a,i));
        }
        for(i=0;i<usize;i++){
            unext[i] = u[i]/3.0 + 2.0/3.0*(u2[i] - cfl*Newscheme(u2,usize,a,i));
        }
        tmp = u;
        u = unext;
        unext = tmp;   
    }
    delete[] u1;
    delete[] u2;
    delete[] unext;
}
void ordertest(double a[]){
    int minN = 6,N;
    double maxdx = 2.*M_PI/minN;
    double dx;
    double err;
    double ux_calc,ux_exact;
    ofstream csvfile;
    csvfile.open("ordertest.csv", ios::out | ios::trunc);
    csvfile << "dx" << ","<<"err" << endl;
    for(int i=1;i<=10;i++){
        dx = maxdx/i;
        N = minN*i+1;
        double *u = new double[N];
        for(int j=0;j<N;j++){
            u[j] = sin(j*dx);
        }
        err = 0.0;
        for(int j=0;j<N;j++){
            ux_calc = Newscheme(u,N,a,j)/dx;
            ux_exact = cos(j*dx);
            err += abs(ux_calc - ux_exact);
        }
        err = err/N;
        csvfile << dx<<","<<err<<endl;
        delete[] u;
    }
    csvfile.close();
}
int main(){
    double a[6]={-0.7500,-0.3333,9.0000,-18.0000,13.0833,-3.0000};
    outputcurve("5order",a);
    ordertest(a);//精度验证
    double w ;

    //开始向时间步推进
    double dt = 0.01;
    int Nx = 21;
    double dx = 2*M_PI/(Nx-1);
    
    double *u = new double[Nx];
    int Nt = 100000;
    for(int i=0;i<Nx;i++){
        u[i] = sin(i*dx);
    }
    updateu(u,Nx,a,dx,dt,Nt);  //计算Nt*dt时刻的u
    
    
    w = golds(Jw,-10.,10.);  //黄金分割法求最佳参数w
    cout << w <<endl;
    awithw(a,w);  //更新a
    for(int i=0;i<6;i++){
        cout << a[i] << "\t" ;
    }
    outputcurve("",a);   //输出耗散色散曲线
    ordertest(a);//精度验证
    double *uf = new double[Nx];
    for(int i=0;i<Nx;i++){
        uf[i] = sin(i*dx);
    }
    updateu(uf,Nx,a,dx,dt,Nt);
    
    ofstream csvfile;
    csvfile.open("calc_exact.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," <<  "5ord" <<","<< "this"<<","<< "Exact" << endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx << "," <<  u[i] <<","<< uf[i]<<","<< sin(i*dx-dt*Nt) << endl;
    }
    csvfile.close();
    delete[] u;
    return 0;
}