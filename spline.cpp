#include "spline.h"
#include <iostream>

int ksi_append(int n, const double *x,double *ksi)
{
    for(int i = 1; i < n; i++)
        ksi[i] = (x[i] + x[i-1])/2;
        
    ksi[0] = x[0] - (x[1] - x[0])/2;
    ksi[n] = x[n-1] + (x[n-1] - x[n-2])/2;

    return 0;
}

int tridiag(int n, double *lowd, double *maind, double *upd, double *b, double *solution)
{
    //upd 0,..,n-2
    //maind 0,...,n-1
    //lowd 0,..,n-2

    for(int i = 0; i < n -1 ; i++)
    {
        double h = lowd[i]/maind[i];
        maind[i+1] -= h*upd[i];
        b[i+1] -= h * b[i];
    }

    solution[n-1] = b[n-1]/maind[n-1];

    for(int i = n - 2; i >= 0;i--)
    {
        solution[i] = (b[i] - solution[i+1]*upd[i])/maind[i];
    }

    return 0;

}
int spline_coeff(int n,double *x, double *y /*array of f(x_i)*/, double *coeff, double *dopvector)
{
    //coeff [3 * (n-1)]
    // Распределяем dopvector как:
    // ksi 0,..,n+1
    // lowd 0,..,n
    // maind 0,...,n+1
    // upd 0,...,n
    //solution = v 0,...,n+1

    double *ksi = dopvector;
    double *lowd = ksi + n + 1;
    double *maind = lowd + n + 1;
    double *upd = maind + n + 1;
    double *b = upd + n + 1;
    double *v = b + n + 1;
    //выделил для всех n+1 для удобства 

    ksi_append(n,x,ksi);

    for (int i = 0; i < n + 1; i++)
    {
        lowd[i] = 0;
        maind[i] = 0;
        upd[i] = 0;
        v[i] = 0;
    }

    //Естественное условие 1 Р_0''(x_0) = 0

    double A = 1 / (x[0] - ksi[0]);
    double B = 1 / (ksi[1] - x[0]);
    maind[0] = A;
    upd[0] = B;
    b[0] = y[0] * (A + B);

    //Естественное условие 2 Р_n''(x_n) = 0

    A = 1 / (x[n-1] - ksi[n-1]);
    B = 1 / (ksi[n] - x[n-1]);
    maind[n] = A;
    upd[n-1] = B;
    b[n] = y[n-1] * (A + B); 


    //Метод параболических сплайнов

    for (int i = 1; i < n; i++)
    {
        double left = 1 / (x[i-1] - ksi[i-1]) - 1 / (ksi[i] - ksi[i-1]);

        double mid = 1 / (ksi[i] - x[i-1]) + 1 / (ksi[i] - ksi[i-1]) + 
        1 / (x[i] - ksi[i]) + 1 / (ksi[i+1] - ksi[i]);

        double right = 1 / (ksi[i+1] - x[i]) - 1 / (ksi[i+1] - ksi[i]);

        double b_i = (1 / (x[i-1] - ksi[i-1]) + 1 / (ksi[i] - x[i-1])) * y[i-1] +
                     (1 / (x[i] - ksi[i]) + 1 / (ksi[i+1] - x[i])) * y[i];

        lowd[i-1] = left;
        maind[i] = mid;
        upd[i] = right;
        b[i] = b_i;
    }

    tridiag(n,lowd,maind,upd,b,v);

    for(int i = 0 ; i < n; i++)
    {
        double denom = x[i] - ksi[i];
        
        double a1 = v[i];
        double a2 = (y[i] - a1) / denom;
        double a3 = ( (v[i+1] - y[i])/(ksi[i+1] - x[i]) - a2  ) / (ksi[i+1] -ksi[i]);

        coeff[3 * i + 0] = a1;
        coeff[3 * i + 1] = a2;
        coeff[3 * i + 2] = a3;
    }

    return 0;
}

double spline_in_point(double t,/*, double a, double b,*/ int n ,const double *x, const double *coeff,double *dopvecor)
{

    double *ksi = dopvecor;

    ksi_append(n,x,ksi);

    int i = find_local_func(t,n,ksi);

    double a1 = coeff[3 * i + 0];
    double a2 = coeff[3 * i + 1];
    double a3 = coeff[3 * i + 2];

    return a1 + a2 * (t - ksi[i]) + a3 * (t-x[i])*(t - ksi[i]);

}


int find_local_func(double t, int n, const double *ksi)
{
    if(t <= ksi[0]) return 0;
    if(t >= ksi[n]) return n - 1;

    int left,mid,right;

    left = 0;
    right = n;

    while (right - left > 1)
    {
        mid = (right + left)/2;
        if(t <= ksi[mid])
            right = mid;
        else
            left = mid;
    }

    return left;

}