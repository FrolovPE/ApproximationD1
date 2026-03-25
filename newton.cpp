#include "newton.h"
#include "stdio.h"

int newton_coeff(int n,double *x, double *y /*array of f(x_i)*/, double *coeff, double *tmp)
{
    for(int i = 0 ; i < n; i++)
    {
        tmp[i] = y[i];
    }

    coeff[0] = tmp[0];

    for(int j = 1; j < n; j++)
    {
        for(int i = 0; i < n -j; i++)
        {
            tmp[i] = (tmp[i+1] - tmp[i])/(x[i+j] - x[i]);
        }
        coeff[j] = tmp[0];
    }
    return 0;
}


double newton_in_point(double t,/*, double a, double b,*/ int n ,double *x, double *coeff,double *tmp)
{
    (void)tmp;
    double res = coeff[n-1];

    for(int i = n - 2; i >= 0; i--)
    {
        res = res * (t - x[i]) + coeff[i];
    }

    return res;
}