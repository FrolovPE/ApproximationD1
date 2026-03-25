#ifndef SPLINE_H
#define SPLINE_H

int spline_coeff(int n,double *x, double *y /*array of f(x_i)*/, double *coeff, double *dopvector);
double spline_in_point(double t,/*, double a, double b,*/ int n ,const double *x,const double *coeff,double *dopvector);
int ksi_append(int n, const double *x,double *ksi);
int tridiag(int n, double *lowd, double *maind, double *upd, double *b, double *solution);
int find_local_func(double t, int n, const double *ksi);


#endif
