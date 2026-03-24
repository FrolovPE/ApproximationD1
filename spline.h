#ifndef SPLINE_H
#define SPLINE_H

int spline_coeff(int n,double *x, double *y /*array of f(x_i)*/, double *coeff, double *tmp);
double spline_in_point(double t,/*, double a, double b,*/ int n ,double *x, double *coeff);

#endif
