#ifndef NEWTON_H
#define NEWTON_H

int newton_coeff(int n,double *x, double *y /*array of f(x_i)*/, double *coeff, double *tmp);
double newton_in_point(double t,/*, double a, double b,*/ int n ,double *x, double *coeff);

#endif
