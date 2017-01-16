#ifndef _SOMEFUNC_H_
#define _SOMEFUNC_H_

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

//structure pour mettre les data dans le stream
typedef struct data_{
  double* h ;
  double* fft_h;
  double* buffer_prec;
  unsigned int L;
  unsigned int M;
} * data ;



int fft(double *x, double *y, const int m);
  int ifft(double *x, double *y, const int m);
int fftr(double *x, double *y, const int m);
int ifftr(double *x, double *y, const int l);
  static int checkm(const int m);
int get_nextpow2(int n);
char *getmem(int leng, unsigned size);
double *dgetmem(int leng);
double get_process_time();



#endif
