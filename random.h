#ifndef INCLUDE_RANDOM
#define INCLUDE_RANDOM

double ran3(long idum=-1);
double gasdev(long idum=-1);
void gauss_corr_variables(double aves [], double** C, int Nvar, int Nsamp, double** rands);
void diagonalize_symm(double** A, int N, double lambda[], double** S);
void tred2(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);
double pythag(double a, double b);

#endif
