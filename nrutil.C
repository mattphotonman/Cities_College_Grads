#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#define NR_END 1
#define FREE_ARG char*
using namespace std;

//Numerical Recipes standard error handler
void nrerror(const char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

//Allocate a double vector with subscript range v[nl..nh]
double *vector(long nl, long nh)
{
  double *v;
  
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

//Free a double vector allocated with vector()
void free_vector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}
