#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

void nrerror(const char error_text[]);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);

#endif
