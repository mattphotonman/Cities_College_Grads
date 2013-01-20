#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "random.h"
#include "nrutil.h"

//Definitions for ran3
//According to Knuth, any large MBIG, and any smaller (but still large) MSEED
//can be substituted for the above values.
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

//Definitions for tqli
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

using namespace std;

//Original NR description:  Returns a uniform random deviate between 0.0 and 1.0.  
//Set idum to any negative value to initialize or reinitialize the sequence.
//My modifications:  I've changed it so that reinitialize isn't possible
//anymore.  I never use this and it makes the code simpler.  It also alleviates
//the problems I had putting long idum=-1 in the header file.  In my version
//you can specify a seed in the argument to ran3 (the default is -1), but
//this only does anything the very first time you call ran3 in your main
//program, which is the only time that things are initialized.
double ran3(long idum)
{
  static int inext, inextp;
  static long ma[56];         //The value 56 (range ma[1..55]) is special and
                              //should not be modified; see Knuth.
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (iff == 0) {   //Initialization
    iff=1;
    mj=labs(MSEED-labs(idum));  //Initialize ma[55] using the seed idum and the
                                 //large number MSEED.
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1; i<=54; i++) {      //Now initialize the rest of the table, in a
                                 //slightly random order, with numbers that
                                 //are not especially random.
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1; k<=4; k++)         //We randomize them by "warming up the
                                 //generator".
      for (i=1; i<=55; i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;      //Prepare indices for our first generated number.
    inextp=31;    //The constant 31 is special; see Knuth.
  }
  //Here is where we start, except on initialization.
  if (++inext == 56) inext=1;    //Increment inext and inextp, wrapping around
                                 //56 to 1.
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];       //Generate a new random number subtractively.
  if (mj < MZ) mj += MBIG;       //Be sure that it is in range.
  ma[inext]=mj;                  //Store it,
  return mj*FAC;                 //and output the derived uniform deviate.
}

//Original NR description: Returns a normally distributed deviate with zero mean 
//and unit variance, using ran3(idum) as the source of uniform deviates. 
//My modifications: Don't allow for reinitialization.  See comments for ran3
//function.  The function gasdev still has an argument idum which it passes
//to ran3, but this is just a seed for ran3 if it's the first time ran3 has been
//called, otherwise it does nothing.
double gasdev(long idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (iset == 0) {
    //We don't have an extra deviate handy, so pick two uniform numbers
    //in the square extending from -1 to +1 in each direction, see if 
    //they are in the unit circle, and if they are not, try again.
    do {
      v1=2.0*ran3(idum)-1.0;
      v2=2.0*ran3(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);	
    fac=sqrt(-2.0*log(rsq)/rsq);
    //Now make the Box-Muller transformation to get two normal deviates.
    //Return one and save the other for next time.
    gset=v1*fac;
    iset=1;  //Set flag.
    return v2*fac;
  } else {
    //We have an extra deviate handy, so unset the flag, and return it.
    iset=0;
    return gset;
  }
}

//Given the means (aves) and covariance matrix (C) for Nvar random variables, produces
//Nsamp samples of the Nvar variables, which are returned to the array rands.
//The different samples are not correlated to each other.
void gauss_corr_variables(double aves [], double** C, int Nvar, int Nsamp, double** rands)
{
  //Assumed array sizes.
  //Inputs:
  //double aves[Nvar]
  //double C[Nvar][Nvar]
  //Output:
  //double rands[Nvar][Nsamp]

  //Diagonalize C.
  double lambda[Nvar];
  double** S= new double* [Nvar];
  for (int i=0; i<Nvar; i++)
    S[i]=new double [Nvar];
  diagonalize_symm(C,Nvar,lambda,S);

  //Check that all the eigenvalues of C are non-negative.
  for (int i=0; i<Nvar; i++)
    if (lambda[i]<0) {
      cout << "gauss_corr_variables: Error, negative eigenvalue of covariance matrix found.  Abort.\n";
      exit(1);
    }
  
  double u[Nvar];
  for (int ii=0; ii<Nsamp; ii++) {
    //For each sample, first generate random variables u[k],
    //where
    //<u[k]u[l]> = lambda[k]*delta[k][l]  (no sum on k)
    //i.e. the u[k] are independent gaussian random variables
    //with mean 0 and standard deviation sqrt(lambda[k]).
    for (int k=0; k<Nvar; k++)
      u[k]=gasdev()*sqrt(lambda[k]);

    //If we define x[i]=sum(S[i][k]*u[k],k=0..Nvar-1), where
    //S is the matrix we found that diagonalizes C, then
    //<x[i]x[j]> = C[i][j]
    //but the x[i] still have 0 mean.  Thus we want to take
    //rands[i][ii]=x[i]+aves[i].
    for (int i=0; i<Nvar; i++) {
      rands[i][ii]=0.0;
      for (int k=0; k<Nvar; k++)
	rands[i][ii]+=S[i][k]*u[k];
      rands[i][ii]+=aves[i];
    }
  }

  
  //Delete all memory allocated using "new" in this function
  for (int i=0; i<Nvar; i++)
    delete [] S[i];
  delete [] S;
  
}

//Returns the eigenvalues and eigenvectors of the symmetric 
//matrix A[0..N-1][0..N-1].
//Eigenvalues are returned in the array lambda[0..N-1] and normalized
//eigenvectors in S[0..N-1][0..N-1].  This function is just a wrapper
//around nr functions tred2 and tqli, that also translates the index
//numbering to nr conventions.
void diagonalize_symm(double** A, int N, double lambda[], double** S)
{
  //Check that A is symmetric.
  for (int i=0;i<N;i++)
    for (int j=i+1;j<N;j++)
      if (A[i][j]!=A[j][i]) {
	cout << "diagonalize_symm:  Error, matrix not symmetric.  Abort.\n";
	exit(1);
      }
  
  //Put matrix A into matrix "a" which has the indexing for nr
  //functions and which will be destroyed.
  double** a=new double* [N+1];
  a[0]=0;
  for (int i=1;i<N+1;i++) {
    a[i]=new double [N+1];
    for (int j=1;j<N+1;j++)
      a[i][j]=A[i-1][j-1];
  }
  
  double d[N+1], e[N+1];
  
  tred2(a,N,d,e);
  tqli(d,e,N,a);
  
  for (int i=0;i<N;i++) {
    lambda[i]=d[i+1];
    for (int j=0;j<N;j++)
      S[i][j]=a[i+1][j+1];
  }
  
  //Have to delete matrix "a" since it was declared locally.
  for (int i=1;i<N+1;i++)
    delete [] a[i];
  delete [] a;
}

//Householder reduction of a real, symmetric matrix a[1..n][1..n]. On
//output, a is replaced by the orthogonal matrix Q effecting the 
//transformation. d[1..n] returns the diagonal elements of the 
//tridiagonal matrix, and e[1..n] the off-diagonal elements, with e[1]=0. 
//Several statements, as noted in comments, can be omitted if only 
//eigenvalues are to be found, in which case a contains no useful 
//information on output. Otherwise they are to be included.
void tred2(double **a, int n, double d[], double e[])
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) { 
    l=i-1;
    h=scale=0.0; 
    if (l > 1) {
      for (k=1;k<=l;k++) 
	scale += fabs(a[i][k]);
      if (scale == 0.0) //Skip transformation.
	e[i]=a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale; //Use scaled a's for transformation.
	  h += a[i][k]*a[i][k]; //Form sigma in h.
	} 
	f=a[i][l]; 
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h)); 
	e[i]=scale*g; 
	h -= f*g;  //Now h is equation (11.2.4). 
	a[i][l]=f-g;	//Store u in the ith row of a. 
	f=0.0; 
	for (j=1;j<=l;j++) { 
	  //Next statement can be omitted if eigenvectors not wanted
	  a[j][i]=a[i][j]/h; //Store u/H in ith column of a.
	  g=0.0; //Form an element of A · u in g.
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k]; 
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k]; 
	  e[j]=g/h; //Form element of p in temporarily unused element of e.
	  f += e[j]*a[i][j];
	} 
	hh=f/(h+h); //Form K, equation (11.2.11).
	for (j=1;j<=l;j++) { //Form q and store in e overwriting p.
	  f=a[i][j]; 
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++) //Reduce a, equation (11.2.13).
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  //Next statement can be omitted if eigenvectors not wanted
  d[1]=0.0;
  e[1]=0.0;
  //Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i];
  for (i=1;i<=n;i++) { //Begin accumulation of transformation matrices.
    l=i-1;
    if (d[i]) { //This block skipped when i=1.
      for (j=1;j<=l;j++) {
	g=0.0;
	for (k=1;k<=l;k++) //Use u and u/H stored in a to form P·Q.
	  g += a[i][k]*a[k][j]; 
	for (k=1;k<=l;k++)
	  a[k][j] -= g*a[k][i];
      }
    } 
    d[i]=a[i][i];  //This statement remains.
    a[i][i]=1.0;  //Reset row and column of a to identity matrix for next iteration.
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

//QL algorithm with implicit shifts, to determine the eigenvalues and 
//eigenvectors of a real, symmetric, tridiagonal matrix, or of a real, 
//symmetric matrix previously reduced by tred2 §11.2. On input, d[1..n] 
//contains the diagonal elements of the tridiagonal matrix. On output, it 
//returns the eigenvalues. The vector e[1..n] inputs the subdiagonal 
//elements of the tridiagonal matrix, with e[1] arbitrary. On output e is 
//destroyed. When finding only the eigenvalues, several lines may be 
//omitted, as noted in the comments. If the eigenvectors of a tridiagonal 
//matrix are desired, the matrix z[1..n][1..n] is input as the identity 
//matrix. If the eigenvectors of a matrix that has been reduced by tred2 
//are required, then z is input as the matrix output by tred2. In either 
//case, the kth column of z returns the normalized eigenvector 
//corresponding to d[k].
void tqli(double d[], double e[], int n, double **z)
{
  int m,l,iter,i,k; 
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];  //Convenient to renumber the elements of e.
  e[n]=0.0; 
  for (l=1;l<=n;l++) {
    iter=0; 
    do {
      for (m=l;m<=n-1;m++) {  //Look for a single small subdiagonal element to split the matrix.
	dd=fabs(d[m])+fabs(d[m+1]); 
	if ((double)(fabs(e[m])+dd) == dd) break;	
      } 
      if (m != l) {
	if (iter++ == 30) nrerror("Too many iterations in tqli"); 
	g=(d[l+1]-d[l])/(2.0*e[l]);  //Form shift. 
	r=pythag(g,1.0); 
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));	 //This is dm - ks. 
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) { //A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
	  f=s*e[i];b=c*e[i]; 
	  e[i+1]=(r=pythag(f,g)); 
	  if (r == 0.0) { //Recover from underflow.
	    d[i+1] -= p; 
	    e[m]=0.0; 
	    break;
	  } 
	  s=f/r; 
	  c=g/r; 
	  g=d[i+1]-p; 
	  r=(d[i]-g)*s+2.0*c*b; 
	  d[i+1]=g+(p=s*r); 
	  g=c*r-b; 
	  //Next loop can be omitted if eigenvectors not wanted
	  for (k=1;k<=n;k++) { //Form eigenvectors.
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f; 
	    z[k][i]=c*z[k][i]-s*f;
	  }
	} 
	if (r == 0.0 && i >= l) continue; 
	d[l] -= p; 
	e[l]=g; 
	e[m]=0.0;
      }
    } while (m != l);
  }
}

//Computes (a^2 + b^2)^(1/2) without destructive underflow or overflow.
double pythag(double a, double b)
{
  double absa,absb; 
  absa=fabs(a); 
  absb=fabs(b); 
  if (absa > absb) 
    return absa*sqrt(1.0+(absb/absa)*(absb/absa)); 
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}
