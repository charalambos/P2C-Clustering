#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include <limits>
#include "Matrix.h"

using namespace std;


#define XX 0
#define YY 1
#define ZZ 2

#define EMAX 0
#define EMID 1
#define EMIN 2
static short eig_sys3d(float a[][3], float d[])	{
	#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
	#define NO 0
	#define YES 1
/* Adopted From Numerical Recipes in C (2nd Edition) pp. 474-475 :
   Householder reduction of a real, symmetric matrix a[0..2][0..2]. On output,
   a is replaced by the orthogonal matrix Q effecting the transformation.
   d[0..2] returns the diagonal elements of the tridiagonal matrix,
   and e[0..2] the off-diagonal elements, with e[0] = 0. Several
   statements, as noted in comments, can be omitted if only eigenvalues are
   to be found, in which case a contains no useful information on output.
   Otherwise they are to be included. */

    static float scale,hh,h,g,f, e[3];

    // i = 2; l = 1;
    h = 0.0;
    scale = (float) (fabs(a[2][0])+fabs(a[2][1]));
    if (scale == 0.0) e[2] = a[2][1];
    else {
      a[2][0] /= scale;
      a[2][1] /= scale;
      h = a[2][0]*a[2][0]+a[2][1]*a[2][1];
      f = a[2][1];
      g = f>=0 ? ((float) -sqrt(h)) : ((float) sqrt(h));
      e[2] = scale*g;
      h -= f*g;
      a[2][1] = f-g;

      a[0][2] = a[2][0]/h;
      e[0] = (a[0][0]*a[2][0]+a[1][0]*a[2][1])/h;
      a[1][2] = a[2][1]/h;
      e[1] = (a[1][0]*a[2][0]+a[1][1]*a[2][1])/h;
      f = e[0]*a[2][0]+e[1]*a[2][1];

      hh = f/(h+h);
      f = a[2][0];
      e[0] = g = e[0]-hh*f;
      a[0][0] -= f*e[0]+g*a[2][0];
      f = a[2][1];
      e[1] = g = e[1]-hh*f;
      a[1][0] -= f*e[0]+g*a[2][0];
      a[1][1] -= f*e[1]+g*a[2][1];
    }
    d[2] = h;

    // i = 1; l = 0;
    e[1] = a[1][0];
    d[1] = 0.0;

    // Next statement can be omitted if eigenvectors not wanted
    // i=0; l=-1;
    d[0] = a[0][0];
    a[0][0] = 1.0;

    // i=1; l=0;
    if (d[1]) a[0][0] -= (a[1][0]*a[0][0])*a[0][1];
    d[1] = a[1][1];
    a[1][1] = 1.0;
    a[0][1]=a[1][0]=0.0;

    // i=2; l=1;
    if (d[2]) {
      g = a[2][0]*a[0][0]+a[2][1]*a[1][0];
      a[0][0] -= g*a[0][2];
      a[1][0] -= g*a[1][2];
      g = a[2][0]*a[0][1]+a[2][1]*a[1][1];
      a[0][1] -= g*a[0][2];
      a[1][1] -= g*a[1][2];
    }
    d[2] = a[2][2];
    a[2][2] = 1.0;
    a[0][2] = a[2][0] = 0.0;
    a[1][2] = a[2][1] = 0.0;

/* Adopted From Numerical Recipes in C (2nd Edition) pp. 480-481 :
   QL algorithm with implicit shifts, to determine the eigenvalues and
   eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
   symmetric matrix previously reduced by tred2.  On input, d[1..n] contains
   the diagonal elements of the tridiagonal matrix. On output, it returns the
   eigenvalues.  The vector e[1..n] inputs the subdiagonal elements of the
   tridiagonal matrix, with e[1] arbitrary.  On output e is destroyed.  When
   finding only the eigenvalues, several lines may be omitted, as noted in the
   comments.  If the eigenvectors of a tridiagonal matrix are desired, the
   matrix a[1..n][1..n] is input as the identity matrix.  If the eigenvectors
   of a matrix that has been reduced by tred2 are required, then z is input as
   the matrix output by tred2. In either case, the kth column of z returns the
   normalized eigenvector corresponding to d[k]. */

    static int m,l,iter,i,k;
    static float s,r,p,dd,c,b;

    e[0] = e[1];
    e[1] = e[2];
    e[2]=0.0;

    for (l=0;l<=2;l++) {
        iter=0;
        do {
            for (m=l;m<=1;m++) {
                dd= (float) (fabs(d[m])+fabs(d[m+1]));
                if (fabs(e[m])+dd == dd) break;
            }
            if (m != l) {
                if (iter++ == 50) {
                    std::cout << "Too many iterations in TQLI" << std::endl;
                    return(NO);
                }
                g=(d[l+1]-d[l])/(2.0f*e[l]);
                r= (float) sqrt((g*g)+1.0);
                g=d[m]-d[l]+e[l]/(g+(float) SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    if (fabs(f) >= fabs(g)) {
                        c=g/f;
                        r= (float) sqrt((c*c)+1.0);
                        e[i+1]=f*r;
                        c *= (s=1.0f/r);
                    } else {
                        s=f/g;
                        r= (float) sqrt((s*s)+1.0);
                        e[i+1]=g*r;
                        s *= (c=1.0f/r);
                    }
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0f*c*b;
                    p=s*r;
                    d[i+1]=g+p;
                    g=c*r-b;
                    // Next loop can be omitted if eigenvectors not wanted
                    for (k=0;k<=2;k++) {
                        f=a[k][i+1];
                        a[k][i+1]=s*a[k][i]+c*f;
                        a[k][i]=c*a[k][i]-s*f;
                    }
                }
                d[l]=d[l]-p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }

	// sort the eigenvalues and eigenvectors
	static float max_val;
	static int max_index;
	for(i=0;i<=1;i++) {
	  max_val = d[i];
	  max_index = i;
	  for(k=i+1;k<=2;k++)
		if (max_val < d[k]) { max_val = d[k]; max_index = k; }
	  if (max_index != i) {
		e[0] = d[i]; d[i] = d[max_index]; d[max_index] = e[0];
		e[0] = a[0][i]; a[0][i] = a[0][max_index]; a[0][max_index] = e[0];
		e[0] = a[1][i]; a[1][i] = a[1][max_index]; a[1][max_index] = e[0];
		e[0] = a[2][i]; a[2][i] = a[2][max_index]; a[2][max_index] = e[0];
	  }
	}

    return(YES);
}

#endif //__ALGEBRA_H__
