
#include "c_headers.h"
/*
routines for linear systems processing via LU factorization
from Handbook of C tools for Scientists and engineers by L. Baker

SUPPORT ROUTINES FOR LU FACTORIZATION

izamax(n,sx,incx)  finds the location of the element
		of greatest absolute value in a vector sx of length n.
		Each incx-th element is examined ( hence, if sx is a
		2-d matrix, may be used to find largest elements in each
		row or column, depending upon whether incx is 1 or n

zaxpy(n,sa,sx,incx,sy,incy)
		performs an elementary row operation sy= sy+sa sx where sx
		and sy	are vectors and sa is a scalar.  Used to subtract
		a scaled row sx from the sy row.  incx,incy are as in isamax.
		Vectors of length n.

zscal(n,sa,sx,incx)   scale a vector  sx= sa sx where a is a
		scalar and sx a vector.

the above are all based upon BLAS  routines

*/

int izamax(int n,COMPLEX *sx, int incx)
{
  int maxi,ix,i;
  double temp,smax;
/*returns 1 less than corresponding FORTRAN version*/
  if (n<=0)return -1;
  if(n==1)return 0;
/* ix=0*/
  maxi=0;
  smax=CABS2(sx[0]);
  ix=incx;/*ix=ix+incx=incx*/
  DFOR(i,2,n)
  {
    temp=CABS2(sx[ix]);
    if (temp>smax)
    {
      smax=temp;
      maxi=i;
/* return ith element as max,NOT subscript a[ix] ix=i*incx*/
    }
    ix+=incx;
  }
  return maxi;
}

void zaxpy (int n, COMPLEX sa, COMPLEX *sx, int incx, COMPLEX *sy, int incy)
{/*sy=sa*sx+sy*/
  int i,iy,ix;
  COMPLEX ctemp;

  if(n<=0)return;
  if(CABS2(sa) == 0.)return;

  iy=ix=0;
  if(incx<0) ix=incx*(1-n);
  if(incy<0) iy=incy*(1-n);
  DOFOR(i,n)
  {
    CMULT(ctemp, sa, sx[ix]);
    sy[iy].x = sy[iy].x + ctemp.x;
    sy[iy].y = sy[iy].y + ctemp.y;
    iy+=incy;
    ix+=incx;
  }
  return;
}

void zscal(int n, COMPLEX sa, COMPLEX *sx, int incx)
{/*scale vector*/
  int i,nincx;
  COMPLEX ctemp;

  if (n<=0) return;
  nincx=incx*n;
  DOV(i,nincx,incx)
  {
    CMULT(ctemp, sx[i], sa);
    sx[i].x = ctemp.x;
    sx[i].y = ctemp.y;
  }
  return;
}