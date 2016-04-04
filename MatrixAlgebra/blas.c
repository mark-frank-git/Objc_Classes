
/* in-line functions for use with 2D arrays:*/

/* various loop constructors */
#define DOFOR(i,to) for(i=0;i<to;i++)
#define DFOR(i,from,to) for(i=from-1;i<to;i++)
#define DOBY(i,from,to,by) for(i=from-1;i<to;i+=by)
#define DOBYY(i,from,to,by) for(i=from;i<t;i+=by)
#define DOBYYY(i,from,to) for(i=from;i<to;i++)
#define DOV(i,to,by) for(i=0;i<to;i+=by)

/* row major order as in C indices run 0 .. n-1 */
#define INDEX(i,j) [j+(i)*number_columns]

/* to index vectors starting with 1 */
#define VECTOR(i) [i-1]


/*
routines for linear systems processing via LU factorization
from Handbook of C tools for Scientists and engineers by L. Baker

SUPPORT ROUTINES FOR LU FACTORIZATION

isamax(n,sx,incx)  finds the location of the element
		of greatest absolute value in a vector sx of length n.
		Each incx-th element is examined ( hence, if sx is a
		2-d matrix, may be used to find largest elements in each 
		row or column, depending upon whether incx is 1 or n. 
		

saxpy(n,sa,sx,incx,sy,incy)	
		performs an elementary row operation sy= sy+sa sx where sx 
		and sy	are vectors and sa is a scalar.  Used to subtract 
		a scaled row sx from the sy row.  incx,incy are as in isamax.
		Vectors of length n.

sdot(n,sx,incx,sy,incy)	
		takes the dot product of 2 vectors sx and sy, of length n.

sswap(n,sx,incx,sy,incy)   
		exchanges two vectors sx and sy.  Used for row exchanges
		which occur during pivoting operation.

sscal(n,sa,sx,incx)   scale a vector  sx= sa sx where a is a 
		scalar and sx a vector.

sasum(n,sx,invx) 	function type float which returns the
		sum of the absolute values of the elements of a vector
*/

int isamax(int n, double *sx, int incx)
{int maxi,ix,i;
    double temp,smax;
/*returns 1 less than corresponding FORTRAN version*/
    if (n<=0)return -1;
    if(n==1)return 0;
/* ix=0*/
    maxi=0;
    smax=abs(sx[0]);
    ix=incx;/*ix=ix+incx=incx*/
    DFOR(i,2,n)
    { temp=abs(sx[ix]);
        if (temp>smax)
        {smax=temp;
            maxi=i;
/* return ith element as max,NOT subscript a[ix] ix=i*incx*/
        }
        ix+=incx;
    }
    return maxi;
}

void saxpy (int n, double sa, double *sx, int incx, double *sy, int incy)
{/*sy=sa*sx+sy*/
    int i,iy,ix;
    if(n<=0)return;
    if(sa==0.)return;

    iy=ix=0;
    if(incx<0) ix=incx*(1-n);
    if(incy<0) iy=incy*(1-n);
    DOFOR(i,n)
    {
        sy[iy]=sy[iy]+sa*sx[ix];
        iy+=incy;
        ix+=incx;
    }
    return;
}

double sdot(int n, double *sx, int incx, double *sy, int incy)
{double stemp;
    int i,ix,iy;
    if(n<=0)return(0.);
    ix=iy=0;
    stemp=0.0;
    if(incx<0) ix=incx*(1-n);
    if(incy<0) iy=incy*(1-n);
    DOFOR(i,n)
    {
        stemp+=sy[iy]*sx[ix];
        iy+=incy;
        ix+=incx;
    }
    return stemp;
}

void sswap(int n, double *sx, int incx, double *sy, int incy)
{
    int ix,iy,i;
    double t;
    if(n<=0)return;
    ix=iy=0;
    if(incx<0) ix=incx*(1-n);
    if(incy<0) iy=incy*(1-n);

    DOFOR(i,n)
    {
        t=sx [ix];
        sx [ix]= sy [iy];
        sy [iy]=t;
        ix+=incx;
        iy+=incy;
    }
    return;
}

void sscal(int n, double sa, double *sx, int incx)
{/*scale vector*/
    int i,nincx;

    if (n<=0) return;
    nincx=incx*n;
    DOV(i,nincx,incx)
    sx[i]=sx[i]*sa;
    return;
}

double sasum(int n, double *sx, int incx)
{/* ssum abs values*/
    int i,nincx;
    double stemp;
    stemp=0.0;
    nincx=n*incx;
    if (n<=0)return 0.0;
    DOV(i,n,nincx) stemp=stemp+abs(sx[i]);
    return (stemp);
}