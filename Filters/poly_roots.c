/************************************************************************
 *																		*
 * This file finds the roots of a polynomial with complex coefficients	*
 *																		*
 * File: poly_roots.c													*
 *																		*
 *  Rev History:														*
 *    8/31/90  - Started												*
 *																		*
 * This is a Jenkins-Traub polynomial root finder from Handbook of		*
 * C tools for scientists and engineers by L. Baker						*
 ************************************************************************/
/**************************
 * Include files:         *
 **************************/
#include "c_headers.h"
#include <stdio.h>   
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
 
/**************************
 * Constants:             *
 **************************/
#define MAXDEG  50


double baker;

/**************************
 * Function prototypes:   *
 **************************/
void   jt(COMPLEX *coeff, int degree, COMPLEX *ans, int scale, int *fail);
void   printc(COMPLEX *z);
double errev(int nn, COMPLEX *q, double ms, double mp, double are, 
                 double mre);
double cauchy(int nn, double *pt, double *q);
void   calct(int *boolean);
int    vrshft(int l3, COMPLEX *z, int *conv);
void   fxshft(int l2, COMPLEX *z, int *conv);
void   noshft( int l1);
void   polyev(int nn, COMPLEX *s, COMPLEX *p, COMPLEX *q, COMPLEX *pv);
void   mcon(double *eta, double *infin, double *smalno, double *base);
double cmod(COMPLEX *x);
void   cdivid(COMPLEX *a, COMPLEX *b, COMPLEX *c);
double scale(int nn, double *pt, double *eta, double *infin, double *smalno,
               double *base);
void nexth(int *boolean);


                  
/*******************************************************************
 *                                                                 *
 * void printc(COMPLEX *z)                                         *
 *******************************************************************/
void printc(COMPLEX *z)
{
  printf("%f %f",z->x,z->y);
}

/*******************************************************************
 *                                                                 *
 *  void jt(COMPLEX *op, int degree, COMPLEX *zero, int scale,     *
 *         int *fail)                                              *
 *                                                                 *
 *******************************************************************/
struct complex p[MAXDEG],h[MAXDEG],qp[MAXDEG],qh[MAXDEG],sh[MAXDEG];
struct complex s,t,pv;
double are,mre,eta,infin,smalno,base,tempa[MAXDEG],tempb[MAXDEG];
double xx,yy,cosr,sinr;
int nn;

void jt(COMPLEX *op, int degree, COMPLEX *zero, int accuracy, int *fail)
{
double bnd,xxx;
int i,idnn2,cntl1,cntl2,conv;
struct complex z;

baker = (double)accuracy;
mcon(&eta,&infin,&smalno,&base);
are=eta;/* factor of 2 my addition*/
mre=2.*sqrt(2.)*eta;
xx=.70710678;
yy = -xx;
cosr = -.069756474;
sinr = .99756405;
*fail=0;
nn=degree+1;

CMPLX(t,-99.,-99.);/* is t initially undefined in fxshft?*/
if (op[0].x==0. && op[0].y==0.)
	{
	*fail=1;
	return;
	}
for(; (op[nn-1].x==0. && op[nn-1].y==0.)&&nn>=0 ;nn--)
	{
	printf(" zero constant term found\n");
	idnn2=degree-nn+1;
	CMPLX(zero[idnn2],0.,0.);
	}
for(i=0;i<nn;i++)
	{
	CLET(p[i],op[i]);
	tempa[i]=cmod(&(p[i]));	
	}
bnd=scale(nn,tempa,&eta,&infin,&smalno,&base);
if(bnd!=1.)
	{
	for(i=0;i<nn;i++)
		{
		CTREAL( (p[i]),(p[i]),bnd);
		}
	}

findzero:
if(nn==1)return;
if(nn<=2)
	{
	cdivid(&p[1],&p[0],&zero[degree-1]);
	CTREAL( zero[degree-1],zero[degree-1],-1.);
	return;
	}
for(i=0;i<nn;i++)
	{
	tempa[i]=cmod(&(p[i]));
	}
bnd=cauchy(nn,tempa,tempb);
for(cntl1=0;cntl1<=1;cntl1++)
	{
	noshft(5);
	for(cntl2=1;cntl2<=9;cntl2++)
		{
		xxx=cosr*xx-sinr*yy;
		yy=sinr*xx+cosr*yy;
		xx=xxx;
		s.x=bnd*xx;
		s.y=bnd*yy;
		fxshft(10*cntl2,&z,&conv);
		if(conv)
			{
			idnn2=degree-nn+1;
			CLET(zero[idnn2],z);
			nn--;
			for(i=0;i<nn;i++)
				{CLET(p[i],qp[i]);
				}
			goto findzero;
			}
		}
	}
*fail=1;
return;
}

/**********************************************************************
 * double errev(int nn, COMPLEX *q, double ms, double mp, double are, *
 *               double mre)                                          *
 *                                                                    *
 **********************************************************************/

double errev(nn,q,ms,mp,are,mre)
int nn;
struct complex *q;
double ms,mp,are,mre;
{
double ans,e,cmod();
int i;
e=cmod(&(q[0])) * mre/(are+mre);
for(i=0;i<nn;i++)
	{
	e=e*ms+cmod(&(q[i]));
	}
ans=e*(are+mre)-mp*mre;
return(ans);
}

void nexth(int *boolean)
{
int n,j;
n=nn-1;
if(*boolean)
	{
	for(j=1;j<n;j++)
		{
		CLET(h[j],qh[j-1]);
		}
	CMPLX(h[0],0.,0.);
	return;
	}	
/*else*/
for(j=1;j<n;j++)
	{
	CMULT(h[j],t,qh[j-1])
	CADD(h[j],h[j],qp[j])
	}
CLET(h[0],qp[0]);
return;
}

/**********************************************************************************
 *                                                                                *
 *  double cauchy(int nn, double *pt, double *q)                                  *
 *                                                                                *
 **********************************************************************************/
double cauchy(int nn, double *pt, double *q)     
{
int n,i,nm;
double x,dx,df, f,xm;
pt[nn-1]*=-1.;
n=nn-1;
nm=n-1;
x= exp(log(-pt[n]))-log(pt[0])/((double)n);
if(pt[nm]!=0.)
   {xm = -pt[n]/pt[nm];
	x=MIN(xm,x);
	}
repeat:
xm=x*.1;
f=pt[0];
for(i=1;i<nn;i++)
	f=f*xm+pt[i];
if(f>0.)
	{
	x=xm;
	goto repeat;
	}
dx=x;
while( abs(dx/x)> .005)
	{
	q[0]=pt[0];
	for (i=1;i<nn;i++)
		q[i]=q[i-1]*x+pt[i];
	f=q[n];
	df=q[0];
	for (i=1;i<n;i++)
		df=df*x+q[i];
	dx=f/df;
		x-=dx;
        }
return(x);
}
 
/**************************************************************************
 *                                                                        *
 *  void calct(int *boolean)                                              *
 *                                                                        *
 **************************************************************************/
void calct(int *boolean)
{
int n,nm1;
struct complex hv;
n=nn-1;   nm1=n-1;
polyev(n,&s,h,qh,&hv);
*boolean=  (cmod(&hv)<= are*10.*cmod(&(h[nm1])) );
if(*boolean)
	{
	CMPLX(t,0.,0.);	
	return;
	}
cdivid(&pv,&hv,&t);
CTREAL(t,t,-1.);
return;
}   


/*****************************************************************************
 *                                                                           *
 * int vrshft(int l3, COMPLEX *z, int *conv)                                 *
 *                                                                           *
 *****************************************************************************/
int vrshft(int l3, COMPLEX *z, int *conv)
{
double mp,ms,omp,relstp,r1,r2,tp;
int boolean,b,j,i;
b=0;
*conv=0;
omp  = relstp = 0.;
CASSN(s,z);
for(i=0;i<l3;i++)
	{
	polyev(nn,&s,p,qp,&pv);
	mp=cmod(&pv);
	ms=cmod(&s);
	if(mp<=(20.*errev(nn,qp,ms,mp,are,mre)) )
		{
		*conv=1;
		CSET(z,s);
		return 0;
		}	
	if(i!=0)
		{
		if( !b && !(mp<omp) && (relstp<.05))
				{
				tp=relstp;
				b=1;
				if(relstp<eta)tp=eta;
				r1=sqrt(tp);
				r2=s.x*(1.+r1)-s.y*r1;
				s.y=s.x*r1+s.y*(1.+r1);
				s.x=r2;
				polyev(nn,&s,p,qp,&pv);
				for(j=0;j<5;j++)
						{
						calct(&boolean);
						nexth(&boolean);	
						}
				omp=infin;
				goto skp;
				}
		if(mp*.1 >omp)return 0;
		}
	omp=mp;
skp:
	calct(&boolean);
	nexth(&boolean);
	calct(&boolean);
	if(!boolean)
		{
		relstp=cmod(&t)/cmod(&s);
		CADD(s,s,t);
		}
	}
return 0;
}
    
/*******************************************************************************
 *                                                                             *
 * void fxshft(int l2, COMPLEX *z, int *conv)                                  *
 *                                                                             *
 *******************************************************************************/
void fxshft(int l2, COMPLEX *z, int *conv)
{
int i,j,n,test,pasd,boolean;
struct complex svs,ot,lou;
n=nn-1;
polyev(nn,&s,p,qp,&pv);
test=1;
pasd=0;
calct(&boolean);
for(j=0;j<l2;j++)
	{
	CLET(ot,t);
	nexth(&boolean);
	calct(&boolean);
	CADD((*z),s,t);
	if( (!boolean)&&test && (j!=(l2-1))  )
			{
			CSUB(lou,t,ot);
			if( cmod(&lou)>= .5*cmod(z) )
					{pasd=0;}
			else if (!pasd) {pasd=1;}
			else
				{
				for(i=0;i<n;i++)
						{
						CLET(sh[i],h[i]);
						}
				CLET(svs,s);
				vrshft(10,z,conv);
				if(*conv)
					{
					return;
					}
				test=0;
				for(i=0;i<n;i++)
						{
						CLET(h[i],sh[i]);
						}
				CLET(s,svs);
				polyev(nn,&s,p,qp,&pv);
				calct(&boolean);
				}
			}
		else
			{
			}
	}
vrshft(10,z,conv);
return;
} 

/*********************************************************************************
 *                                                                               *
 *  void noshft( int l1)                                                         *
 *                                                                               *
 *********************************************************************************/
void noshft( int l1)
{
int n,nm1,i,j,jj;
double xni;
n=nn-1;
nm1=n-1;
for(i=0;i<n;i++)
	{
	xni=((double)(n-i))/((double) n);
	CTREAL(h[i],p[i],xni);
	}
for(jj=0;jj<l1;jj++)
	{
	if(cmod(&(h[nm1]))> eta*10.*cmod(&(p[nm1])))
		{
		cdivid(&(p[n]),&(h[nm1]),&t);
		CTREAL(t,t,-1.);
		for(i=1;i<=nm1;i++)
			{
			j=nn-i-1;
		/*	t1=h[j-1].x;
			t2=h[j-1].y;
			h[j].x=t.x*t1-t.y*t2;
			h[j].y=t.x*t2+t.y*t1;*/
			CMULT(h[j],t,h[j-1]);
			CADD(h[j],h[j],p[j]);
			}
		CLET(h[0],p[0]);
		}
	else
		{
		for(i=1;i<=nm1;i++)
			{
			j=nn-i-1;
			CLET(h[j],h[j-1]);
			}
		CMPLX(h[0],0.,0.);
		}
	}
return;
}

/************************************************************************
 *                                                                      *
 * void polyev(int nn, COMPLEX *s, COMPLEX *p, COMPLEX *q, COMPLEX *pv) *
 *                                                                      *
 ************************************************************************/
void polyev(int nn, COMPLEX *s, COMPLEX *p, COMPLEX *q, COMPLEX *pv)
{/* nested polynomial evaluation*/
int i;
struct complex temp;
CLET(q[0],p[0]);
CLET((*pv),q[0]);
for(i=1;i<nn;i++)
	{
	CMULT(temp,(*pv),(*s));
	CADD((*pv),temp,(p[i]));
	CASSN(q[i],pv);
	}
return;
}

/************************************************************************
 *                                                                      *
 *  void mcon(double *eta, double *infin, double *smalno, double *base) *
 *                                                                      *
 ************************************************************************/
void mcon(double *eta, double *infin, double *smalno, double *base) 
{
/* will compute eta*/
for(*eta=1.; (1.+(*eta))>1.;(*eta)*=.5);
*eta *= baker;


*smalno=(1.e-30);
/* set by hand as some compilers can't handle overflows/underflows*/
*infin=(1.e30);
*base=2.;

return;
}

/***************************************************************************
 *                                                                         *
 *  double cmod(COMPLEX *x)                                                *
 *                                                                         *
 ***************************************************************************/
double cmod(COMPLEX *x)    
{/* robust but expensive version of cabs()*/
double ans,rpart,ipart,aux;
rpart= x->x; ipart= x->y;
rpart= abs(rpart);
ipart= abs(ipart);
if(rpart>ipart)
	{
	aux= ipart/rpart;
	ans=rpart*sqrt(1.+ aux*aux);
	return(ans);
	}
else if (rpart<ipart)
	{
	aux=rpart/ipart;
	ans=ipart*sqrt(1.+aux*aux);
	return(ans);
	}
/*else*/
return (rpart*sqrt(2.));
}

/******************************************************************************
 *                                                                            *
 *  void cdivid(COMPLEX *a, COMPLEX *b, COMPLEX *c)                           *
 *                                                                            *
 ******************************************************************************/
void cdivid(COMPLEX *a, COMPLEX *b, COMPLEX *c)
{/* robust c=a/b*/
double r,d,dummy,infin;
if( b->x==0. && b->y==0.)
	{
	mcon(&dummy,&infin,&dummy,&dummy);
	c->x=infin;
	c->y=infin;
	return;
	}
if (abs(b->x)>=abs(b->y))
	{
	r=b->y/b->x;
	d=1./(b->x + r*b->y);
	c->x=(a->y *r + a->x)*d;
	c->y=(-a->x *r +a->y)*d;
	}
else
	{
	r=b->x/b->y;
	d=1./(b->y + r*b->x);
	c->x=(a->x *r + a->y)*d;
	c->y=(a->y *r - a->x)*d;
	}
return;
}


/*********************************************************************************
 *                                                                               *
 *  double scale(int nn, double *pt, double *eta, double *infin, double *smalno, *
 *               double *base)                                                   *
 *                                                                               *
 *********************************************************************************/
double scale(int nn, double *pt, double *eta, double *infin, double *smalno,
             double *base)
{
int i;
double scal,maxx,minn,hi,lo,x,sc;
hi=sqrt(*infin);
lo = *smalno/(*eta);
maxx=0.;
minn = *infin;
for(i=0;i<nn;i++)
	{
	x=pt[i];
	maxx= MAX(x,maxx);
	if(x!=0.)minn=MIN(x,minn);
	}
scal=1.;
if(minn>=lo && maxx<=hi)
	{
	return(scal);
	}
x=lo/minn;
if(x>1.)
	{
	sc=x;
	if( *infin/sc >maxx) sc=1.;
	}
else
	{
	sc=1./(sqrt(maxx)*sqrt(minn));
	}
printf(" before log,pow %e\n",sc);
i= (int) (log(sc)/log(*base)+.5);
scal= pow(*base, ((double) i));
return(scal);
}
