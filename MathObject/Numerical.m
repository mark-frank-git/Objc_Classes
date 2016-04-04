/************************************************************************
 * This subclass of object implements certain numerical math methods    *
 *                                                                      *
 * File:Numerical.m                                                     *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/26/92  - Started                                              *
 *  2. 11/23/93  - Add numerical integration                            *
 *  3. 10/11/95  - Check and return error in zeroIn:                    *
 ************************************************************************/
 
#import "Numerical.h"
#import <stdlib.h>
#import <stdio.h>

#define ABS(a)      ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)   ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)    ((a)>(b)?(a):(b))
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )


   
@implementation Numerical

- init
{
  double tol1;
  [super init];

/*******************
 * Find machine eps*
 *******************/
  machineEps = 1.;
  tol1 = 1.1;
  while(tol1 > 1.)
  {
    machineEps /= 2.;
    tol1 = 1. + machineEps;
  }
  return self;
}

   
/*###############################*
 * Find a zero in [ax, bx] of    *
 * a function supplied by the    *
 * sender.                       *
 *###############################*/
- (double) zeroIn:(double)ax    to:(double)bx    withTol:(double)tol
             from: (id)sender   error:(BOOL *)error
{
  int    begin_step, converged, bisection;
  int    sign_fa, sign_fb;
  double a, b, c, d, e, fa, fb, fc, tol1;
  double xm, p, q, r, s;

  if(![ sender respondsTo:sel_getUid("functionToFindRoot:") ])
  {
    printf("Can't call functionToFindRoot in zeroIn\n");
    return 0.;
  }
  
/********************
 * Initialization   *
 ********************/
  c = d = e = fc = 0.;
  a  = ax;
  b  = bx;
  fa = [sender functionToFindRoot:a];
  fb = [sender functionToFindRoot:b];

/********************
 * Check for error: *
 ********************/
  sign_fa   = SGN(fa);
  sign_fb   = SGN(fb);
  if(sign_fa == sign_fb)
  {
    *error  = YES;
    return  bx;
  }
  *error    = NO;
  
/*******************
 * Begin step:     *
 *******************/
  converged  = 0;
  begin_step = 1;
  while(!converged)
  {
    if(begin_step)
    {
      c  = a;
      fc = fa;
      d  = b - a;
      e  = d;
    }
    if( ABS(fc) < ABS(fb) )
    {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
  
/*******************
 * Convergence test*
 *******************/
  tol1 = 2.*machineEps*ABS(b) + 0.5*tol;
  xm   = 0.5*(c-b);
  if( (ABS(xm)<=tol1) || (fb == 0.) )
  {
    converged = 1;
    break;
  }

/********************
 * Bisection        *
 * necessary?       *
 ********************/
  bisection = 1;
  if( (ABS(e)>=tol1) && (ABS(fa)>ABS(fb)) )
  {
  
/********************
 * quadratic interp *
 * possible?        *
 ********************/
     if( a == c)
     {
/********************
 * linear interp.   *
 ********************/
        s = fb/fa;
        p = 2.*xm*s;
        q = 1. - s;
     }
     else
     {
/*******************
 * Inverse quadrat.*
 * interp.         *
 *******************/
       q = fa/fc;
       r = fb/fc;
       s = fb/fa;
       p = s*(2.*xm*q*(q-r) - (b-a)*(r-1.));
       q = (q-1.)*(r-1.)*(s-1.);
     }
/*******************
 * Adjust signs:   *
 *******************/
     q = (p > 0.) ? -q : q;
     p = ABS(p);
/*******************
 * Is interpolation*
 * acceptable?     *
 *******************/
     if(  ( 2.*p < (2.*xm*q - ABS(tol1*q)) )  &&
          ( p < ABS(0.5*e*q) )                 )
     {
       e = d;
       d = p/q;
       bisection = 0;
     }
     else
       bisection = 1;
   }
/*******************
 * Bisection       *
 *******************/
   if(bisection)
   {
     d = xm;
     e = d;
   }
/*******************
 * Complete step   *
 *******************/
     a  = b;
     fa = fb;
     if( ABS(d) > tol1)
       b += d;
     if( ABS(d) <= tol1)
       b += SIGN(tol1, xm);
     fb = [sender functionToFindRoot:b];
     if(fb*fc/ABS(fc) > 0.)
       begin_step = 1;
     else
       begin_step = 0;
   }
   return b;
}

#define NUMSEG  100
/*##################################*
  This a Simpson's integrator from  *
 * Quinn Curtis                 *
 *##################################*/
- (double) integrateFrom:(double) ax to:(double)bx from:(id)sender
{
  int    numseg;
  double xl, xh;
  double strtpnt, endpnt, area1;
  double area2, area, segwidth;

  if(![ sender respondsTo:sel_getUid("functionToIntegrate:") ])
  {
    printf("Can't call functionToIntegrate in simpsons:\n");
    return 0.;
  }

   xl = ax;
   xh = bx;
   numseg = NUMSEG;
   segwidth = (xh - xl) / numseg;
   if(segwidth<=machineEps) 
    return 0.;
   endpnt = xh;
   area1 = 0.0;
   area2 = 0.0;
   area = 0.0;
   if ( (numseg%2) != 0 ) {
      area1 = 3.0 / 8.0 * segwidth * 
             ([sender functionToIntegrate:(endpnt - 3.0 * segwidth)] +
              3.0 * [sender functionToIntegrate:(endpnt - 2.0 * segwidth)] +
              3.0 * [sender functionToIntegrate:(endpnt - segwidth)] + 
                    [sender functionToIntegrate:(endpnt)]
             );
      endpnt = endpnt - 3.0 * segwidth;
   }
   else
      area1 = 0.0;

   if ( numseg != 3 ) {
      strtpnt = xl;
      while ( (strtpnt < endpnt - segwidth) )
      {
         area2 += 1.0 / 3.0 * (segwidth * 
                 ([sender functionToIntegrate:(strtpnt)] + 
                  4.0 *[sender functionToIntegrate:(strtpnt + segwidth)] + 
                        [sender functionToIntegrate:(strtpnt + 2.0 * segwidth)]
                 ));
         strtpnt += 2.0 * segwidth;
     }
   }
   area = area1 + area2;
   return area;
}

/*######################################*
 * This a Simpson's integrator from     *
 * Quinn Curtis.  It facilitates        *
 * two dimensional integration by       *
 * calling 'outsideFunctionToIntegrate' *
 * in sender.  This function must do    * 
 * the inner integration.               *
 *######################################*/
- (double) twoDimIntegrateFrom:(double) ax to:(double)bx from:(id)sender
{
  int    numseg;
  double xl, xh;
  double strtpnt, endpnt, area1;
  double area2, area, segwidth;

  if(![ sender respondsTo:sel_getUid("outsideFunctionToIntegrate:") ])
  {
    printf("Can't call functionToIntegrate in simpsons:\n");
    return 0.;
  }

   xl = ax;
   xh = bx;
   numseg = NUMSEG;
   segwidth = (xh - xl) / numseg;
   if(segwidth<=machineEps) 
    return 0.;
   endpnt = xh;
   area1 = 0.0;
   area2 = 0.0;
   area = 0.0;
   if ( (numseg%2) != 0 ) {
      area1 = 3.0 / 8.0 * segwidth * 
             ([sender outsideFunctionToIntegrate:(endpnt - 3.0 * segwidth)] +
              3.0 * [sender outsideFunctionToIntegrate:(endpnt - 2.0 * segwidth)] +
              3.0 * [sender outsideFunctionToIntegrate:(endpnt - segwidth)] + 
                    [sender outsideFunctionToIntegrate:(endpnt)]
             );
      endpnt = endpnt - 3.0 * segwidth;
   }
   else
      area1 = 0.0;

   if ( numseg != 3 ) {
      strtpnt = xl;
      while ( (strtpnt < endpnt - segwidth) )
      {
         area2 += 1.0 / 3.0 * (segwidth * 
                 ([sender outsideFunctionToIntegrate:(strtpnt)] + 
                  4.0 *[sender outsideFunctionToIntegrate:(strtpnt + segwidth)] + 
                        [sender outsideFunctionToIntegrate:(strtpnt + 2.0 * segwidth)]
                 ));
         strtpnt += 2.0 * segwidth;
     }
   }
   area = area1 + area2;
   return area;
}
      
@end