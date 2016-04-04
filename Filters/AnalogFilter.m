/************************************************************************
 * This subclass of object implements an analog filter object           *
 *                                                                      *
 * File:AnalogFilter.                                                   *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/12/92  - Started                                              *
 *  2. 12/22/93  - Added MAX_ORDER                                      *
 *  3. 01/03/95  - remove complex_response array in                     *
 *                 analogResponseAtFrequencies.                         *
 ************************************************************************/

#import "AnalogFilter.h"
#import "Polynomial.h"
#import "Numerical.h"
#import <stdlib.h>
#import <math.h>

#define EPSILON     1.                          /* for chebyshev */
#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )
#define MAX_ORDER   10                          /* max prototype order */

@implementation AnalogFilter

- init
{
  [super init];
  
// init instance variables:
  filterPassType    = LOW_PASS;
  filterStructureType   = TRANSFER_FUNCTION;
  analogType        = BUTTERWORTH;
  responseType      = MAGNITUDE;
  filterOrder       = 2;
  passBandGain      = 1.;
  realAxisPoles     = NO;
  numberPoles       = numberZeros = filterOrder;
  analogPoles = analogZeros = NULL;
  aPolyObject = [[Polynomial alloc] init];
  bPolyObject = [[Polynomial alloc] init];
  
  return self;
}


- initWithRealCoeffs:(double *)aCoeffs b:(double *)bCoeffs order:(int)order
{
  [self init];
  filterOrder = order;
  numberPoles  = numberZeros = filterOrder;
  
  [self initPolesZeros];
  [aPolyObject setAPoly:aCoeffs order:order];
  [bPolyObject setAPoly:bCoeffs order:order];

  return self;
}
  
- initWithPassType:(int)type fo:(double)centerFreq fc:(double)cutoffFreq
            order:(int)order
{
  [self init];

  filterPassType    = type;
  filterOrder       = order;
  numberPoles       = numberZeros = filterOrder;
  fo                = centerFreq;
  fc                = cutoffFreq;
  [self initPolesZeros];
  
  return self;
}

- initPolesZeros
{
  if(analogPoles != NULL)
  {
    free(analogPoles);
    free(analogZeros);
  }
  analogPoles   = (COMPLEX *) malloc(2*numberPoles*sizeof(COMPLEX));
  analogZeros   = (COMPLEX *) malloc(2*numberPoles*sizeof(COMPLEX));
  return self;
}

/*##############################*
 * Set Parameters:      *
 *##############################*/
- setPassType:(int)type {filterPassType = type; return self;}
- setFilterStructureType:(int)type {filterStructureType = type; return self;}
- setAnalogType:(int)type {analogType = type; return self;}
- setResponseType:(int)type {responseType = type; return self;}
- setFilterOrder:(int)order
{
  filterOrder = order;
  numberPoles = numberZeros = order;
  [self initPolesZeros];
  return self;
}
- setFilterFo: (double)centerFreq fc:(double)cutoff;
{
  fo = centerFreq;
  fc = cutoff;
  return self;
}

- setFilterACoeff:(double *)aCoeffs bCoeff:(double *)bCoeffs;
{
  
  if(aPolyObject == nil)
  {
    aPolyObject = [[Polynomial alloc] initWithAPoly:aCoeffs order:filterOrder];
    bPolyObject = [[Polynomial alloc] initWithAPoly:bCoeffs order:filterOrder];
  }
  else
  {
    [aPolyObject setAPoly:aCoeffs order:filterOrder];
    [bPolyObject setAPoly:bCoeffs order:filterOrder];
  }
  
  return self;
}

/*##############################*
 * Get parameters:      *
 *##############################*/
- (COMPLEX *)getAnalogPoles { [self findPolesZeros]; return analogPoles; }
- (COMPLEX *)getAnalogZeros { [self findPolesZeros]; return analogZeros; }
- (int) getFilterOrder      {return filterOrder; }
- (int) getNumberOfPoles    {return numberPoles; }
- (int) getNumberOfZeros    {return numberZeros; }
- (double) fo               {return fo; }
- (double) fc               {return fc; }


/*##############################*
 * Find the filter's poles and  *
 * zeros:                       *
 *##############################*/
- findPolesZeros
{
  double  wc, wo;

/***************************
 * Conversions:            *
 ***************************/
  wc = TWOPI*fc;
  wo = TWOPI*fo;
  
  [self findAnalogPolesZerosFor:wc omegaO:wo bandwidth:(wc+wc)];

  return self;
                             
}

/*##############################*
 * Find the filter's poles and  *
 * zeros, also find passband    *
 * gain.                        *
 *##############################*/
- findAnalogPolesZerosFor:(double)wc omegaO:(double)wo bandwidth:(double)bw
{
  int i, n;
  COMPLEX carg, response;
  
  realAxisPoles = NO;
  filterOrder = MIN(MAX_ORDER, filterOrder);
  n = filterOrder;
  if(n<=0)
    analogPoles[0].x=analogPoles[0].y=analogZeros[0].x=analogZeros[0].y = 0.;
  else if( (wc==0.)  )              /* degenerate case     */
    for(i=0; i<n; i++)
      analogPoles[0].x=analogPoles[0].y=analogZeros[0].x=analogZeros[0].y = 0.;
  else
  {
    switch(analogType)              /* LPF wc = 1 */
    {
      case BUTTERWORTH:
      default:
        [self butterPrototype];     
        break;
      case CHEBYSHEV:
        [self chebyPrototype];
        break;
      case BESSEL:
        [self besselPrototype];
        break;
    }
    switch(filterPassType)          /* convert to wc filter */
    {
      case LOW_PASS:
      default:
        [self lowPassToLowPass:wc];
        CMPLX(carg, 0., 0.);
        response = [self poleZeroResponseAt:carg poles:analogPoles zeros:analogZeros];
        break;
      case HIGH_PASS:
        [self lowPassToHighPass:wc];
        CMPLX(carg, 0., 4.*wc);
        response = [self poleZeroResponseAt:carg poles:analogPoles zeros:analogZeros];
        break;
      case BAND_PASS:
        if(wo == 0.)
          [self lowPassToLowPass:wc];
        else
          [self lowPassToBandPass:wo bw:bw];
          CMPLX(carg, 0., wo);
          response = [self poleZeroResponseAt:carg poles:analogPoles zeros:analogZeros];
          break;
      case BAND_STOP:
        if(wo==0.)
          [self lowPassToHighPass:wc];
        else
          [self lowPassToBandStop:wo bw:bw];
        CMPLX(carg, 0., 0.);
       response = [self poleZeroResponseAt:carg poles:analogPoles zeros:analogZeros];
       break;
    }
    passBandGain = sqrt(response.x*response.x + response.y*response.y);
    passBandGain = MAX(1.e-40, passBandGain);
  }
  return self;
}


/*######################################################*
 * This method converts poles zeros to xfer fn form:    *
 *                                                      *
 *######################################################*/
- transferFromPoles:(COMPLEX *)poles andZeros:(COMPLEX *)zeros
{
  int    k, n, number_conj_poles;
  double real_pole, imag_pole, real_zero, imag_zero, wo;
  double c[1], d[3];

  wo = TWOPI*fo;
  n = numberPoles;

/***************************
 * Find the denominator    *
 * from the poles:         *
 ***************************/
  c[0] = 1.;
  [aPolyObject setAPoly:c order:0];
  
  number_conj_poles = n/2;
  if(realAxisPoles)
    number_conj_poles--;
  for(k=0; k<number_conj_poles; k++)
  {                         
    real_pole = poles[k].x;
    imag_pole = poles[k].y;
    d[2]      = 1.;
    d[1]      = -(real_pole + real_pole);
    d[0]      = real_pole*real_pole + imag_pole*imag_pole;
    [aPolyObject setBPoly:d order:2];
    [aPolyObject multABToC];
    [aPolyObject copyCToA];
  }   
/*********************************
 * If n is odd, or real axis     *
 * poles from band pass or band  *
 * stop transformations, mult    *
 * by these poles:               *
 *********************************/
  if( n%2 )
  {
    d[1] = 1.;
    d[0] = -poles[n/2].x;
    [aPolyObject setBPoly:d order:1];
    [aPolyObject multABToC];
    [aPolyObject copyCToA];
  }
  else if(realAxisPoles)
  {
    for(k=0; k<2; k++)
    {
      d[1] = 1.;
      d[0] = -poles[n/2-1].x;
      if(k == 0)
        d[0] = -poles[n-1].x;
      [aPolyObject setBPoly:d order:1];
      [aPolyObject multABToC];
      [aPolyObject copyCToA];
    }
  }
                                     
/***************************
 * Find the numerator      *
 * from the zeros:         *
 ***************************/
  n = numberZeros;
  c[0] = 1.;
  [bPolyObject setAPoly:c order:0];
  
  for(k=0; k<n/2; k++)
  {                         
    real_zero = zeros[k].x;
    imag_zero = zeros[k].y;
    d[2]      = 1.;
    if(filterPassType == BAND_PASS) /* bandpass has all real axis */
    {
      d[1]  = 0.;
      d[0]  = -real_zero*real_zero;
    }
    else
    {
      d[1]      = -(real_zero + real_zero);
      d[0]      = real_zero*real_zero + imag_zero*imag_zero;
    }
    [bPolyObject setBPoly:d order:2];
    [bPolyObject multABToC];
    [bPolyObject copyCToA];
  }   
/*********************************
 * If n is odd, mult by the real *
 * axis zero:                    *
 *********************************/
  if( n%2 )
  {
    d[1] = 1.;
    d[0] = -zeros[n/2].x;
    [bPolyObject setBPoly:d order:1];
    [bPolyObject multABToC];
    [bPolyObject copyCToA];
  }
  
  return self;
}

/*##############################*
 * This routine calculates the  *
 * poles for a Butterworth filter*
 * with cut-off frequency of wc *
 * = 1 rad/s and a pass band    *
 * ripple of 3 dB.  See Budak,  *
 *  Passive and Active Network  *
 * Analysis Synthesis, p. 506   *
 *##############################*/
- butterPrototype
{
  int   i, k, n;
  double arg;

/***************************
 * Find poles for wc = 1   *
 ***************************/
  n = filterOrder;
  for(k=0; k<n/2; k++)
  {                         
    arg       = PI*(1.+k+k)/2./n;
    analogPoles[k].x = -sin(arg);
    analogPoles[k].y = cos(arg);
  }
/*********************************
 * If n is odd, add pole at s = -1*
 *********************************/
  if( n%2 )
  {
    analogPoles[k].x = -1.;
    analogPoles[k++].y = 0.;
  }
/********************************
 * Copy upper half to lower *
 ********************************/
  i = 0;
  for( ; k<n; k++)
  {
    analogPoles[k].x = analogPoles[i].x;
    analogPoles[k].y = -analogPoles[i++].y;
  }
  return self;
}

/*##############################*
 * This routine calculates the  *
 * poles for a Chebyshev filter *
 * with cut-off frequency of wc *
 * = 1 rad/s and a pass band    *
 * ripple of 3 dB.  See Budak,  *
 *  Passive and Active Network  *
 * Analysis Synthesis, p. 515   *
 *##############################*/
- chebyPrototype
{
  int   i, k, n;
  double arg, sinh_value, cosh_value;

/***************************
 * Find poles for wc = 1   *
 ***************************/
  n = filterOrder;
  sinh_value = sinh(asinh(1./EPSILON)/n);
  cosh_value = cosh(asinh(1./EPSILON)/n);
  for(k=0; k<n/2; k++)
  {                         
    arg       = PI*(1.+k+k)/2./n;
    analogPoles[k].x = -sin(arg)*sinh_value;
    analogPoles[k].y = cos(arg)*cosh_value;
  }
/*********************************
 * If n is odd, add real axis pole*
 *********************************/
  if( n%2 )
  {
    analogPoles[k].x = -sinh_value;
    analogPoles[k++].y = 0.;
  }
/********************************
 * Copy upper half to lower *
 ********************************/
  i = 0;
  for( ; k<n; k++)
  {
    analogPoles[k].x = analogPoles[i].x;
    analogPoles[k].y = -analogPoles[i++].y;
  }
  return self;
}

/*##############################*
 * This routine calculates the  *
 * poles for a Bessel filter    *
 * with cut-off frequency of wc *
 * = 1 rad/s.       See Budak,  *
 *  Passive and Active Network  *
 * Analysis Synthesis, p. 525   *
 *##############################*/
- besselPrototype
{
  int   i, j, n;
  double *a, omega_c;
  COMPLEX *roots, ctemp, response;
  id    polynomial;
  id    root_finder;
  BOOL  error_flag;

/***************************
 * Find poles for wc = 1   *
 ***************************/
  n = filterOrder;
  a = (double *)malloc((n+1)*sizeof(double));
  
  if(n>0)
  {
    a[0] = a[1] = 1.;
    for(i=0; i<n; i++)
      a[i+1] = a[i]*2.*(n-i)/(n+n-i)/(i+1.);
    polynomial = [[Polynomial alloc] initWithAPoly:a order:n];
    if( (roots = [polynomial rootsOfA]) != NULL)
    {
        for(i=0; i<n; i++)
          analogPoles[i] = roots[i];
    }
    [polynomial free];
  }
  free(a);
//
// We need to normalize frequency
// so that response passes through 3 dB point
// at w == 1
//
  j = numberZeros;
  numberZeros = 0;
  CMPLX(ctemp, 0., 0.);
  response = [self poleZeroResponseAt:ctemp poles:analogPoles
                                  zeros:analogZeros];
  passBandGain = sqrt(response.x*response.x + response.y*response.y);
  passBandGain = MAX(1.e-40, passBandGain);
  root_finder = [[Numerical alloc] init];
  omega_c = [root_finder zeroIn:0.5 to:5. withTol:1.e-5 from:self error:&error_flag];
  omega_c = MAX(omega_c, 1.e-5);
  omega_c = 1./omega_c;
  numberZeros = j;
  [root_finder free];
  for(i=0; i<numberPoles; i++)
  {
    CTREAL(ctemp, analogPoles[i], omega_c);
    analogPoles[i] = ctemp;
  }
    
  
  return self;
}

/*##############################*
 * This routine moves the poles *
 * of a filter having a cut off *
 * frequency of  1 rad/s to wc  *
 * rad/s                        *
 *##############################*/
- lowPassToLowPass: (double)wc
{
  int n;

  numberZeros = 0;
  n = numberPoles = filterOrder;
/***************************
 * Multiply poles by wc    *
 ***************************/
  while(n--)
  {
    analogPoles[n].x *= wc;
    analogPoles[n].y *= wc;
  }
  
  return self;
}

/*##############################*
 * This routine moves the poles *
 * of a filter having a cut off *
 * frequency of 1 rad/s to wc   *
 * rad/s and transforms them to *
 * high pass poles      *
 *##############################*/
- lowPassToHighPass: (double)wc
{
  int n;
  double c_norm;
  COMPLEX ctemp, divisor;

  n = numberZeros = numberPoles = filterOrder;
/***************************
 * pole -> wc/pole         *
 * zeros: s**n             *
 ***************************/
  divisor.x = wc;
  divisor.y = 0.;
  while(n--)
  {
    if(CNORM(analogPoles[n])>0.)
    {
      CDIV(ctemp, divisor, analogPoles[n]);
      analogPoles[n].x = ctemp.x;
      analogPoles[n].y = ctemp.y;
    }
    analogZeros[n].x = analogZeros[n].y = 0.;
  }
  
  return self;  
}

/*##############################*
 * This routine moves the poles *
 * of a filter having a cut off *
 * frequency of 1 rad/s to wc   *
 * rad/s and transforms them to *
 * band pass poles      *
 *##############################*/
- lowPassToBandPass:(double)wo bw:(double)bw
{
  int     i, j, k, n;
  double  wo_lp, q_lp, delta, q_bp, wo_bp[2];
  double  real_part, imag_part, temp1, temp2;
  COMPLEX *new_poles;
 
/***************************
 * Find # of conjugate low *
 * pass poles:             *
 ***************************/
  numberZeros = n = filterOrder;
  k = filterOrder/2;
  numberPoles = n+n;
  real_part = imag_part = temp1 = temp2 = 0.;
  
/***************************
 * Note: Each low pass pole*
 * produces two bp poles   *
 ***************************/
  new_poles = (COMPLEX *)malloc(n*sizeof(COMPLEX));

/***************************
 * Find bandpass poles     *
 * as in p. 580, Budak     *
 ***************************/
  i = 0;
  while(k--)
  {
    wo_lp = sqrt(analogPoles[k].x*analogPoles[k].x +
                 analogPoles[k].y*analogPoles[k].y);
    q_lp  = -wo_lp/2./analogPoles[k].x;
    delta = bw*wo_lp/wo;
    q_bp  = sqrt((1+4./delta/delta)*(1+4./delta/delta) 
                  -4./delta/delta/q_lp/q_lp);
    q_bp  = q_lp*sqrt(1+4./delta/delta + q_bp)/SQRT2;
    wo_bp[0] = 0.5*wo*(delta*q_bp/q_lp - 
                   sqrt(delta*delta*q_bp*q_bp/q_lp/q_lp -4.) );
    wo_bp[1] = 0.5*wo*(delta*q_bp/q_lp + 
                   sqrt(delta*delta*q_bp*q_bp/q_lp/q_lp -4.) );
    for(j=0; j<2; j++)
    {
      real_part = -wo_bp[j]/2./q_bp;
      imag_part = sqrt(wo_bp[j]*wo_bp[j] - real_part*real_part);
      new_poles[i].x = real_part;
      new_poles[i].y = imag_part;
      i++;
    }
  }
/**********************
 * Real axis pole:    *
 **********************/
  if(n%2)
  {
    real_part = bw*analogPoles[n/2].x;
    temp1     = 4.*wo*wo;
    temp2     = real_part*real_part;
    if(temp1>=temp2)
    {
        imag_part = sqrt(temp1 - temp2)/2.;
        new_poles[i].x = real_part/2.;
    }
    else
    {
        realAxisPoles = YES;
        imag_part = 0.;
    new_poles[i].x = (real_part+sqrt(temp2-temp1))/2.;
    }
    new_poles[i].y = imag_part;
  }
/*********************
 * Copy over poles   *
 * zeros: s**(n/2)   *
 *********************/
  for(i=0; i<n; i++)
  {
    analogPoles[i].x = new_poles[i].x;
    analogPoles[i].y = new_poles[i].y;
    analogZeros[i].x = analogZeros[i].y = 0.;
  }
  for( ; i<numberPoles; i++)        /* conjugate poles */
  {
    analogPoles[i].x = analogPoles[i-n].x;
    analogPoles[i].y = -analogPoles[i-n].y;
  }
  if(realAxisPoles)
    analogPoles[numberPoles-1].x = (real_part - sqrt(temp2-temp1))/2.;
  free(new_poles);
  return self;
}

/*##############################*
 * This routine moves the poles *
 * of a filter having a cut off *
 * frequency of 1 rad/s to wc   *
 * rad/s and transforms them to *
 * band stop poles      *
 *##############################*/
- lowPassToBandStop: (double)wo bw:(double)bw
{
  int     i, j, k, n;
  double  wo_lp, q_lp, delta, q_bs, wo_bs[2];
  double  real_part, imag_part, temp1, temp2;
  COMPLEX *new_poles;
 
/***************************
 * Find # of conjugate low *
 * pass poles:             *
 ***************************/
  n = filterOrder;
  k = filterOrder/2;
  numberZeros  = numberPoles = n+n;
  real_part = imag_part = temp1 = temp2 = 0.;
  
/***************************
 * Note: Each low pass pole*
 * produces two bs poles   *
 ***************************/
  new_poles = (COMPLEX *)malloc(n*sizeof(COMPLEX));

/***************************
 * Find bandstop poles     *
 * as in p. 630, Budak     *
 ***************************/
  i = 0;
  while(k--)
  {
    wo_lp = sqrt(analogPoles[k].x*analogPoles[k].x +
                 analogPoles[k].y*analogPoles[k].y);
    q_lp  = -wo_lp/2./analogPoles[k].x;
    delta = bw/wo_lp/wo;
    q_bs  = sqrt((1+4./delta/delta)*(1+4./delta/delta) 
                  -4./delta/delta/q_lp/q_lp);
    q_bs  = q_lp*sqrt(1+4./delta/delta + q_bs)/SQRT2;
    wo_bs[0] = 0.5*wo*(delta*q_bs/q_lp - 
                   sqrt(delta*delta*q_bs*q_bs/q_lp/q_lp -4.) );
    wo_bs[1] = 0.5*wo*(delta*q_bs/q_lp + 
                   sqrt(delta*delta*q_bs*q_bs/q_lp/q_lp -4.) );
    for(j=0; j<2; j++)
    {
      real_part = -wo_bs[j]/2./q_bs;
      imag_part = sqrt(wo_bs[j]*wo_bs[j] - real_part*real_part);
      new_poles[i].x = real_part;
      new_poles[i].y = imag_part;
      i++;
    }
  }
/**********************
 * Real axis pole:    *
 **********************/
  if(n%2)
  {
    if(analogPoles[n/2].x != 0)
    {
      real_part = bw/analogPoles[n/2].x;
      temp1     = 4.*wo*wo;
      temp2     = real_part*real_part;
      if(temp1>=temp2)
      {
        imag_part = sqrt(temp1 - temp2)/2.;
        new_poles[i].x = real_part/2.;
      }
      else
      {
        realAxisPoles = YES;
        imag_part = 0.;
    new_poles[i].x = (real_part+sqrt(temp2-temp1))/2.;
      }
    }
    new_poles[i].y = imag_part;
  }
/**********************
 * Copy over poles,   *
 * zeros:             *
 * (s-jwo)(s+jwo)     *
 **********************/
  for(i=0; i<n; i++)
  {
    analogPoles[i].x = new_poles[i].x;
    analogPoles[i].y = new_poles[i].y;
    analogZeros[i].x = 0.;
    analogZeros[i].y = -wo;
  }
  for( ; i<numberPoles; i++)        /* conjugate poles */
  {
    analogPoles[i].x = analogPoles[i-n].x;
    analogPoles[i].y = -analogPoles[i-n].y;
    analogZeros[i].x = 0.;
    analogZeros[i].y = wo;
  }
  if(realAxisPoles)
    analogPoles[numberPoles-1].x = (real_part - sqrt(temp2-temp1))/2.;
  free(new_poles);
  return self;
}

/*##############################*
 * This routine calculates the  *
 * response of the analog   *
 * filter.          *
 *##############################*/
- (float *)analogResponseAtFrequencies:(float *)omega
                          numberFreq:(int)numberPts
{
  int     n;
  static  int number_pts_store = 0;
  static  float   *response = NULL;
  COMPLEX complex_response, carg;

  [self findPolesZeros];
  if(filterStructureType == TRANSFER_FUNCTION)
    [self transferFromPoles:analogPoles andZeros:analogZeros];
/***************************
 * Get space for output    *
 * arrays:                 *
 ***************************/  
  if(numberPts > number_pts_store)
  {
    if(response != NULL)
      free(response);
    response         = (float *) malloc(numberPts*sizeof(float));
    number_pts_store = numberPts;
  }

  if(numberPts>1)
    deltaOmega  = (omega[numberPts-1] - omega[0])/numberPts;
  deltaOmega    = MAX(1e-30, deltaOmega);
  thetaOld      = 0.01;
  subAngle      = 0.;
/****************************
 * Loop over omega:         *
 ****************************/
  for(n=0; n<numberPts; n++)
  {                  
                                          
    CMPLX(carg,0.,omega[n])
    if(filterStructureType == POLE_ZERO)
       complex_response = [self poleZeroResponseAt:carg poles:analogPoles
                                        zeros:analogZeros];
    else 
       complex_response = [self transferResponseAt:carg];
    
    response[n] = [self outputResponseFor:&complex_response at:omega[n]];
  }
  return response;

}

/*##############################*
 * Find complex pole zero   *
 * response:            *
 *##############################*/
- (COMPLEX) poleZeroResponseAt: (COMPLEX)point poles:(COMPLEX *)poles
                         zeros:(COMPLEX *)zeros
{
  int     i;
  double  c_norm;
  COMPLEX num,den,ctemp; 
  COMPLEX cprod,response;

/****************************
 * Find numerator and de-   *
 * nominator for pole zero: *
 ****************************/
  CMPLX(num,1.,0.);
  for(i=0; i<numberZeros; i++)
  {
    CSUB(ctemp,point,zeros[i]);
    CMULT(cprod,num,ctemp);
    num = cprod;
  }
  CMPLX(den,1.,0.);
  for(i=0; i<numberPoles; i++)
  {
    CSUB(ctemp,point,poles[i]);
    CMULT(cprod,den,ctemp);
    den = cprod;  
  }
  CDIV(response, num, den);
  
  return response;
}
    

/*##############################*
 * Find complex xfer function   *
 * response:            *
 *##############################*/
- (COMPLEX) transferResponseAt:(COMPLEX)point
{
  double  c_norm;
  COMPLEX num, den, response;

/****************************
 * Find numerator value     *
 * use Horner's method:     *
 ****************************/
  num = [bPolyObject evaluateAAtComplexPoint:point];

/****************************
 * Find denominator value   *
 * use Horner's method:     *
 ****************************/ 
  den = [aPolyObject evaluateAAtComplexPoint:point];

/***************************
 * Divide numerator by     *
 * denominator:            *
 ***************************/
  CDIV(response, num, den);
  
  return response;
}

/*##############################*
 * Convert complex response to  *
 * form desired by user:    *
 *##############################*/
- (float)outputResponseFor:(COMPLEX *)complexResponse at:(float)omega
{
  float response;
  double temp, theta_new;

  switch(responseType)
  {
    case MAGNITUDE:
    default:
      response = sqrt(complexResponse->x*complexResponse->x +
                       complexResponse->y*complexResponse->y);
      response /= passBandGain;
      break;
    case DB_MAGNITUDE:
      temp = complexResponse->x*complexResponse->x +
                       complexResponse->y*complexResponse->y;
      temp /= passBandGain*passBandGain;
      if(temp > 0.)
        response = 10.*log10(temp);
      else
        response = 0.;
      break;
    case PHASE:
      response = atan2(complexResponse->y, complexResponse->x);
      if((thetaOld<=0.) && (response>0.))       /* Pi crossing*/
        subAngle += TWOPI;
      thetaOld = response;
      response = response - subAngle;
      break;
    case PHASE_DELAY:
      response = atan2(complexResponse->y, complexResponse->x);
      if((thetaOld<=0.) && (response>0.))       /* Pi crossing*/
        subAngle += TWOPI;
      thetaOld = response;
      response = response - subAngle;
      if(omega > 0.)
        response /= -omega;
      break;
    case GROUP_DELAY:
      theta_new = atan2(complexResponse->y, complexResponse->x);
      if(SGN01(theta_new*thetaOld))         /* Pi crossing */
        response = (thetaOld - theta_new)/deltaOmega;
      else
        response = (thetaOld - theta_new + TWOPI)/deltaOmega;
      thetaOld = theta_new;
      break;
    case ONE_MINUS_MAG:
      response = sqrt((passBandGain-complexResponse->x)*(passBandGain-complexResponse->x) +
                       complexResponse->y*complexResponse->y);
      response /= passBandGain;
      break;
  }
  return response;
}

/*#####################################*
 * Used in besselPrototype to find  *
 * 3 dB point               *
 *######################################*/
- (double)functionToFindRoot:(double)x
{
  COMPLEX ctemp, response;
  double  output;
  
  CMPLX(ctemp, 0., x);
  response = 
         [self poleZeroResponseAt:ctemp poles:analogPoles zeros:analogZeros];
  output = response.x*response.x + response.y*response.y;
  output /= passBandGain*passBandGain;
  output -= 0.5;
 
  return output;
}

@end