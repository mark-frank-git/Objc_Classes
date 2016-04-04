/************************************************************************
 * This subclass of analogFilter implements a digital filtering object  *
 *                                                                      *
 * File:DigitalFilter.m                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/12/92  - Started                                              *
 *  2. 01/28/93  - Divide output by passBandGain                        *
 *  3. 12/22/93  - Fix in findTransferFunction to account for FIR.      *
 ************************************************************************/
#import "DigitalFilter.h"
#import "Polynomial.h"
#import <stdlib.h>
#import <math.h>

#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )


@implementation DigitalFilter

- init
{
  [super init];
  
// init instance variables:
  filterType        = DIGITAL_FILTER;
  windowType        = HAMMING;
  doubleShift       = NULL;
  complexShift      = NULL;
  digitalPoles = digitalZeros   = NULL;
  interpolateFactor = decimateFactor = 1;
  
  return self;
}


- initWithPassType:(int)type fo:(double)centerFreq fc:(double)cutoffFreq
                fs:(double)samplingFreq order:(int)order
{
  [self init];

  filterPassType    = type;
  filterOrder       = order;
  numberPoles       = numberZeros = filterOrder;
  fo                = centerFreq;
  fc                = cutoffFreq;
  fs                = samplingFreq;
  
  [self initPolesZeros];
  
  return self;
}

- initPolesZeros
{
  int i;
  
  [super initPolesZeros];
  if(digitalPoles != NULL)
  {
    free(digitalPoles);
    free(digitalZeros);
    free(doubleShift);
    free(complexShift);
  }
  digitalPoles  = (COMPLEX *) malloc(2*numberPoles*sizeof(COMPLEX));
  digitalZeros  = (COMPLEX *) malloc(2*numberPoles*sizeof(COMPLEX));
  doubleShift   = (double *) malloc(2*filterOrder*sizeof(double));
  complexShift  = (COMPLEX *) malloc(2*filterOrder*sizeof(COMPLEX));
  for(i=0; i<2*filterOrder; i++)
    doubleShift[i] = complexShift[i].x = complexShift[i].y = 0.;
  return self;
}

/*##############################*
 * Filter the input array and   *
 * overwrite it with output.    *
 * It uses the direct form II   *
 * network, p. 151, Oppenheim & *
 * Schafer.         *
 *##############################*/
- filterFloatArray:(float *)x numberPts:(int)numberPts
{
   double yout, xin, *a, *b;
   int    i, j, i_minus_1, n;

   if(aPolyObject == nil)
     return self;
   
   a = [aPolyObject getAPoly];
   b = [bPolyObject getAPoly];
   if( (n=numberPoles) == 0)
   {
     for(j=0; j<numberPts; j++)
       x[j] *= b[0]/passBandGain;
     return self;
   }
      
   for(j=0; j<numberPts; j++)
   {
     xin  = x[j] - a[n-1]*doubleShift[0];
     yout = b[n-1]*doubleShift[0];
     for(i=n; i>1 ;i--)
     {
       i_minus_1 = i-1;
       yout  += b[n-i]*doubleShift[i_minus_1];
       xin   -= a[n-i]*doubleShift[i_minus_1];
       doubleShift[i_minus_1] = doubleShift[i-2];
     }
     doubleShift[0] = xin;
     yout    += b[n]*xin;
     x[j]     = (float )yout/passBandGain;
   }
   return self;
}

/*##############################*
 * Filter the input data point  *
 * It uses the direct form II   *
 * network, p. 151, Oppenheim & *
 * Schafer.         *
 *##############################*/
- (double) filterFloatData:(float)xin
{
   double yout, *a, *b;
   int    i, i_minus_1, n;

   if(aPolyObject == nil)
     return 0.;
   
   a = [aPolyObject getAPoly];
   b = [bPolyObject getAPoly];
   if( (n=numberPoles) == 0)
     return (xin*b[0]/passBandGain);
      
   yout = b[n-1]*doubleShift[0];
   for(i=n; i>1 ;i--)
   {
       i_minus_1 = i-1;
       yout  += b[n-i]*doubleShift[i_minus_1];
       xin   -= a[n-i]*doubleShift[i_minus_1];
       doubleShift[i_minus_1] = doubleShift[i-2];
   }
   xin -= a[n-1]*doubleShift[0];
   doubleShift[0] = xin;
   yout    += b[n]*xin;
   yout    /= passBandGain;
   
   return yout;
}

/*##############################*
 * Filter the complex data pt.  *
 * It uses the direct form II   *
 * network, p. 151, Oppenheim & *
 * Schafer.                     *
 *##############################*/
- (FCOMPLEX) filterComplexData:(FCOMPLEX )xin
{
   FCOMPLEX yout, temp1, temp2;
   int    i, i_minus_1, n;
   double *a, *b, dtemp;

   if(aPolyObject == nil)
     return xin;
   
   a = [aPolyObject getAPoly];
   b = [bPolyObject getAPoly];

   n = numberPoles;
   if(n==0)
   {
      CTREAL(yout, xin, b[0]);
      return yout;
   }
   CTREAL(yout,complexShift[0],b[n-1]);
   for(i=n; i>1 ;i--)
   {
     i_minus_1 = i-1;
     CTREAL(temp1,complexShift[i_minus_1],b[n-i]);
     CADD(yout, yout, temp1);
     CTREAL(temp2,complexShift[i_minus_1],a[n-i]);
     CSUB(xin, xin, temp2);
     CLET(complexShift[i_minus_1], complexShift[i-2]);
   }
   CTREAL(temp2,complexShift[0],a[n-1]);
   CSUB(xin, xin, temp2);
   CTREAL(temp2,xin,b[n]);
   CADD(yout, yout, temp2);
   CLET(complexShift[0], xin);
   dtemp = 1./passBandGain;
   CTREAL(temp2,yout,dtemp);        /* normalize gain   */
   CLET(yout, temp2);
   return yout;
}

/*##############################*
 * Filter the complex data pt.  *
 * It uses the direct form II   *
 * network, p. 151, Oppenheim & *
 * Schafer.                     *
 *##############################*/
- filterComplexArray:(FCOMPLEX *)x numberPts:(int)numberPts
{
   FCOMPLEX yout, temp1, temp2, xin;
   int    i, j, i_minus_1, n;
   double *a, *b, dtemp;

   if(aPolyObject == nil)
     return self;
   
   a = [aPolyObject getAPoly];
   b = [bPolyObject getAPoly];

   n = numberPoles;
   if(n==0)
   {
     for(j=0; j<numberPts; j++)
     {
       CTREAL(yout, x[j], b[0]/passBandGain);
       x[j] = yout;
     }
   }
   for(j=0; j<numberPts; j++)
   {
     CTREAL(temp2,complexShift[0],a[n-1]);
     CSUB(xin, x[j], temp2);
     CTREAL(yout,complexShift[0],b[n-1]);
     for(i=n; i>1 ;i--)
     {
       i_minus_1 = i-1;
       CTREAL(temp1,complexShift[i_minus_1],b[n-i]);
       CADD(yout, yout, temp1);
       CTREAL(temp2,complexShift[i_minus_1],a[n-i]);
       CSUB(xin, xin, temp2);
       CLET(complexShift[i_minus_1], complexShift[i-2]);
     }
     CLET(complexShift[0], xin);
     CTREAL(temp2,complexShift[0],b[n]);
     CADD(yout, yout, temp2);
     dtemp = 1./passBandGain;
     CTREAL(temp2,yout,dtemp);      /* normalize gain   */
     CLET(x[j], temp2);
  }
  return self;
}

#define MAX_WARP  11            /* interp, decimate < MAX_WARP*/
/*##############################*
 * Find the interpolation and   *
 * decimation factors for a     *
 * given sampling ratio     *
 *##############################*/
- findWarpFactorsFor:(double)ratio
{
  int    min_interpolate, min_decimate;
  double min_delta, delta;

/*************************
 * Now find warp factor: *
 *************************/                                              
  min_delta = 100.;
  min_interpolate = min_decimate = 1;
  for(decimateFactor=1; decimateFactor<MAX_WARP; decimateFactor++)
  {                                  
    delta = [self checkWarpFor:ratio];
    if(delta < min_delta)
    {                        
      min_delta = delta;
      min_interpolate = interpolateFactor;
      min_decimate    = decimateFactor;
    }
  }
  interpolateFactor = min_interpolate;
  decimateFactor    = min_decimate;

  return self;
}

/*##############################*
 * Check the interpolation and  *
 * decimation factors versus the*
 * input ratio, return delta.   *
 *##############################*/
- (double)checkWarpFor:(double)ratio
{
   int i;
   double float_decimate, float_interpolate;                                 
   double new_delta, min_delta;
                   
   float_decimate   = (double) decimateFactor;    
   float_interpolate    = 1.; 
   min_delta = 100.;
   for(i=1; i<=MAX_WARP; i++)
   {     
     new_delta = fabs(float_interpolate/float_decimate - ratio);
     if(new_delta < min_delta)
     {
       interpolateFactor = i;
       min_delta = new_delta;
     }
     float_interpolate++;
   }
   return min_delta;
}

/*##############################*
 * Interpolate and input    *
 * array by inserting zeros *
 * return the interpolated array*
 *##############################*/
- (float *)interpolateInput:(float *)input numberPts:(int)numberPts
{
  int   i, j, k, output_pts;
  static int old_number_pts = 0;
  static float *output = NULL;
  float  temp;
  
  output_pts = numberPts*interpolateFactor;
  if(old_number_pts < output_pts)
  {
    if(output != NULL)
      free(output);
    output = (float *) malloc(output_pts*sizeof(float));
    old_number_pts = output_pts;
  }

  k = 0;
  temp = (float)interpolateFactor;
  for(i=0; i<numberPts; i++)
  {
    output[k++] = temp*input[i];
    for(j=1; j<interpolateFactor; j++)
      output[k++] = 0.;
  }
  return output;
}


/*##############################*
 * Decimate an input array  *
 * by selecting the ith input   *
 * return the decimated array.  *
 *##############################*/
- (float *)decimateInput:(float *)input numberPts:(int)numberPts
{
  int   i, k, output_pts;
  static int old_number_pts = 0;
  static float *output = NULL;
  
  if(decimateFactor)
  {
    output_pts = numberPts/decimateFactor;
    if(old_number_pts < output_pts)
    {
      if(output != NULL)
        free(output);
      output = (float *) malloc(output_pts*sizeof(float));
      old_number_pts = output_pts;
    }
    k = 0;
    for(i=0; i<output_pts; i++)
    {
      output[i] = input[k];
      k += decimateFactor;
    }
  }
  return output;
}

/*##############################*
 * Set Parameters:      *
 *##############################*/
- setDigitalFilterType:(int)type {filterType = type; return self;}
- setWindowingType:(int)type {windowType = type; return self;}
- setFilterFo: (double)centerFreq fc:(double)cutoff fs:(double)samplingFreq;
{
  fo = centerFreq;
  fc = cutoff;
  fs = samplingFreq;
  return self;
}
- zeroOutTaps
{
  int i;
  for(i=0; i<2*filterOrder; i++)
    doubleShift[i] = complexShift[i].x = complexShift[i].y = 0.;
  return self;
}

/*##############################*
 * Get parameters:      *
 *##############################*/
- (int) getFilterType       {return filterType; }
- (int) getInterpolateFactor    {return interpolateFactor; }
- (int) getDecimateFactor   {return decimateFactor; }

- (COMPLEX *)getPoles
{ 
  [self findPolesZeros];
  if(filterType == DIGITAL_FILTER)
    return digitalPoles;
  else
    return analogPoles;
}
- (COMPLEX *)getZeros
{ 
  [self findPolesZeros];
  if(filterType == DIGITAL_FILTER)
    return digitalZeros;
  else
    return analogZeros;
}

/*##############################*
 * Find the filter's poles and  *
 * zeros:                       *
 *##############################*/
- findPolesZeros
{
  int     i;
  double  wc, wo, delta_t;
  double  wc_warp, wo_warp, w1_warp, w2_warp;
  double  real_pole, imag_pole;
  COMPLEX carg, response, exponent;

  if(filterType == ANALOG_FILTER)
  {
    [super findPolesZeros];
    return self;
  }
/***************************
 * Conversions:            *
 ***************************/
  wc = TWOPI*fc;
  wo = TWOPI*fo;
  fs = MAX(1.e-30, fs);
  delta_t = 1./fs;
  wo_warp = wo/fs;              /* pre- warp center freq*/
  wo_warp = 2.*tan(wo_warp/2.)/delta_t;
  wc_warp = wc/fs;              /* pre- warp cut off    */
  wc_warp = 2.*tan(wc_warp/2.)/delta_t;
  w1_warp = (wo-wc)/fs;             /* pre-warp bandwidth   */
  if(w1_warp>0.)
    w1_warp = 2.*tan(w1_warp/2.)/delta_t;
  w2_warp = (wo+wc)/fs;             /* pre-warp bandwidth   */
  if(w2_warp>0.)
    w2_warp = 2.*tan(w2_warp/2.)/delta_t;

/***************************
 * Find Analog poles and   *
 * zeros, then convert to  *
 * z domain using bilinear *
 * tranformation.          *
 ***************************/
  [self findAnalogPolesZerosFor:wc_warp omegaO:wo_warp
                                     bandwidth:(w2_warp-w1_warp)];
  for(i=0; i<numberPoles; i++)
    [self sToz:analogPoles[i].x: analogPoles[i].y: &digitalPoles[i].x:
           &digitalPoles[i].y];


/*********************************
 * Now find zeros:               *
 *    LPF = (z+1)**m             *
 *    HPF = (z-1)**m             *
 *    BPF = [(z-1)(z+1)]**m      *
 *    BSF = [(s-jwo)(s+jwo)]**m  *
 *********************************/
  numberZeros = numberPoles;
  switch(filterPassType)
  {
    case LOW_PASS:
    default:
      for(i=0; i<numberPoles; i++)
      {
        digitalZeros[i].x = -1.;
        digitalZeros[i].y =  0.;
      }
      CMPLX(exponent,0.,0.)
      CEXP(carg,exponent)
      response = [self poleZeroResponseAt:carg poles:digitalPoles
                                  zeros:digitalZeros];
      break;
    case HIGH_PASS:
      for(i=0; i<numberPoles; i++)
      {
        digitalZeros[i].x = 1.;
        digitalZeros[i].y =  0.;
      }
      CMPLX(exponent,0.,PI)
      CEXP(carg,exponent)
      response = [self poleZeroResponseAt:carg poles:digitalPoles
                                  zeros:digitalZeros];
      break;
    case BAND_PASS:
      for(i=0; i<numberPoles; i++)
      {
        if(i%2)
          digitalZeros[i].x =  1.;
        else
          digitalZeros[i].x = -1.;
        digitalZeros[i].y =  0.;
      }
      [self sToz:0.: wo_warp: &carg.x: &carg.y];
      response = [self poleZeroResponseAt:carg poles:digitalPoles
                                  zeros:digitalZeros];
      break;
    case BAND_STOP:
      [self sToz:0.: wo_warp: &real_pole: &imag_pole];
      for(i=0; i<numberPoles; i++)
      {
        if(i<numberPoles/2)
          digitalZeros[i].y =  imag_pole;
        else
          digitalZeros[i].y = -imag_pole;
        digitalZeros[i].x =  real_pole;
      }
      CMPLX(exponent,0.,0.)
      CEXP(carg,exponent)
      response = [self poleZeroResponseAt:carg poles:digitalPoles
                                  zeros:digitalZeros];
      break;
  }
  passBandGain = sqrt(response.x*response.x + response.y*response.y);
  passBandGain = MAX(1.e-40, passBandGain);
  
  return self;
                             
}

/***************************************************************************
 *                                                                         *
 *  void s_to_z(double s_real, double s_imag, double t,                    *
 *                                 double *z_real, double *z_imag)         *
 *                                                                         *
 *    double input variables                                               *
 *    -------------------                                                  *
 *    s_real = real part of s plane pole                                   *
 *    s_imag = imag part of s plane pole                                   *
 *    t      = sampling time in seconds                                    *
 *                                                                         *
 *    double  output variables                                             *
 *    ----------------------                                               *
 *    z_real = real part of z plane pole                                   *
 *    z_imag = imag part of z plane pole                                   *
 *                                                                         *
 *                                                                         *
 *  Convert s plane pole to z plane pole using Bilinear transformation     *
 *  s -> 2(z-1)/(z+1)/T
 *                                                                         *
 ***************************************************************************/
- sToz: (double)sReal: (double)sImag: (double *)zReal: (double *)zImag
{
  double  t, c_norm;
  COMPLEX p_t_2, cnum, cden;
  
  if(fs>0.)
    t = 1/fs;
  else
    t = 1.;
  p_t_2.x = sReal*t/2.;
  p_t_2.y = sImag*t/2.;
  
  cnum.x = 1. + p_t_2.x;
  cnum.y = p_t_2.y;
  
  cden.x = 1. - p_t_2.x;
  cden.y = -p_t_2.y;
   
  CDIV(p_t_2, cnum, cden);
  
  *zReal = p_t_2.x;
  *zImag = p_t_2.y;
  
  return self;
}

/*##############################*
 * This routine finds the   *
 * coefficients of the low pass *
 * filter found by multiplying  *
 * the window function by a sinc*
 * pulse whose spectrum is an   *
 * ideal filter.                *
 *##############################*/
- transferForFIRWindow
{
  int    i, j, alpha, n; 
  double wo, wc, w1, w2;
  double *a, *b;
  
  if(aPolyObject == nil)
     return self;
   
  n = filterOrder+1;
  a = (double *)malloc(n*sizeof(double));
  b = (double *)malloc(n*sizeof(double));
  [self windowFunction: n];
  alpha  = (n-1)/2;
  wo     = TWOPI*fo/fs;
  wc     = TWOPI*fc/fs;
  w1     = wo - wc;
  w2     = wo + wc;
  j      = filterOrder;
  passBandGain = 1.;
  if(filterStructureType == FIR_RAISED_COS)
  {
    [self raisedCosineCoefficients:b number:n];
    for(i=0; i<n; i++)
    {
      a[i] = 0.;
      b[i] *= window[i];
    }
  }
  else
  {
   for(i=0; i<n; i++)
   {
    a[i] = 0.;
    switch(filterPassType)
    {
      case LOW_PASS:
      default:
        if(i != alpha)
          b[j] = sin(wc*(i-alpha))*window[i]/PI/(i-alpha);
        else
          b[j] = wc*window[i]/PI;
        break;
      case HIGH_PASS:
        if(i != alpha)
          b[j] = (sin(PI*(i-alpha)) - sin(wc*(i-alpha)))
                      *window[i]/PI/(i-alpha);
        else
          b[j] = (PI-wc)*window[i]/PI;
        break;
      case BAND_PASS:
        if(i != alpha)
          b[j] = (sin(w2*(i-alpha)) - sin(w1*(i-alpha)))
                      *window[i]/PI/(i-alpha);
        else
          b[j] = (w2-w1)*window[i]/PI;
        break;
      case BAND_STOP:
        if(i != alpha)
          b[j] = (sin(w1*(i-alpha)) + sin(PI*(i-alpha)) - 
                   sin(w2*(i-alpha))) *window[i]/PI/(i-alpha);
        else
          b[j] = (TWOPI - 2*(w2-w1))*window[i]/TWOPI;
        break;
    }
    j--;
   }
  }
  a[filterOrder] = 1.;

  [aPolyObject setAPoly:a order:filterOrder];
  [bPolyObject setAPoly:b order:filterOrder];
  free(a);
  free(b);

  return self;
}

/*##############################*
 * Find the transfer function   *
 * from the poles and zeros:    *
 *##############################*/
- findTransferFunction
{
  if(filterType == DIGITAL_FILTER)
  {
    if( (filterStructureType==FIR_WINDOW) || (filterStructureType==FIR_RAISED_COS) )
      [self transferForFIRWindow];
    else
      [self transferFromPoles:digitalPoles andZeros:digitalZeros];
  }
  else
    [self transferFromPoles:analogPoles  andZeros:analogZeros];

  return self;
}


/*##############################*
 * This routine calculates the  *
 * response of the filter   *
 *##############################*/
- (float *)filterResponseAtFrequencies:(float *)omega
                          numberFreq:(int)numberPts
{
  if(filterType == DIGITAL_FILTER)
    return [self digitalResponseAtFrequencies:omega numberFreq:numberPts];
  else
    return [self analogResponseAtFrequencies:omega numberFreq:numberPts];
}

/*##############################*
 * This routine calculates the  *
 * response of the digital  *
 * filter.          *
 *##############################*/
- (float *)digitalResponseAtFrequencies:(float *)omega
                          numberFreq:(int)numberPts
{
  int     n;
  static  int number_pts_store = 0;
  static  float   *response = NULL;
  COMPLEX *complex_response, carg, exponent;

  switch(filterStructureType)
  {
    case TRANSFER_FUNCTION:
    case POLE_ZERO:
    default:
      [self findPolesZeros];
      if(filterStructureType == TRANSFER_FUNCTION)
       [self transferFromPoles:digitalPoles andZeros:digitalZeros];
      break;
    case FIR_WINDOW:
    case FIR_RAISED_COS:
      [self transferForFIRWindow];
      break;
  }
  if(numberPts>1)
    deltaOmega  = (omega[numberPts-1] - omega[0])/numberPts;
  deltaOmega    = MAX(1e-30, deltaOmega);
    

/***************************
 * Get space for output    *
 * arrays:                 *
 ***************************/  
  complex_response = (COMPLEX *) malloc(numberPts*sizeof(COMPLEX));
  if(numberPts > number_pts_store)
  {
    if(response != NULL)
      free(response);
    response         = (float *) malloc(numberPts*sizeof(float));
    number_pts_store = numberPts;
  }

/****************************
 * Loop over omega:         *
 ****************************/
  for(n=0; n<numberPts; n++)
  {                  
                                          
    CMPLX(exponent,0.,omega[n])
    CEXP(carg,exponent)
    if(filterStructureType == POLE_ZERO)
      complex_response[n] = [self poleZeroResponseAt:carg poles:digitalPoles
                                        zeros:digitalZeros];
    else 
      complex_response[n] = [self transferResponseAt:carg];
    
    response[n] = [self outputResponseFor:&complex_response[n] at:fs*omega[n]];
  }
  return response;

}


/*##############################*
 * This routine calculates the  *
 * windowing function:      *
 *##############################*/
- (float *)windowFunction: (int)length
{
  int   i;
  double arg;
  
  if(window != NULL)
    free(window);
  window = (float *)malloc(length*sizeof(float));

  switch(windowType)
  {
    case RECTANGULAR:
      for(i=0; i<length; i++)
        window[i] = 1.;
      break;
    case HAMMING:
      arg = TWOPI/(length-1.);
      for(i=0; i<length; i++)
        window[i] = 0.54 - 0.46*cos(arg*i);
      break;
    case HANNING:
      arg = TWOPI/(length-1.);
      for(i=0; i<length; i++)
        window[i] = 0.50 - 0.50*cos(arg*i);
      break;
    case TRIANGULAR:
      arg = 2.0/(length-1.);
      for(i=0; i<=(length-1)/2; i++)
        window[i] = i*arg;
      for(   ; i<length; i++)
        window[i] = 2. - i*arg;
      break;
    case BLACKMAN:
      arg = TWOPI/(length-1.);
      for(i=0; i<length; i++)
        window[i] = 0.42 - 0.50*cos(arg*i) + 0.08*cos((arg+arg)*i);
      break;
    case BLACKMAN_HARRIS:
      arg = TWOPI/(length-1.);
      for(i=0; i<length; i++)
        window[i] = 0.35875 - 0.48829*cos(arg*i) +
                      0.14128*cos((arg+arg)*i) - 0.01168*cos(3.*arg*i);
      break;
    default:      
      for(i=0; i<length; i++)
        window[i] = 1.;
      break;
  }
  
  return window;
}

#define ROLLOFF 0.4
/*##################################*
 * Find the impulse response of     *
 * the 40% sq root raised cosine    *
 * filter.  See Couch p. 166.  We   *
 * need to do an IFFT, since we     *
 * need the square root of the      *
 * filter (1/2 in transmitter, 1/2  *
 * in receiver).                    *
 * Mult by exp(-jwalpha), see Oppen-*
 * heim and Schafer, p. 244.        *
 *##################################*/
- raisedCosineCoefficients:(double *)coeff number:(int)number
{
  int i, n, half_length, alpha;
  float delta_f, f1, f_delta, f, b;
  float *imag_part, *real_part;

  if(number<2)
    return self;
  alpha     = (number-1)/2;
  imag_part = (float *)malloc(number*sizeof(float));
  real_part = (float *)malloc(number*sizeof(float));
  delta_f   = fs/(float)number;
  f1        = (1.-ROLLOFF)*fc;
  b         = (1.+ROLLOFF)*fc;
  f_delta   = ROLLOFF*fc;
  f         = 0.;
  half_length = number/2+1;
  if(number%2)
    i           = half_length - 1;
  else
    i           = half_length - 2;
  for(n=0; n<half_length; n++)
  {
    if(f<f1)
      real_part[n] = 1.;
    else if( (f>f1) && (f<b) )
      real_part[n] = sqrt(0.5*(1.+cos(PI*(f-f1)/2./f_delta)));
    else if( (f>b) && (f<fs/2.) )
      real_part[n] = 0.0;
    f += delta_f;
  }
  for(n=half_length; n<number; n++)
    real_part[n] = real_part[i--];
  for(n=0; n<number; n++)
    imag_part[n] = 0.;
/* Now, do inverse FFT */
  fft(real_part, imag_part, number, INVERSE_FFT);
/* Do, a time delay, for causal filter */
  i = 0;
  for(n=half_length; n<number; n++)
    coeff[i++] = (double)real_part[n];
  for(n=0; n<half_length; n++)
    coeff[i++] = (double)real_part[n];
  free(real_part);
  free(imag_part);
  return self;
}

@end