/****************************************************************************
 * This subclass of object emulates fading channel characteristics.         *
 *                                                                          *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/Fading.m                      *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 11/23/93 - Started                                                   *
 ***************************************************************************/

#import "Fading.h"
#import "Numerical.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

#define CORRECTION_FACTOR   5.57    /* See Hata and Nagatsu: "Mobile location ..." */
#define TOLERANCE           0.0001
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


@implementation Fading

/*###############################*
 * Initialization:               *
 *###############################*/
- initWithType:(int)type doppler:(float)doppler
{
  
  [super init];

// 
// Initialize instance variables:
//
  integrator        = [[Numerical alloc] init];
  fadingType        = type;
  dopplerFrequency  = doppler;
  correctionFactor  = CORRECTION_FACTOR;

  return self;
}

/*##############################*
 * These methods set parameters *
 *##############################*/
- setDoppler:(float)doppler
{
  dopplerFrequency = doppler;
  return self;
}

- setFading:(int)type
{
  fadingType = type;
  return self;
}

- setCorrection:(float)factor
{
  correctionFactor = factor;
  return self;
}

/*##############################*
 * These methods get parameters *
 *##############################*/
- (int) fadingType          {return fadingType;}
- (float) dopplerFrequency  {return dopplerFrequency;}
- (float) correctionFactor  {return correctionFactor;}

/*###############################*
 * Return the variance of the   *
 * sample mean of the fading    *
 * signal strength.             *
 * See Hata and Nagatsu.        *
 *###############################*/
- (float) meanVarianceForTime:(float)time
{
  float variance;

  detectionTime = time;
  if(detectionTime>0.)
  {
    variance = [integrator integrateFrom:-detectionTime to:detectionTime withTol:TOLERANCE
                from:self];
    variance *= correctionFactor*correctionFactor/detectionTime/detectionTime;
  }
  else
    variance = correctionFactor;
  return variance;
}

/*###############################*
 * Return the integrand evalu-  *
 * ated at the point, x.        *
 * See Hata and Nagatsu.        *
 *###############################*/
- (double)functionToIntegrate:(double)x
{
  double bessel_arg, function_value;

  bessel_arg = TWOPI*dopplerFrequency*x;
  function_value = j0_bessel(bessel_arg);
  function_value *= function_value*(detectionTime-ABS(x));

////////  function_value = cos(x);

  return function_value;
}

- (double)functionToIntegrate2:(double)x
{
  return sin(x);
}

- (double) integrateSinFrom:(double)x1 to:(double)x2
{
  detectionTime = x2;
  return correctionFactor*correctionFactor*[integrator integrateFrom:x1 to:x2 withTol:TOLERANCE
                from:self]/detectionTime/detectionTime;
}



@end
