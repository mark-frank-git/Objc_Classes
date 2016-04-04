/****************************************************************************
 * This subclass of object calculates power measurement statistics.         *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/PowerMeasurement.m            *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 02/15/94 - Started                                                   *
 *  2. 08/21/95 - Updated to only consider fading effects.                  *
 ***************************************************************************/

#import "PowerMeasurement.h"
#import "Numerical.h"
#import "Probability.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

#define NUMBER_STD_DEVIATIONS   4.
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


@implementation PowerMeasurement

- init
{
  [super init];
  integrator        = [[Numerical alloc] init];
  probability[0]    = [ [Probability alloc] init];
  probability[1]    = [ [Probability alloc] init];
  hysteresis        = 2.;
  
  return self;
}

/*###############################*
 * Initialization:               *
 *###############################*/
- initWithHysteresis:(float)h
{
  
  [self init];

// 
// Initialize instance variables:
//
  hysteresis        = h;

  return self;
}

/*##############################*
 * These methods set parameters *
 *##############################*/
- setHysteresis:(float)h
{
  hysteresis = h;
  return self;
}

- setBeam0Mean:(float)mean andVariance:(float)var
{
  double vars[1];

  means[0]  = mean;
  vars[0]   = var;
  stdDev[0] = sqrt(var);
  [probability[0] setMean:means];
  [probability[0] setVariance:vars];
  return self;
}

- setBeam1Mean:(float)mean andVariance:(float)var
{
  double vars[1];

  means[1]  = mean;
  vars[0]   = var;
  stdDev[1] = sqrt(var);
  [probability[1] setMean:&means[1]];
  [probability[1] setVariance:vars];
  return self;
}
 
/*##############################*
 * These methods get parameters *
 *##############################*/
- (float)   hysteresis          {return hysteresis;}


#define NUMBER_STDEVS   10.                 /* Number of std deviations for limiting cases */
/*##################################*
 * Return the handoff probability   *
 * where power0 = old beam,         *
 * power1 = candidate beam.         *
 *##################################*/
- (float) handoffProbability
{
  float handoff_prob;

// Check limiting cases first:
  if( ((means[0]-NUMBER_STDEVS*stdDev[0])-(means[1]+NUMBER_STDEVS*stdDev[1])) > hysteresis )
    handoff_prob = 0.;
  else if ( ((means[1]-NUMBER_STDEVS*stdDev[1])-(means[0]+NUMBER_STDEVS*stdDev[0])) > hysteresis )
    handoff_prob = 1.;
  else
  {
    if(means[0]>means[1])
    {
      xMin = means[1] - NUMBER_STD_DEVIATIONS*stdDev[0];
      xMax = means[0] + NUMBER_STD_DEVIATIONS*stdDev[1];
    }
    else
    {
      xMin = means[0] - NUMBER_STD_DEVIATIONS*stdDev[1];
      xMax = means[1] + NUMBER_STD_DEVIATIONS*stdDev[0];
    }
    handoff_prob = [integrator integrateFrom:xMin to:xMax from:self];
  }
  return handoff_prob;
}

/*##################################*
 * Return the the inside integrand  *
 * evaluated at the point, x.       *
 *###############################*/
- (double)functionToIntegrate:(double)x
{
  double f1, f2;
  f1 =  qx((x+hysteresis-means[1])/stdDev[1]);
  f2 = [probability[0] densityFunctionAtX:&x];
  return f1*f2;
}


@end
