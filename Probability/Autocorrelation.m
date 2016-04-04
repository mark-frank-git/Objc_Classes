/****************************************************************************
 * This subclass of object returns sampled autocorrelation functions for    *
 * various signals.                                                         *
 *                                                                          *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/Autocorrelation.h             *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 01/26/94 - Started                                                   *
 ***************************************************************************/

#import "Autocorrelation.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define MIN_TIME    1.e-30


@implementation Autocorrelation

/*###############################*
 * Initialization:               *
 *###############################*/
- initWithType:(int)type bandwidth:(float)bw sampleTime:(float)time symbolTime:(float)sTime
{
  
  [super init];

// 
// Initialize instance variables:
//
  signalType        = type;
  bandwidth         = bw;
  sampleTime        = MAX(MIN_TIME, time);
  symbolTime        = sTime;
  piBTs             = PI*bandwidth*sampleTime;

  return self;
}

/*##############################*
 * These methods set parameters *
 *##############################*/
- setSignalType:(int)type
{
  signalType = type;
  return self;
}

- setBandwidth:(float)bw
{
  bandwidth = bw;
  piBTs     = PI*bandwidth*sampleTime;
  return self;
}

- setSampleTime:(float)time
{
  sampleTime = time;
  sampleTime = MAX(MIN_TIME, time);
  piBTs      = PI*bandwidth*sampleTime;
  return self;
}

- setSymbolTime:(float)time
{
  symbolTime = time;
  symbolTime = MAX(MIN_TIME, time);
  return self;
}

/*##############################*
 * These methods get parameters *
 *##############################*/
- (int)   signalType        {return signalType;}
- (float) bandwidth         {return bandwidth;}
- (float) sampleTime        {return sampleTime;}
- (float) symbolTime        {return symbolTime;}

/*###############################*
 * Return the autocorrelation   *
 * according to the signal type *
 *###############################*/
- (float) autocorrelationAt:(int)index
{
  float time, auto_corr, temp;

  index = ABS(index);
  auto_corr = 1.;
  switch(signalType)
  {
    case IDEAL_FILTERED_NOISE:
      if(index == 0)
        break;
      temp = piBTs*(float)index;
      auto_corr = sin(temp)/temp;
      break;
    case ONE_POLE_FILTERED_NOISE:
      auto_corr = PI2*exp(-piBTs*(float)index);
      break;
    case RECTANGULAR_WAVE:
    default:
      time = index*sampleTime;
      if(time>symbolTime)
        return 0.;
      auto_corr = 1. - time/symbolTime;
      break;
  }
  return auto_corr;
}

@end
