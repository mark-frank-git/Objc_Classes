/****************************************************************************
 * This subclass of object returns sampled autocorrelation functions for    *
 * various signals.                                                         *
 *                                                                          *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/Autocorrelation.h             *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 01/26/94 - Started                                                   *
 *                                                                          *
 * NOTE: For the noise autocorrelation, the NoB factor is not included here.*
 ***************************************************************************/

#import <objc/Object.h>

#define IDEAL_FILTERED_NOISE    0               /* ideal filtered noise                 */
#define ONE_POLE_FILTERED_NOISE 1               /* 1 pole filtered noise                */
#define RECTANGULAR_WAVE        2               /* Rectangular pulse shape binary data */



@interface Autocorrelation:Object
{
  int   signalType;
  float sampleTime;                         /* Time between samples             */
  float symbolTime;                         /* Symbol duration                  */
  float bandwidth;                          /* filter 3 dB 2 sided bw in Hz     */
  float piBTs;                              /* Pi * B * sampleTime              */
}


/*********************************
 * Initialization:               *
 *********************************/
- initWithType:(int) type bandwidth:(float)bw sampleTime:(float)time symbolTime:(float)sTime;

/*******************************
 * These methods set parameters:*
 *******************************/
- setSignalType:(int)type;
- setBandwidth: (float)bw;
- setSampleTime: (float)time;
- setSymbolTime: (float)time;

/*******************************
 * These methods get parameters*
 *******************************/
- (int)   signalType;
- (float) bandwidth;
- (float) sampleTime;
- (float) symbolTime;

/*******************************
 * These methods get outputs:  *
 *******************************/
- (float) autocorrelationAt:(int)index;

@end