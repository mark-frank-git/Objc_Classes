/****************************************************************************
 * This subclass of object emulates fading channel characteristics.         *
 *                                                                          *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/Fading.h                      *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 11/23/93 - Started                                                   *
 ***************************************************************************/

#import <objc/Object.h>

#define RAYLEIGH                0
#define RICIAN                  1
#define LOG_NORMAL_FADING       2



@interface Fading:Object
{
  id    integrator;
  int   fadingType;
  float dopplerFrequency;                   /* doppler shift frequency          */
  float correctionFactor;                   /* Correction from Hata's paper     */
  float detectionTime;                      /* Integration time                 */
}


/*********************************
 * Initialization:               *
 *********************************/
- initWithType:(int) type doppler:(float)doppler;

/*******************************
 * These methods set parameters:*
 *******************************/
- setDoppler:(float) doppler;
- setFading: (int)type;
- setCorrection:(float)factor;

/*******************************
 * These methods get parameters*
 *******************************/
- (int) fadingType;
- (float) dopplerFrequency;
- (float) correctionFactor;

/*******************************
 * These methods get outputs:  *
 *******************************/
- (float) meanVarianceForTime:(float)detectionTime;
- (double)functionToIntegrate:(double)x;
- (double)integrateSinFrom:(double)x1 to:(double)x2;

@end