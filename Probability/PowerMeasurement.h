/****************************************************************************
 * This subclass of object calculates power measurement statistics.         *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/PowerMeasurement.h            *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 02/15/94 - Started                                                   *
 *  2. 08/21/95 - Updated to only consider fading effects.                  *
 ***************************************************************************/

#import <objc/Object.h>


@interface PowerMeasurement:Object
{
  id    integrator;
  id    probability[2];
  float hysteresis;                         /* Handoff hysteresis level         */
  double xMin, xMax;                        /* Integration limits               */
  double means[2];                          /* Means out of integrators         */
  double stdDev[2];                         /* Standard deviations of beam 0,1  */
}

/*********************************
 * Initialization:               *
 *********************************/
- init;
- initWithHysteresis:(float)h;

/*******************************
 * These methods set parameters:*
 *******************************/
- setBeam0Mean:     (float)mean andVariance:(float) variance;
- setBeam1Mean:     (float)mean andVariance:(float) variance;
- setHysteresis:    (float)h;

/*******************************
 * These methods get parameters*
 *******************************/
- (float) hysteresis;

/*******************************
 * These methods get outputs:  *
 *******************************/
- (float) handoffProbability;

/*******************************
 * Miscellaneous                *
 ********************************/
- (double)functionToIntegrate:(double)x;

@end