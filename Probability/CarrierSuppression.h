/****************************************************************************
 * This subclass of object calculates QAM carrier suppression statistics.   *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/CarrierSuppression.h          *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 07/11/97 - Started                                                   *
 ***************************************************************************/

#import <objc/Object.h>


@interface CarrierSuppression:Object
{
  id        integrator;                     // Numerical integrator

  int       numberSamples;                  // N

  BOOL      integrateVariance;              // YES = integrate to find variance
  BOOL      integratePDF;                   // YES = integrate to find PDF
  BOOL      integrateMeanPDF;               // YES = integrate to find mean of PDF

  double    noiseVariance;                  // sigma^2
  double    signalPower;                    // signal power
  double    noncentralVariance;             // noncentral chi-square variance
  double    iqOffset;                       // s^2
  double    xMin, xMax;                     // Integration limits
  double    machineEps;                     // machine smallest double
}

/*********************************
 * Initialization:               *
 *********************************/
- init;

/*******************************
 * Local methods:               *
 *******************************/
- findNoncentralVariance;

/*******************************
 * These methods set parameters:*
 *******************************/
- setNumberSamples: (int) samples;
- setNoiseVariance: (double) variance;
- setIQOffset:      (double) sSquared;
- setSignalPower:   (double) power;

/*******************************
 * These methods get parameters*
 *******************************/

/*******************************
 * These methods get outputs:  *
 *******************************/
- (double) getNoncentralVariance;
- (double) outputVariance;
- (double) outputVarianceFrom:(double)xLo to:(double)xHi;
- (double) outputExpectedValue;
- (double) outputExpectedValueFrom:(double)xLo to:(double)xH;
- (double) outputMeanOfPDFFrom:(double)xLo to:(double)xH;
- (double) pdfFrom:(double)xLo to:(double)xHi;

/*******************************
 * Miscellaneous                *
 ********************************/
- (double)pdf:(double)x;
- (double)pdfCS:(double)x;
- (double)functionToIntegrate:(double)x;

@end