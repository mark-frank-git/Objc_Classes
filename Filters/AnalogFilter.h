/************************************************************************
 * This subclass of object implements an analog filter object       *
 *                                                                  *
 * File:AnalogFilter.h                                              *
 *                                                                  *
 * The filter is stored in two forms:                               *
 *                                                                  *
 *          b[n] + b[n-1]z**-1 + ... + b[0]s**-(n)                  *
 *    H(s) = ------------------------------------                   *
 *          1    + a[n-1]z**-1 + ... + a[0]s**-(n)                  *
 *                                                                  *
 *           (s-zero[0]) * (s-zero[1]) ... (s-zero[n_zero])         *
 *    H(s) = ----------------------------------------------         *
 *           (s-pole[0]) * (s-pole[1]) ... (s-pole[n_pole])         *
 *                                                                  *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                  *
 * Revision history:                                                *
 *  1. 08/12/92  - Started                                          *
 ************************************************************************/
#import <objc/Object.h>
#import "c_headers.h"

#define PI          3.141592654                 /* PI                   */
#define SQRT2       1.414213562                 /* sqrt(2)              */
#define TWOPI       6.283185307                 /* 2*PI                 */  

#define LOW_PASS    0                           /* filterPassType       */
#define HIGH_PASS   1
#define BAND_PASS   2
#define BAND_STOP   3

#define TRANSFER_FUNCTION   0                   /* filterStructureType  */
#define POLE_ZERO           1
#define FIR_WINDOW          2
#define FIR_RAISED_COS      3                   /* 40% sqrt raised cos  */

#define BUTTERWORTH         0                   /* analogType           */
#define CHEBYSHEV           1
#define BESSEL              2

#define MAGNITUDE           0                   /* response type        */
#define DB_MAGNITUDE        1
#define PHASE               2
#define PHASE_DELAY         3
#define GROUP_DELAY         4
#define POLE_ZERO_PLOT      5
#define SIGNAL_RESPONSE     6
#define SIGNAL_FFT          7
#define ONE_MINUS_MAG       8                   /* |1-H(w)|             */

@interface AnalogFilter:Object
{
  id    aPolyObject, bPolyObject;   /* polynomials for xfer fnction */
  int   filterPassType;             /* Type of pass band        */
  int   filterStructureType;        /* transfer, pole-zero, etc.    */
  int   analogType;                 /* Butter, cheby, etc.      */
  int   responseType;               /* magnitude, phase, etc.   */
  int   filterOrder;                /* Order of the filter      */
  int   numberPoles, numberZeros;   /* # of poles not always = order*/
  BOOL  realAxisPoles;              /* real poles from bp, bs?  */
  double fo, fc;                    /* center, cutoff       */
  double thetaOld, subAngle, deltaOmega;/* Used in phase calculations   */
  double passBandGain;              /* Divisor for unity pass gain  */

  COMPLEX *analogPoles, *analogZeros;   /* Complex poles and zeros  */
}

/**********************
 * Initialization:    *
 **********************/
- init;
- initWithRealCoeffs:(double *)aCoeffs b:(double *)bCoeffs order:(int)order;
- initWithPassType:(int)type fo:(double)centerFreq fc:(double)cutoffFreq
            order:(int)order;
- initPolesZeros;
  

/**********************
 * Running:           *
 **********************/

/**********************
 * Set parameters:    *
 **********************/
- setPassType: (int)type;
- setFilterStructureType: (int)type;
- setAnalogType: (int)type;
- setResponseType: (int)type;
- setFilterOrder: (int)order;
- setFilterFo: (double)centerFreq fc:(double)cutoff;
- setFilterACoeff: (double *)a bCoeff:(double *)b;

/**********************
 * Get parameters:    *
 **********************/
- (COMPLEX *)getAnalogPoles;
- (COMPLEX *)getAnalogZeros;
- (int) getFilterOrder;
- (int) getNumberOfPoles;
- (int) getNumberOfZeros;
- (double) fo;
- (double) fc;

/**********************
 * Converting:        *
 **********************/
- findPolesZeros;
- findAnalogPolesZerosFor:(double)wc omegaO:(double)wo bandwidth:(double)bw;
- transferFromPoles:(COMPLEX *)poles andZeros:(COMPLEX *)zeros;
- butterPrototype;
- chebyPrototype;
- besselPrototype;
- lowPassToLowPass: (double)wc;
- lowPassToHighPass: (double)wc;
- lowPassToBandPass: (double)wo bw:(double)bw;
- lowPassToBandStop: (double)wo bw:(double)bw;

/**********************
 * Getting Outputs:   *
 **********************/
- (float *)analogResponseAtFrequencies:(float *)frequencies
                          numberFreq:(int)number;
- (COMPLEX) poleZeroResponseAt: (COMPLEX)point poles:(COMPLEX *)poles
                         zeros:(COMPLEX *)zeros;
- (COMPLEX) transferResponseAt:(COMPLEX)point;
- (float)outputResponseFor:(COMPLEX *)complexResponse at:(float)omega;

/**********************
 * Miscellaneous:     *
 **********************/
- (double)functionToFindRoot:(double)x;
 
 
@end