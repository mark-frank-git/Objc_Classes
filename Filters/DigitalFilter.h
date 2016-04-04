/************************************************************************
 * This subclass of analogFilter implements a digital filtering object  *
 *                                                                      *
 * File:DigitalFilter.h                                                 *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *          b[n] + b[n-1]z**-1 + ... + b[0]z**-(n)                      *
 *    H(z) = ------------------------------------                       *
 *          1    + a[n-1]z**-1 + ... + a[0]z**-(n)                      *
 *                                                                      *
 *           (z-zero[0]) * (z-zero[1]) ... (z-zero[n_zero])             *
 *    H(z) = ----------------------------------------------             *
 *           (z-pole[0]) * (z-pole[1]) ... (z-pole[n_pole])             *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/12/92  - Started                                              *
 *  2. 01/28/93  - Add interpolate, decimate function                   *
 *  3. 09/12/93  - Add -zeroOutTaps                                     *
 ************************************************************************/
#import "AnalogFilter.h"

#define DIGITAL_FILTER  0                       /* filterType       */
#define ANALOG_FILTER   1

#define HAMMING     0                           /* windowType       */
#define HANNING     1
#define RECTANGULAR 2
#define TRIANGULAR  3
#define BLACKMAN    4
#define BLACKMAN_HARRIS 5

@interface DigitalFilter:AnalogFilter
{
  int   filterType;                             /* digital or analog        */
  int   windowType;                             /* Window function      */
  int   interpolateFactor, decimateFactor;
  float *window;                                /* Window function      */
  double fs;                                    /* sampling frequency       */
  double *doubleShift;                          /* Real shift register      */
  COMPLEX *complexShift;                        /* complex shift register   */
  COMPLEX *digitalPoles, *digitalZeros;         /* Complex poles and zeros  */
}

/**********************
 * Initialization:    *
 **********************/
- init;
- initWithPassType:(int)type fo:(double)centerFreq  fc:(double)cutoffFreq
            fs:(double)samplingFreq order:(int)order;
        
/**********************
 * Running:           *
 **********************/
- filterFloatArray: (float *)input numberPts: (int)numberPts;
- (double) filterFloatData: (float)input;
- (FCOMPLEX) filterComplexData: (FCOMPLEX)input;
- filterComplexArray:(FCOMPLEX *)x numberPts:(int)numberPts;
- findWarpFactorsFor:(double)ratio;
- (double)checkWarpFor:(double)ratio;
- (float *)interpolateInput:(float *)input numberPts:(int)numberPts;
- (float *)decimateInput:(float *)input numberPts:(int)numberPts;

/**********************
 * Set parameters:    *
 **********************/
- setDigitalFilterType: (int)type;
- setWindowingType: (int)type;
- setFilterFo: (double)centerFreq fc:(double)cutoff fs: (double)samplingFreq;
- zeroOutTaps;

/**********************
 * Get parameters:    *
 **********************/
- (int)getFilterType;
- (int)getInterpolateFactor;
- (int)getDecimateFactor;
- (COMPLEX *)getPoles;
- (COMPLEX *)getZeros;

/**********************
 * Converting:        *
 **********************/
- findPolesZeros;
- sToz: (double)sReal: (double)sImag: (double *)zReal: (double *)zImag;
- transferForFIRWindow;
- findTransferFunction;

/**********************
 * Getting Outputs:   *
 **********************/
- (float *)filterResponseAtFrequencies:(float *)frequencies
                           numberFreq:(int)number;
- (float *)digitalResponseAtFrequencies:(float *)frequencies
                           numberFreq:(int)number;
- (float *)windowFunction:(int)length;
- raisedCosineCoefficients:(double *)coeff number:(int)number;

@end