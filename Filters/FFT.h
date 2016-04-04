/************************************************************************
 * This subclass of object implements an object for computing fast      *
 * Fourier transforms, as well as DFTs and DCTs                         *
 *                                                                      *
 * File: FFT.h                                                          *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 09/03/93  - Started                                              *
 ************************************************************************/

#import <objc/Object.h>

#define COMPLEX_MAGNITUDE_DB    0               /* Complex to real conversions  */
#define COMPLEX_MAGNITUDE   1
#define COMPLEX_PHASE       2
#define COMPLEX_REAL_PART   3
#define COMPLEX_IMAG_PART   4

#define MEAN_FROM_SCM       0               /* Calc mean index from spectral center mass    */
#define MEAN_FROM_THRESHOLD 1               /* Calc from threshold dB down                  */

#define DB_THRESHOLD        15.             /* Find mean index from this threshold          */

@interface FFT:Object
{
  id    digitalFilter;                      /* Used for windowing           */
  int   numberPoints, oldSpectrumPoints;    /* Number of FFT points         */
  int   oldAutoPoints;
  int   meanIndexType;                      /* Method to calc mean index    */
  float *outputSpectrum;                    /* results of FFT               */
  float *autocorrelation;                   /* Autocorrelation function     */
  float *frequencyPoints;                   /* frequencies of FFT samples   */
  float *twoSidedFrequency;                 /* Negative and positive freqs  */
  float samplingFreq;                       /* Sampling frequency           */
  float twoSidedBandwidth;                  /* Calc'd in meanIndexOf        */
}

/**********************
 * Initialization:    *
 **********************/
- init;
- initWithSamplingFreq:(float)frequency;
- allocateMemory;
- freeMemory;
- free;

/**********************
 * Set parameters:    *
 **********************/
- setSamplingFrequency:(float)frequency;
- setWindowType:(int)type;
- setMeanIndexType:(int)type;

/**********************
 * Get parameters:    *
 **********************/
- (float *)frequencyPoints;
- (float  )deltaFrequency;
- (float  )twoSidedBandwidth;
- (float *)twoSidedFrequency;

/**********************
 * Calculating:       *
 **********************/
- (float *)powerSpectrum:(float *)input numberPoints:(int)number;
- (float *)fftMagnitude:(float *)input numberPoints:(int)number;
- (float *)cosineTransform:(float *)input numberPoints:(int)number;
- (float *)psdFromAuto:(float *)input numberPoints:(int)number;
- (float *)avgPowerSpectrum:(float *)realInput imag:(float *)imagInput
               numberPoints:(int)number numberAvgs:(int)avgs;
- (float *)autoFromPSD:(float *)inputSpectrum numberPoints:(int)number type:(int)type
                 scale:(float *)scale;
- (float *)autocorrelationOf:(float *)realInput imag:(float *)imagInput
          numberPoints:(int)number windowLength:(int)length lagLength:(int)lag type:(int)type
          scale:(float *)scale;
- (float  )meanIndexOf:(float *)inputSpectrum size:(int)points;
- shiftRight:(float *)inputArray by:(int)shift size:(int)points;
- dBSpectrum:(float *)inputSpectrum size:(int)points;
- normalizeAuto:(float *)autoc size:(int)points by:(float)scale;
@end