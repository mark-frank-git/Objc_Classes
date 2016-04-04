/************************************************************************
 * This object is a subclass of Object.  It provides processing of      *
 * speech signals.  Some of the code originated from Bill Kushner of    *
 * CSRL.                                                                *
 *                                                                      *
 * File: /User/frank/Objc_Classes/ProcessedSound/SpeechProcessor.h      *
 *                                                                      *
 * The speech is processed on a frame-by-frame basis.  The following    *
 * types of data are calculated:                                        *
 *   1. spectral data: e.g., cepstral, filter bank outputs, etc.        *
 *   2. zero crossing data                                              *
 *   3. peak to peak amplitude                                          *
 *   4. frame energy                                                    *
 *   5. speech/noise classification                                     *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/20/93 - Converted from perceptual.c                           *
 *  2. 06/17/93 - Added lpcMagnitudeResponse                            *
 *  3. 07/12/93 - Add spectral subtraction instance variables.          *
 *  4. 09/21/93 - Add oldCepstralLength, oldFrameLength                 *
 ************************************************************************/

/*************************************
 * Constants:                        *
 *************************************/
#define MAX_FEATURES    6

#import <objc/Object.h>

@interface SpeechProcessor:Object
{
  int   spectralType;               /* Type of spectral data            */
  int   windowType, lifterType;     /* Types of window, lifter          */
  int   lpcOrder;                   /* Order of LPC filter              */
  int   cepstralLength;             /* # of cepstral coefficients       */
  int   frameLength;                /* Length of speech frame           */
  int   oldCepstralLength, oldFrameLength;
  int   overlapLength;              /* Length of frame overlap          */
  int   powerOf2Length;             /* nearest power of 2 frame lngth   */
  int   numberFilterBanks;          /* # of perceptual filter banks     */
  int   numberFrames;               /* # of frames in input speech      */
  int   spectralLength;             /* Length of spectral data/frame    */
  int   *speechOrNoiseFromSpectral; /* array, 1 = speech, 0=noise       */
  int   *speechOrNoiseFromWave;     /* array, 1 = speech, 0=noise       */
  BOOL  perceptualProcessing;       /* YES = perceptual filtering       */
  BOOL  spectralSubtraction;        /* YES = perform subtraction        */
  float samplingFrequency;          /* Sampling rate                    */
  float foldingFrequency;           /* Nyquist frequency                */
  float spectralFloor, noiseFactor; /* Used in spectral subtraction     */
  float *window;                    /* Speech window function           */
  float *lifter;                    /* Cepstral lifter                  */
  float *spectralData;              /* Speech spectral data             */
  float *ptpData;                   /* Speech PTP data                  */
  float *zeroCrossData;             /* Speech zero crossing data        */
  float *energyData;                /* Speech log energy data           */
  float *filterBanks, *loudnessWeighting, *centerFrequencies;
                                    /* Used in perceptual processing    */
}

/*********************************
 * Initialization:               *
 *********************************/
- initWithOrder:(int)order cepstralLength:(int)clength
               frameLength:(int)flength overlapLength:(int)olength
               numberBanks:(int)numberBanks samplingFrequency:(float)frequency;
- initializeSpectralLength;
- initializeWindows;
- initializeFilterBanks;

/*******************************
 * These methods return         *
 * instance variables:          *
 *******************************/
- (int)spectralType;
- (int)windowType;
- (int)lifterType;
- (int)lpcOrder;
- (int)cepstralLength;
- (int)frameLength;
- (int)overlapLength;
- (int)powerOf2Length;
- (int)numberFilterBanks;
- (int)numberFrames;
- (int)spectralLength;
- (int *)speechOrNoise;
- (BOOL)perceptualProcessing;
- (BOOL)spectralSubtraction;
- (float *)windowData;
- (float *)lifterData;
- (float *)spectralData;
- (float *)ptpData;
- (float *)zeroCrossData;
- (float *)energyData;
- (float)samplingFrequency;
- (float)spectralFloor;
- (float)noiseFactor;
- (char *)spectralName;
- (char *)lifterName;
- (char *)windowName;
- (char *) frequencyAnalysisName;

/********************************
 * These methods set instance   *
 * variables:                   *
 *******************************/
- setSpectralType:      (int)type;
- setWindow:            (int)type;
- setLifter:            (int)type;
- setLPCFilterOrder:    (int)order;
- setCepLength:         (int)length;
- setSpeechFrameLength: (int)length;
- setOverlapLength:     (int)length;
- setNumberFilterBanks: (int)number;
- setPerceptualProcessing:(BOOL)flag;
- setSpectralSubtraction:(BOOL)flag;
- setSpectralFloor:     (float)floor;
- setNoiseFactor:       (float)factor;
- setSamplingFrequency: (float)frequency;


/*********************************
 * Speech Processing methods:    *
 *********************************/
- (float *)calculateSpectralData:(float *)speechInput start:(int)start
           length:(int)length;
- (float *)calculatePTPData:(float *)speechInput start:(int)start
           length:(int)length;
- (float *)calculateZeroCrossData:(float *)speechInput start:(int)start
           length:(int)length;
- (float *)calculateEnergyData:(float *)speechInput start:(int)start
           length:(int)length;
- (int *)  calculateSpeechOrNoiseFromEnergy:(float *)energyInput numberFrames:(int)frames;

/************************************
 * Local methods:                   *
 ************************************/
- (float *)windowSpeechData:(float *)speech length:(int)length;
- (float *)powerSpectrumOf:(float *)speech length:(int)length
               frameEnergy:(float *)frameEnergy spectralCenter:(float *)scm;
- (float *)powerSpectrumFromAuto:(float *)autocorr length:(int)length;
- (float *)filterBankOutputs:(float *)powerSpectrum length:(int)length;
- (float *)autoFromFilterBanks:(float *)filterOutputs length:(int)length;
- (float *)autoFromSpeech:(float *)speech length:(int)length
                                    lagLength:(int)lagLength;
- (float *)lpcFromAuto:(float *)autoCorrelation order:(int)order;
- (float *)cepstralFromLPC:(float *)lpcCoefficients lpcOrder:(int)order 
                    length:(int)outputLength;
- (float *)cepstralFromAuto:(float *)autoCorrelation lpcOrder:(int)order
                    length:(int)outputLength;
- (float *)cepstralFromSpectrum:(float *)spectrum length:(int)outputLength;
- (float *)lpcMagnitudeResponse:(float *)lpcCoeffs lpcOrder:(int)order
                    digitalFrequencies:(float *)omega numberPoints:(int)number;
- (int)speechOrNoise:(float)frameEnergy frame:(int)frame;
- subtractNoiseFrom:(float *)spectrum noise:(float *)noiseSpectrum
          speechFlag:(int)speech length:(int)length
         frameEnergy:(float *)frameEnergy;


@end
