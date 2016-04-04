/************************************************************************
 * This object is a subclass of the Sound object.  In addition to the   *
 * to the normal sound methods, it provides for processing of the sound *
 *                                                                      *
 * File: /User/frank/Objc_Classes/ProcessedSound/ProcessedSound.h       *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 06/17/92 - Started                                               *
 *  2. 11/17/92 - Add numberBytes for output data                       *
 *              - Add intScaleFactor.                                   *
 *  3. 12/15/92 - Add getWindowData, getLifterData                      *
 *  4. 02/04/93 - Add extendEndpointsBy:                                *
 *  5. 02/26/93 - Add getIntScaleFactor                                 *
 *  6. 04/16/93 - Move savePrcFile to NoisySound                        *
 *  7. 04/26/93 - Add SpeechProcessor object                            *
 *                Delete strobe stuff.                                  *
 *  8. 05/28/93 - Combined maxs and mins into arrays.                   *
 *  9. 09/27/93 - Move findFeatureVectors inside calculateFloat...      *
 ************************************************************************/

/*************************************
 * Constants:                        *
 *************************************/
#define SUB_VECTORS 5       /* Spectral, z/c, ptp, engy, s/n    */

//
// inputSoundType definitions:
//
#define SPEECH_FROM_SND     0                           /* From a .snd file             */
#define SPEECH_FROM_PRC     1                           /* From a .prc file             */
#define SPEECH_FROM_WAV     2                           /* From a NIST .wav file        */
#define SPEECH_FROM_MIKE    3                           /* From NeXT microphone         */
#define NOISE_SND           4                           /* open a noise file            */

#define STRING_LEN      255

#import <soundkit/soundkit.h>

@interface ProcessedSound:Sound
{
  id  digitalFilter;                                    /* used for preemphasis         */
  id  soundWindow;                                      /* Displays the sound           */
  id  speechProcessor;                                  /* Does the processing          */
  int inputSoundType;                                   /* input sound origination      */
  int numberBytes;                                      /* # of storage bytes           */
  int numberSoundSamples;                               /* From sound header file       */
  int numberPRCFrames;                                  /* # of frames in PRC file      */
  int subvectorCount[SUB_VECTORS];                      /* Size of each subvector       */
  int subvectorMin[SUB_VECTORS];                        /* Min value of each subvect    */
  int subvectorMax[SUB_VECTORS];                        /* Max value of each subvect    */
  int selectionStart, selectionLength;                  /* Sound selection range        */
  int featureVectorLength;                              /* feature vector array lng     */
  int minPRC[SUB_VECTORS];                              /* mins of subvectors fr .prc   */
  int maxPRC[SUB_VECTORS];                              /* maxs of subvectors fr .prc   */
  float minVector[SUB_VECTORS];                         /* mins of subvectors fr .snd   */
  float maxVector[SUB_VECTORS];                         /* maxs of subvectors fr .snd   */
  float preemphasisCoefficient;                         /* Sound preemphasis            */
  unsigned char  *rawSoundData;                         /* Input from sound file        */
  char  *charFeatureVectors;                            /* For saving to .prc files     */
  int   *intFeatureVectors;                             /* vectors for int calcs        */
  float *floatFeatureVectors;                           /* Speech feature vectors       */
  float *soundData;                                     /* Pre-emphasized sound         */
    
  BOOL floatVectorsNeeded;                              /* YES to calculate floats      */
  BOOL *forcedRange, *linearMapping;                    /* For feature vector calc.     */
  double outputSamplingRate;                            /* New sampling rate            */
  double subvectorScale[SUB_VECTORS];                   /* scale factor for each sub    */
  double subvectorOffset[SUB_VECTORS];                  /* and offset for each subvec   */
  double intScaleFactor;                                /* float -> int scale factor    */
  struct header_t *nistHeader;                          /* NIST header for file         */
  char  utteranceId[STRING_LEN];
  char  inputSpeechFile[STRING_LEN];
}

/*********************************
 * Initialization:               *
 *********************************/
- initWithSoundWindow: (id)newSoundWindow;

/*******************************
 * These methods respond to    *
 * other objects request for   *
 * information:                *
 *******************************/
- (id)speechProcessor;
- (int)getInputSoundType;
- (int)getVectorLength;
- (int)getNumberBytes;
- (int)getNumberFrames;
- (int)getSelectionLength;
- (int)getSelectionStart;
- (int)getFrameLength;
- (int)getOverlapLength;
- (int)getCepstralLength;
- (int)getLPCOrder;
- (int)getFilterBanks;
- (int) getIntScaleFactor;
- (int *)getIntFeatureVector;
- (char *)getUtteranceId;
- (float) getMinVector:(int)index;
- (float) getMaxVector:(int)index;
- (float) getSamplingFrequency;
- (float) getPreemphasisCoefficient;
- (float *)getSoundData;
- (float *)getFloatFeatureVectors;
- (char *) windowName;
- (char *) spectralName;
- (char *) lifterName;
- (char *) frequencyAnalysisName;

/*******************************
 * These methods set instance  *
 * variables.                  *
 *******************************/
- setMinVector:     (int)index to:(float)value;
- setMaxVector:     (int)index to:(float)value;
- setNumberBytes:   (int)number;
- setPreemphasis:   (float)coefficient;
- setFrequencyType: (int)type;
- setLPCOrder:      (int)order;
- setLiftering:     (int)type;
- setWindowType:    (int)type;
- setProcessType:   (int)type;
- setFilterBanks:   (int)number;
- setFrameLength:   (int)length;
- setOverlapLength: (int)length;
- setCepstralLength:    (int)length;
- setInputSoundType:    (int)type;
- setFloatVectorsNeeded: (BOOL)flag;

/*********************************
 * These methods open and read   *
 * data files:                   *
 *********************************/
- (BOOL) readSndFile:(const char *)fileName;
- (BOOL) readMikeFile:(const char *)fileName;
- (BOOL) readPrcFile:(const char *)fileName;
- setDataFromSoundStructure;

/*********************************
 * These methods get info from   *
 * NIST header:                  *
 *********************************/
- getNISTHeaderData;
- getSubVectorMinMax: (int) vectorNumber: (int *)min: (int *)max;
- getNISTInteger: (char *)name: (int *)data;
- getNISTFloat: (char *)name: (float *)data;
- getNISTChar: (char *)name: (char *)data;

/*******************************
 * These methods adjust the    *
 * selected speech range.      *
 *******************************/
- findEndpoints;
- resetEndpoints;
- extendEndpointsBy:(int)samples;
- (BOOL) getSelectionRange;
- setSelectionRange: (int)start: (int)length;

/*********************************
 * Signal Processing methods:    *
 *********************************/
- (BOOL)calculateFloatFeatureVectors;
- (float *) rawSoundToFloat;
- floatVectorsToInt;
- floatVectorsToChar;
- intVectorsToChar;

@end
