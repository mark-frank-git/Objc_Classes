/************************************************************************
 * This object is a subclass of the NISTSound object. It provides for   *
 *  adding noise to a sound object.  The data from noise.snd is in:     *
 *  *floatNoise.  The data from signal.snd is in: *soundData.           *
 *  (Integer data for both is elsewhere).                               *
 *                                                                      *
 * File path: /User/frank/Objc_Classes/ProcessedSound                   *
 * File name: NoisySound.h                                              *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/13/93  - Started with skeleton by MF                          *
 *  2. 04/16/93  - First version by SG                                  *
 *  3. 09/27/93  - Add addNoise                                         *
 ************************************************************************/
#import "NISTSound.h"

#define NOISELESS_CLASSIFIER    0               /* classify prior to adding noise       */
#define NOISY_CLASSIFIER        1               /* classify after adding noise          */

@interface NoisySound:NISTSound
{
  id    noise;                                  /* Raw noise sound                      */
  id    playSound;                              /* Sound object for playing             */
  int   noiseSamples;                           /* # of noise samples                   */
  int   noiseOffset;                            /* Offset into noise data for SNR       */
  int   hiEnergyFrameCt;                        /* nbr frames counted in signal energy  */
  int   speechNoiseClassifierType;              /* Type of classification, see above    */
  int   *cleanClassifier;                       /* Speech or noise from clean data.     */
  BOOL  addNoise;                               /* Yes = add noise to input sound.      */
  float *floatNoise;                            /* float point of above                 */
  double snr;                                   /* signal-to-noise ratio                */
  double noiseEnergy;                           /* ave. energy per sample               */
  double signalEnergy;                          /* ave. energy per non-silence sample   */
  double threshold;                             /* ave. energy per signal sample necessary to 
                                                    include frame's energy in signal energy */
  char  noiseFile[255];                         /* name of noise file                   */
}

/********************************
 * Initialization & free:       *
 * - initialize instance vars   *
 * and call super's init.       *
 * - free floatNoise            *
 ********************************/
- initWithSoundWindow: (id)newSoundWindow;
- free;

/********************************
 * These methods set or         *
 * calculate instance variables *
 *******************************/
- setAddNoise:  (BOOL)flag;
- setSpeechNoiseClassifier: (int)type;
- setThreshold: (float)threshold;
- setSignalToNoiseRatio: (float)signalToNoise;
- setNoiseOffset:(int)offset;
- setSpectralSubtraction:(BOOL)flag;
- setSpectralFloor:(float)floor;
- setNoiseFactor:(float)factor;
- (int)calculateEnergies;

/*******************************
 * These methods return         *
 * instance variables           *
 *******************************/
- (int)noiseSamples;
- (int)noiseOffset;
- (float)snr;
- (char *)noiseFile;
- (float)threshold;
- (float)noiseEnergy;
- (float)signalEnergy;
- (int)hiEnergyFrameCt;
- (float)spectralFloor;
- (float)noiseFactor;

/****************************
 * These methods process    *
 * the sound and noise      * 
 ****************************/
- (BOOL) addNoiseToSound;
- playNoisySound;
- (BOOL) calculateFloatFeatureVectors;              /* Override super's method  */

/*********************************
 * These methods read/write to  *
 * disk:                        *
 *********************************/
- (BOOL) readSndFile:(const char *)fileName;        /* Override super's method  */
- (BOOL) readMikeFile:(const char *)fileName;       /* Override super's method  */
- (BOOL) readWavFile:(const char *)fileName;        /* Override super's method  */
- (int) readNoiseFile:(const char *)fileName;
- (BOOL) savePrcFile: (char *)processedFile;
- (BOOL) writeProcessedHeader: (FILE *)filePtr: (int)numberFrames:
                              (int) vectorSize;
- writeProcessedData: (FILE *)filePtr: (int)outputCount:
                              (char *)outputData;

@end
