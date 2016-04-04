/***************************************************************************
 * This object is a subclass of the processedSound object. It provides      *
 * additional processing of NIST .wav files.                                *
 *                                                                          *
 * File: /User/frank/Objc_Classes/ProcessedSound/NISTSound.h                *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 07/02/92  - Started                                                  *
 *  2. 09/18/92  - Add DigitalFilter object                                 *
 *  3. 01/27/93  - Add playConvertedSound and saveSndFile: methods.         *
 *  4. 04/26/93  - Delete -findStrobes with addition of speechProcessor     *
 *  5. 06/03/93  - Remove fc                                                *
 ***************************************************************************/

#import "ProcessedSound.h"

/**************************************
 * The following define the radio     *
 * buttons for outputSamplingType     *
 **************************************/
#define SAMPLE_CODEC    0                                   /* 8012 Hz */
#define SAMPLE_LOW      1                                   /* 22050 Hz  */
#define SAMPLE_HIGH     2                                   /* 44100 Hz  */
#define SAMPLE_SAME     3                                   /* output =  input  */

#define UCHAR       unsigned char

/**************************************
 * The following define the types     *
 * of low pass filters:               *
 **************************************/
#define BUTTERWORTH     0
#define HAMMING_SINC    1

#define DEFAULT_POLES   50

#define MAX_PHONEMES    100                             /* Max phonemes in a .wav file  */
#define MAX_PHONE_LENGTH 4                              /* phoneme descriptor length    */
#define PHONE_FILE      ".phn"                          /* extension of phoneme file    */
#define WAVE_FILE       ".wav"                          /* extension of .wav file   */

@interface NISTSound:ProcessedSound
{
  id    phonemeMap;                                     /* Converts phonemes to integers*/
  id    convertedSound;                                 /* Converted for play       */
  int   filterType;                                     /* LPF for sampling conversion  */
  int   numberPoles;                                    /* LPF poles            */
  int   outputSamplingType;                             /* Type of output sampling  */
  int   numberPhonemes;
  int   phonemeStartSamples[MAX_PHONEMES];              /* start of phonemes    */
  int   phonemeMiddleSamples[MAX_PHONEMES];             /* middle of phonemes   */
  int   phonemeEndSamples[MAX_PHONEMES];                /* end of phonemes  */
  char  phonemeList[MAX_PHONEMES][MAX_PHONE_LENGTH];    /* list of phones   */
  float offset;                                         /* LPF cutoff freq parameter    */
  BOOL  phonemeLabelsAvailable;                         /* phoneme file in directory?   */
}

/********************************
 * Initialization:              *
 ********************************/
- initWithSoundWindow: (id)newSoundWindow;

/********************************
 * These methods get info from  *
 * the NIST header:             *
 ********************************/
- (int) NISTdataFormat;
- (int) NISTdataSize;
- (int) NISTsamplingRate;
- (int) NISTchannelCount;
- (int) NISTinfoSize;
- (int) NISTnumberSamples;
- (char *) NISTsoundInfo;
- (int) NISTsampleMin;
- (int) NISTsampleMax;

/********************************
 * These methods get info from  *
 * the converted sound object   *
 ********************************/
- (int) convertInfoSize;
- (char *) convertSoundInfo: (float)sampleFreq;
- (int) soundSampleMin;
- (int) soundSampleMax;

/********************************
 * These methods set instance   *
 * variables:                   *
 ********************************/
- setFilterType: (int)type;
- setNumberPoles: (int)poles;
- setOutputSamplingType: (int)type;
- setOffset: (float)newOffset;

/*******************************
 * Higher level sound procesing*
 *******************************/
- convertSound;
- playConvertedSound;

/*********************************
 * These methods open and read   *
 * data files:                   *
 *********************************/
- (BOOL) readWavFile:(const char *)fileName;
- (BOOL) saveSndFile:(const char *)fileName;

/*******************************
 * These methods adjust the    *
 * selected speech range.      *
 *******************************/
- findEndpoints;

/*******************************
 * These methods process the   *
 * sound                       *
 *******************************/
- swapBytes: (char *)data: (int)no_points:
                            (int)data_format;
- convert: (UCHAR *)data length:(int)no_points format:(int)data_format
                         fromFloat:(float *)fdata;
- (int) dataWidth: (int)soundFormat;


/********************************
 * Dummy prototype              *
 ********************************/
- setSamplingFrequency: (float)frequency;

@end
