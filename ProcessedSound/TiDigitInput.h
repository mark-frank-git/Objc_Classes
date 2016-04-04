/************************************************************************
 * This subclass of object inputs processed speech files from the TI    *
 * digit data base.  It also contains a tapped delay lines for cal-     *
 * culating delta cepstral, delta energy values.                        *
 *                                                                      *
 * File: TiDigitInput.h                                                 *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 02/19/93 - Started                                               *
 *  2. 09/17/93 - Dynamically allocate speechFiles.                     *
 *              - Add capability of reading in list of winning neurons  *
 ************************************************************************/

/*************************************
 * Constants:                        *
 *************************************/
#define MAX_FILE_NAME       256

#define CEPSTRAL_DATA       0                           /* Cepstral data        */
#define ENERGY_DATA         1                           /* Energy data          */

#define CEPSTRAL_OFFSET     0                           /* index into speech data vector*/
#define ENERGY_OFFSET       9                           /* index into speech data vector*/

#define CEPSTRAL_SIZE       7                           /* Cepstral vector size     */
#define ENERGY_SIZE         1                           /* Energy vector size       */

#import <objc/Object.h>


@interface TiDigitInput:Object
{
  id  processedSound;                                   /* Reads in speech frames       */
  id  cepstralDelayLine;                                /* tapped delay line for cep    */
  id  energyDelayLine;                                  /* delay line for energy        */

  short ***winningNeurons;                              /* List of winning neurons      */
  short *framesPerFile;                                 /* # of speech frames in each file  */
  int numberCodebooks;                                  /* number of output vectors     */
  int *deltas, maxDelta;                                /* delta frame sizes            */
  int *dataTypes;                                       /* cepstral or energy           */
  int *inputSize;                                       /* Speech vector size           */
  int outputSize;                                       /* Train vector size            */
  int frameCounter;                                     /* Current speech frame         */
  int fileCounter;                                      /* Current speech file          */
  int trainFileCounter;                                 /* Counts files for training    */
  int numberFiles;                                      /* Total number of speech files */
  int numberFrames;                                     /* Frames in speech file        */
  int speechVectorSize;                                 /* Size of the speech vector    */
  int numberStates, numberModels;                       /* For HMM-type training        */
  int stateIncrement;                                   /* state increment HMM training */
  int inputModel;                                       /* The current speech model     */
  int trainValue;                                       /* Value for supervised trning  */
  int speechScale;                                      /* Scale for feature map        */
  int *speechData;                                      /* From the speech file         */
  int oldIndex;                                         /* Used in generating train data*/
  int **trainingData;                                   /* For training fm output wgts  */
  BOOL  winnerListAvailable;                            /* YES if winning neurons read in   */
  float **floatTrainingData;                            /* For floating point neural nets*/
  char **speechFiles;                                   /* processed speech file names  */
}

/********************************
 * Initialization:              *
 ********************************/
- initWithDeltas:(int *)newDeltas number:(int)number types:(int *)types
          models:(int)models states:(int)states;
- allocateArrays;
- free;
- freeArrays;
- initializeData;
- setNewDeltas:(int *)newDeltas types:(int *)types;
- (BOOL)openNewSpeechFile;

/********************************
 * These methods get info from  *
 * this object:                 *
 ********************************/
- (int) numberFrames;
- (int) inputSize:(int)number;
- (int) numberFiles;
- (int) fileCounter;
- (int) trainFileCounter;
- (BOOL) speechDataAvailable;
- (BOOL) winnerListAvailable;
- (int) modelNumber;
- (int) modelState;
- (short **) winningNeuronsForFile:(int)fileNumber;
- (int) lastFileCounter;
- (char *) lastFile;
- (BOOL) currentNoiseFrame;
- (int **) trainingOutput;
- (int **) oldTrainingOutput;
- (int **) testOutput;
- (float **) trainingFloatOutput;
- (float **) oldTrainingFloatOutput;
- (float **) testFloatOutput;


/********************************
 * These methods set instance   *
 * variables:                   *
 ********************************/
- setLabelScale:(int)scale;
- resetTrainFileCounter;
- setTrainingOutputsFor:(int)modelNumber state:(int)state;
 
/*******************************
 * Action methods              *
 *******************************/
- (BOOL)setSpeechFileList: (char *)fileName;
- (BOOL)setWinnersFromFile:(char *)fileName;    
- (BOOL)stepData;                                       /* step in new data */
- (BOOL)incrementFrame;


/********************************
 * Local methods                *
 ********************************/
- (int)getUtteranceDigit: (char *)string;

@end
