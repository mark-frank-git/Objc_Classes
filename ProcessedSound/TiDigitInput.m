/************************************************************************
 * This subclass of object inputs processed speech files from the TI    *
 * digit data base.  It also contains a tapped delay lines for cal-     *
 * culating delta cepstral, delta energy values.                        *
 *                                                                      *
 * File: TiDigitInput.m                                                 *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 02/19/93 - Started                                               *
 *  2. 09/17/93 - Dynamically allocate speechFiles.                     *
 ************************************************************************/

/*************************************
 * Constants:                        *
 *************************************/
#define PRINT_ON        1
#define TRAIN_VALUE     100
#define MAX_TRY         10                  /* Max # of speech files to */
                                            /* try before giving up */
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))


#import "TiDigitInput.h"
#import "ProcessedSound.h"
#import "DelayLine.h"
#import "c_headers.h"
#import <string.h>
#import <math.h>
#import <ctype.h>



@implementation TiDigitInput


- init { return [super init];}

/*###################################*
 * Initialize instance variables     *
 *###################################*/
- initWithDeltas:(int *)newDeltas number:(int)number types:(int *)types
          models:(int)models states:(int)states
{
  int i;

  [self init];
// 
// Instantiate objects
//
  processedSound    = [[ProcessedSound alloc] initWithSoundWindow:nil];
  [processedSound setFloatVectorsNeeded:NO];
  numberCodebooks   = number;
  numberStates      = MAX(1,states);
  numberModels      = models;
  outputSize        = numberStates*numberModels;
  oldIndex          = 0;
  speechScale       = 1;
  deltas            = (int *)malloc(numberCodebooks*sizeof(int));
  dataTypes         = (int *)malloc(numberCodebooks*sizeof(int));
  inputSize         = (int *)malloc(numberCodebooks*sizeof(int));
  maxDelta          = 0;
  for(i=0; i<numberCodebooks; i++)
  {
    deltas[i]       = newDeltas[i];
    maxDelta        = MAX(maxDelta, deltas[i]);
    dataTypes[i]    = types[i];
    if( (dataTypes[i] = types[i]) == CEPSTRAL_DATA)
      inputSize[i]  = CEPSTRAL_SIZE;
    else
      inputSize[i]  = ENERGY_SIZE;
  }

//
// Initialize instance variables:
//
  trainValue            = TRAIN_VALUE;
  winnerListAvailable   = NO;
  fileCounter           = numberFiles  = trainFileCounter = 0;
  numberFrames          = frameCounter = 0;
  speechData            = NULL;
  
// Allocate arrays
  [self allocateArrays];

  return self;
}

/*##############################*
 * Allocate arrays and objects: *
 *##############################*/
- allocateArrays
{
  int i, j, training_size;

  trainingData = (int **)malloc(numberCodebooks*sizeof(int *));
  floatTrainingData = (float **)malloc(numberCodebooks*sizeof(float *));
  for(i=0; i<numberCodebooks; i++)
  {
    training_size   = inputSize[i] + outputSize;
    trainingData[i] = (int *)malloc(training_size*sizeof(int));
    floatTrainingData[i] = (float *)malloc(training_size*sizeof(float));
    for(j=0; j<training_size; j++)
    {
      trainingData[i][j] = 0;
      floatTrainingData[i][j] = 0.;
    }
  }
  cepstralDelayLine = [[DelayLine alloc] initDelayLineTaps:2*maxDelta+1 
                              width:CEPSTRAL_SIZE dataType:INT_TYPE];
  energyDelayLine   = [[DelayLine alloc] initDelayLineTaps:2*maxDelta+1 
                              width:ENERGY_SIZE dataType:INT_TYPE];
  return self;
}

/*##############################*
 * Free up arrays and objects:  *
 *##############################*/
- freeArrays
{
  int i;

  for(i=0; i<numberCodebooks; i++)
  {
    free(trainingData[i]);
    free(floatTrainingData[i]);
  }
  for(i=0; i<numberFiles; i++)
    free(speechFiles[i]);
  if(numberFiles)
    free(speechFiles);
  free(trainingData);
  free(floatTrainingData);
  [cepstralDelayLine free];
  [energyDelayLine   free];
  return self;
}

/*##################################*
 * Free the instantiated objects:   *
 *##################################*/
- free
{

  [self freeArrays];
  free(deltas);
  free(dataTypes);
  free(inputSize);
  [processedSound free];

  return self;
}

/*##################################*
 * Re-start the data at the first   *
 * speech data file. Don't open the *
 * file.                            *
 *##################################*/
- initializeData
{
  frameCounter  = fileCounter = trainFileCounter = 0;
  numberFrames  = 0;
  return self;
}

/*###################################*
 * Create the plot window and        *
 * initialize instance variables     *
 *###################################*/
- setNewDeltas:(int *)newDeltas types:(int *)types
{
  int i;

  [self freeArrays];
  maxDelta     = 0;
  for(i=0; i<numberCodebooks; i++)
  {
    deltas[i]    = newDeltas[i];
    maxDelta    = MAX(maxDelta, deltas[i]);
    dataTypes[i] = types[i];
    if( (dataTypes[i] = types[i]) == CEPSTRAL_DATA)
      inputSize[i]  = CEPSTRAL_SIZE;
    else
      inputSize[i]  = ENERGY_SIZE;
  }
  [self allocateArrays];
  
  return self;
}

/*##################################*
 * Open a new speech file       *
 *##################################*/
- (BOOL)openNewSpeechFile
{
  int i;
  char *utterance_id;
  float temp;

  if( (numberFiles==0)||(![processedSound readPrcFile:speechFiles[fileCounter]]) )
  {
    printf("Can't open speech file %s\n", speechFiles[fileCounter]);
    if(++fileCounter == numberFiles)
      fileCounter = 0;
    return NO;
  }
  if(++fileCounter == numberFiles)
    fileCounter = 0;
  trainFileCounter++;
  speechData        = [processedSound getIntFeatureVector];
  speechVectorSize  = [processedSound getVectorLength];
  numberFrames      = [processedSound getNumberFrames];
  temp              = [processedSound getIntScaleFactor];
  if(temp>0.)
    speechScale    = ROUND(trainValue/temp);
  frameCounter     = 0;
  stateIncrement   = numberFrames/numberStates;
// Find the digit:
  utterance_id     = [processedSound getUtteranceId];
  inputModel       = [self getUtteranceDigit:utterance_id];
  
// Fill up the delay line
  for(i=0; i<(maxDelta+1); i++)
  {
    [cepstralDelayLine writeDelay:&speechData[CEPSTRAL_OFFSET]];
    [energyDelayLine   writeDelay:&speechData[ENERGY_OFFSET]];
  }
  for(i=0; i<maxDelta; i++)
  {
    [cepstralDelayLine writeDelay:&speechData[CEPSTRAL_OFFSET]];
    [energyDelayLine   writeDelay:&speechData[ENERGY_OFFSET]];
    speechData += speechVectorSize;
  }
  
  return YES;
}

/*#################################*
 * The following methods return    *
 * instance variables.         *
 *#################################*/
- (int) numberFrames            { return numberFrames;}
- (int) inputSize:(int)number   { return inputSize[number];}
- (int) numberFiles             { return numberFiles;}
- (int) fileCounter             { return fileCounter;}
- (int) trainFileCounter        { return trainFileCounter; }
- (BOOL)speechDataAvailable     {  return ((speechData != NULL));}
- (BOOL)winnerListAvailable     { return winnerListAvailable;}
- (int) modelNumber             { return inputModel; }
- (int) modelState              { return frameCounter/stateIncrement; }
- (short **)winningNeuronsForFile:(int)fileNumber { return winningNeurons[fileNumber];}
- (int) lastFileCounter
{
  int count;
  count = fileCounter - 1;
  if(count < 0)
    count += numberFiles;
  return count;
}
- (char *)lastFile      
{
  int count;
  count = fileCounter - 1;
  if(count < 0)
    count += numberFiles;
  return speechFiles[count];
}
/*#################################*
 * Check if the current frame is   *
 * speech or noise.  Return YES    *
 * if it is noise.                 *
 * NOTE: In implementing this      *
 * method we assume stepData has   *
 * just been called, so we need to *
 * decrement the speech pointer.   *
 * We also assume, that the speech *
 * or noise data is the last       *
 * element of the speechData vector*
 *#################################*/
- (BOOL) currentNoiseFrame
{
  char subvector_name[255];
  int  *speech_ptr;
  [processedSound getNISTChar:"name_subvector_4":subvector_name];
  if( strncmp(subvector_name, "Speech_Or_Noise",15) == 0)
  {
     speech_ptr = speechData;
     speech_ptr--;
     if( *speech_ptr > 0.5)
        return NO;      /* classified as speech */
     else
        return YES;     /* classified as noise */
  }
  printf("Can't find Speech_Or_Noise data in speech file\n");
  return NO;
}

/*##############################*
 * Return the training data:    *
 *##############################*/
- (int **)trainingOutput
{
  int i, j, state, train_index, tap;
  int *speech_data, *delay_output;
  id  delay_line;

// Find the supervised training data:
  state = frameCounter/stateIncrement;
  state = MIN(state, numberStates-1);
  train_index = inputModel*numberStates + state;
  for(i=0; i<numberCodebooks; i++)
  {
    if(dataTypes[i] == CEPSTRAL_DATA)
      delay_line = cepstralDelayLine;
    else
      delay_line = energyDelayLine;
    tap = maxDelta-deltas[i];
    speech_data = [delay_line readTap: tap];
    for(j=0; j<inputSize[i]; j++)
       trainingData[i][j] = speechScale*speech_data[j];
    if(deltas[i]>0)
    {
      delay_output = [delay_line readTap: maxDelta+deltas[i]];
      for(j=0; j<inputSize[i]; j++)
        trainingData[i][j] -= speechScale*delay_output[j];
    }
    trainingData[i][inputSize[i]+oldIndex]    = 0;
    trainingData[i][inputSize[i]+train_index] = trainValue;
  }
  oldIndex = train_index;
  return trainingData;
}

/*##############################*
 * Return the last training data:*
 *##############################*/
- (int **)oldTrainingOutput
{
  return trainingData;
}

/*##############################*
 * Return the test data:    *
 *##############################*/
- (int **)testOutput
{
  int i, j, tap;
  int *speech_data, *delay_output;
  id  delay_line;

  for(i=0; i<numberCodebooks; i++)
  {
    if(dataTypes[i] == CEPSTRAL_DATA)
      delay_line = cepstralDelayLine;
    else
      delay_line = energyDelayLine;
    tap = maxDelta-deltas[i];
    speech_data = [delay_line readTap: tap];
    for(j=0; j<inputSize[i]; j++)
       trainingData[i][j] = speechScale*speech_data[j];
    if(deltas[i]>0)
    {
      delay_output = [delay_line readTap: maxDelta+deltas[i]];
      for(j=0; j<inputSize[i]; j++)
        trainingData[i][j] -= speechScale*delay_output[j];
    }
  }
  return trainingData;
}

/*##############################*
 * Return the training data:    *
 *##############################*/
- (float **)trainingFloatOutput
{
  int i, j, state, train_index, tap;
  int *speech_data, *delay_output;
  id  delay_line;

// Find the supervised training data:
  state = frameCounter/stateIncrement;
  state = MIN(state, numberStates-1);
  train_index = inputModel*numberStates + state;
  for(i=0; i<numberCodebooks; i++)
  {
    if(dataTypes[i] == CEPSTRAL_DATA)
      delay_line = cepstralDelayLine;
    else
      delay_line = energyDelayLine;
    tap = maxDelta-deltas[i];
    speech_data = [delay_line readTap: tap];
    for(j=0; j<inputSize[i]; j++)
       floatTrainingData[i][j] = speechScale*speech_data[j];
    if(deltas[i]>0)
    {
      delay_output = [delay_line readTap: maxDelta+deltas[i]];
      for(j=0; j<inputSize[i]; j++)
        floatTrainingData[i][j] -= speechScale*delay_output[j];
    }
    floatTrainingData[i][inputSize[i]+oldIndex]    = 0;
    floatTrainingData[i][inputSize[i]+train_index] = trainValue;
  }
  oldIndex = train_index;
  return floatTrainingData;
}

/*##############################*
 * Return the last training data:*
 *##############################*/
- (float **)oldTrainingFloatOutput
{
  return floatTrainingData;
}

/*##############################*
 * Return the test data:    *
 *##############################*/
- (float **)testFloatOutput
{
  int i, j, tap;
  int *speech_data, *delay_output;
  id  delay_line;

  for(i=0; i<numberCodebooks; i++)
  {
    if(dataTypes[i] == CEPSTRAL_DATA)
      delay_line = cepstralDelayLine;
    else
      delay_line = energyDelayLine;
    tap = maxDelta-deltas[i];
    speech_data = [delay_line readTap: tap];
    for(j=0; j<inputSize[i]; j++)
       floatTrainingData[i][j] = speechScale*speech_data[j];
    if(deltas[i]>0)
    {
      delay_output = [delay_line readTap: maxDelta+deltas[i]];
      for(j=0; j<inputSize[i]; j++)
        floatTrainingData[i][j] -= speechScale*delay_output[j];
    }
  }
  return floatTrainingData;
}

/*################################*
 * Set the training scale factor  *
 *################################*/
- setLabelScale:(int)scale
{
  trainValue = scale;
  return self;
}

/*################################*
 * Reset the training file counter*
 *################################*/
- resetTrainFileCounter
{
  trainFileCounter = 0;
  return self;
}

/*################################*
 * Set the training data array for*
 * supervised training:       *
 *################################*/
- setTrainingOutputsFor:(int)modelNumber state:(int)state
{
  int i, train_index;

  state = MIN(state, numberStates-1);
  state = MAX(state, 0);
  train_index = modelNumber*numberStates + state;
  for(i=0; i<numberCodebooks; i++)
  {
    trainingData[i][inputSize[i]+oldIndex]         = 0;
    trainingData[i][inputSize[i]+train_index]      = trainValue;
    floatTrainingData[i][inputSize[i]+oldIndex]    = 0;
    floatTrainingData[i][inputSize[i]+train_index] = trainValue;
  }
  oldIndex = train_index;

  return self;
}

/*#################################*
 * Open a new list of speech files *
 * Return NO on error.             *
 *#################################*/
- (BOOL)setSpeechFileList: (char *)fileName
{
  int i;
  FILE *file_list;
  char file_name[MAX_FILE_NAME];

  if( (file_list = fopen(fileName, "r")) == NULL)
  {
    printf("Can't open file list: %s\n", fileName);
    return NO;
  }
  if(numberFiles > 0)
  {
    for(i=0; i<numberFiles; i++)
      free(speechFiles[i]);
    free(speechFiles);
  }
// How many files in the list?
  for(numberFiles=0; ; numberFiles++)
    if(fscanf(file_list, "%s", file_name) == EOF)
      break;
  fclose(file_list);
  speechFiles = (char **)malloc(numberFiles*sizeof(char *));
  for(i=0; i<numberFiles; i++)
    speechFiles[i] = (char *)malloc(MAX_FILE_NAME*sizeof(char));

// Read in the files
  file_list = fopen(fileName, "r");
  for(i=0; i<numberFiles; i++)
    if(fscanf(file_list, "%s", speechFiles[i]) == EOF)
      break;
  fclose(file_list);
  
  fileCounter = 0;
  return [self openNewSpeechFile];
}

/*#################################*
 * Open a list of winning neurons  *
 *#################################*/
- (BOOL)setWinnersFromFile: (char *)fileName
{
  int i, j, k, itemp;
  static int old_number_files = 0;
  FILE *winner_file;
  id   processed_sound;

  if( (winner_file = fopen(fileName, "r")) == NULL)
  {
    printf("Can't open file list: %s\n", fileName);
    return NO;
  }
  if(numberFiles == 0)
  {
    printf("File list hasn't been opened\n");
    return NO;
  }
    
  if(winnerListAvailable)
  {
    for(i=0; i<old_number_files; i++)
    {
      for(j=0; j<numberCodebooks; j++)
        free(winningNeurons[i][j]);
      free(winningNeurons[i]);
    }
    free(winningNeurons);
    free(framesPerFile);
    old_number_files = numberFiles;
  }
  winnerListAvailable = YES;

// Allocate the arrays
  processed_sound   = [[ProcessedSound alloc] initWithSoundWindow:nil];
  [processed_sound setFloatVectorsNeeded:NO];
  framesPerFile     = (short *)malloc(numberFiles*sizeof(short));
  winningNeurons    = (short ***)malloc(numberFiles*sizeof(short **));
  for(i=0; i<numberFiles; i++)
  {
    if(![processed_sound readPrcFile:speechFiles[i]])
      return NO;
    framesPerFile[i]  = [processed_sound getNumberFrames];
    winningNeurons[i] = (short **)malloc(numberCodebooks*sizeof(short *));
    for(j=0; j<numberCodebooks; j++)
      winningNeurons[i][j] = (short *)malloc(framesPerFile[i]*sizeof(short));
  }
  [processed_sound free];

// Read in the winning neurons
  for(i=0; i<numberFiles; i++)
  {
    for(j=0; j<numberCodebooks; j++)
      for(k=0; k<framesPerFile[i]; k++)
      {
        fscanf(winner_file, "%d", &itemp);
        winningNeurons[i][j][k] = itemp;
      }
  }
  fclose(winner_file);
  return YES;
}

/*################################*
 * Step in new input data, return *
 * YES if step OK, return NO, if  *
 * no speech data.                *
 * Notes:                         *
 *  1. We have to flush out the   *
 *     delay line at the end.     *
 *################################*/
- (BOOL)stepData
{
  if(speechData == NULL)
    return NO;
  [cepstralDelayLine writeDelay:&speechData[CEPSTRAL_OFFSET]];
  [energyDelayLine   writeDelay:&speechData[ENERGY_OFFSET]];
  if(frameCounter < numberFrames - maxDelta - 1)
    speechData += speechVectorSize;
  return YES;
}

/*################################*
 * Increment the frame counter.   *
 * Open a new file if we are at   *
 * the end of the old one.        *
 *################################*/
- (BOOL)incrementFrame
{
  int i;

  if(++frameCounter >= numberFrames)
  {
    for(i=0; i<MAX_TRY; i++)
      if([self openNewSpeechFile])
        break;
    if(i==MAX_TRY)
      return NO;
  }
  return YES;
}

/*################################*
 * Convert the utterance id to a  *
 * digit.                 *
 *################################*/
- (int)getUtteranceDigit: (char *)string
{
  char  *digit_ptr, c;
  int   digit, len;
  
  len = strlen(string);
  digit_ptr = string;
  digit = 0;
  while( (digit_ptr - string)<len )
  {
    digit_ptr = index(digit_ptr,'_');
    digit_ptr++;
    c = *digit_ptr;
    if(c == 'o')
    {
      digit = 10;
      break;
    }
    else if(c == 'z')
    {
      digit = 0;
      break;
    }
    digit = c - '1' + 1;
    if( (digit>0) && (digit<10) )
      break;
  }
  return digit;
}


@end