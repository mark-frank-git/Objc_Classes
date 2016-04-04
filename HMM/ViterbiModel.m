/************************************************************************
 * This subclass of HMMModel allows for training of the model.          *
 *                                                                      *
 *                                                                      *
 * File: ViterbiModel.m                                                 *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/09/92  - Started                                              *
 *  2. 05/27/92  - Use continuous transition probabilities              *
 *  3. 06/01/92  - Free output arrays on next call                      *
 *  4. 12/09/92  - Udates for training on TIMIT phonemes                *
 *  5. 01/15/92  - Enhancements from Ed Srenger's code.                 *
 *               - remove reestimation code, see OldViterbiModel if     *
 *                 needed                                               *
 *  6. 02/12/93  - Add -addSkipState                                    *
 *  7. 02/18/93  - Add weights: to initWithLength:...                   *
 *                 Correct RecognitionStepWithProbs for mult codebooks  *
 *  8. 04/15/93  - Change order of writing/reading probabilities to file*
 *  9. 05/11/93  - Make stateProbability 2-D for score variance metric. *
 * 10. 07/30/93  - Deleted score variance stuff.                        *
 * 11. 08/06/93  - Add -setModelLength:                                 *
 * 12. 08/24/93  - Abstracted out HMMModel                              *
 ***********************************************************************/

#import "ViterbiModel.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>
#import <string.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

@implementation ViterbiModel

/*######################################*
 * Initialize instance variables        *
 *######################################*/
- initWithLength:(int)length name:(char *)name states:(int *)states
                 codebooks:(int)numberCodebooks weights:(float *)weights
         zone:(NXZone *)zone
{

  int i;

  [super init];
  myZone                    = zone;
  modelLength               = length;
  codebooks                 = numberCodebooks;
  strncpy(modelDescription, name, MAX_NAME);
  maxScore                  = 0.;
  codebookDimensions        = NULL;
  backTrace                 = NULL;
  useCodebookIndices        = NO;
  oldState                  = 0;
  maxFrames                 = 0;
  newTrainingFile           = YES;
  [self allocateArrays];
  if(states != NULL)
   for(i=0; i<modelLength; i++)
    modelStates[i] = states[i];
  if(weights != NULL)
   for(i=0; i<codebooks; i++)
    codebookWeights[i] = weights[i];
  
  return self;
}
      
/*#######################################*
 * Initialize the codebook dimensions,   *
 * and accordingly set up the observation*
 * probabilities arrays:                 *
 *#######################################*/
- initCodebookDimensions: (int *)dimensions withWeights:(float *)weights
{
  int i, j;

  [super initCodebookDimensions:dimensions withWeights:weights];
  observationsCount  = (int ***)NXZoneMalloc(myZone,codebooks*sizeof(int **));
  for(i=0; i<codebooks; i++)
  {
    observationsCount[i]  = (int **)NXZoneMalloc(myZone,
                                    modelLength*sizeof(int *));
    for(j=0; j<modelLength; j++)
      observationsCount[i][j] = (int *)NXZoneMalloc(myZone,
                                 codebookDimensions[i]*sizeof(int));
  }
  return self;
}

/*#######################################*
 * Allocate the training, and model     *
 * probability arrays based on model    *
 * length:                              *
 *#######################################*/
- allocateArrays
{
  int i;

  [super allocateArrays];
// 2 dimensional arrays:
  transitionsCount  = (int **)NXZoneMalloc(myZone,modelLength*sizeof(int *));
  for(i=0; i<modelLength; i++)
    transitionsCount[i] = (int *)NXZoneMalloc(myZone,modelLength*sizeof(int ));

// 1 dimensional arrays
  initialStateCount = (int *)NXZoneMalloc(myZone, modelLength*sizeof(int));
  finalStateCount   = (int *)NXZoneMalloc(myZone, modelLength*sizeof(int));

  return self;
}

/*######################################*
 * Free up alloc'd arrays               *
 *######################################*/
- freeArrays
{
  int i, j;

  [super freeArrays];
// 2 dimensional arrays:
  for(i=0; i<modelLength; i++)
    NXZoneFree(myZone, transitionsCount[i]);
  NXZoneFree(myZone, transitionsCount);

// 1 dimensional arrays
  NXZoneFree(myZone, initialStateCount);
  NXZoneFree(myZone, finalStateCount);

// 3-d arrays
  if(codebookDimensions != NULL)
  {
    for(i=0; i<codebooks; i++)
    {
      for(j=0; j<modelLength; j++)
        NXZoneFree(myZone, observationsCount[i][j]);
      NXZoneFree(myZone, observationsCount[i]);
    }
    NXZoneFree(myZone, observationsCount);
  }
  return self;
}
      
/*#######################################*
 * Initialize for a new training run:    *
 *#######################################*/
- initTrainingCounts
{
  int  i, j, k;
  
  for( i=0; i<modelLength; i++ )
  {
    initialStateCount[i] = finalStateCount[i] = 0;
    for(j=0; j<modelLength; j++)
      transitionsCount[i][j] = 0;
  }
// Are we using codebookIndices
  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
    {
      for(j=0; j<modelLength; j++)
      {
        for(k=0; k<codebookDimensions[i]; k++)
           observationsCount[i][j][k] = 0;
      }
    }
  }

  return self;
}

/*#######################################*
 * Set the flag indicating the next frame*
 * is from a new training file.          *
 *#######################################*/
- initTrainingFile
{
  finalStateCount[oldState]++;
  newTrainingFile = YES;
  return self;
}

/*######################################*
 * Update the counts for a step during  *
 * initial training:            *
 *######################################*/
- initialTrainingStep:(int)state indices:(short *)indices
{
  int i;

  if(newTrainingFile)
  {
    newTrainingFile = NO;
    initialStateCount[state]++;
  }
  else
    transitionsCount[oldState][state]++;
  for(i=0; i<codebooks; i++)
    observationsCount[i][state][indices[i]]++;

// For final state count
  oldState = state;

  return self;
}

/*######################################*
 * Update the counts for a training step*
 *######################################*/
- initialTrainingStep: (int)state finalFrame:(BOOL)final
{
  if(newTrainingFile)
  {
    newTrainingFile = NO;
    initialStateCount[state]++;
  }
  else
    transitionsCount[oldState][state]++;
  oldState = state;
  if(final)
    finalStateCount[state]++;
  return self;
}

/*######################################*
 * Add skip counts to transitionsCount  *
 * for initial training:                *
 *######################################*/
- addSkipState
{
  int i, j, row_sum;
  float temp;

  temp = SKIP_PROBABILITY/(1. - SKIP_PROBABILITY);
  for(i=0; i<(modelLength-2); i++)
  {
    row_sum = 0;
    for(j=i; j<modelLength; j++)
      row_sum += transitionsCount[i][j];
    transitionsCount[i][i+2] = ROUND((float)row_sum*temp);
  }

  return self;
}

/*######################################*
 * Calculate the initial transition and *
 * other probabilities from the train-  *
 * ing counts.                          *
 *######################################*/
- updateInitialModel
{
  int i, j, k;

  for(i=0; i<modelLength; i++)
  {
    logPiInitial[i] = (float)initialStateCount[i];
    logPiFinal[i]   = (float)finalStateCount[i];
    for(j=0; j<modelLength; j++)
       logA[i][j] = (float)transitionsCount[i][j];
  }
  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
    {
      for(j=0; j<modelLength; j++)
      {
        for(k=0; k<codebookDimensions[i]; k++)
          logB[i][j][k] = (float)observationsCount[i][j][k];
      }
    }
  }
  [self normalizeProbabilities];
  return self;
}

/*#######################################*
 * Normalize transition, initial, and   *
 * final state probabilities to add to 1*
 * set zero probs -> MIN_PROBABILITY.   *
 * Take log of probabilities after  *
 * normalization.           *
 *#######################################*/
- normalizeProbabilities
{
  int    i,j,k;
  float init_norm, final_norm, row_norm;
  float init_scale, final_scale, row_scale;
  
  init_scale = final_scale = 1.;
  init_norm  = final_norm = 0.;
  for(i=0; i<modelLength; i++)
  {
    if(logPiInitial[i]<MIN_PROBABILITY)
       init_scale  -= MIN_PROBABILITY;
    if(logPiFinal[i]<MIN_PROBABILITY)
       final_scale -= MIN_PROBABILITY;
    init_norm  += logPiInitial[i];
    final_norm += logPiFinal[i];
    row_norm    = 0.;
    row_scale   = 1.;
    for(j=0; j<modelLength; j++)
    {
      if(logA[i][j]<MIN_PROBABILITY)
        row_scale -= MIN_PROBABILITY;
      row_norm += logA[i][j];
    }
    if(row_norm > 0.)
      for(j=0; j<modelLength; j++)
      {
        logA[i][j] *= row_scale/row_norm;
        if(logA[i][j] > MIN_PROBABILITY)
          logA[i][j] = log(logA[i][j]);
        else
          logA[i][j] = LOG_MIN;
      }
  }
  for(i=0; i<modelLength; i++)
  {
    if(init_norm>0.)
    {
      logPiInitial[i] *= init_scale/init_norm;
      if(logPiInitial[i] > MIN_PROBABILITY)
        logPiInitial[i] = log(logPiInitial[i]);
      else
        logPiInitial[i] = LOG_MIN;
    }
    if(final_norm>0.)
    {
      logPiFinal[i]   *= final_scale/final_norm;
      if(logPiFinal[i] > MIN_PROBABILITY)
        logPiFinal[i] = log(logPiFinal[i]);
      else
        logPiFinal[i] = LOG_MIN;
    }
  }
  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
    {
      for(j=0; j<modelLength; j++)
      {
        row_norm  = 0.;
        row_scale = 1.;
        for(k=0; k<codebookDimensions[i]; k++)
        {
          if(logB[i][j][k]<MIN_PROBABILITY)
            row_scale -= MIN_PROBABILITY;
          row_norm  += logB[i][j][k];
        }
        if(row_norm > 0.)
          for(k=0; k<codebookDimensions[i]; k++)
          {
             logB[i][j][k] *= row_scale/row_norm;
             if(logB[i][j][k] > MIN_PROBABILITY)
               logB[i][j][k] = log(logB[i][j][k]);
             else
               logB[i][j][k] = LOG_MIN;
          }
      }
    }
  }

          
  return self;
}


@end

