/***********************************************************************
 * This subclass of object implements a model used for HMM speech       *
 * recognition.  For training the model, see the subclass,              *
 * ViterbiModel.                                                        *
 *                                                                      *
 * File: HMMModel.m                                                     *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 08/24/93  - Abstracted from ViterbiModel.                        *
 ***********************************************************************/

#import "HMMModel.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>
#import <string.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

@implementation HMMModel

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

  useCodebookIndices = YES;
  codebookDimensions = (int *)NXZoneMalloc(myZone,codebooks*sizeof(int));
  logB               = (float ***)NXZoneMalloc(myZone, codebooks*sizeof(float **));
  for(i=0; i<codebooks; i++)
  {
    codebookDimensions[i] = dimensions[i];
    codebookWeights[i]    = weights[i];
    logB[i]               = (float **)NXZoneMalloc(myZone,
                                    modelLength*sizeof(float *));
    for(j=0; j<modelLength; j++)
      logB[i][j]              = (float *)NXZoneMalloc(myZone,
                                 codebookDimensions[i]*sizeof(float));
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

// 2 dimensional arrays:
  logA  = (float **)NXZoneMalloc(myZone,modelLength*sizeof(float *));
  for(i=0; i<modelLength; i++)
    logA[i]   = (float *)NXZoneMalloc(myZone, modelLength*sizeof(float));
  for(i=0; i<MAX_FRAMES; i++)
    psi[i] = (int *)NXZoneMalloc(myZone, modelLength*sizeof(int ));

// 1 dimensional arrays
  modelStates       = (int *)NXZoneMalloc(myZone, modelLength*sizeof(int));
  
  logPiInitial    = (float *)NXZoneMalloc(myZone, modelLength*sizeof(float));
  logPiFinal      = (float *)NXZoneMalloc(myZone, modelLength*sizeof(float));
  stateProbability= (float *)NXZoneMalloc(myZone, modelLength*sizeof(float));
  maxProbability  = (float *)NXZoneMalloc(myZone, modelLength*sizeof(float));
  codebookWeights = (float *)NXZoneMalloc(myZone, codebooks*sizeof(float));

  return self;
}

/*#######################################*
 * Free up the allocated arrays and my-  *
 * self.                                 *
 *#######################################*/
- free
{
  [self freeArrays];
  return [super free];
}

/*######################################*
 * Free up alloc'd arrays               *
 *######################################*/
- freeArrays
{
  int i, j;

// 2 dimensional arrays:
  for(i=0; i<modelLength; i++)
    NXZoneFree(myZone, logA[i]);
  for(i=0; i<MAX_FRAMES; i++)
    NXZoneFree(myZone, psi[i]);
  NXZoneFree(myZone, logA);
  NXZoneFree(myZone, psi);

// 1 dimensional arrays
  NXZoneFree(myZone, modelStates);
  NXZoneFree(myZone, logPiInitial);
  NXZoneFree(myZone, logPiFinal);
  NXZoneFree(myZone, stateProbability);
  NXZoneFree(myZone, maxProbability);
  NXZoneFree(myZone, codebookWeights);

// 3-d arrays
  if(codebookDimensions != NULL)
  {
    NXZoneFree(myZone, codebookDimensions);
    for(i=0; i<codebooks; i++)
    {
      for(j=0; j<modelLength; j++)
        NXZoneFree(myZone, logB[i][j]);
      NXZoneFree(myZone, logB[i]);
    }
    NXZoneFree(myZone, logB);
  }
  return self;
}

/*#######################################*
 * Initialize for a new search:          *
 *#######################################*/
- initNewSearch
{
  int  i;
  
  for( i=0; i<modelLength; i++ )
    stateProbability[i] = logPiInitial[i];
  return self;
}

/*#######################################*
 * Step the model:                      *
 * Assume a left->right model.          *
 *#######################################*/
- recognitionStepWithProbs:(double **)logProbabilities count:(int)frameCount
{
  int     i,j, max_transition;
  double  transition_probability;

  if(frameCount>=MAX_FRAMES)
  {
    printf("Too many frames in ViterbiModel step\n");
    return self;
  }
    
/************************
 * Find max probability *
 * score over all in-   *
 * coming transitions   *
 ************************/
  if(frameCount == 0)
  {
    for(i=0; i<modelLength; i++)
      for(j=0; j<codebooks; j++)
        stateProbability[i] += codebookWeights[j]*
                               logProbabilities[j][modelStates[i]];
  }
  else
  {
    for(i=0; i<modelLength; i++ )
    {
      max_transition    = 0;
      maxProbability[i] = logA[0][i] + stateProbability[0];
      for(j=1; j<=i; j++)
      {
        transition_probability = logA[j][i]+stateProbability[j];
        if(transition_probability > maxProbability[i])
        {
          maxProbability[i] = transition_probability;
          max_transition    = j;
        }
      }
      psi[frameCount-1][i]  = max_transition;   /* save for backtrace */
    }
/************************
 * Save max probability *
 * for each state:      *
 ************************/
    for(i=0; i<modelLength; i++)
    {
      stateProbability[i]    = maxProbability[i];
      for(j=0; j<codebooks; j++)
        stateProbability[i] += codebookWeights[j]*
                               logProbabilities[j][modelStates[i]];
    }
  }

  return self;
}

/*#######################################*
 * Step the model with codebook indices *
 * Assume a left->right model.      *
 *#######################################*/
- recognitionStepWithIndices:(short *)codebookIndices count:(int)frameCount
{
  int     i,j, max_transition;
  double  transition_probability;

  if(frameCount>=MAX_FRAMES)
  {
    printf("Too many frames in ViterbiModel step\n");
    return self;
  }
    
/************************
 * Find max probability *
 * score over all in-   *
 * coming transitions   *
 ************************/
  if(frameCount == 0)
  {
    for(i=0; i<modelLength; i++)
      for(j=0; j<codebooks; j++)
        stateProbability[i] += codebookWeights[j]*
                                logB[j][i][codebookIndices[j]];
  }
  else
  {
    for(i=0; i<modelLength; i++ )
    {
      max_transition    = 0;
      maxProbability[i] = logA[0][i] + stateProbability[0];
      for(j=1; j<=i; j++)
      {
        transition_probability = logA[j][i]+stateProbability[j];
        if(transition_probability > maxProbability[i])
        {
          maxProbability[i] = transition_probability;
          max_transition    = j;
        }
      }
      psi[frameCount-1][i]  = max_transition;   /* save for backtrace */
    }
/************************
 * Save max probability *
 * for each state:      *
 ************************/
    for(i=0; i<modelLength; i++)
    {
      stateProbability[i] = maxProbability[i];
      for(j=0; j<codebooks; j++)
        stateProbability[i] += codebookWeights[j]*
                               logB[j][i][codebookIndices[j]];
    }
  }

  return self;
}

/*######################################*
 * Return the model's backtrace:    *
 *######################################*/
- (int *)getBackTrace:(int)frameCount 
{
  int    i;
  int    model_end;
  
  if(frameCount > maxFrames)
  {
    maxFrames = frameCount;
    if(backTrace != NULL)
      free(backTrace);
    backTrace = (int *)malloc((frameCount+1)*sizeof(int));
  }
  if(frameCount >= MAX_FRAMES)
  {
    printf("Too many frames in HMMModel getBackTrace\n");
    return NULL;
  }
/******************************
 * Construct the path for the *
 * model with the best score  *
 * from psi.                  *
 ******************************/
  backTrace[frameCount-1] = model_end = bestModelEnd;
  for(i=frameCount-2; i>=0; i-- )
    backTrace[i] = model_end = psi[i][model_end];
  
  return backTrace;
}


/*#######################################*
 * Return Some of our instance variables:*
 *#######################################*/
- (char *)  getModelDescription         {  return modelDescription;}
- (int)     getModelLength              {  return modelLength;}
- (double)  getLastScore                {  return maxScore;}
- (float *) getTransitionProbabilities  {  return &logA[0][0];}
- (float *) getInitialProbabilities     {  return &logPiInitial[0];}
- (float *) getFinalProbabilities       {  return &logPiFinal[0];}
- (float *) stateProbability            { return stateProbability;}

/*#######################################*
 * Return the current score:             *
 *#######################################*/
- (double) getModelScore:(int)frameCount
{
  int    j;
  double new_probability;
  
  maxScore    = stateProbability[0] + logPiFinal[0];
  bestModelEnd = 0;
  for(j=1; j<modelLength; j++ )
  {
    new_probability = stateProbability[j] + logPiFinal[j];
    if(new_probability>maxScore)
    {
      maxScore    = new_probability;
      bestModelEnd = j;
    }
  }
  return maxScore;
}

/*######################################*
 * Set a new model length:              *
 *######################################*/
- setModelLength:(int)length states:(int *)states
{
  int i;

  modelLength = length;
  [self freeArrays];
  [self allocateArrays];
  if(states != NULL)
   for(i=0; i<modelLength; i++)
    modelStates[i] = states[i];
  
  return self;
}

/*######################################*
 * Save probabilities to disk:      *
 *######################################*/
- (BOOL)saveProbabilitiesToFile:(FILE *)fp;
{
  int i, j;

// Initial state probabilities:
  fwrite(logPiInitial, sizeof(float), modelLength, fp);

// Final state probabilities:
  fwrite(logPiFinal, sizeof(float), modelLength, fp);

// Transition probabilities:
  for(i=0; i<modelLength; i++)
    fwrite(logA[i], sizeof(float), modelLength, fp);

// Output probabilities:
  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
      for(j=0; j<modelLength; j++)
         fwrite(logB[i][j], sizeof(float), codebookDimensions[i], fp);
  }

  return YES;
}

/*######################################*
 * Read probabilities from disk:    *
 *######################################*/
- (BOOL)readProbabilitiesFromFile:(FILE *)fp;
{
  int i, j;

// Initial state probabilities:
  fread(logPiInitial, sizeof(float), modelLength, fp);

// Final state probabilities:
  fread(logPiFinal, sizeof(float), modelLength, fp);

// Transition probabilities:
  for(i=0; i<modelLength; i++)
    fread(logA[i], sizeof(float), modelLength, fp);

// Output probabilities:
  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
      for(j=0; j<modelLength; j++)
         fread(logB[i][j], sizeof(float), codebookDimensions[i],fp);
  }

  return YES;
}

/*#######################################*
 * Print model probabilities:         *
 *#######################################*/
- printModel:(FILE *)fp
{
  int i, j, k;

  if(useCodebookIndices)
  {
    for(i=0; i<codebooks; i++)
    {
      fprintf(fp,"bij for codebook %d\n", i);
      for(j=0; j<modelLength; j++)
      {
        for(k=0; k<codebookDimensions[i]; k++)
          fprintf(fp,"%6.3e ", exp(logB[i][j][k]));
        fprintf(fp,"\n\n");
      }
    }
  }
  fprintf(fp,"aij:\n");
  for(i=0; i<modelLength; i++)
  {
    for(j=0; j<modelLength; j++)
      fprintf(fp,"%5.3f ", exp(logA[i][j]));
    fprintf(fp,"\n");
  }
  fprintf(fp,"\nPiI = ");
  for(i=0; i<modelLength; i++)
    fprintf(fp,"%f ", exp(logPiInitial[i]));
  fprintf(fp,"\n");
  fprintf(fp,"\nPiF = ");
  for(i=0; i<modelLength; i++)
    fprintf(fp,"%f ", exp(logPiFinal[i]));
  fprintf(fp,"\n");
  return self;
}


@end

