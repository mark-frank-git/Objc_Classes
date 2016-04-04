/***********************************************************************
 * This subclass of object implements an HMM recognition algorithm.     *
 * It uses the ViterbiModel object to solve the problem.  Future        *
 * versions could use either ViterbiModel or BaumWelchModel.            *
 * Note that it can operate in two modes depending on the output of the *
 * vector quantizer: 1) The VQ outputs codebook indices, and 2) The     *
 * VQ outputs class probabilties.  Choose stepWithInputProbabilities:   *
 * or stepWithCodeIndices: as appropriate.                              *
 * Also set dimensions = NULL in initialization if VQ outputs class     *
 * probabilities.                                                       *
 *                                                                      *
 * File: HMM.m                                                          *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/09/92  - Started                                              *
 *  2. 04/24/92  - Went to integer probability inputs                   *
 *  3. 06/02/92  - Eliminated debug                                     *
 *  4. 06/04/92  - Added forcedViterbi                                  *
 *  5. 12/08/92  - Changes going to TIMIT                               *
 *  6. 01/15/92  - Enhancements from Ed Srenger's code.                 *
 *  7. 02/12/93  - Add -addSkipState                                    *
 *  8. 02/18/93  - send weights to initWithLength:... method            *
 *  9. 03/26/93  - Add recognitionStepWithFloatProbs                    *
 * 10. 04/02/93  - Add (float **)probabilityScale                       *
 * 14. 08/06/93  - Add stateCounter, -initModels, scaleType, etc.       *
 ***********************************************************************/

#import "HMM.h"
#import "ViterbiModel.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define DEFAULT_SCALE   0.001       /******* Feature map scale  ***/

@implementation HMM

/*######################################*
 * Initialize instance variables        *
 *######################################*/
- initWithNumberModels:(int)models numberStates:(int)states
         numberCodebooks:(int)numberCodebooks
         codebookDimensions:(int *)dimensions codebookWeights:(float*)weights
         zone:(NXZone *)zone;
{  
  int i, j, k;
  char model_name[10];
  int  *model_states;

  [super init];

/******************************************
 * Initialize Models:                     *
 ******************************************/
  modelLength   = states;
  codebooks     = numberCodebooks;
  model_states  = (int *)malloc(states*sizeof(int));
  myZone        = zone;
  numberModels  = models;
  wordModels    = (id *)NXZoneMalloc(myZone,numberModels*sizeof(id));
  k = 0;
  for(i=0; i<numberModels; i++ )
  {
    for(j=0; j<states; j++)
      model_states[j] = k++;
    sprintf(model_name,"model_%d",i);
    wordModels[i] = [[ViterbiModel allocFromZone:myZone] 
       initWithLength:states name:model_name states:model_states 
       codebooks:codebooks weights:weights zone:myZone];
    if(dimensions != NULL)
       [wordModels[i] initCodebookDimensions:dimensions withWeights:weights];
  }
  free(model_states);

  [self allocateArrays];
  [self initTrainingCounts];
// 
// Other initialization:
//
  bestModel     = 0;
  forcedModel   = 0;
  bestScore     = 0.;
  forcedViterbi = NO;
 
  return self;
}

/*######################################*
 * Allocate arrays:                     *
 *######################################*/
- allocateArrays
{
  int i, j;

  numberInputProbabilities  = numberModels*modelLength;
  stateCounter          = (int *)NXZoneMalloc(myZone, numberInputProbabilities*sizeof(int));
  stateProbabilities    = (float *)NXZoneMalloc(myZone, numberInputProbabilities*sizeof(float));
  probabilityScale      = (float **)NXZoneMalloc(myZone,codebooks*sizeof(float ));
  logProbabilities      = (double **)NXZoneMalloc(myZone,codebooks*sizeof(double**));
  for(i=0; i<codebooks; i++)
  {
    logProbabilities[i] = (double *)NXZoneMalloc(myZone,
                                       numberInputProbabilities*sizeof(double));
    probabilityScale[i] = (float *)NXZoneMalloc(myZone,
                                       numberInputProbabilities*sizeof(double));
    for(j=0; j<numberInputProbabilities; j++)
      probabilityScale[i][j]     = DEFAULT_SCALE;
  }
  return self;
}

/*######################################*
 * Free up arrays and models            *
 *######################################*/
- free
{
  int i;
  [self freeArrays];
  for(i=0; i<numberModels; i++)
    [wordModels[i] free];
  NXZoneFree(myZone, wordModels);

  return [super free];
}

/*######################################*
 * Free up allocated arrays:            *
 *######################################*/
- freeArrays
{
  int i;

  NXZoneFree(myZone, stateCounter);
  NXZoneFree(myZone, stateProbabilities);
  for(i=0; i<codebooks; i++)
  {
    NXZoneFree(myZone, probabilityScale[i]);
    NXZoneFree(myZone, logProbabilities[i]);
  }
  NXZoneFree(myZone, probabilityScale);
  NXZoneFree(myZone, logProbabilities);
  return self;
}
/*#######################################*
 * Initialize for a new recognition     *
 * search:                              *
 *#######################################*/
- initNewSearch
{
  int i;
  
  frameCount     = 0;       /* number of frames so far in word */
  if(forcedViterbi)
    [wordModels[forcedModel] initNewSearch];
  else
    for( i=0; i<numberModels; i++ ) 
      [wordModels[i] initNewSearch];
  return self;
}

/*#######################################*
 * Reset the models to start training:   *
 *#######################################*/
- initTrainingCounts
{
  int i;

  frameCount = 0;
  for(i=0; i<numberModels; i++)
    [wordModels[i] initTrainingCounts];
  for(i=0; i<numberInputProbabilities; i++)
    stateCounter[i] = 0;
  return self;
}

/*######################################*
 * Reset the model's counters           *
 * to start training:                   *
 *######################################*/
- initTrainingCountsFor:(int)modelNumber
{
  frameCount = 0;
  [wordModels[modelNumber] initTrainingCounts];
  return self;
}

/*######################################*
 * Tell the model we got a new training *
 * file.                                *
 *######################################*/
- initTrainingFileFor:(int)modelNumber
{
  [wordModels[modelNumber] initTrainingFile];
  return self;
}

/*######################################*
 * Update with new input for recognition*
 *######################################*/
- recognitionStepWithProbs: (int **)inputProbabilities
{
  int    i, j;
  double probability;

  for(j=0; j<codebooks; j++)
  {
    for(i=0; i<numberInputProbabilities; i++)
    {
     if((probability =(double)inputProbabilities[j][i]*probabilityScale[j][i])
                     > MIN_PROBABILITY )
        logProbabilities[j][i] = log( probability );
     else
        logProbabilities[j][i] = LOG_MIN;
    }
  }

  if(forcedViterbi)
    [wordModels[forcedModel] recognitionStepWithProbs:logProbabilities
                                                      count:frameCount];
  else
   for( i=0; i<numberModels; i++ )             /* cycle viterbi */
    [wordModels[i] recognitionStepWithProbs:logProbabilities count:frameCount];
  frameCount++;
    
  return self;   
}

/*######################################*
 * Update with new input for recognition*
 *######################################*/
- recognitionStepWithFloatProbs: (float **)inputProbabilities
{
  int    i, j;
  double probability;

  for(j=0; j<codebooks; j++)
  {
    for(i=0; i<numberInputProbabilities; i++)
    {
     if((probability =(double)inputProbabilities[j][i]*probabilityScale[j][i])
                     > MIN_PROBABILITY )
        logProbabilities[j][i] = log( probability );
     else
        logProbabilities[j][i] = LOG_MIN;
    }
  }

  if(forcedViterbi)
    [wordModels[forcedModel] recognitionStepWithProbs:logProbabilities
                                                      count:frameCount];
  else
   for( i=0; i<numberModels; i++ )             /* cycle viterbi */
    [wordModels[i] recognitionStepWithProbs:logProbabilities count:frameCount];
  frameCount++;
    
  return self;   
}

/*######################################*
 * Update with new input for recog-     *
 * nition/training.                     *
 *######################################*/
- recognitionStepWithIndices: (short *)indices
{
  int    i;

  if(forcedViterbi)
    [wordModels[forcedModel] recognitionStepWithIndices:indices
                                                      count:frameCount];
  else
   for( i=0; i<numberModels; i++ )             /* cycle viterbi */
    [wordModels[i] recognitionStepWithIndices:indices count:frameCount];
  frameCount++;
    
  return self;   
}

/*######################################*
 * Update with new input for initiali-  *
 * zation training                      *
 *######################################*/
- initialTrainingStepFor:(int)modelNumber state:(int)state
            codebookIndices: (short *)indices
{
  stateCounter[modelNumber*modelLength + state]++;
  [wordModels[modelNumber] initialTrainingStep:state indices:indices];
  return self;
}

/*######################################*
 * Update with new input for initiali-  *
 * zation training                      *
 *######################################*/
- initialTrainingStepFor:(int)modelNumber state:(int)state
              finalFrame:(BOOL)flag
{
  stateCounter[modelNumber*modelLength + state]++;
  [wordModels[modelNumber] initialTrainingStep:state finalFrame:flag];
  return self;
}

/*######################################*
 * Set new bestModel, for use in forced *
 * mode of operation.                   *
 *######################################*/
- setForcedModel: (int)model
{
  forcedModel = MAX(model,0);
  forcedModel = MIN(forcedModel, numberModels-1);
  return self;
}

/*######################################*
 * Set forced flag:                     *
 *######################################*/
- setForced:(BOOL)flag
{
  forcedViterbi = flag;
  return self;
}

/*######################################*
 * Load in an external probability      *
 * scale factor.                        *
 *######################################*/
- setProbabilityScale:(float *)scale codebook:(int)number
{
  int i;
  for(i=0; i<numberInputProbabilities; i++)
    probabilityScale[number][i] = scale[i];
  return self;
}

/*######################################*
 * Set a new number of states per model *
 *######################################*/
- setNumberModelStates:(int)states
{
  int i, j, k;
  int  *model_states;

  modelLength   = states;
  model_states  = (int *)malloc(modelLength*sizeof(int));
  k = 0;
  for(i=0; i<numberModels; i++ )
  {
    for(j=0; j<states; j++)
      model_states[j] = k++;
    [wordModels[i] setModelLength:modelLength states:model_states];
  }
  free(model_states);
  [self freeArrays];
  [self allocateArrays];

  return self;
}

/*#######################################*
 * Return the current best model from    *
 *#######################################*/
- (int) getBestModel
{
  int i;
  double model_score;

/********************************
 * Compare the scores of all    *
 * the models and select the    *
 * model with best score.  If   *
 * forced, then use bestModel   *
 * which should have been pre-  *
 * viously set by the sender.   *
 ********************************/
  if(forcedViterbi)
  {
    bestModel = forcedModel;
    bestScore = [wordModels[bestModel] getModelScore:(frameCount-1)];
  }
  else
  {
    bestModel = 0;
    bestScore = [wordModels[0] getModelScore:(frameCount-1)];
    for( i=1; i<numberModels; i++ )
    {
      model_score = [wordModels[i] getModelScore:(frameCount-1)];
      if(model_score > bestScore)
      {
        bestScore = model_score;
        bestModel  = i;
      }
    }
  }

  return bestModel;
}

/*#######################################*
 * Return the backtrace of the best model*
 * In an array with array[0] = last,     *
 * array[length] = -1.                   *
 *#######################################*/
- (int *) getBestBackTrace
{
  return [wordModels[bestModel] getBackTrace:frameCount];
}

/*######################################*
 * return a pointer to a given model:   *
 *######################################*/
- getModel: (int)model
{
  if( (model<numberModels) && (model>=0) )
    return wordModels[model];
  else
    return nil;
}

/*######################################*
 * Return instance variables:           *
 *######################################*/
- (double) getBestScore     {return bestScore;}
- (int)getFrameCount        {return frameCount;}
- (BOOL)isForced            {return forcedViterbi;}
- (int)getNumberCodebooks   {return codebooks;}
- (int)getNumberStates      {return modelLength;}
- (int)getNumberModels      {return numberModels;}
- (int)numberInputProbabilities {return numberInputProbabilities;}
- (float *)probabilityScaleForCodebook:(int)number {return probabilityScale[number];}

/*#####################################*
 * Calculate the state probabilities   *
 *#####################################*/
- (float *)stateProbabilities
{
  int i, total_counts;
  float  temp;

  total_counts = 0;
  for(i=0; i<numberInputProbabilities; i++)
    total_counts += stateCounter[i];
  temp = total_counts;
  if(total_counts > 0)
  {
   for(i=0; i<numberInputProbabilities; i++)
     if(stateCounter[i] > 0)
       stateProbabilities[i] = stateCounter[i]/temp;
  }

  return stateProbabilities;
}

/*#######################################*
 * Update the model probabilities        *
 * after training:                       *
 *#######################################*/
- updateInitialModels
{
  int i;
  
  for( i=0; i<numberModels; i++ )   
    [wordModels[i] updateInitialModel];
  return self;
}

/*#######################################*
 * Update the model probabilities        *
 * after training:                       *
 *#######################################*/
- updateInitialModel:(int)modelNumber
{
  [wordModels[modelNumber] updateInitialModel];
  return self;
}

/*######################################*
 * Add skip counts to transitionsCount  *
 * for initial training:        *
 *######################################*/
- addSkipStateFor:(int)modelNumber
{
  [wordModels[modelNumber] addSkipState];
  return self;
}

/*######################################*
 *Save out a model's  probabilities *
 *######################################*/
- (BOOL) saveProbabilitiesToFile:(char *)fileName model:(int)modelNumber
{
  BOOL return_code;
  FILE *fp;
  if( (fp = fopen(fileName, "w")) == NULL)
    return NO;
  return_code = [wordModels[modelNumber] saveProbabilitiesToFile:fp];
  fclose(fp);
  return return_code;
}

/*######################################*
 * Read in a model's  probabilities *
 *######################################*/
- (BOOL) readProbabilitiesFromFile:(char *)fileName model:(int)modelNumber
{
  BOOL return_code;
  FILE *fp;
  if( (fp = fopen(fileName, "r")) == NULL)
    return NO;
  return_code = [wordModels[modelNumber] readProbabilitiesFromFile:fp];
  fclose(fp);
  return return_code;
}

/*######################################*
 *Save out a model's  probabilities *
 *######################################*/
- (BOOL) saveProbabilityScaleToFile:(char *)fileName
{
  int i, j, k, n;
  FILE *fp;

  if( (fp = fopen(fileName, "w")) == NULL)
    return NO;
  for(n=0; n<codebooks; n++)
  {
    k = 0;
    for(i=0; i<numberModels; i++)
    {
      for(j=0; j<modelLength; j++)
        fprintf(fp,"%7.5g ", probabilityScale[n][k++]);
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  return YES;
}

/*######################################*
 * Read in a model's  probabilities *
 *######################################*/
- (BOOL) readProbabilityScaleFromFile:(char *)fileName
{
  int i, n;
  FILE *fp;
  if( (fp = fopen(fileName, "r")) == NULL)
    return NO;
  for(n=0; n<codebooks; n++)
    for(i=0; i<numberInputProbabilities; i++)
      fscanf(fp,"%f", &probabilityScale[n][i]);
  fclose(fp);
  return YES;
}

/*######################################*
 * Print model probabilities        *
 *######################################*/
- (BOOL)printAllModelsTo:(char *)fileName
{
  int i;
  FILE *fp;
  if( (fp = fopen(fileName, "w")) == NULL)
    return NO;
  
  for( i=0; i<numberModels; i++ )
  {
    [wordModels[i] printModel:fp];
  }
  fclose(fp);
  return YES;
}

/*######################################*
 * Print model probabilities        *
 *######################################*/
- (BOOL)printModel:(int)modelNumber to:(char *)fileName
{
  FILE *fp;
  if( (fp = fopen(fileName, "w")) == NULL)
    return NO;
  [wordModels[modelNumber] printModel:fp];
  fclose(fp);

  return YES;
}


@end