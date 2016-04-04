/************************************************************************
 * This subclass of object implements a feed forward self-organizing    *
 * feature map.                                                         *
 *                                                                      *
 * File:IntFeatureMap.m                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/27/92  - Started                                              *
 *  2. 04/24/92  - Modified FeatureMapEngine -> IntFeatureMap           *
 *  3. 04/28/92  - Added updateNeuronInput/OutputWeights:               *
 *  4. 05/04/92  - Use rounding on weight updates                       *
 *  5. 05/07/92  - Add separate eta for output weights                  *
 *  6. 06/01/92  - Free output arrays on next call                      *
 *               - Went to arrays for mapNeurons and inputNeurons       *
 *  7. 11/24/92  - Combine updateWinner into applyInput.                *
 *               - Add readWeightsFile, saveWeightsFile.                *
 *  8. 12/03/92  - malloc() -> NXZoneMalloc()                           *
 *               - Remove inputNeurons                                  *
 *                 numberInputNeurons -> numberWeightsPerNeuron         *
 *                 numberInputs -> numberInputWeights                   *
 *                 numberOutputs -> numberOutputWeights                 *
 *  9.  03/01/93  - Add normalizeOutputWeightsOfNeuron, *outputWeights, *
 *                  *deltaWeights.                                      *
 * 10.  03/12/93  - Add weightUpdateType                                *
 * 11.  04/02/93  - Add probabilityScale                                *
 * 12.  04/08/93  - Fixed getNearestNeurons for DOT_PRODUCT, use        *
 *                  runnerUp for smoothing output with numberSmooth=2   *
 * 13.  05/17/93  - Use rounding in statistical learning                *
 * 14.  09/22/93  - Changes to make FltFeatureMap a subclass of this    *
 *                - Add -setWinner                                      *
 ************************************************************************/
 
#import "IntFeatureMap.h"
#import "Random.h"
#import "IntNeuron.h"
#import "c_headers.h"
#import <stdlib.h>
#import <stdio.h>
#import <string.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

int     *neuron_output, *nearest_index=NULL;    /* list of closest neurons   */ 
void bubble_sort(int n, int (*comp)(), void (*swap)());
void index_swap(int i, int j);
int ascend_comp(int i, int j);
int descend_comp(int i, int j);
/************************************************************************
 *                                                                      *
 *  int ascend_comp(int i, int j)                                       *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    i = index of 1st element                                          *
 *    j = index of 2nd element                                          *
 *                                                                      *
 ************************************************************************/
int ascend_comp(int i, int j)
{
  return(
  neuron_output[nearest_index[i]]-neuron_output[nearest_index[j]] < 0 ? -1 : 
  neuron_output[nearest_index[i]]-neuron_output[nearest_index[j]] ? 1 : 0);
}
/************************************************************************
 *                                                                      *
 *  int descend_comp(int i, int j)                                      *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    i = index of 1st element                                          *
 *    j = index of 2nd element                                          *
 *                                                                      *
 ************************************************************************/
int descend_comp(int i, int j)
{
  return(
  neuron_output[nearest_index[i]]-neuron_output[nearest_index[j]] < 0 ? 1 : 
  neuron_output[nearest_index[i]]-neuron_output[nearest_index[j]] ? -1 : 0);
}
/************************************************************************
 *                                                                      *
 *  int index_swap(i, j)                                                *
 *                                                                      *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    i = index of 1st element                                          *
 *    j = index of 2nd element                                          *
 *                                                                      *
 ************************************************************************/
void index_swap(int i, int j)
{
   int temp;
   temp             = nearest_index[i];
   nearest_index[i] = nearest_index[j];
   nearest_index[j] = temp;
   
}

   
@implementation IntFeatureMap

- init
{
   [super init];
   
   mapNeurons   = NULL;
   mapInputs    = NULL;
   myZone       = NULL;
   deltaWeights = outputWeights = NULL;
   
   neighborhoodSize = DEFAULT_SIZE;
   etaInput = etaOutput = DEFAULT_ETA;
   neuronType       = EUCLIDEAN_DISTANCE;
   feedforward      = YES;
   learningEnabled  = YES;
   smoothOutput     = NO;
   weightUpdateType = TRUNCATED;
   numberSmooth     = DEFAULT_SMOOTH;
   winnerNumber     = runnerUp = 0;
   synapseScaleFactor   = DEFAULT_SCALE;
   
   return self;
}

/*######################*
 * Initialization:  *
 *######################*/
- initWithInputs:(int)totalInputs outputs:(int)outputs
            rows:(int)numberRows columns:(int)numberCols
{
   
   [self init];
   
   myZone         = NXDefaultMallocZone();
   randomObject       = [[Random allocFromZone:myZone] initSeeds: SEED1: 
                                                           SEED2: SEED3];
   rows               = numberRows;
   columns            = numberCols;
   numberOutputWeights= outputs;
   numberInputWeights = totalInputs - outputs;
   numberWeightsPerNeuron = totalInputs;
   numberMapNeurons       = rows*columns;
   
   [self initConnections];
   return self;
}

/*######################*
 * Initialization:  *
 *######################*/
- initWithInputs:(int)totalInputs outputs:(int)outputs
            rows:(int)numberRows columns:(int)numberCols zone:(NXZone*)zone
{
   
   [self init];
   
   myZone             = zone;
   randomObject       = [[Random allocFromZone:myZone] initSeeds: SEED1: 
                                                           SEED2: SEED3];
   rows               = numberRows;
   columns            = numberCols;
   numberOutputWeights= outputs;
   numberInputWeights = totalInputs - outputs;
   numberWeightsPerNeuron = totalInputs;
   numberMapNeurons       = rows*columns;
   
   [self initConnections];
   return self;
}

/*##############################*
 * Initialize Connections:  *
 *##############################*/
- initConnections
{
  int i;

//
// Allocate work arrays
//
  deltaWeights = (int *)NXZoneMalloc(myZone,
                        numberWeightsPerNeuron*sizeof(int));
  outputWeights = (int *)NXZoneMalloc(myZone,
                        numberWeightsPerNeuron*sizeof(int));
  probabilityScale = (float *)NXZoneMalloc(myZone,
                        numberOutputWeights*sizeof(float));
//
// create the nodes
//
  if(mapNeurons!=NULL)
    NXZoneFree(myZone,mapNeurons);
  mapNeurons  = (id *)NXZoneMalloc(myZone,
                              (numberMapNeurons*sizeof(IntNeuron)));
  for(i=0; i<numberMapNeurons; i++)
      mapNeurons[i] = [[IntNeuron allocFromZone:myZone]
         initWithNumberInputs:numberWeightsPerNeuron zone:myZone random:randomObject];
      
  return self;
}

/*##############################*
 * Free up Connections:     *
 *##############################*/
- freeConnections
{
  int i;

  NXZoneFree(myZone,deltaWeights);
  NXZoneFree(myZone,outputWeights);
  NXZoneFree(myZone,probabilityScale);

  for(i=0; i<numberMapNeurons; i++)
  {
      [mapNeurons[i] freeStorage];
      [mapNeurons[i] free];
  }
  return self;
}

/*##############################*
 * Initialize the feature   *
 * map weights.         *
 *##############################*/
- initializeWeights
{
  int j;
  
 [randomObject setSeeds: SEED1: SEED2: SEED3];
  for(j=0; j<numberMapNeurons; j++)
    [mapNeurons[j] initializeWeights];
  return self;
}

/*##############################*
 * Update input neurons with new*
 * input, and feed to map. Also,*
 * find the new winning neuron  *
 *##############################*/
- applyInput:(void *)input
{
  int   i;
  int    max_activity, min_activity, activity;
   
  mapInputs = (int *)input;
  winnerNumber = 0;
  if(feedforward)       /* use only input weights */
  {
    max_activity = min_activity = 
       [mapNeurons[0] stepStart:0 length:numberInputWeights input:mapInputs];
    if(neuronType == EUCLIDEAN_DISTANCE)
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] stepStart:0 length:numberInputWeights 
                                      input:mapInputs];
        if(activity < min_activity)
        {
          min_activity = activity;
          winnerNumber = i;
        }
      }
    }
    else
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] stepStart:0 length:numberInputWeights
                                      input:mapInputs];
        if(activity > max_activity)
        {
          max_activity = activity;
          winnerNumber = i;
        }
      }
    }
  }
  else              /* use input and output weights */
  {
    max_activity = min_activity = [mapNeurons[0] step:mapInputs];
    if(neuronType == EUCLIDEAN_DISTANCE)
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] step:mapInputs];
        if(activity < min_activity)
        {
          min_activity = activity;
          winnerNumber = i;
        }
      }
    }
    else
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] step:mapInputs];
        if(activity > max_activity)
        {
          max_activity = activity;
          winnerNumber = i;
        }
      }
    }
  }
   return self;
}

/*##############################*
 * Update input neurons with new*
 * input, but don't feed to map.*
 *##############################*/
- setInput:(void *)input
{   
  mapInputs = (int *)input;
  return self;
}

/*##############################*
 * Update the stored weights    *
 * in the feature map:          *
 *##############################*/
- updateWeights
{
   int   i;
   int  *winner_neighbors;
   
   //
   // First find the winner neighbors:
   //
   winner_neighbors     = [self winnerNeighbors];
   i = 0;
   if(learningEnabled)
   {
     while(YES)
     {
       [self updateWeightsOfNeuron: winner_neighbors[i]];
       if(winner_neighbors[++i] == -1)
         break;
     }
   }
   
   return self;
}

/*##############################*
 * Find the current winning     *
 * (i.e., highest activity      *
 * neuron.  This is now per-    *
 * formed in applyInput, and is *
 * included here for compati-   *
 * bility.                      *
 *##############################*/
- updateWinner { return self;}

/*##############################*
 * Find the runner up neuron    *
 * (i.e., 2nd highest activity  *
 * neuron;                      *
 *##############################*/
- updateRunnerUp
{
  int    i;
  int    max_activity, min_activity, activity;
   
  if(winnerNumber != 0)
  {
    runnerUp = 0;
    max_activity = min_activity = [mapNeurons[0] lastOutput];
    i = 1;
  }
  else
  {
    runnerUp = 1;
    max_activity = min_activity = [mapNeurons[1] lastOutput];
    i = 2;
  }
  if(neuronType == EUCLIDEAN_DISTANCE)
  {
    for(; i<numberMapNeurons; i++)
    {
      activity = [mapNeurons[i] lastOutput];
      if( (activity<min_activity) && (i!=winnerNumber) )
      {
        min_activity = activity;
        runnerUp = i;
      }
    }
  }
  else
  {
    for(; i<numberMapNeurons; i++)
    {
      activity = [mapNeurons[i] lastOutput];
      if( (activity>min_activity) && (i!=winnerNumber) )
      {
        max_activity = activity;
        runnerUp = i;
      }
    }
  }
  return self;
}

/*##############################*
 * Update the weights of        *
 * a single neuron              *
 *##############################*/
- updateWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  int    *old_weights;
  double delta;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = etaInput*(mapInputs[i] - old_weights[i]);
    switch(weightUpdateType)
    {
      case TRUNCATED:
      default:
        deltaWeights[i] = (int)delta;
        break;
      case ROUNDED:
        deltaWeights[i] = ROUND(delta);
        break;
      case MINIMUM_STEP:
        if( (deltaWeights[i] = (int)delta) == 0 )
          deltaWeights[i] = SGN(delta);
        break;
      case STATISTICAL:
        if( (deltaWeights[i] = (int)delta) == 0 )
        {
          if([randomObject percent] < fabs(delta))
            deltaWeights[i] = SGN(delta);
          else
            deltaWeights[i] = 0;
        }
        break;
      case INCREMENT:
        deltaWeights[i] = etaInput*mapInputs[i];
        break;
    }
  }
  for(; i<numberWeightsPerNeuron; i++)
  {
    delta        = etaOutput*(mapInputs[i] - old_weights[i]);
    switch(weightUpdateType)
    {
      case TRUNCATED:
      default:
        deltaWeights[i] = (int)delta;
        break;
      case ROUNDED:
        deltaWeights[i] = ROUND(delta);
        break;
      case MINIMUM_STEP:
        if( (deltaWeights[i] = (int)delta) == 0 )
          deltaWeights[i] = SGN(delta);
        break;
      case STATISTICAL:
        if( (deltaWeights[i] = (int)delta) == 0 )
        {
          if([randomObject percent] < fabs(delta))
            deltaWeights[i] = SGN(delta);
          else
            deltaWeights[i] = 0;
        }
        break;
      case INCREMENT:
        deltaWeights[i] = etaOutput*mapInputs[i];
        break;
    }
  }
  [mapNeurons[neuronNumber] changeAllWeightsBy: deltaWeights];

  return self;
}

/*##############################*
 * Update the input weights     *
 * of a single neuron:          *
 *##############################*/
- updateInputWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  int    *old_weights;
  double delta;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = etaInput*(mapInputs[i] - old_weights[i]);
    switch(weightUpdateType)
    {
      case TRUNCATED:
      default:
        deltaWeights[i] = (int)delta;
        break;
      case ROUNDED:
        deltaWeights[i] = ROUND(delta);
        break;
      case MINIMUM_STEP:
        if( (deltaWeights[i] = (int)delta) == 0 )
          deltaWeights[i] = SGN(delta);
        break;
      case STATISTICAL:
        if( (deltaWeights[i] = (int)delta) == 0 )
        {
          if([randomObject percent] < fabs(delta))
            deltaWeights[i] = SGN(delta);
          else
            deltaWeights[i] = 0;
        }
        break;
      case INCREMENT:
        deltaWeights[i] = etaInput*mapInputs[i];
        break;
    }
  }
  for(; i<numberWeightsPerNeuron; i++)
    deltaWeights[i] = 0;
  [mapNeurons[neuronNumber] changeAllWeightsBy: deltaWeights];
  
  return self;
}

/*##############################*
 * Update the output weights    *
 * of a single neuron:          *
 *##############################*/
- updateOutputWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  int    *old_weights;
  double delta;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
    deltaWeights[i] = 0;
  for(; i<numberWeightsPerNeuron; i++)
  {
    delta = etaOutput*(mapInputs[i] - old_weights[i]);
    switch(weightUpdateType)
    {
      case TRUNCATED:
      default:
        deltaWeights[i] = (int)delta;
        break;
      case ROUNDED:
        deltaWeights[i] = ROUND(delta);
        break;
      case MINIMUM_STEP:
        if( (deltaWeights[i] = (int)delta) == 0 )
          deltaWeights[i] = SGN(delta);
        break;
      case STATISTICAL:
        if( (deltaWeights[i] = (int)delta) == 0 )
        {
          if([randomObject percent] < fabs(delta))
            deltaWeights[i] = SGN(delta);
          else
            deltaWeights[i] = 0;
        }
        break;
      case INCREMENT:
        deltaWeights[i] = etaOutput*mapInputs[i];
        break;
    }
  }
  [mapNeurons[neuronNumber] changeAllWeightsBy: deltaWeights];

  return self;
}

/*##############################*
 * Normalize the output weights *
 * of all the neurons such that *
 * the sum over all the output  *
 * weights for each neuron is   *
 * equal to synapseScaleFactor. *
 * This is used in MAP criterion*
 *##############################*/
- normalizeMAPOutputWeights
{
  int    j;

// Loop over all the map neurons:
  for(j=0; j<numberMapNeurons; j++)
    [self normalizeOutputWeightsOfNeuron:j];

  return self;
}

/*##############################*
 * Normalize the output weights *
 * of the given neuron such that*
 * the sum is equal to synapse- *
 * ScaleFactor:                 *
 *##############################*/
- normalizeOutputWeightsOfNeuron:(int)neuronNumber
{
  int    i, weight_sum;
  int    *old_weights;
  double scale;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  weight_sum = 0;
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
  {
    old_weights[i] = MAX(0, old_weights[i]);
    weight_sum    += old_weights[i];
  }
  if(weight_sum>0)
    scale = synapseScaleFactor/weight_sum;
  else
    scale = 0.;
  for(i=0; i<numberInputWeights; i++)
    deltaWeights[i] = old_weights[i];
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
    deltaWeights[i] = ROUND((scale*old_weights[i]));
  [mapNeurons[neuronNumber] setAllWeightsTo: deltaWeights];

  return self;
}

/*##############################*
 * Normalize the output weights *
 * of all the neurons such that *
 * the sum over all the neurons *
 * for a particular state is    *
 * equal to synapseScaleFactor. *
 * This is used in ML criterion *
 *##############################*/
- normalizeMLOutputWeights
{
  int i, j; 
  int weight_sum, new_weight;
  double scale = 1.;

// Loop over the output (probability) weights:
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
  {
    weight_sum = 0;
    for(j=0; j<numberMapNeurons; j++)
      weight_sum += MAX(0, [mapNeurons[j] getWeight:i]);
    if(weight_sum>0)
      scale = synapseScaleFactor/weight_sum;
    for(j=0; j<numberMapNeurons; j++)
    {
      new_weight = ROUND(scale*MAX(0., [mapNeurons[j] getWeight:i]));
      [mapNeurons[j] setWeight:i to:new_weight];
    }
  }

  return self;
}

/*##############################*
 * Return a list of the         *
 * winner's neighbors:          *
 *##############################*/
- (int *)winnerNeighbors
{
  int  i, j, k, n, index, size;
  int  winner_row, winner_col;
  int  min_row, max_row, min_col, max_col;
  static int  *neighbors = NULL;
  static int  old_size   = 0;
  
  if(columns>0)
  {
    winner_row = winnerNumber/columns;
    winner_col = winnerNumber%columns;
  }
  else
    winner_row = winner_col = 0;
  
  min_row = MAX(winner_row-neighborhoodSize, 0);
  max_row = MIN(winner_row+neighborhoodSize+1, rows);
  min_col = MAX(winner_col-neighborhoodSize, 0);
  max_col = MIN(winner_col+neighborhoodSize+1, columns);
  
  size = neighborhoodSize+neighborhoodSize + 1;
  if(size>old_size)
  {
    old_size = size;
    if(neighbors != NULL)
      NXZoneFree(myZone,neighbors);
    neighbors = (int *)NXZoneMalloc(myZone, ((size*size+1)*sizeof(int)));
  }
  
  k = 0;
  n = min_row*columns + min_col;
  for(i=min_row; i<max_row; i++)
  {
    index = n;
    for(j=min_col; j<max_col; j++)
      neighbors[k++] = index++;
    n += columns;
  }
  neighbors[k] = -1;
  
  return neighbors;
}

/*##############################*
 * Set feedforward flag         *
 *##############################*/
- setFeedforward: (BOOL)flag
{
  feedforward = flag;
  return self;
}

/*##############################*
 * Set new learning             *
 * parameters:                  *
 *##############################*/
- setLearningInput: (double)inputEta output:(double)outputEta
{
  etaInput  = inputEta;
  etaOutput = outputEta;
  learningEnabled = (etaInput>ETA_EPS || etaOutput>ETA_EPS);
  return self;
}

/*##############################*
 * Set new neuron type:         *
 *##############################*/
- setNeuronType: (int)type
{
  int i;
  
  neuronType = type;
// update the map neurons:
  for(i=0; i<numberMapNeurons; i++)
    [mapNeurons[i] setNeuronType:type];

  return self;
}

/*##############################*
 * Set number of neurons for    *
 * smoothing:                   *
 *##############################*/
- setNumberSmooth: (int)number
{
  numberSmooth = MIN(number, MAX_SMOOTH);
  return self;
}

/*##############################*
 * Set new neighborhood         *
 * size:                        *
 *##############################*/
- setRadius: (int) newSize
{
  neighborhoodSize = newSize;
  return self;
}

/*##############################*
 * Set smoothOutput flag        *
 *##############################*/
- setSmoothOutput: (BOOL)flag
{
  smoothOutput = flag;
  return self;
}

/*##############################*
 * Set a new number of inputs   *
 *##############################*/
- setNumberWeightsPerNeuron: (int)number
{
  numberWeightsPerNeuron = number;
  numberInputWeights    = numberWeightsPerNeuron - numberOutputWeights;
  return self;
}

/*##############################*
 * Set a new number of outputs  *
 *##############################*/
- setReverseFlowOutputs: (int)number
{
  numberOutputWeights   = number;
  numberInputWeights    = numberWeightsPerNeuron - numberOutputWeights;
  return self;
}

/*##############################*
 * Set a new number of rows     *
 * and columns in the map.      *
 *##############################*/
- setNumberRows:(int)numberRows columns:(int)numberCols
{
  rows      = numberRows;
  columns   = numberCols;
  numberMapNeurons = rows*columns;
  return self;
}

/*##############################*
 * Set the synapse scale factor *
 *##############################*/
- setScaleFactor:(double)scale
{
  synapseScaleFactor = scale;
  return self;
}

/*##############################*
 * Set the learning type:       *
 *##############################*/
- setWeightUpdateType:(int)type
{
  weightUpdateType = type;
  return self;
}

/*##############################*
 * Set a new winning neuron #   *
 *##############################*/
- setWinner:(int)newWinner
{
  winnerNumber = newWinner;
  return self;
}

/*##############################*
 * Get a list of the            *
 * neuron #'s closest to        *
 * the last input.  I.e.,       *
 * winner = list[0]             *
 *##############################*/
- (int *)getNearestNeurons
{
  int i;
  static int old_number=0;
  
  if( old_number<numberMapNeurons )
  {
    old_number = numberMapNeurons;
    if(nearest_index!=NULL)
    {
      free(nearest_index);
      free(neuron_output);
    }
    neuron_output  = (int *)NXZoneMalloc(myZone,
                                   (numberMapNeurons*sizeof(int)));
    nearest_index  = (int *)NXZoneMalloc(myZone,
                                   (numberMapNeurons*sizeof(int)));
  }
  for(i=0; i<numberMapNeurons; i++)
  {
    neuron_output[i] = [mapNeurons[i] lastOutput];
    nearest_index[i] = i;
  }
  if(neuronType == EUCLIDEAN_DISTANCE)
    bubble_sort(numberMapNeurons, ascend_comp, index_swap);
  else
    bubble_sort(numberMapNeurons, descend_comp, index_swap);
  return nearest_index;
}

/*##############################*
 * Return instance variables:   *
 *##############################*/
- (int) getRadius           {return neighborhoodSize;}
- (int) numberMapNeurons    {return numberMapNeurons;}

/*##############################*
 * Get class of input vector:   *
 *##############################*/
- (int)getInputClass
{
  int i,j;
  int max_class, max_value, new_value;
//
// Find the class of the input vector by looking at 'class' neurons
//
  max_class = 0;
  max_value = mapInputs[numberInputWeights];
  j = 1;
  for(i=numberInputWeights+1; i<numberWeightsPerNeuron; i++)
  {
    new_value = mapInputs[i];
    if(new_value > max_value)
    {
      max_class = j;
      max_value = new_value;
    }
    j++;
  }
  return max_class;
}

/*##############################*
 * Get learning parameters      *
 *##############################*/
- getLearningInput:(double *)inputEta output:(double *)outputEta
{
  *inputEta    = etaInput;
  *outputEta   = etaOutput;
  return self;
}

/*##############################*
 * Get neuron type:             *
 *##############################*/
- (int) getNeuronType;
{
  return neuronType;
}

/*##############################*
 * Get a neuron:                *
 *##############################*/
- getNeuron: (int)neuronNumber
{
  return mapNeurons[neuronNumber];
}

/*##############################*
 * Return the reverse flow      *
 * output by getting the output *
 * weights of the winning       *
 * neuron:                      *
 *##############################*/
- (void *)getReverseFlow
{
  int    i, n, start;
  static int *distances=NULL;
  int   *nearest_neurons, runner_up[2];
  double  sum, temp, smooth_minus_one, total_distance;
  
  if(distances==NULL)
    distances = (int *)NXZoneMalloc(myZone,(MAX_SMOOTH*sizeof(int)));
  if(smoothOutput && numberSmooth>1)
  {
    if(numberSmooth == 2)
    {
      [self updateRunnerUp];
      runner_up[0] = winnerNumber;
      runner_up[1] = runnerUp;
      nearest_neurons = runner_up;
    }
    else
      nearest_neurons   = [self getNearestNeurons];
    total_distance  = 0.;
    for(n=0; n<numberSmooth; n++)
    {
      distances[n]  = [mapNeurons[nearest_neurons[n]] lastOutput];
      total_distance    += distances[n];
    }
    total_distance = MAX(1., total_distance);
    start = numberInputWeights;
    smooth_minus_one = numberSmooth-1;
    for(i=0; i<numberOutputWeights; i++)
    {
      sum = 0;
      for(n=0; n<numberSmooth; n++)
      {
        temp = [mapNeurons[nearest_neurons[n]]
                      getWeight:start];
        sum += temp * (total_distance - smooth_minus_one*distances[n]);
      }
      outputWeights[i] = sum/total_distance;
      start++;
    }
  }
  else
  {
    start = numberInputWeights;
    for(i=0; i<numberOutputWeights; i++)
      outputWeights[i] = [mapNeurons[winnerNumber] getWeight:start++];
  }
    
  return outputWeights;
}

/*##############################*
 * Return the winning           *
 * neuron number:               *
 *##############################*/
- (int)getWinner
{
  return winnerNumber;
}

/*##############################*
 * Return the synaptic          *
 * weights of the winning       *
 * neuron:                      *
 *##############################*/
- (void *)getWinnerWeightsStart:(int) start length:(int)length
{
  int    i;
  id     winner_neuron;
  
  winner_neuron = mapNeurons[winnerNumber];
  for(i=0; i<length; i++)
    outputWeights[i] = [winner_neuron getWeight:start++];
  
  return outputWeights;
}


/*##############################*
 * Return the synaptic          *
 * weights of the given         *
 * neuron:                      *
 *##############################*/
- (void *)getNeuronWeightsStart:(int)start length:(int)length
                        neuron:(int)neuronNumber
{
  int    i;
  id     neuron;
  
  neuron  = mapNeurons[neuronNumber];
  for(i=0; i<length; i++)
    outputWeights[i] = [neuron getWeight:start++];
  
  return outputWeights;
}

/*##############################*
 * Return the feedforward       *
 * flag:                        *
 *##############################*/
- (BOOL)isFeedforward
{
  return feedforward;
}

/*##############################*
 * Return the synapse scale     *
 * factor:                      *
 *##############################*/
- (double)getScaleFactor
{
  return synapseScaleFactor;
}

/*##############################*
 * Calculate the scale factors  *
 * such that the sum over all   *
 * the neurons for a particular *
 * state is equal to 1.         *
 *##############################*/
- (float *)probabilityScale
{
  int i, j, k, weight_sum;

// Loop over the output (probability) weights:
  k = 0;
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
  {
    weight_sum = 0;
    for(j=0; j<numberMapNeurons; j++)
      weight_sum += MAX(0, [mapNeurons[j] getWeight:i]);
    if(weight_sum>0)
      probabilityScale[k] = 1./weight_sum;
    k++;
  }

  return probabilityScale;
}


#define INPUT_NAME "DELAY1"
#define OUTPUT_NAME "CONT1"
/*###############################*
 * Read the stored weights into  *
 * the feature map:              *
 *###############################*/
-(BOOL) readWeightsFile:(char *)fileName
{
  char c;
  int  string_index, new_index, weight_index, neuron_number;
  FILE *weights_file;
  char line[MAX_LINE];
  float temp;
  
  if( (weights_file = fopen(fileName, "r")) == NULL)
    return NO;
 
  weight_index = neuron_number = 0; 
  while(fgets(line, MAX_LINE, weights_file) != NULL)
  {
    if(sindex(line, INPUT_NAME) >= 0)
    {
      string_index = 0;
      while( (new_index = sindex(&line[string_index], INPUT_NAME)) >= 0)
      {
        string_index += new_index;
        c = line[string_index];
        while( (c != '+') && (c != '-') )
          c = line[++string_index];
        sscanf(&line[string_index], "%f", &temp);
        if(weight_index>=numberWeightsPerNeuron)
        {
          printf("Too many input weights in readWeightsFile\n");
          return NO;
        }
        deltaWeights[weight_index++] = ROUND(synapseScaleFactor*temp);
      }
    }
    if(sindex(line, OUTPUT_NAME) >= 0)
    {
      string_index = 0;
      while( (new_index = sindex(&line[string_index], OUTPUT_NAME)) >= 0)
      {
        string_index += new_index;
        c = line[string_index];
        while( (c != '+') && (c != '-') )
          c = line[++string_index];
        sscanf(&line[string_index], "%f", &temp);
        if(weight_index>=numberWeightsPerNeuron)
        {
          printf("Too many output weights in readWeightsFile\n");
          return NO;
        }
        deltaWeights[weight_index++] = ROUND((synapseScaleFactor*temp));
      }
    }
    if(weight_index >= numberWeightsPerNeuron)
    {
      weight_index = 0;
      [mapNeurons[neuron_number++] setAllWeightsTo:deltaWeights];
      if(neuron_number>numberMapNeurons)
        break;
    }
  }
  fclose(weights_file);
  return YES;
}


#define LINE_LENGTH 20
/*###############################*
 * Save the map's weights:       *
 *###############################*/
- (BOOL)saveWeightsFile: (char *)weightFile
{
  int   i, j;
  int   *weights;
  char  weight_string[20];
  float temp;
  FILE  *file_write;

  if( (file_write = fopen(weightFile,"w")) == NULL)
    return NO;

  for(i=0; i<numberMapNeurons; i++)
  {
    fprintf(file_write,"<%3d %d; ",i+1, numberWeightsPerNeuron);
    weights = [mapNeurons[i] getAllWeights];
    for(j=0; j<numberWeightsPerNeuron; j++)
    {
      if(j<numberInputWeights)
      {
        fprintf(file_write,INPUT_NAME);
        fprintf(file_write,":%d",j+1);
      }
      else
      {
        fprintf(file_write,OUTPUT_NAME);
        fprintf(file_write,":%d",j+1-numberInputWeights);
      }
      temp = (double)weights[j]/synapseScaleFactor;
      sprintf(weight_string,"%4.3f", temp);
      if(rindex(weight_string,'-')==NULL)
        fprintf(file_write,"+");
      fprintf(file_write,"%4.3f ",temp);
      if(!(j%LINE_LENGTH))
        fprintf(file_write,"\n");
    }
    fprintf(file_write,">\n");
  }
  fclose(file_write);
  return YES;
}

      
@end