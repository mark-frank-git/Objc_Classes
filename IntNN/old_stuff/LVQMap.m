/************************************************************************
 * This subclass of IntFeatureMap implements the LVQ2 supervised learn- *
 * ing algorithm.  Normally, the feature map should be initially    *
 * trained using self organization, and then tuned using LVQ.       *
 *                                  *
 * File:LVQMap.h                            *
 *                                  *
 * Revision history:                            *
 *  1. 05/21/92  - Started                      *
 *  2. 06/25/92  - If feedforward don't update output weights in LVQ.   *
 *  3. 06/30/92  - Add LVQ_MOD                                          *
 *  4. 08/03/92  - Add MOD_PLASTICITY                   *
 *  5. 11/15/92  - Add LBG algorithm                    *
 *  6. 11/18/92  - Add getDistortion                    *
 *  7. 12/03/92  - malloc() -> NXZoneMalloc()                           *
 *               - Remove inputNeurons                                  *
 *                 numberInputNeurons -> numberWeightsPerNeuron         *
 *                 numberInputs -> numberInputWeights                   *
 *                 numberOutputs -> numberOutputWeights                 *
 ************************************************************************/
 
#import "LVQMap.h"
#import "IntNeuron.h"
#import <stdio.h>
#import <stdlib.h>
#import <math.h>

#define DEFAULT_CONVERGENCE 0.001
#define DEFAULT_EPSILON     0.35
#define DELTA .1   // default delta value for update..RangeErr

#define MIN(a,b)    ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)    ( ((a)>(b)) ? (a) : (b) )
#define SGN(a)      ( ((a)>=0.) ? 1 : -1 )
#define SQ(a)       ( (a)*(a) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )


@implementation LVQMap

- init
{
   [super init];
   
   learningType = SELF_ORGANIZING;
   epsilon  = DEFAULT_EPSILON;      /* LVQ window size */
   convergenceEpsilon   = DEFAULT_CONVERGENCE;
   converged        = NO;
   lbgDistortion    = 0.;

   return self;
}

/*##############################*
 * Initialize Connections:  *
 *##############################*/
- initConnections
{
  int i;

  [super initConnections];
  
//
// malloc the training arrays:
//
  converged = NO;
  oldDistortion = 1.e30;

  meanCount = (int *)NXZoneMalloc(myZone,numberMapNeurons*sizeof(int));
  means     = (double **)NXZoneMalloc(myZone,
                                        numberMapNeurons*sizeof(double *));
  classCount    = (int **)NXZoneMalloc(myZone,numberMapNeurons*sizeof(int *));
  localPlasticities = (double *)NXZoneMalloc(myZone,
                                             numberMapNeurons*sizeof(double));
  for(i=0; i<numberMapNeurons; i++)
  {
    means[i] = (double *)NXZoneMalloc(myZone,
                                      numberInputWeights*sizeof(double));
    classCount[i] = (int *)NXZoneMalloc(myZone,
                                      numberOutputWeights*sizeof(int));
    localPlasticities[i] = 1.;
  }

  [self initTrainingCounts];
  return self;
}

/*##############################*
 * Initialize Connections:  *
 *##############################*/
- initTrainingCounts
{
  int i, j;
  
  lbgDistortion = 0.;
  for(i=0; i<numberMapNeurons; i++)
  {
    meanCount[i] = 0;
    for(j=0; j<numberInputWeights; j++)
      means[i][j] = 0.;
    for(j=0; j<numberOutputWeights; j++)
      classCount[i][j] = 0;
  }
  
  return self;
}

/*##############################*
 * Initialize the feature   *
 * map weights.         *
 *##############################*/
- initializeWeights: sender
{
  int i;
  int *new_weight;
  
  if(learningType != LBG)
    [super initializeWeights];
  else
  {
    converged = NO;
    oldDistortion = 1.e30;
    if([sender respondsTo:sel_getUid("trainingData")])
    {
      for(i=0; i<numberMapNeurons; i++)
      {
        new_weight = [sender trainingData];
    [mapNeurons[i] setAllWeightsTo:new_weight];
      }
    }
  }
  return self;
}

/*##############################*
 * Free up Connections:     *
 *##############################*/
- freeConnections
{
  int i;
  [super freeConnections];
  NXZoneFree(myZone, meanCount);
  for(i=0; i<numberMapNeurons; i++)
  {
    NXZoneFree(myZone, classCount[i]);
    NXZoneFree(myZone, means[i]);
  }
  NXZoneFree(myZone, classCount);
  NXZoneFree(myZone, means);
  NXZoneFree(myZone, localPlasticities);

  return self;
}


/*##########################*
 * Update the winning neuron    *
 * and the counters:        *
 *##########################*/
- updateWinner
{
  int i, max_class;

  [super updateWinner];
  if(learningType == LBG)
  {
    lbgDistortion  += [mapNeurons[winnerNumber] lastOutput];
    meanCount[winnerNumber]++;
    max_class = [self getInputClass];
    classCount[winnerNumber][max_class]++;
    for(i=0; i<numberInputWeights; i++)
      means[winnerNumber][i] += mapInputs[i];
  }
    
  return self;
}
  

/*##########################*
 * Update the stored weights    *
 * in the feature map:      *
 *##########################*/
- updateWeights
{
  switch(learningType)
  {
    case SELF_ORGANIZING:
    default:
      [super updateWeights];
      break;
    case LVQ_LEARNING:
      [self updateWeightsUsingLVQ];
      break;
    case RANGE_ERR_LEARNING:
      [self updateWeightsUsingRangeErr];
      break;
    case LVQ_MOD_LEARNING:
      [self updateWeightsUsingLVQMod];
      break;
    case LBG:
      [self updateWeightsUsingLBG];
      break;
    case MOD_PLASTICITY:
      [self updateWeightsUsingModulatedPlasticity];
      break;
   }
      
   return self;
}

/*##########################*
 * Update the stored weights    *
 * using LVQ2 learning:     *
 *##########################*/
- updateWeightsUsingLVQ
{
  int   max_class;
  int   di, dj;

//
// Find the class of the input vector by looking at 'class' neurons
//
  max_class = [self getInputClass];

//
// We need to check three conditions before modifying the weights
// of the top two neurons:
//  i.  Winner is not in class of input
//  ii. Runner up is in class of input
//  iii. Input is in window between winner and runner up
//
  if(max_class != [self getClassOfNeuron:winnerNumber])
  {
    [self updateRunnerUp];
    if(max_class == [self getClassOfNeuron:runnerUp])
    {
      di = [mapNeurons[winnerNumber] lastOutput];
      dj = [mapNeurons[runnerUp] lastOutput];
      if(di>0 && dj>0)
      {
        if( MIN( ((double)di/dj), ((double)dj/di) ) >  (1.-epsilon) )
    {
      if(feedforward)
      {
        [self negativeUpdateInputWeightsOfNeuron:winnerNumber];
        [self updateInputWeightsOfNeuron:runnerUp];
      }
      else
      {
        [self negativeUpdateWeightsOfNeuron:winnerNumber];
        [self updateWeightsOfNeuron:runnerUp];
      }
    }
      }
    }
  }
  return self;
}

#define NUMBER_WINNERS  5
/*##########################*
 * Update the stored weights    *
 * using LVQ2 learning, with    *
 * custom modifications.    *
 *##########################*/
- updateWeightsUsingLVQMod
{
  int   i, n, max_class;
  int   *top_neurons;
  int   di, dj;
  
  top_neurons = [self getNearestNeurons];
  
//
// Find the class of the input vector by looking at 'class' neurons
//
  max_class = [self getInputClass];

//
// We need to check three conditions before modifying the weights
// of the top two neurons:
//    i.  Winner is not in class of input
//    ii. Runner up is in class of input
//    iii. Input is in window between winner and runner up
// Also for modified: check these two:
//    i. Winner and runner up in class of input
//
  if(max_class != [self getClassOfNeuron:top_neurons[0]])
  {
    if(max_class == [self getClassOfNeuron:top_neurons[1]])
    {
      di = [mapNeurons[top_neurons[0]] lastOutput];
      dj = [mapNeurons[top_neurons[1]] lastOutput];
      if(di>0 && dj>0)
      {
        if( MIN( ((double)di/dj), ((double)dj/di) ) >  (1.-epsilon) )
    {
      if(feedforward)
      {
        [self negativeUpdateInputWeightsOfNeuron:top_neurons[0]];
        [self updateInputWeightsOfNeuron:top_neurons[1]];
      }
      else
      {
        [self negativeUpdateWeightsOfNeuron:top_neurons[0]];
        [self updateWeightsOfNeuron:top_neurons[1]];
      }
    }
      }
    }
  }
  else
  {
    for(i=1; i<NUMBER_WINNERS; i++)
    {
      if([self getClassOfNeuron:top_neurons[i]] != max_class)
        break;
    }
    if(i==NUMBER_WINNERS)
    {
      for( ; i<numberMapNeurons; i++)
      {
        if([self getClassOfNeuron:top_neurons[i]] != max_class)
          break;
      }
      if(i<numberMapNeurons)
      {
        n = NUMBER_WINNERS-1;
        if(feedforward)
          [self moveInputWeightsOfNeuron:top_neurons[n] toward:top_neurons[i]];
        else
          [self moveWeightsOfNeuron:top_neurons[n] toward:top_neurons[i]];
      }
    }
  }
  return self;
}


/*##############################*
 * This is a method to imple-   *
 * a different wgt tuning   *
 * using range error learning.  *
 *##############################*/
- updateWeightsUsingRangeErr
{
  int       i, j, *old_weights, *delta_weights, *winner_neighbors;
  double    rate=0.,rateI, rateO;
  id        best_neuron, map_neuron;
  
//
// Winning neuron has been calculated 
// earlier with respect to newest data
//
  
  best_neuron = mapNeurons[winnerNumber];

//
// Only update if error in domain > delta, avoids /0
//

  if( [best_neuron lastOutput] > DELTA * synapseScaleFactor )
  { 

//
// Calculate rate = eta*(range error)/(domain error)
//

    for( i=numberInputWeights; i<numberWeightsPerNeuron; i++ )
      rate += SQ([best_neuron getWeight:i]
                  - mapInputs[i]);
    rate /= [best_neuron lastOutput];
    rate = MIN( sqrt(rate),1.);
    rateI = etaInput * rate;
    rateO = etaOutput * rate;
    delta_weights = (int *)malloc(numberWeightsPerNeuron*sizeof(int));

//
// Iterate over neighbors
//
    winner_neighbors    = [self winnerNeighbors];
    i = 0;
    while(YES)
    {
      map_neuron = mapNeurons[winner_neighbors[i]];
      old_weights = [map_neuron getAllWeights];

//
// Iterate over weights
//

      for( j=0; j<numberInputWeights; j++ )
        delta_weights[j] = (int) rateI *(mapInputs[j] - old_weights[j]);
      for( j=numberInputWeights; j<numberWeightsPerNeuron; j++ )
        delta_weights[j] = (int) rateO *(mapInputs[j] - old_weights[j]);
      [map_neuron changeAllWeightsBy: delta_weights];
      if(winner_neighbors[++i] == -1)
        break;
    }
    free(delta_weights);
  }
  return self;
}
  
/*##########################*
 * Update the input weights of  *
 * the neurons in the map   *
 * these will just be the means *
 * of the training inputs which *
 * were in each neurons Voronoi *
 * cell.            *
 *##########################*/
- updateWeightsUsingLBG
{
  int   i, j;
  int   *old_weights;
  static int old_number = 0;
  static int *new_weights;
  long  training_count;
  double temp;
  
  if(old_number!=numberWeightsPerNeuron)
  {
    if(new_weights!=NULL)
      NXZoneFree(myZone, new_weights);
    new_weights = (int *)NXZoneMalloc(myZone, numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  
  if(!converged)        /* Train input weights */
  {
    training_count = 0;
    for(i=0; i<numberMapNeurons; i++)
    {
      training_count += meanCount[i];
      old_weights = [mapNeurons[i] getAllWeights];
      if(meanCount[i]>0)
        temp    = 1./meanCount[i];
      else
        temp      = 10000.;
      for(j=0; j<numberInputWeights; j++)
        new_weights[j] = (int)(means[i][j]*temp);
      for(j=numberInputWeights; j<(numberWeightsPerNeuron); j++)
        new_weights[j] = old_weights[j];
      [mapNeurons[i] setAllWeightsTo:new_weights];
    }
    if(lbgDistortion > 0.)
    {
      lbgDistortion /= (double)training_count;
      if( (oldDistortion - lbgDistortion)/lbgDistortion < convergenceEpsilon)
        converged = YES;
      else
        converged = NO;
      oldDistortion = lbgDistortion;
    }
    if(converged)
      printf("LBG converged\n");
  }
  else              /* Train output weights */
  {
    for(i=0; i<numberMapNeurons; i++)
    {
      old_weights = [mapNeurons[i] getAllWeights];
      if(meanCount[i]>0)
        temp    = 1./meanCount[i];
      else
        temp      = 10000.;
      for(j=0; j<numberInputWeights; j++)
        new_weights[j] = old_weights[j];
      for(j=numberInputWeights; j<(numberWeightsPerNeuron); j++)
         new_weights[j] = synapseScaleFactor*
                      ROUND(classCount[i][j-numberInputWeights]*temp);
      [mapNeurons[i] setAllWeightsTo:new_weights];
    }
  }
  [self initTrainingCounts];
  
  return self;
}

#define PLASTICITY_SCALE    0.1
/*##########################*
 * Update the stored weights    *
 * using SOM with local     *
 * plasticities modulated by    *
 * class errors:        *
 *##########################*/
- updateWeightsUsingModulatedPlasticity
{
  int   max_class;

//
// Find the class of the input vector by looking at 'class' neurons
//
  max_class = [self getInputClass];

//
// If the winner's class == input class, then decrease the plasticity
// otherwise, increase the plasticity.
//
  if(max_class != [self getClassOfNeuron:winnerNumber])
    localPlasticities[winnerNumber] += PLASTICITY_SCALE*
                                      localPlasticities[winnerNumber];
  else
    localPlasticities[winnerNumber] -= PLASTICITY_SCALE*
                                      localPlasticities[winnerNumber];
  [self updateWeightsOfNeuron:winnerNumber 
                      withEtaScale:localPlasticities[winnerNumber] ];

  return self;
}
      

/*##########################*
 * Update the weights of    *
 * a single neuron in the   *
 * opposite direction of    *
 * input:           *
 *##########################*/
- negativeUpdateWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  id     map_neuron;
  double delta;
  
  map_neuron    = mapNeurons[neuronNumber];
  if(old_number!=numberWeightsPerNeuron)
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone, delta_weights);
    delta_weights = (int *)NXZoneMalloc(myZone,numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [map_neuron getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
   delta    = etaInput*(old_weights[i] - mapInputs[i]);
   if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
  {
   delta    = etaOutput*(old_weights[i] - mapInputs[i]);
   if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  [map_neuron changeAllWeightsBy: delta_weights];
  free(delta_weights);

  return self;
}

/*##############################*
 * Update the input weights *
 * of a single neuron:      *
 *##############################*/
- negativeUpdateInputWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  id     map_neuron;
  double delta;
  
  map_neuron    = mapNeurons[neuronNumber];
  if(old_number!=numberWeightsPerNeuron)
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone, delta_weights);
    delta_weights = (int *)NXZoneMalloc(myZone,numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [map_neuron getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = etaInput*(old_weights[i] - mapInputs[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
    delta_weights[i] = 0;
  [map_neuron changeAllWeightsBy: delta_weights];
  
  return self;
}

/*##############################*
 * Update the weights of    *
 * a single neuron      *
 *##############################*/
- updateWeightsOfNeuron: (int) neuronNumber withEtaScale: (double)etaScale
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  double delta, eta_input, eta_output;
  
  eta_input  = etaScale*etaInput;
  eta_output = etaScale*etaOutput;
  if(old_number!=numberWeightsPerNeuron)
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone, delta_weights);
    delta_weights = (int *)NXZoneMalloc(myZone,numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = eta_input*(mapInputs[i] - old_weights[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
  {
    delta        = eta_output*(mapInputs[i] - old_weights[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  [mapNeurons[neuronNumber] changeAllWeightsBy: delta_weights];

  return self;
}



/*##############################*
 * Update the weights of    *
 * a single neuron to be more   *
 * like another neuron.     *
 *##############################*/
- moveWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  int    *goal_weights;
  id     map_neuron;
  double delta;
  
  map_neuron    = mapNeurons[neuronNumber];
  if(old_number!=numberWeightsPerNeuron)
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone, delta_weights);
    delta_weights = (int *)NXZoneMalloc(myZone,numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [map_neuron getAllWeights];
  goal_weights = [mapNeurons[toNeuron] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = etaInput*(goal_weights[i] - old_weights[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
  {
    delta        = etaOutput*(goal_weights[i] - old_weights[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  [map_neuron changeAllWeightsBy: delta_weights];

  return self;
}


/*##############################*
 * Update the weights of    *
 * a single neuron to be more   *
 * like another neuron.     *
 *##############################*/
- moveInputWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  int    *goal_weights;
  id     map_neuron;
  double delta;
  
  map_neuron    = mapNeurons[neuronNumber];
  if(old_number!=numberWeightsPerNeuron)
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone, delta_weights);
    delta_weights = (int *)NXZoneMalloc(myZone,numberWeightsPerNeuron*sizeof(int));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [map_neuron getAllWeights];
  goal_weights = [mapNeurons[toNeuron] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = etaInput*(goal_weights[i] - old_weights[i]);
    if( (delta_weights[i] = (int)delta) == 0 )
      delta_weights[i] = SGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
    delta_weights[i] = 0;

  [map_neuron changeAllWeightsBy: delta_weights];

  return self;
}

/*########################*
 * Set convergence epsilon:*
 *########################*/
- setConvergenceEpsilon: (int)newEpsilon
{
  convergenceEpsilon    = newEpsilon;
  return self;
}

/*########################*
 * Set converged:   *
 *########################*/
- setConverged: (BOOL)flag
{
  converged = flag;
  return self;
}

/*########################*
 * Set epsilon:         *
 *########################*/
- setEpsilon: (double)newEpsilon
{
  epsilon   = newEpsilon;
  return self;
}

/*########################*
 * Set learning type        *
 *########################*/
- setLearningType: (int)type
{
  learningType = type;
  return self;
}

/*########################*
 * K means converged?:  *
 *########################*/
- (BOOL) isConverged
{
  return converged;
}

/*########################*
 * Get epsilon:         *
 *########################*/
- (double) getEpsilon
{
  return epsilon;
}

/*########################*
 * Get Class of neuron by   *
 * finding output synapse   *
 * having largest value:    *
 *########################*/
- (int) getClassOfNeuron:(int)neuronNumber
{
  int i, start, max_class, max_value, new_value;
  id  neuron;
  
  neuron    = mapNeurons[neuronNumber];
  start     = numberInputWeights;
  max_class = 0;
  max_value = [neuron getWeight:start++];
  for(i=1; i<numberOutputWeights; i++)
  {
    new_value = [neuron getWeight:start++];
    if(new_value > max_value)
    {
      max_value = new_value;
      max_class = i;
    }
  }

  return max_class;
}

/*########################*
 * Get Learning Type:   *
 *########################*/
- (int) getLearningType;
{
  return learningType;
}

/*########################*
 * Get convergence epsilon:*
 *########################*/
- (int) getConvergenceEpsilon
{
  return convergenceEpsilon;
}
/*########################*
 * Return the average     *
 * distortion from LBG:   *
 *########################*/
- (double)getDistortion
{
  return oldDistortion;
}

/*########################*
 * Dummy:       *
 *########################*/
- (int *)trainingData
{
  return NULL;
}

@end