/************************************************************************
 * This subclass of object implements a feed forward self-organizing    *
 * feature map.                                                         *
 *                                                                      *
 * File:FltFeatureMap.m                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/24/93  - Modified from IntFeatureMap                          *
 *  2. 08/10/93  - Add probabilityScale                                 *
 *  3. 08/16/93  - Add normalizeMAPOutputWeights, etc.                  *
 *  4. 09/22/93  - Made subclass of IntFeatureMap.                      *
 ************************************************************************/
 
#import "FltFeatureMap.h"
#import "Random.h"
#import "FltNeuron.h"
#import "c_headers.h"
#import <stdlib.h>
#import <stdio.h>
#import <string.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

float *flt_neuron_output;
int   *closest_index=NULL;  /* list of closest neurons   */ 
void bubble_sort(int n, int (*comp)(), void (*swap)());
void flt_index_swap(int i, int j);
int flt_ascend_comp(int i, int j);
int flt_descend_comp(int i, int j);
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
int flt_ascend_comp(int i, int j)
{
  return(
  flt_neuron_output[closest_index[i]]-flt_neuron_output[closest_index[j]] < 0. ? -1 : 
  (flt_neuron_output[closest_index[i]]-flt_neuron_output[closest_index[j]]>0.) ? 1 : 0);
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
int flt_descend_comp(int i, int j)
{
  return(
  flt_neuron_output[closest_index[i]]-flt_neuron_output[closest_index[j]] < 0. ? 1 : 
  (flt_neuron_output[closest_index[i]]-flt_neuron_output[closest_index[j]]>0.) ? -1 : 0);
}

/************************************************************************
 *                                                                      *
 *  int flt_index_swap(i, j)                                            *
 *                                                                      *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    i = index of 1st element                                          *
 *    j = index of 2nd element                                          *
 *                                                                      *
 ************************************************************************/
void flt_index_swap(int i, int j)
{
   int temp;
   temp             = closest_index[i];
   closest_index[i] = closest_index[j];
   closest_index[j] = temp;
   
}
   
@implementation FltFeatureMap

/*##############################*
 * Initialize connections,      *
 * ride super's method.         *
 *##############################*/
- initConnections
{
  int i;

//
// Allocate work arrays
//
  floatDeltaWeights  = (float *)NXZoneMalloc(myZone,
                        numberWeightsPerNeuron*sizeof(float));
  floatOutputWeights = (float *)NXZoneMalloc(myZone,
                        numberWeightsPerNeuron*sizeof(float));
  probabilityScale = (float *)NXZoneMalloc(myZone,
                        numberOutputWeights*sizeof(float));
//
// create the nodes
//
  if(mapNeurons!=NULL)
    NXZoneFree(myZone,mapNeurons);
  mapNeurons  = (id *)NXZoneMalloc(myZone,
                              (numberMapNeurons*sizeof(FltNeuron)));
  for(i=0; i<numberMapNeurons; i++)
      mapNeurons[i] = [[FltNeuron allocFromZone:myZone]
         initWithNumberInputs:numberWeightsPerNeuron zone:myZone random:randomObject];
      
  return self;
}

/*##############################*
 * Update input neurons with new*
 * input, and feed to map. Also,*
 * find the new winning neuron  *
 * Override super's method.     *
 *##############################*/
- applyInput:(void *)input
{
  int   i;
  float max_activity, min_activity, activity;
   
  floatMapInputs = (float *)input;
  winnerNumber = 0;
  if(feedforward)       /* use only input weights */
  {
    max_activity = min_activity = 
       [mapNeurons[0] stepStart:0 length:numberInputWeights input:floatMapInputs];
    if(neuronType == EUCLIDEAN_DISTANCE)
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] stepStart:0 length:numberInputWeights 
                                      input:floatMapInputs];
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
                                      input:floatMapInputs];
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
    max_activity = min_activity = [mapNeurons[0] step:floatMapInputs];
    if(neuronType == EUCLIDEAN_DISTANCE)
    {
      for(i=1; i<numberMapNeurons; i++)
      {
        activity = [mapNeurons[i] step:floatMapInputs];
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
        activity = [mapNeurons[i] step:floatMapInputs];
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
 * Override super's method.     *
 *##############################*/
- setInput:(void *)input
{   
  floatMapInputs = (float *)input;
  return self;
}


/*##############################*
 * Find the runner up neuron    *
 * (i.e., 2nd highest activity  *
 * neuron;                      *
 * Override super's method.     *
 *##############################*/
- updateRunnerUp
{
  int    i;
  float  max_activity, min_activity, activity;
   
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
 * a single neuron.             *
 * Override super's method.     *
 *##############################*/
- updateWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  float  *old_weights;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
    floatDeltaWeights[i] = etaInput*(floatMapInputs[i] - old_weights[i]);
  for(; i<numberWeightsPerNeuron; i++)
    floatDeltaWeights[i] = etaOutput*(floatMapInputs[i] - old_weights[i]);
  [mapNeurons[neuronNumber] changeAllWeightsBy: floatDeltaWeights];

  return self;
}

/*##############################*
 * Update the input weights     *
 * of a single neuron:          *
 * Override super's method.     *
 *##############################*/
- updateInputWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  float  *old_weights;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
    floatDeltaWeights[i] = etaInput*(floatMapInputs[i] - old_weights[i]);
  for(; i<numberWeightsPerNeuron; i++)
    floatDeltaWeights[i] = 0.;
  [mapNeurons[neuronNumber] changeAllWeightsBy: floatDeltaWeights];
  
  return self;
}

/*##############################*
 * Update the output weights    *
 * of a single neuron:          *
 * Override super's method.     *
 *##############################*/
- updateOutputWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  float  *old_weights;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
    floatDeltaWeights[i] = 0.;
  for(; i<numberWeightsPerNeuron; i++)
    floatDeltaWeights[i] = etaOutput*(floatMapInputs[i] - old_weights[i]);
  [mapNeurons[neuronNumber] changeAllWeightsBy: floatDeltaWeights];

  return self;
}

/*##############################*
 * Normalize the output weights *
 * of the given neuron such that*
 * the sum is equal to synapse- *
 * ScaleFactor:                 *
 * Override super's method.     *
 *##############################*/
- normalizeOutputWeightsOfNeuron:(int)neuronNumber
{
  int    i;
  float  weight_sum, *old_weights, scale;
  
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  weight_sum = 0.;
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
  {
    old_weights[i] = MAX(0., old_weights[i]);
    weight_sum    += old_weights[i];
  }
  if(weight_sum>0.)
    scale = synapseScaleFactor/weight_sum;
  else
    scale = 0.;
  for(i=0; i<numberInputWeights; i++)
    floatDeltaWeights[i] = old_weights[i];
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
    floatDeltaWeights[i] = scale*old_weights[i];
  [mapNeurons[neuronNumber] setAllWeightsTo: floatDeltaWeights];

  return self;
}

/*##############################*
 * Normalize the output weights *
 * of all the neurons such that *
 * the sum over all the neurons *
 * for a particular state is    *
 * equal to synapseScaleFactor. *
 * This is used in ML criterion *
 * Override super's method.     *
 *##############################*/
- normalizeMLOutputWeights
{
  int i, j; 
  float weight_sum, new_weight;
  float scale = 1.;

// Loop over the output (probability) weights:
  for(i=numberInputWeights; i<numberWeightsPerNeuron; i++)
  {
    weight_sum = 0.;
    for(j=0; j<numberMapNeurons; j++)
      weight_sum += MAX(0., [mapNeurons[j] getWeight:i]);
    if(weight_sum>0.)
      scale = synapseScaleFactor/weight_sum;
    for(j=0; j<numberMapNeurons; j++)
    {
      new_weight = scale*MAX(0., [mapNeurons[j] getWeight:i]);
      [mapNeurons[j] setWeight:i to:new_weight];
    }
  }

  return self;
}

/*##############################*
 * Get a list of the            *
 * neuron #'s closest to        *
 * the last input.  I.e.,       *
 * winner = list[0]             *
 * Override super's method.     *
 *##############################*/
- (int *)getNearestNeurons
{
  int i;
  static int old_number=0;
  
  if( old_number<numberMapNeurons )
  {
    old_number = numberMapNeurons;
    if(closest_index!=NULL)
    {
      free(closest_index);
      free(flt_neuron_output);
    }
    flt_neuron_output  = (float *)NXZoneMalloc(myZone,
                                   (numberMapNeurons*sizeof(float)));
    closest_index  = (int *)NXZoneMalloc(myZone,
                                   (numberMapNeurons*sizeof(int)));
  }
  for(i=0; i<numberMapNeurons; i++)
  {
    flt_neuron_output[i] = [mapNeurons[i] lastOutput];
    closest_index[i] = i;
  }
  if(neuronType == EUCLIDEAN_DISTANCE)
    bubble_sort(numberMapNeurons, flt_ascend_comp, flt_index_swap);
  else
    bubble_sort(numberMapNeurons, flt_descend_comp, flt_index_swap);
  return closest_index;
}

/*##############################*
 * Get class of input vector:   *
 * Override super's method.     *
 *##############################*/
- (int)getInputClass
{
  int i,j;
  float max_class, max_value, new_value;
//
// Find the class of the input vector by looking at 'class' neurons
//
  max_class = 0;
  max_value = floatMapInputs[numberInputWeights];
  j = 1;
  for(i=numberInputWeights+1; i<numberWeightsPerNeuron; i++)
  {
    new_value = floatMapInputs[i];
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
 * Return the reverse flow      *
 * output by getting the output *
 * weights of the winning       *
 * neuron:                      *
 * Override super's method.     *
 *##############################*/
- (void *)getReverseFlow
{
  int    i, n, start;
  int   *nearest_neurons, runner_up[2];
  static float *distances=NULL;
  double  sum, temp, smooth_minus_one, total_distance;
  
  if(distances==NULL)
    distances = (float *)NXZoneMalloc(myZone,(MAX_SMOOTH*sizeof(float)));
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
      floatOutputWeights[i] = sum/total_distance;
      start++;
    }
  }
  else
  {
    start = numberInputWeights;
    for(i=0; i<numberOutputWeights; i++)
      floatOutputWeights[i] = [mapNeurons[winnerNumber] getWeight:start++];
  }
    
  return floatOutputWeights;
}

/*##############################*
 * Return the synaptic          *
 * weights of the winning       *
 * neuron:                      *
 * Override super's method.     *
 *##############################*/
- (void *)getWinnerWeightsStart:(int) start length:(int)length
{
  int    i;
  id     winner_neuron;
  
  winner_neuron = mapNeurons[winnerNumber];
  for(i=0; i<length; i++)
    floatOutputWeights[i] = [winner_neuron getWeight:start++];
  
  return floatOutputWeights;
}


/*##############################*
 * Return the synaptic          *
 * weights of the given         *
 * neuron:                      *
 * Override super's method.     *
 *##############################*/
- (void *)getNeuronWeightsStart:(int)start length:(int)length
                        neuron:(int)neuronNumber
{
  int    i;
  id     neuron;
  
  neuron  = mapNeurons[neuronNumber];
  for(i=0; i<length; i++)
    floatOutputWeights[i] = [neuron getWeight:start++];
  
  return floatOutputWeights;
}

/*##############################*
 * Calculate the scale factors  *
 * such that the sum over all   *
 * the neurons for a particular *
 * state is equal to 1.         *
 * Override super's method.     *
 *##############################*/
- (float *)probabilityScale
{
  int i, j, k;
  float weight_sum;

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
 * Override super's method.     *
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
    floatDeltaWeights[weight_index++] = synapseScaleFactor*temp;
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
      printf("Too many input weights in readWeightsFile\n");
      return NO;
    }
    floatDeltaWeights[weight_index++] = synapseScaleFactor*temp;
      }
    }
    if(weight_index >= numberWeightsPerNeuron)
    {
      weight_index = 0;
      [mapNeurons[neuron_number++] setAllWeightsTo:floatDeltaWeights];
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
 * Override super's method.     *
 *###############################*/
- (BOOL)saveWeightsFile: (char *)weightFile
{
  int   i, j;
  char  weight_string[20];
  float temp;
  float *weights;
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
      fprintf(file_write,"%4.3f ", temp);
      if(!(j%LINE_LENGTH))
        fprintf(file_write,"\n");
    }
    fprintf(file_write,">\n");
  }
  fclose(file_write);
  return YES;
}

      
@end