/************************************************************************
 * This subclass of object implements a neuron object.                  *
 *                                                                      *
 * File:FltNeuron.m                                                     *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/24/93  - Modified from IntNeuron                              *
 *  2. 09/22/93  - Make 'random' object a static variable to method     *
 ************************************************************************/
#import "FltNeuron.h"
#import "Random.h"          // Random number generator
#import <objc/HashTable.h>
#import "math.h"
#import "stdlib.h"
#import "stdio.h"

#define RANDOM_WEIGHT   5   // Random weights [-5,5]
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define INT_SCALE   100


//----------------------------------------------------------

@implementation FltNeuron

- setType:(int)type { nodeType = type; return self; }
- (int)getType { return nodeType; }
- setNeuronType:(int)type {neuronType = type; return self; }
- (int)getNeuronType { return neuronType; }
- setRandom:theRandom { random = theRandom; return self; }
- setRandomMin:(int)min max:(int)max 
         {randomMin=min; randomMax=max; return self;}

//----------------------------------------------------------

- (float)activation:(float)net
{
   float output;
   double temp;

   switch (nodeType) {
   case Binary :
      temp = (double)net/INT_SCALE; 
      output = (temp > 0.5) ? INT_SCALE : 0;
      break;
   case Sigmoid :
      temp = (double)net/INT_SCALE; 
      output = INT_SCALE/(1.0+exp(-temp));
      break;
   case Sign : 
      temp = (double)net/INT_SCALE; 
      output = (temp > 0.0) ? INT_SCALE : -INT_SCALE;
      break;
   case Tanh :
      temp = (double)net/INT_SCALE; 
      output = INT_SCALE*tanh(temp);
      break;
   case None:
   default:
      output = net;
      break;
   }
   
   return output;
}

//----------------------------------------------------------

- initWithNumberInputs: (int)number
{
   [super init];
   numberInputs    = number;
   myZone          = NXDefaultMallocZone();
   weights         = (float *)NXZoneMalloc(myZone, number*sizeof(float));
   lastOutput      = 0.;
   nodeType        = None;            // default node type
   neuronType      = EUCLIDEAN_DISTANCE;  // default step calculation
   randomMin       = -RANDOM_WEIGHT;
   randomMax       = RANDOM_WEIGHT;
   [self initializeWeights];
    
   return self;
}

//----------------------------------------------------------

- initWithNumberInputs: (int)number zone:(NXZone *)zone random:inputRandom
{
   [super init];
   numberInputs    = number;
   myZone          = zone;
   random          = inputRandom;
   weights         = (float *)NXZoneMalloc(myZone, number*sizeof(float));
   lastOutput      = 0.;
   nodeType        = None;            // default node type
   neuronType      = EUCLIDEAN_DISTANCE;  // default step calculation
   randomMin       = -RANDOM_WEIGHT;
   randomMax       = RANDOM_WEIGHT;
   [self initializeWeights];
    
   return self;
}

//-----------------------------------------------------------

- initializeWeights
{
  int        i;
  
  if(random == nil)
   random          = [[Random allocFromZone:myZone] init];
        
  for(i=0; i<numberInputs; i++)
     weights[i] = (float)[random randMin:randomMin max:randomMax];
  return self;
}

- freeStorage
{
  NXZoneFree(myZone, weights);
  return self;
}

//-----------------------------------------------------------

- (float)step:(float *)inputs
// update the output value based on our inputs
{
   int    i;
   float  output;
   float  temp = 0.;
   double norm_input, norm_weight;
   
   switch(neuronType)
   {
     case DOT_PRODUCT:
       for(i=0; i<numberInputs; i++)
         temp += weights[i]*inputs[i];
       break;
     case NORMALIZED_DOT_PRODUCT:
       norm_input = norm_weight = 0.;
       for(i=0; i<numberInputs; i++)
       {
         output       = inputs[i];
         temp        += weights[i]*output;
         norm_input  += output*output;
         norm_weight += weights[i]*weights[i];
       }
       temp /= sqrt(norm_input*norm_weight);
       break;
     case EUCLIDEAN_DISTANCE:
       for(i=0; i<numberInputs; i++)
       {
         output       = inputs[i] - weights[i];
         temp        += output*output;
       }
       break;
     default:
       printf("Unknown neuron type\n");
       break;
   }

   lastOutput = [self activation:temp];
   return lastOutput;
}

//-----------------------------------------------------------

- (float)stepStart:(int)start length:(int)length input:(float *)inputs
// update the output value based on our inputs
{
   int    i, end;
   float  output;
   float  temp=0.;
   double norm_input, norm_weight;
   
   end = start+length;
   switch(neuronType)
   {
     case DOT_PRODUCT:
       for(i=start; i<end; i++)
         temp += weights[i]*inputs[i];
       break;
     case NORMALIZED_DOT_PRODUCT:
       norm_input = norm_weight = 0.;
       for(i=start; i<end; i++)
       {
         output       = inputs[i];
         temp        += weights[i]*output;
         norm_input  += output*output;
         norm_weight += weights[i]*weights[i];
       }
       temp /= sqrt(norm_input*norm_weight);
       break;
     case EUCLIDEAN_DISTANCE:
       for(i=start; i<end; i++)
       {
         output       = inputs[i] - weights[i];
         temp        += output*output;
       }
       break;
     default:
       printf("Unknown neuron type\n");
       break;
   }

   lastOutput = [self activation:temp];
   return lastOutput;
}

//-----------------------------------------------------------

- (float)lastOutput
{
   return lastOutput;
}

//-----------------------------------------------------------


//-----------------------------------------------------------
- (float *)getAllWeights
{
   return weights;
}

//-----------------------------------------------------------
- (float )getWeight:(int)index
{
   return weights[index];
}


//-----------------------------------------------------------
- setAllWeightsTo: (float *)newWeights
{
   int        i;
   
   for(i=0; i<numberInputs; i++)
     weights[i] = newWeights[i];

   return self;
}


//-----------------------------------------------------------
- setWeight:(int)weightIndex to:(float)newWeight
{
   weights[weightIndex] = newWeight;
   return self;
}


//-----------------------------------------------------------

- setOutput:(float)output
{
   lastOutput = output;
   
   return self;
}
//-----------------------------------------------------------


//-----------------------------------------------------------

- changeAllWeightsBy: (float *)delta
{
   int        i;
   
   for(i=0; i<numberInputs; i++)
     weights[i] += delta[i];
   return self;
}

@end
