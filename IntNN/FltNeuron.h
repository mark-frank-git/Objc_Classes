/************************************************************************
 * This subclass of object implements a neuron object.                  *
 *                                                                      *
 * File:FltNeuron.h                                                     *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/24/93  - Modified from IntNeuron                              *
 *  2. 09/22/93  - Make 'random' object a static variable to method     *
 ************************************************************************/
#import <objc/Object.h>


@interface FltNeuron:Object

#define Binary              0                   // Node Types
#define Sigmoid             1
#define Sign                2
#define Tanh                3
#define None                4

#define DOT_PRODUCT         0                   // Neuron types
#define NORMALIZED_DOT_PRODUCT  1
#define EUCLIDEAN_DISTANCE  2
#define MANHATTAN_DISTANCE  3


{
  float *weights;                               // Input weights
  float lastOutput;                             // output value from last time step
    
  id    random;                                 // Random instance
  int   nodeType;                               // the type of this node
  int   neuronType;                             // Type of step calculation
  int   numberInputs;                           // Number of input connections to this neuron
  int   randomMin, randomMax;                   // Random weight ranges

  NXZone *myZone;                               // Malloc zone
}
/***************************
 * Initialization:         *
 ***************************/
- initWithNumberInputs:(int)number;
- initWithNumberInputs:(int)number zone:(NXZone *)zone random:(id)inputRandom;
- initializeWeights;                            // Randomize weights
- freeStorage;


/***************************
 * Running:                *
 ***************************/
- (float)step:(float *)input;                   // read inputs and generate output 
- (float) stepStart:(int)start length:(int)length input:(float *)inputs;
- (float)activation:(float)net;

/***************************
 * Setting parameters:     *
 ***************************/
- setAllWeightsTo:(float *)weight;
- setWeight:(int)weightIndex to:(float)newWeight;
- setOutput:(float)output;                      // set the output value manually (used for inputs)
- setType:(int)type;                            // set the type of node - determines the activation
                                                // function
- setNeuronType: (int)type;                     // set the type of neuron calculation
- setRandom:theRandom;                          // set random number generator pointer
- setRandomMin:(int)min max:(int)max;           // Range on random weights
- changeAllWeightsBy: (float *)delta;

/***************************
 * Getting parameters:     *
 ***************************/
- (float)lastOutput;                            // returns the last output value of the Neuron
- (float *)getAllWeights;                       // returns a pointer to all input weights
- (float) getWeight:(int)index;
- (int)getType;                                 // returns the node type for this Neuron
- (int)getNeuronType;                           // returns the neuron calculation type
@end