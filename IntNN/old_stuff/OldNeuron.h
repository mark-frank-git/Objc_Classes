/************************************************************************
 *  Neural Network Classes for the NeXT Computer                    *
 *  Written by: Ralph Zazula                                        *
 *  University of Arizona                                           *
 *  zazula@bonehead.tucson.az.us (NeXT Mail)                        *
 *                                                                      *
 * Revisions:                                                           *
 *   1.1  01/02/92  12:42:30  zazula                                    *
 *   1.11 04/01/92  added Sigmoid type = None   (frank)                 *
 *        04/06/92  added numberInputs to speed up weight calculations. *
 *        04/16/92  added setAllWeightsTo:                              *
 *        04/23/92  added stepStart:length:                             *
 *   1.2  04/24/92  Changed to new class: IntNeuron (Integer calcs.)    *
 *                  deleted things I wasn't using (temperature, etc.)   *
 *                                                                      *
 ************************************************************************/
#import <objc/Object.h>
#import <objc/List.h>
#import <objc/Storage.h>
#import "Random.h"          // Random number generator

typedef struct connected_list{
    id      source;
    int weight;
    struct  connected_list *next;
} connection;

@interface IntNeuron:Object

#define Binary 0            // Node Types
#define Sigmoid 1
#define Sign 2
#define Tanh 3
#define None 4

#define DOT_PRODUCT     0       // Neuron types
#define NORMALIZED_DOT_PRODUCT  1
#define EUCLIDEAN_DISTANCE  2


{
    id  inputs;     // list of Neuron's that are inputs to this Neuron
    int lastOutput; // output value from last time step
    
    id  random;     // Random instance
    int nodeType;   // the type of this node
    int     neuronType; // Type of step calculation
    int numberInputs;// Number of input connections to this neuron
    int     randomMin, randomMax; // Random weight ranges
    int *intOutputs;
    connection *head;   // the head of the linked list
    connection *tail;   // the tail of the linked list
}
/***************************
 * Initialization:         *
 ***************************/
- init;         // initialization method
- initializeWeights;    // Randomize weights
- freeStorage;

/***************************
 * Creating connections:   *
 ***************************/
- connect:sender;   // add sender to the list of inputs to this Neuron
- connect:sender withWeight:(int)weight;

/***************************
 * Running:                *
 ***************************/
- step;         // read inputs and generate output 
- stepStart:(int)start length:(int)length;
- (int)activation:(int)net;

/***************************
 * Setting parameters:     *
 ***************************/
- setAllWeightsTo:(int *)weight;
- setWeightFor:source to:(int)weight;
            // sets the weight for the given input Neuron
- setOutput:(int)output;// set the output value manually (used for inputs)
- setType:(int)type;// set the type of node - determines the activation
            // function
- setNeuronType: (int)type; // set the type of neuron calculation
- setRandom:theRandom;  // set random number generator pointer
- setRandomMin:(int)min max:(int)max; // Range on random weights
- changeWeightFor:source by:(int)delta;
- changeAllWeightsBy: (int *)delta;

/***************************
 * Getting parameters:     *
 ***************************/
- inputs;       // returns a pointer to the List of inputs  
- (int)lastOutput;  // returns the last output value of the Neuron
- (int)getWeightFor:source;         
            // returns the weight for the given input Neuron
- (int *)getAllWeights; // returns a pointer to all input weights
- (int)getType;     // returns the node type for this Neuron
- (int)getNeuronType;   // returns the neuron calculation type
@end