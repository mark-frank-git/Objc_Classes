/************************************************************************
 * This subclass of object implements a feed forward self-organizing    *
 * feature map.                                                         *
 *                                                                      *
 * File:FltFeatureMap.h                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/24/93  - Modified from IntFeatureMap                          *
 *  2. 08/10/93  - Add probabilityScale                                 *
 *  3. 08/16/93  - Add normalizeMAPOutputWeights, etc.                  *
 *  4. 09/22/93  - Made subclass of IntFeatureMap.                      *
 ************************************************************************/
#import "IntFeatureMap.h"

@interface FltFeatureMap:IntFeatureMap

{
  float *floatMapInputs;                            /* input activities array   */
  float *floatOutputWeights, *floatDeltaWeights;    /* Arrays for I/O and updates   */
}

/**********************
 * Initialization:    *
 **********************/
- initConnections;
        
/**********************
 * Running:           *
 **********************/
- applyInput: (void *)input;
- setInput:   (void *)input;
- updateRunnerUp;
- updateWeightsOfNeuron:      (int)neuronNumber;
- updateInputWeightsOfNeuron: (int)neuronNumber;
- updateOutputWeightsOfNeuron:(int)neuronNumber;
- normalizeOutputWeightsOfNeuron:(int)neuronNumber;
- normalizeMLOutputWeights;

/**********************
 * Set parameters:    *
 **********************/

/**********************
 * Get parameters:    *
 **********************/
- (int *)getNearestNeurons;
- (void *)getReverseFlow;
- (void *)getWinnerWeightsStart:(int)start length:(int)length;
- (void *)getNeuronWeightsStart:(int)start length:(int)length
                neuron:(int)neuronNumber;
- (float *)probabilityScale;


/***********************
 * Reading, writing    *
 * weights files.      *
 ***********************/
- (BOOL)readWeightsFile: (char *)fileName;
- (BOOL)saveWeightsFile: (char *)weightFile;


@end