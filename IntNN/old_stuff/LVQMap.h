/************************************************************************
 * This subclass of IntFeatureMap implements the LVQ2 supervised learn- *
 * ing algorithm.  Normally, the feature map should be initially    *
 * trained using self organization, and then tuned using LVQ.       *
 *                                  *
 * File:LVQMap.h                            *
 *                                  *
 * Revision history:                            *
 *  1. 05/21/92  - Started                      *
 *  3. 06/30/92  - Add LVQ_MOD                      *
 *  3. 07/17/92  - Add K Means algorithm                *
 *  4. 08/03/92  - Add MOD_PLASTICITY                   *
 *  5. 11/15/92  - Add LBG algorithm                    *
 *  6. 11/18/92  - Add getDistortion                    *
 ************************************************************************/

#import "IntFeatureMap.h"

// Learning types
#define SELF_ORGANIZING     0
#define LVQ_LEARNING        1
#define RANGE_ERR_LEARNING  2
#define LVQ_MOD_LEARNING    3       /* modified LVQ */
#define LBG         4       /* LBG clustering   */
#define MOD_PLASTICITY      5       /* modulated plasticity */


@interface LVQMap:IntFeatureMap

{
  int   learningType;       /* SOM, LVQ, Range err, etc.*/
  int   *meanCount;     /* counter for training     */
  int   **classCount;       /* class counter for training   */
  BOOL  converged;      /* training has converged   */
  double *localPlasticities;    /* Plasticities local to each neuron */
  double **means;       /* mean values          */
  double convergenceEpsilon;    /* Convergence criterion    */
  double oldDistortion;     /* used in LBG          */
  double epsilon;       /* defines learning window */
  double lbgDistortion;     /* Convergence distortion   */
}

/**********************
 * Initialization:    *
 **********************/
- init;
- initConnections;
- initTrainingCounts;
- freeConnections;
        
/**********************
 * Training:          *
 **********************/
- updateWinner;
- updateWeights;
- updateWeightsUsingLVQ;
- updateWeightsUsingLVQMod;
- updateWeightsUsingRangeErr;
- updateWeightsUsingLBG;
- updateWeightsUsingModulatedPlasticity;
- negativeUpdateWeightsOfNeuron: (int) neuronNumber;
- negativeUpdateInputWeightsOfNeuron: (int) neuronNumber;
- updateWeightsOfNeuron: (int) neuronNumber withEtaScale: (double)etaScale;
- moveWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron;
- moveInputWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron;

/**********************
 * Set parameters:    *
 **********************/
- setConvergenceEpsilon: (int)newEpsilon;
- setConverged: (BOOL)flag;
- setEpsilon: (double)newEpsilon;
- setLearningType: (int)type;

/**********************
 * Get parameters:    *
 **********************/
- (BOOL)isConverged;
- (double)getEpsilon;
- (int)getClassOfNeuron:(int)neuronNumber;
- (int)getLearningType;
- (int)getConvergenceEpsilon;
- (double)getDistortion;


/********************
 * Dummy:           *
 ********************/
- (int *) trainingData;

@end