/************************************************************************
 * This subclass of IntFeatureMap implements special feature map-type   *
 * algorithms other than the self-organized learning.                   *
 *                                                                      *
 * File:SpecialMaps.h                                                   *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/08/93  - LVQMap -> SpecialMaps                                *
 ************************************************************************/

#import "IntFeatureMap.h"

// Learning types
#define SELF_ORGANIZING     0
#define LVQ_LEARNING        1
#define RANGE_ERR_LEARNING  2
#define LVQ_MOD_LEARNING    3       /* modified LVQ */
#define LBG         4       /* LBG clustering   */
#define MOD_PLASTICITY      5       /* modulated plasticity */
#define CYCLONNE        6       /* Simulate the hardware */


@interface SpecialMaps:IntFeatureMap

{
  int   learningType;       /* SOM, LVQ, Range err, etc.    */
  int   plasticitySelect;
  int   *meanCount;         /* counter for training     */
  int   **classCount;       /* class counter for training   */
  BOOL  converged;          /* LBG training has converged   */
  double *localPlasticities;    /* Plasticities local to each neuron */
  double **means;               /* mean values          */
  double convergenceEpsilon;    /* Convergence criterion    */
  double oldDistortion;         /* used in LBG          */
  double epsilon;               /* defines learning window */
  double lbgDistortion;        /* Convergence distortion   */
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
- updateWeightsUsingCyclonne;
- negativeUpdateWeightsOfNeuron: (int) neuronNumber;
- negativeUpdateInputWeightsOfNeuron: (int) neuronNumber;
- updateWeightsOfNeuron: (int) neuronNumber withEtaScale: (double)etaScale;
- updateWeightsOfCyclonneNeuron: (int)neuronNumber;
- moveWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron;
- moveInputWeightsOfNeuron: (int) neuronNumber toward: (int)toNeuron;

/**********************
 * Set parameters:    *
 **********************/
- setConvergenceEpsilon: (int)newEpsilon;
- setConverged: (BOOL)flag;
- setEpsilon: (double)newEpsilon;
- setLearningType: (int)type;
- setPlasticitySelect: (int)index;

/**********************
 * Get parameters:    *
 **********************/
- (BOOL)isConverged;
- (double)getEpsilon;
- (int)getLearningType;
- (int)getConvergenceEpsilon;
- (double)getDistortion;
- (int)getClassOfNeuron:(int)neuronNumber;


/********************
 * Dummy:           *
 ********************/
- (int *) trainingData;

@end