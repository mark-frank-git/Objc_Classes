/************************************************************************
 * This subclass of object implements a feed forward self-organizing    *
 * feature map.                                                         *
 *                                                                      *
 * File:IntFeatureMap.h                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/27/92  - Started                                              *
 *  2. 04/24/92  - Modified FeatureMapEngine -> IntFeatureMap           *
 *  3. 04/28/92  - Added updateNeuronInput/OutputWeights:               *
 *  4. 05/04/92  - Use rounding on weight updates                       *
 *  5. 05/07/92  - Add separate eta for output weights                  *
 *  6. 06/01/92  - Free output arrays on next call                      *
 *               - Went to arrays for mapNeurons and inputNeurons       *
 *  7. 07/01/92  - Add runnerUp                                         *
 *  8. 09/04/92  - Add learningEnabled                                  *
 *  9. 11/24/92  - Combine updateWinner into applyInput.                *
 *               - Add readWeightsFile, saveWeightsFile.                *
 *               - FEATURE_MAP_SCALE -> synapseScaleFactor              *
 * 10. 12/03/92  - Remove inputNeurons                                  *
 *                 numberInputNeurons -> numberWeightsPerNeuron         *
 *                 numberInputs -> numberInputWeights                   *
 *                 numberOutputs -> numberOutputWeights                 *
 * 11.  03/01/93  - Add normalizeOutputWeightsOfNeuron, *outputWeights, *
 *                  *deltaWeights.                                      *
 *                - Add -setInput                                       *
 * 12.  03/12/93  - Add weightUpdateType                                *
 * 13.  04/02/93  - Add probabilityScale                                *
 * 14.  08/16/93  - Add normalizeMAPOutputWeights, etc.                 *
 * 15.  09/22/93  - Changes to make FltFeatureMap a subclass of this    *
 *                - Add -setWinner                                      *
 ************************************************************************/
#import <objc/Object.h>
#import <objc/zone.h>

#define TRUNCATED       0                           /* weight update types      */
#define ROUNDED         1
#define MINIMUM_STEP    2                           /* This is in the hardware  */
#define STATISTICAL     3
#define INCREMENT       4

#define DEFAULT_SIZE    1
#define DEFAULT_ETA     0.1
#define DEFAULT_SMOOTH  2                           /* # of outputs to smooth   */
#define MAX_SMOOTH      20                          /* Maximum # of outputs     */
#define DEFAULT_SCALE   1000.

#define ETA_EPS         0.00001
#define SEED1           23
#define SEED2           77
#define SEED3           113
#define MAX_LINE        1024                        /* for reading in weights file  */

@interface IntFeatureMap:Object

{
  id    *mapNeurons;                                /* List of neurons in map   */
  id    randomObject;                               /* Random number generator  */
    
  int   neuronType;
  int   rows, columns;                              /* # of rows, cols in map   */
  int   numberMapNeurons;                           /* # of neurons in map      */
  int   numberWeightsPerNeuron;                     /* # of neurons in input lyer   */
  int   numberInputWeights, numberOutputWeights;    /* input vs output synapses */
  int   neighborhoodSize;                           /* radius of neighborhoods  */
  int   winnerNumber;                               /* Current winning neuron   */
  int   runnerUp;                                   /* runner up to winner      */
  int   numberSmooth;                               /* # of neurons for smoothing   */
  int   weightUpdateType;                           /* Type of delta wgt calculat.  */
  BOOL  feedforward;                                /* feedforward calc. flag   */
  BOOL  smoothOutput;                               /* smooth reverse flow output   */
  BOOL  learningEnabled;                            /* learning turned on       */

  int   *mapInputs;                                 /* input activities array   */
  int   *outputWeights, *deltaWeights;              /* Arrays for I/O and updates   */
  float *probabilityScale;                          /* For scaling output probs */

  NXZone  *myZone;                                  /* zone for mallocs     */
    
  double  etaInput, etaOutput;                      /* learning rates       */
  double  synapseScaleFactor;                       /* scaling for I/O      */
}

/**********************
 * Initialization:    *
 **********************/
- init;
- initWithInputs:(int)totalInputs outputs:(int)numberOutputs
            rows:(int)numberRows columns:(int)numberCols;
- initWithInputs:(int)totalInputs outputs:(int)numberOutputs
            rows:(int)numberRows columns:(int)numberCols zone:(NXZone*)zone;
- initConnections;
- freeConnections;
- initializeWeights;
        
/**********************
 * Running:           *
 **********************/
- applyInput: (void *)input;
- setInput:   (void *)input;
- updateWeights;
- updateWinner;
- updateRunnerUp;
- updateWeightsOfNeuron:      (int)neuronNumber;
- updateInputWeightsOfNeuron: (int)neuronNumber;
- updateOutputWeightsOfNeuron:(int)neuronNumber;
- normalizeMAPOutputWeights;
- normalizeOutputWeightsOfNeuron:(int)neuronNumber;
- normalizeMLOutputWeights;
- (int *)winnerNeighbors;

/**********************
 * Set parameters:    *
 **********************/
- setFeedforward: (BOOL)flag;
- setLearningInput:(double)newInputEta output:(double)newOutputEta;
- setNeuronType: (int) type;
- setNumberSmooth: (int) number;
- setRadius: (int) newSize;
- setSmoothOutput: (BOOL)flag;
- setNumberWeightsPerNeuron:(int)number;
- setReverseFlowOutputs:(int)number;
- setNumberRows:(int)numberRows columns:(int)numberCols;
- setScaleFactor:(double)scale;
- setWeightUpdateType: (int)type;
- setWinner:(int)newWinner;

/**********************
 * Get parameters:    *
 **********************/
- (int)getInputClass;
- getLearningInput:(double *)inputEta output:(double *)outputEta;
- (int *)getNearestNeurons;
- getNeuron: (int)neuronNumber;
- (int) getNeuronType;
- (int) getRadius;
- (int) numberMapNeurons;
- (void *)getReverseFlow;
- (int) getWinner;
- (void *)getWinnerWeightsStart:(int)start length:(int)length;
- (void *)getNeuronWeightsStart:(int)start length:(int)length
                neuron:(int)neuronNumber;
- (BOOL)isFeedforward;
- (double)getScaleFactor;
- (float *)probabilityScale;


/***********************
 * Reading, writing    *
 * weights files.      *
 ***********************/
- (BOOL)readWeightsFile:(char *)fileName;
- (BOOL)saveWeightsFile: (char *)weightFile;


@end