/***********************************************************************
 * This subclass of object implements an HMM recognition algorithm.     *
 * File: HMM.h                                                          *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/09/92  - Started                                              *
 *  2. 04/24/92  - Went to integer probability inputs                   *
 *  3. 06/02/92  - Eliminated debug                                     *
 *  4. 06/04/92  - Added forcedViterbi                                  *
 *  5. 12/08/92  - Changes going to TIMIT                               *
 *  6. 01/15/92  - Enhancements from Ed Srenger's code.                 *
 *  7. 02/12/93  - Add -addSkipState                                    *
 *  8. 02/18/93  - Two dimensional input probabilities                  *
 *  9. 03/17/93  - Add probabilityScale                                 *
 * 10. 03/26/93  - Add recognitionStepWithFloatProbs                    *
 * 11. 04/02/93  - Add (float **)probabilityScale                       *
 * 12. 05/11/93  - Add getMinVarianceModel                              *
 * 13. 07/30/93  - Delete get MinVarianceModel                          *
 * 14. 08/06/93  - Add stateCounter, -initModels, scaleType, etc.       *
 ***********************************************************************/


#import <objc/Object.h>
#import <stdio.h>

/****************************************
 * Constants:                           *
 ****************************************/

@interface HMM:Object
{
  id    *wordModels;                        /* The word models              */
  int   numberModels;                       /* dimension of above           */
  int   numberInputProbabilities;           /* dimension of state probs     */
  int   frameCount;                         /* # of frames in current word  */
  int   bestModel, forcedModel;             /* Current best model           */
  int   codebooks;                          /* Number of codebooks          */
  int   modelLength;                        /* Length of model, # of states */
  int   *stateCounter;                      /* Counter for state scaling    */

  float *stateProbabilities;                /* Derived from state counter   */
  BOOL  forcedViterbi;                      /* use input best model         */
  float **probabilityScale;                 /* scale for feature map outputs*/
  double bestScore;                         /* best model's score           */
  double **logProbabilities;                /* log of input probabilities   */
  NXZone *myZone;                           /* malloc zone                  */
}

/*****************************************
 * Initialize the local  instance        *
 * variables.                            *
 *****************************************/
- initWithNumberModels:(int)models numberStates:(int)states
         numberCodebooks:(int)codebooks codebookDimensions:(int *)dimensions
         codebookWeights:(float *)weights zone:(NXZone *)zone;
- allocateArrays;
- free;
- freeArrays;

/*****************************************
 * Initialization before training:       *
 *****************************************/
- initNewSearch;
- initTrainingCounts;
- initTrainingCountsFor:(int)modelNumber;
- initTrainingFileFor:(int)modelNumber;

/*****************************************
 * Running                               *
 *****************************************/
- recognitionStepWithProbs: (int **)inputProbabilities;
- recognitionStepWithFloatProbs: (float **)inputProbabilities;
- recognitionStepWithIndices: (short *)indices; /* Use during training  */
                                                /* or recognition   */
- initialTrainingStepFor:(int)modelNumber state:(int)state
            codebookIndices: (short *)indices;  /* Use during initialization*/
- initialTrainingStepFor:(int)modelNumber state:(int)state
              finalFrame:(BOOL)flag;

/*****************************************
 * Setting parameters:                   *
 *****************************************/
- setForcedModel:       (int)model;
- setForced:            (BOOL)flag;
- setProbabilityScale:  (float *)scale codebook:(int)number;
- setNumberModelStates: (int)states;

/*****************************************
 * Getting parameters:                   *
 *****************************************/
- (int) getBestModel;
- (int *) getBestBackTrace;
- (double) getBestScore;
- getModel: (int)model;
- (int) getFrameCount;
- (int) getNumberCodebooks;
- (int) getNumberStates;
- (int) getNumberModels;
- (int) numberInputProbabilities;
- (BOOL) isForced;
- (float *)probabilityScaleForCodebook:(int)number;
- (float *)stateProbabilities;

/*****************************************
 * Updating the models:                  *
 *****************************************/
- updateInitialModels;
- updateInitialModel:(int)modelNumber;
- addSkipStateFor:(int)modelNumber;

/****************************************
 * Saving/Reading to file:              *
 ****************************************/
- (BOOL)saveProbabilitiesToFile:(char *)fileName model:(int)modelNumber;
- (BOOL)readProbabilitiesFromFile:(char *)fileName model:(int)modelNumber;
- (BOOL)saveProbabilityScaleToFile:(char *)fileName;
- (BOOL)readProbabilityScaleFromFile:(char *)fileName;

/*****************************************
 * Printing the models:                  *
 *****************************************/
- (BOOL)printAllModelsTo:(char *)fileName;
- (BOOL)printModel:(int)modelNumber to:(char *)fileName;

@end
