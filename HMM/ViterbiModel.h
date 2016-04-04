/***********************************************************************
 * This subclass of HMMModel allows for training of the model.          *
 *                                                                      *
 * File: ViterbiModel.h                                                 *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/09/92  - Started                                              *
 *  2. 05/27/92  - Use continuous transition probabilities              *
 *  3. 12/09/92  - Udates for training on TIMIT phonemes                *
 *  4. 01/15/92  - Enhancements from Ed Srenger's code.                 *
 *               - remove reestimation code, see OldViterbiModel if     *
 *                 needed                                               *
 *  5. 02/12/93  - Add -addSkipState                                    *
 *  6. 02/18/93  - Add weights: to initWithLength:...                   *
 *  7. 05/11/93  - Make stateProbability 2-D for score variance metric. *
 *  8. 07/30/93  - Deleted score variance stuff.                        *
 *  9. 08/06/93  - Add -setModelLength:                                 *
 * 10. 08/24/93  - Abstracted out HMMModel                              *
 ***********************************************************************/
#import <HMMModel.h>


@interface ViterbiModel:HMMModel
{
  int   **transitionsCount;                     /* aij Training counters                */
  int   ***observationsCount;                   /* bij Training counter                 */
  int   *initialStateCount;                     /* piInitial Training counters          */
  int   *finalStateCount;           
  BOOL  newTrainingFile;                        /* YES = next input = initial state     */
}

/*****************************************
 * Initialize the local  instance        *
 * variables.                            *
 *****************************************/
- initWithLength:(int)length name:(char *)name states:(int *)states
                 codebooks:(int)numberCodebooks weights:(float *)weights
                 zone:(NXZone *)zone;
- initCodebookDimensions: (int *)dimensions withWeights:(float *)weights;

/*****************************************
 * Alloc and free memory:                *
 *****************************************/
- allocateArrays;
- freeArrays;

/*****************************************
 * Initialization before training:       *
 *****************************************/
- initTrainingCounts;
- initTrainingFile;

/*****************************************
 * Processing training data:             *
 *****************************************/
- initialTrainingStep:(int)state indices:(short *)indices;
- initialTrainingStep:(int)state finalFrame:(BOOL)final;

/*****************************************
 * Updating the model:                   *
 *****************************************/
- addSkipState;
- updateInitialModel;
- normalizeProbabilities;

/*****************************************
 * Getting Parameters/outputs            *
 *****************************************/

/*****************************************
 * Setting Parameters                    *
 *****************************************/

@end