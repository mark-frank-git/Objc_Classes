/***********************************************************************
 * This subclass of object implements a model used for HMM speech       *
 * recognition.  For training the model, see the subclass,              *
 * ViterbiModel.                                                        *
 *                                                                      *
 * File: HMMModel.h                                                     *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 08/24/93  - Abstracted from ViterbiModel.                        *
 ***********************************************************************/
#import <objc/Object.h>
#import <stdio.h>

/*  Constants:   */
#define MAX_NAME        64                      /* Model name               */
#define MAX_FRAMES      1024                    /* max # of input frames    */
#define MIN_PROBABILITY 1.e-5
#define LOG_MIN         -11.51293               /* Log of above             */
#define SKIP_PROBABILITY 2.e-3                  /* initial skip probability */


@interface HMMModel:Object
{
  int   modelLength;                            /* length of Markov chain               */
  int   *modelStates;                           /* The states for this model            */
  int   *backTrace;                             /* Viterbi back trace                   */
  int   bestModelEnd;                           /* best ending for this model           */
  int   codebooks;                              /* number of codebooks                  */
  int   *psi[MAX_FRAMES];                       /* Viterbi backtrace                    */
  int   *codebookDimensions;                    /* codebook sizes                       */
  int   oldState;                               /* last state in initial trning         */
  int   maxFrames;                              /* Max # of frames seen so far          */
  float *codebookWeights;                       /* relative wgts given to cdbks         */
  float *logPiInitial;                          /* Log of initial state                 */
  float *logPiFinal;                            /* Log of final state probs             */
  float **logA;                                 /* log of transition probs              */
  float ***logB;                                /* log of observations probs            */
  float *stateProbability;                      /* Current probability score            */
  float *maxProbability;                        /* Used in calculating the above        */
  double maxScore;                              /* Max score for this model             */
  BOOL  useCodebookIndices;                     /* YES = codebook indices               */
  char  modelDescription[MAX_NAME];             /* word name                            */
  NXZone *myZone;                               /* Malloc zone                          */
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
- free;
- freeArrays;

/****************************************
 * Initialization before new search     *
 ****************************************/
-initNewSearch;

/*****************************************
 * Processing new input data:            *
 *****************************************/
- recognitionStepWithProbs:(double **)logProbabilities count:(int)frameCount;
- recognitionStepWithIndices:(short *)indices count:(int)frameCount;

/*****************************************
 * Getting Parameters/outputs            *
 *****************************************/
- (int *)   getBackTrace:(int)strobeCount;
- (char *)  getModelDescription;
- (int)     getModelLength;
- (double)  getLastScore;
- (float *) getTransitionProbabilities;
- (float *) getInitialProbabilities;
- (float *) getFinalProbabilities;
- (float *) stateProbability;
- (double)  getModelScore:(int)frameCount;

/*****************************************
 * Setting Parameters                    *
 *****************************************/
- setModelLength:(int)length states:(int *)states;

/*****************************************
 * Reading/Writing to disk              *
 *****************************************/
- (BOOL)saveProbabilitiesToFile:(FILE *)fp;
- (BOOL)readProbabilitiesFromFile: (FILE *)fp;

/*****************************************
 * Printing the model:                   *
 *****************************************/
- printModel:(FILE *)fp;

@end