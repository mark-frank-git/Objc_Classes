/***************************************************************************
 * This subclass of object performs a conversion between NIST phoneme      *
 * classes, and an integer number.                                         *
 *                                                                         *
 * File: /User/frank/Objc_Classes/ProcessedSound/NISTPhonemeMap.h          *
 *                                                                         *
 * Revision History:                                                       *
 *  1. 07/03/92  - Started                                                 *
 *  2. 10/30/92  - Added functionality for checking phoneme classes.       *
 ***************************************************************************/

#import <objc/Object.h>

#define PHONEME_LIST_LENGTH 62  /* Total # of phonemes */
#define MAX_PHONE_LENGTH    4   /* phoneme descriptor length    */
#define MAX_PHONEME_COUNT   100 /* max phonemes in a file */

/***** Phoneme classes  ********/
#define ALL_CLASSES 0
#define CLOSURE     1
#define STOP        2
#define AFFRICATE   3
#define FRICATIVE   4
#define NASAL       5
#define SEMI_GLIDE  6
#define VOWEL       7
#define NON_PHONEME 8

@interface NISTPhonemeMap:Object
{
  id    hashTable;          /* Conversion look up table */
  int   desiredClass;           /* phoneme class of interest    */
  int   validClassStart, validClassEnd; /* Used for checking valid class*/
  int   phonemeStates;          /* states/phoneme in HMM    */
  int   validFrames[MAX_PHONEME_COUNT]; /* valid frame endpoints    */
  int   *frameLabels;           /* phoneme labels/frame this is */
                    /* indexed at 0 for training the*/
                    /* output wgts of feature map   */
}

/********************************
 * Initialization:      *
 ********************************/
- init;

/********************************
 * These methods get info from  *
 * this object:                 *
 ********************************/
- (int) getPhonemeIndex:(char *)phoneme;
- (int) numberPhonemes;
- (int) numberValidFrames;
- (int) phonemeStates;
- (int) getFrameLabel:  (int)frame;
- (BOOL) frameValid:    (int)frame;
- (BOOL) isNoisePhoneme: (char *)phoneme;

/*******************************
 * These methods set instance   *
 * variables:           *
 ********************************/
- setDesiredPhonemeClass:(int)class;
- setValidFramesFrom: (int *)phonemeLabels length:(int)length;
- setPhonemeStates: (int)states;


@end
