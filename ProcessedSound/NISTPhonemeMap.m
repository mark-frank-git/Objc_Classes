/***************************************************************************
 * This subclass of object performs a conversion between NIST phoneme     *
 * classes, and an integer number.                                        *
 *                                                                        *
 * File: /User/frank/Objc_Classes/ProcessedSound/NISTPhonemeMap.m         *
 *                                                                        *
 * Revision History:                                                      *
 *  1. 07/03/92  - Started                                                *
 ***************************************************************************/

#import "NISTPhonemeMap.h"
#import <objc/HashTable.h>
#import <stdio.h>
#import <stdlib.h>
#import <strings.h>

char *phoneme_list[PHONEME_LIST_LENGTH] = 
  {"bcl", "dcl", "gcl", "pcl", "tcl", "kcl",            /*closures*/
   "b", "d", "g", "p", "t", "k", "dx", "q",         /*stops   */
   "jh", "ch",                          /*affricates*/
   "s", "sh", "z", "zh", "f", "th", "v", "dh",          /*fricatives*/
   "m", "n", "ng", "em", "en", "eng", "nx",         /*nasals*/
   "l", "r", "w", "y", "hh", "hv", "el",            /*semi+glides*/
   "iy", "ih", "eh", "ey", "ae", "aa", "aw", "ay", "ah", "ao",
   "oy", "ow", "uh", "uw", "ux", "er", "ax", "ix", "axr", "ax-h", /*vowels*/
   "pau", "epi", "h#", "1"};                    /*others*/

/****!!!!!!!!!!!!!! These must match the above !!!!!!!!!!!!!!!****/
#define CLOSURE_START       0
#define STOP_START      6
#define AFFRICATE_START     14
#define FRICATIVE_START     16
#define NASAL_START     24
#define SEMI_GLIDE_START    31
#define VOWEL_START     38
#define NON_PHONEME_START   58
#define NOISE_PHONEME       60  

#define PHONEME_STATES      3
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

@implementation NISTPhonemeMap

/*****************************************
 * Create the NISTPhonemeMap and    *
 * initialize instance variables:   *
 *****************************************/
- init
{
  int i;
  
  [super init];
// Initialize Hash Table
  hashTable = [[HashTable alloc] initKeyDesc: "*" valueDesc:"i" 
                            capacity:PHONEME_LIST_LENGTH];

  for(i=0; i<PHONEME_LIST_LENGTH; i++)
    [hashTable insertKey:phoneme_list[i] value:(void *)(i+1)];  

  desiredClass  = ALL_CLASSES;
  phonemeStates = PHONEME_STATES;
// set valid frames to entire length;
  validFrames[0] = 0;
  validFrames[1] = 32000;
  validFrames[2] = -1;

  return self;
}

/*######################################*
 * Return the value of the given    *
 * phoneme.             *
 *######################################*/
- (int)getPhonemeIndex:(char *)phoneme
{
  return (int)[hashTable valueForKey:phoneme];
}

/*######################################*
 * Return the # of phonemes in the valid*
 * class.               *
 *######################################*/
- (int)numberPhonemes
{
  return (validClassEnd - validClassStart);
}

/*######################################*
 * Return the # of valid frames in the  *
 * current file:            *
 *######################################*/
- (int)numberValidFrames
{
  int i, number;

  number = 0;;
  for(i=0; i<MAX_PHONEME_COUNT; i+=2)
  {
    if(validFrames[i] == -1)
      break;
    else
      number += validFrames[i+1] - validFrames[i];
  }
  return number;
}

/*######################################*
 * Return the # of states/phoneme:      *
 *######################################*/
- (int)phonemeStates { return phonemeStates; }

/*######################################*
 * Return the phoneme frame label for   *
 * training the output wgts of the  *
 * feature map.             *
 * NOTE: no check if frameLabels==NULL  *
 *######################################*/
- (int)getFrameLabel:(int)frame
{
  return (frameLabels[MAX(frame,0)]);
}

/*######################################*
 * Check if input frame is valid:       *
 *######################################*/
- (BOOL) frameValid: (int)frame
{
  int i;
  BOOL valid;

  valid = NO;
  for(i=0; i<MAX_PHONEME_COUNT; i+=2)
  {
    if(validFrames[i] == -1)
      break;
    if( (frame>=validFrames[i]) && (frame<=validFrames[i+1]) )
    {
      valid = YES;
      break;
    }
  }
  return valid;
}

/*######################################*
 * Check if input phoneme is noise:     *
 *######################################*/
- (BOOL) isNoisePhoneme: (char *)phoneme
{
  if(strcmp(phoneme, phoneme_list[NOISE_PHONEME]) == 0)
    return YES;
  else
    return NO;
}
  
/*######################################*
 * Set the phoneme class desired by the *
 * sender.                              *
 *######################################*/
- setDesiredPhonemeClass:(int)class
{
  desiredClass = class;
  switch(desiredClass)
  {
    default:
    case ALL_CLASSES:
      validClassStart = 0;
      validClassEnd   = PHONEME_LIST_LENGTH;
      break;
    case CLOSURE:
      validClassStart = CLOSURE_START;
      validClassEnd   = STOP_START;
      break;
    case STOP:
      validClassStart = STOP_START;
      validClassEnd   = AFFRICATE_START;
      break;
    case AFFRICATE:
      validClassStart = AFFRICATE_START;
      validClassEnd   = FRICATIVE_START;
      break;
    case FRICATIVE:
      validClassStart = FRICATIVE_START;
      validClassEnd   = NASAL_START;
      break;
    case NASAL:
      validClassStart = NASAL_START;
      validClassEnd   = SEMI_GLIDE_START;
      break;
    case SEMI_GLIDE:
      validClassStart = SEMI_GLIDE_START;
      validClassEnd   = VOWEL_START;
      break;
    case VOWEL:
      validClassStart = VOWEL_START;
      validClassEnd   = NON_PHONEME_START;
      break;
    case NON_PHONEME:
      validClassStart = NON_PHONEME_START;
      validClassEnd   = PHONEME_LIST_LENGTH;
      break;
  }
  return self;
}

/*######################################*
 * Set the valid frame regions from the *
 * input labels.  The label array is of *
 * form: 0, 0, 0, 3, 0, 0, -1, 0, 0,    *
 * where the numbers>0  are the labels, *
 * at mid phoneme, and the -1s indicate *
 * the end of the phoneme.              *
 * The validFrames array is of the form:*
 *   {start1, end1, start2, end2, -1}   *
 *######################################*/
- setValidFramesFrom: (int *)phonemeLabels length:(int)length
{
  BOOL valid_frame;
  int i, j, k, n, start;
  int label, frame_length, state_length;
  static int old_length = 0;
  static int *index = 0;

  if(length>old_length)
  {
    old_length = length;
    if(index == NULL)
    {
      free(index);
      free(frameLabels);
    }
    index = (int *)malloc(length*sizeof(int));
    frameLabels = (int *)malloc(length*sizeof(int));
  }

// Find the valid frames:
  i = j = 0;
  while(i<length)
  {
    start = i;
    valid_frame = NO;
    for( ;i<length; i++)
    {
      if(phonemeLabels[i] == -1)
      {
        i++;
        break;
      }
      if(phonemeLabels[i] > 0)
      {
        index[j] = phonemeLabels[i] - 1;
        if( (index[j]>=validClassStart) && (index[j]<validClassEnd) )
        {
          index[j]   -= validClassStart;
          valid_frame = YES;
        }
      }
    }
    if(valid_frame)
    {
      validFrames[j++] = start;
      validFrames[j++] = i-1;
    }
  }
  validFrames[j] = -1;

// Now, label the frames
  for(i=0; i<length; i++)
    frameLabels[i] = -1;
  i = 0;
  while(validFrames[i]>0)
  {
    label = phonemeStates*index[i];
    k = validFrames[i++];
    frame_length = validFrames[i] - k + 1;
    state_length = frame_length/phonemeStates;
    for(j=0; j<phonemeStates; j++)
    {
      for(n=0; n<state_length; n++)
        frameLabels[k++] = label + j;
    }
    for( ;k<=validFrames[i]; k++)
      frameLabels[k] = label + phonemeStates - 1;
    i++;
  }
     
  return self;
}


/*######################################*
 * Set the number of states per phoneme *
 * for the HMM.             *
 *######################################*/
- setPhonemeStates: (int)states
{
  phonemeStates = states;
  return self;
}

@end
