/***********************************************************************
 * This subclass of object implements a delay line using a ring buffer *
 * File: DelayLine.h                                                   *
 *                                                                     *
 * Revision History:                                                   *
 *  1. 03/30/92  - Started                                             *
 *  2. 04/29/92  - Added zeroOutTaps                                   *
 *  3. 02/22/93  - Add delayOutput                                     *
 ***********************************************************************/

#import <objc/Object.h>

/****************************************
 * Constants: data types                *
 ****************************************/
#define CHAR_TYPE   0
#define INT_TYPE    1
#define FLOAT_TYPE  2
#define DOUBLE_TYPE 4


@interface DelayLine:Object
{
  void  **delayLine;
  int   delayTaps, delayWidth;          /* cols, rows of delay line */
  int   bufferPointer;
  int   dataType;                       /* 0 = char, 1 = int, etc. */
  void  *delayOutput;                   /* delay line output arrays*/
}

/*****************************************
 * Initialize the local  instance        *
 * variables.                            *
 *****************************************/
- free;
- init;
- initDelayLineTaps:(int)taps width:(int)width dataType:(int)type;
- zeroOutTaps;

/*****************************************
 * read and write to buffer              *
 *****************************************/
- writeDelay: (void *) input;
- (void *)readTap:  (int) tap;
- (void *)readAllTaps;
- (void *)readAllRows;

@end
