/***********************************************************************
 * This subclass of object implements a Huffman code generator.        *
 * File: Huffman.h                                                     *
 *                                                                     *
 * Revision History:                                                   *
 *  1. 10/05/92  - Started                                             *
 ***********************************************************************/

#import <objc/Object.h>

/****************************************
 * Constants: data types                *
 ****************************************/
#define MAX_CODE    30          /* max length */
#define MAX_SYMBOL  10          /* max symbol name */
#define MAX_LEVEL   10          /* 10-ary symbols max */


@interface Huffman:Object
{
  id    antecedents[MAX_LEVEL];         /* predecessors */
  int   level;                  /* binary, ternary, etc.*/
  double probability;               /* probability of this word */
  char  codeString[MAX_CODE];
  char  symbolName[MAX_SYMBOL];
}

/*****************************************
 * Initialize the local  instance        *
 * variables.                            *
 *****************************************/
- initWithProbability:(double)initProbability;
- initWithAntecedents:(id *)initAntecedents andProb:(double)newProbability;

/*****************************************
 * Set instance variables:               *
 *****************************************/
- setAntecedents: (id *)newAntecedents;
- setLevel: (int)newLevel;
- setProbability: (double)newProbability;
- setCodeString: (char *)newCodeString;
- setSymbolName: (char *)newSymbolName;

/*****************************************
 * Get instance variables:               *
 *****************************************/
- (int)getLevel;
- (double)getProbability;
- (char *)getCodeString;
- (int )getCodeLength;
- (char *)getSymbolName;
- (id *)getAntecedents;

/*****************************************
 * Checking instance variables:          *
 *****************************************/
- (BOOL)isSymbol:(char *)inputSymbol;
- (BOOL)isCode:  (char *)inputCode;

@end
