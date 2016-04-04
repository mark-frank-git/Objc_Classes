/***********************************************************************
 * This subclass of object implements a Huffman code generator.        *
 * File: Huffman.m                                                     *
 *                                                                     *
 * Revision History:                                                   *
 *  1. 10/05/92  - Started                                             *
 ***********************************************************************/


#import "Huffman.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>
#import <strings.h>


#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )


@implementation Huffman

- init
{
  [super init];
  probability = 0.;
  antecedents[0] = antecedents[1] = nil;
  codeString[0] = '\0';
  level = 2;                /* binary code */
  return self;
}

/*****************************************
 * Initialize the instance variables:    *
 *****************************************/
- initWithProbability:(double)newProbability
{
  [self init];
  probability = newProbability;
  
  return self;
}


- initWithAntecedents:(id *)initAntecedents andProb:(double)newProbability
{
  [self init];
  probability = newProbability;
  antecedents[0] = initAntecedents[0];
  antecedents[1] = initAntecedents[1];
  return self;
}

/*****************************************
 * Setting instance variables:           *
 *****************************************/
- setAntecedents: (id *)newAntecedents
{
  int i;
  for(i=0; i<level; i++)
    antecedents[i] = newAntecedents[i];
  return self;
}

- setLevel: (int)newLevel
{
  level = newLevel;
  return self;
}

- setProbability: (double)newProbability
{
  probability = newProbability;
  return self;
}

- setCodeString: (char *)newCodeString
{
  strcpy(codeString, newCodeString);
  return self;
}

- setSymbolName: (char *)newSymbolName
{
  strcpy(symbolName, newSymbolName);
  return self;
}

/*****************************************
 * Getting instance variables:           *
 *****************************************/
- (int)getLevel {return level;}
- (double)getProbability {return probability;}
- (char *)getCodeString {return codeString;}
- (int )getCodeLength {return strlen(codeString);}
- (char *)getSymbolName {return symbolName;}
- (id *)getAntecedents {return antecedents;}



/*****************************************
 * Check if input symbol is equal to     *
 * local symbol name:                    *
 *****************************************/
- (BOOL)isSymbol:(char *)inputSymbol
{
  if(strcmp(inputSymbol, symbolName) == 0)
    return YES;
  return NO;
}


/*****************************************
 * Check if input code is equal to       *
 * local code string:                    *
 *****************************************/
- (BOOL)isCode:(char *)inputCode
{
  if(strcmp(inputCode, codeString) == 0)
    return YES;
  return NO;
}

@end

