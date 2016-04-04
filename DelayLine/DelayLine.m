/***********************************************************************
 * This subclass of object implements a delay line using a ring buffer *
 * File: DelayLine.m                                                   *
 *                                                                     *
 * Revision History:                                                   *
 *  1. 03/30/92  - Started                                             *
 *  2. 06/01/92  - Free output arrays on next call                     *
 *  3. 02/22/93  - Add delayOutput                                     *
 ***********************************************************************/

#import "DelayLine.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>


#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )


@implementation DelayLine

/****************************************
 * Free up the pointers:                *
 ****************************************/
- free
{
  int i;
  
  if(delayLine != NULL)
  {
    for(i=0; i<delayWidth; i++)
      free(delayLine[i]);

    free(delayLine);
  }
  free(delayOutput);
  [super free];
  return self;
}

/****************************************
 * Initialize instance variables        *
 ****************************************/
- init
{
  [super init];
  delayLine     = NULL;
  bufferPointer = delayTaps = 0;
  return self;
}

/*****************************************
 * Initialize the buffers                *
 *****************************************/
- initDelayLineTaps:(int)taps width:(int)width dataType:(int)type
{
  int     i, j, size;
  char    **char_ptr;
  int     **int_ptr;
  float   **float_ptr;
  double  **double_ptr;
  
  [self init];
  delayTaps   = MAX(taps, 1);
  delayWidth  = MAX(width, 1);
  dataType    = type;
  size        = delayTaps*delayWidth;
  
  switch(dataType)
  {
    case CHAR_TYPE:
      delayOutput = (void *)calloc(size, sizeof(char));
      char_ptr = (char **)calloc(delayWidth, sizeof(char *));
      for(i=0; i<delayWidth; i++)
      {
        char_ptr[i] = (char *)calloc(delayTaps, sizeof(char));
    for(j=0; j<delayTaps; j++)
      char_ptr[i][j] = 0;
      }
      delayLine = (void **)char_ptr;
      break;
    case INT_TYPE:
      delayOutput = (void *)calloc(size, sizeof(int));
      int_ptr = (int **)calloc(delayWidth, sizeof(int *));
      for(i=0; i<delayWidth; i++)
      {
        int_ptr[i] = (int *)calloc(delayTaps, sizeof(int));
    for(j=0; j<delayTaps; j++)
      int_ptr[i][j] = 0.;
      }
      delayLine = (void **)int_ptr;
      break;
    case FLOAT_TYPE:
      delayOutput = (void *)calloc(size, sizeof(float));
      float_ptr = (float **)calloc(delayWidth, sizeof(float *));
      for(i=0; i<delayWidth; i++)
      {
        float_ptr[i] = (float *)calloc(delayTaps, sizeof(float));
    for(j=0; j<delayTaps; j++)
      float_ptr[i][j] = 0.;
      }
      delayLine = (void **)float_ptr;
      break;
    case DOUBLE_TYPE:
      delayOutput = (void *)calloc(size, sizeof(double));
      double_ptr = (double **)calloc(delayWidth, sizeof(double *));
      for(i=0; i<delayWidth; i++)
      {
        double_ptr[i] = (double *)calloc(delayTaps, sizeof(double));
    for(j=0; j<delayTaps; j++)
      double_ptr[i][j] = 0.;
      }
      delayLine = (void **)double_ptr;
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  
  return self;
}

/*****************************************
 * Zero out the delay line:              *
 *****************************************/
- zeroOutTaps
{
  int     i, j;
  char    **char_ptr;
  int     **int_ptr;
  float   **float_ptr;
  double  **double_ptr;
  
  switch(dataType)
  {
    case CHAR_TYPE:
      char_ptr = (char **)delayLine;
      for(i=0; i<delayWidth; i++)
      {
    for(j=0; j<delayTaps; j++)
      char_ptr[i][j] = 0;
      }
      break;
    case INT_TYPE:
      int_ptr = (int **)delayLine;
      for(i=0; i<delayWidth; i++)
      {
    for(j=0; j<delayTaps; j++)
      int_ptr[i][j] = 0.;
      }
      break;
    case FLOAT_TYPE:
      float_ptr = (float **)delayLine;
      for(i=0; i<delayWidth; i++)
      {
    for(j=0; j<delayTaps; j++)
      float_ptr[i][j] = 0.;
      }
      delayLine = (void **)float_ptr;
      break;
    case DOUBLE_TYPE:
      double_ptr = (double **)delayLine;
      for(i=0; i<delayWidth; i++)
      {
    for(j=0; j<delayTaps; j++)
      double_ptr[i][j] = 0.;
      }
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  
  return self;
}

/*****************************************
 * New input to delay line               *
 *****************************************/
- writeDelay: (void *)input
{
  int     i;
  char    **char_delay_ptr,   *char_input_ptr;
  int     **int_delay_ptr,    *int_input_ptr;
  float   **float_delay_ptr,  *float_input_ptr;
  double  **double_delay_ptr, *double_input_ptr;

  switch(dataType)
  {
    case CHAR_TYPE:
      char_input_ptr = (char *)input;
      char_delay_ptr = (char **)delayLine;
      for(i=0; i<delayWidth; i++)
        char_delay_ptr[i][bufferPointer] = char_input_ptr[i];
      break;
    case INT_TYPE:
      int_input_ptr = (int *)input;
      int_delay_ptr = (int **)delayLine;
      for(i=0; i<delayWidth; i++)
        int_delay_ptr[i][bufferPointer] = int_input_ptr[i];
      break;
    case FLOAT_TYPE:
      float_input_ptr = (float *)input;
      float_delay_ptr = (float **)delayLine;
      for(i=0; i<delayWidth; i++)
        float_delay_ptr[i][bufferPointer] = float_input_ptr[i];
      break;
    case DOUBLE_TYPE:
      double_input_ptr = (double *)input;
      double_delay_ptr = (double **)delayLine;
      for(i=0; i<delayWidth; i++)
        double_delay_ptr[i][bufferPointer] = double_input_ptr[i];
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  bufferPointer++;
  bufferPointer %= delayTaps;
  
  return self;
}

/*****************************************
 * Read from delay line, get a tap       *
 * at a time:                            *
 * NOTE: tap 0 = input just written      *
 *****************************************/
- (void *)readTap: (int)tap
{
  int     i;
  int     buffer_pointer;
  char    **char_delay_ptr,   *char_output_ptr;
  int     **int_delay_ptr,    *int_output_ptr;
  float   **float_delay_ptr,  *float_output_ptr;
  double  **double_delay_ptr, *double_output_ptr;

  buffer_pointer = bufferPointer-tap-1;
  while(buffer_pointer<0)
    buffer_pointer += delayTaps;
  buffer_pointer %= delayTaps;
  
  switch(dataType)
  {
    case CHAR_TYPE:
      char_output_ptr = (char *)delayOutput;
      char_delay_ptr = (char **)delayLine;
      for(i=0; i<delayWidth; i++)
        char_output_ptr[i] = char_delay_ptr[i][buffer_pointer];
      break;
    case INT_TYPE:
      int_output_ptr = (int *)delayOutput;
      int_delay_ptr = (int **)delayLine;
      for(i=0; i<delayWidth; i++)
        int_output_ptr[i] = int_delay_ptr[i][buffer_pointer];
      break;
    case FLOAT_TYPE:
      float_output_ptr = (float *)delayOutput;
      float_delay_ptr = (float **)delayLine;
      for(i=0; i<delayWidth; i++)
        float_output_ptr[i] = float_delay_ptr[i][buffer_pointer];
      break;
    case DOUBLE_TYPE:
      double_output_ptr = (double *)delayOutput;
      double_delay_ptr = (double **)delayLine;
      for(i=0; i<delayWidth; i++)
        double_output_ptr[i] = double_delay_ptr[i][buffer_pointer];
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  
  return delayOutput;
}

/*****************************************
 * Read from delay line, get all rows a  *
 * row at a time, start at the oldest in-*
 * put, and move forward:                *
 *****************************************/
- (void *)readAllRows
{
  int     i, j, k;
  int     buffer_pointer, size;
  char    **char_delay_ptr,   *char_output_ptr;
  int     **int_delay_ptr,    *int_output_ptr;
  float   **float_delay_ptr,  *float_output_ptr;
  double  **double_delay_ptr, *double_output_ptr;

  buffer_pointer = bufferPointer;
  size = delayWidth*delayTaps;
  switch(dataType)
  {
    case CHAR_TYPE:
      char_output_ptr = (char *)delayOutput;
      char_delay_ptr = (char **)delayLine;
      for(j=0; j<delayTaps; j++)
      {
        k = j;
        for(i=0; i<delayWidth; i++)
    {
          char_output_ptr[k] = char_delay_ptr[i][buffer_pointer];
      k += delayTaps;
    }
        buffer_pointer++;
    buffer_pointer %= delayTaps;
      }
      break;
    case INT_TYPE:
      int_output_ptr = (int *)delayOutput;
      int_delay_ptr = (int **)delayLine;
      for(j=0; j<delayTaps; j++)
      {
        k = j;
        for(i=0; i<delayWidth; i++)
    {
          int_output_ptr[k] = int_delay_ptr[i][buffer_pointer];
      k += delayTaps;
    }
        buffer_pointer++;
    buffer_pointer %= delayTaps;
      }
      break;
    case FLOAT_TYPE:
      float_output_ptr = (float *)delayOutput;
      float_delay_ptr = (float **)delayLine;
      for(j=0; j<delayTaps; j++)
      {
        k = j;
        for(i=0; i<delayWidth; i++)
    {
          float_output_ptr[k] = float_delay_ptr[i][buffer_pointer];
      k += delayTaps;
    }
        buffer_pointer++;
    buffer_pointer %= delayTaps;
      }
      break;
    case DOUBLE_TYPE:
      double_output_ptr = (double *)delayOutput;
      double_delay_ptr = (double **)delayLine;
      for(j=0; j<delayTaps; j++)
      {
        k = j;
        for(i=0; i<delayWidth; i++)
    {
          double_output_ptr[k] = double_delay_ptr[i][buffer_pointer];
      k += delayTaps;
    }
        buffer_pointer++;
    buffer_pointer %= delayTaps;
      }
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  
  return delayOutput;
}



/*****************************************
 * Read from delay line, get all taps a  *
 * tap at a time, start at the input     *
 * tap, and move forward:                *
 *****************************************/
- (void *)readAllTaps
{
  int     i, j, k;
  int     buffer_pointer, size;
  char    **char_delay_ptr,   *char_output_ptr;
  int     **int_delay_ptr,    *int_output_ptr;
  float   **float_delay_ptr,  *float_output_ptr;
  double  **double_delay_ptr, *double_output_ptr;

  buffer_pointer = bufferPointer-1;
  while(buffer_pointer<0)
    buffer_pointer += delayTaps;
  buffer_pointer %= delayTaps;
  size = delayWidth*delayTaps;
  switch(dataType)
  {
    case CHAR_TYPE:
      char_output_ptr = (char *)delayOutput;
      char_delay_ptr = (char **)delayLine;
      k = 0;
      for(j=0; j<delayTaps; j++)
      {
        for(i=0; i<delayWidth; i++)
          char_output_ptr[k++] = char_delay_ptr[i][buffer_pointer];
        buffer_pointer--;
        while(buffer_pointer<0)
          buffer_pointer += delayTaps;
    buffer_pointer %= delayTaps;
      }
      break;
    case INT_TYPE:
      int_output_ptr = (int *)delayOutput;
      int_delay_ptr = (int **)delayLine;
      k = 0;
      for(j=0; j<delayTaps; j++)
      {
        for(i=0; i<delayWidth; i++)
          int_output_ptr[k++] = int_delay_ptr[i][buffer_pointer];
        buffer_pointer--;
        while(buffer_pointer<0)
          buffer_pointer += delayTaps;
    buffer_pointer %= delayTaps;
      }
      break;
    case FLOAT_TYPE:
      float_output_ptr = (float *)delayOutput;
      float_delay_ptr = (float **)delayLine;
      k = 0;
      for(j=0; j<delayTaps; j++)
      {
        for(i=0; i<delayWidth; i++)
          float_output_ptr[k++] = float_delay_ptr[i][buffer_pointer];
        buffer_pointer--;
        while(buffer_pointer<0)
          buffer_pointer += delayTaps;
    buffer_pointer %= delayTaps;
      }
      break;
    case DOUBLE_TYPE:
      double_output_ptr = (double *)delayOutput;
      double_delay_ptr = (double **)delayLine;
      k = 0;
      for(j=0; j<delayTaps; j++)
      {
        for(i=0; i<delayWidth; i++)
          double_output_ptr[k++] = double_delay_ptr[i][buffer_pointer];
        buffer_pointer--;
        while(buffer_pointer<0)
          buffer_pointer += delayTaps;
    buffer_pointer %= delayTaps;
      }
      break;
    default:
      printf("Unknown data type in DelayLine\n");
      break;;
  }
  
  return delayOutput;
}


@end

