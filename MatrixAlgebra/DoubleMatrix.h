/************************************************************************
 * This subclass of object implements a double precision matrix.        *
 * The matrix is represented as a double array, with the elements       *
 * increasing by rows.                                                  *
 *                                                                      *
 * File:DoubleMatrix.h                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/93  - Started                                              *
 ************************************************************************/
#import <objc/Object.h>

#define INPUT_MATRIX    0       /* defines for printing output  */
#define LU_MATRIX   1
#define INVERSE_MATRIX  2
#define RESULT_MATRIX   3

@interface DoubleMatrix:Object

{
  int      rows, columns;          /* # of rows, cols in matrix    */
  int      resultantRows, resultantCols;   /* # of rows, cols in result    */
  int     *pivot;                 /* pivot for LU decomp      */
  double  *matrix;              /* matrix elements      */
  double  *luMatrix;            /* LU decomposition     */
  double  *inverseMatrix;       /* Inverse from LU decomp   */
  double  *resultantMatrix;     /* Result of x, + operations    */
  double  *xSolution;           /* Solution to Ax = b       */
}

/**********************
 * Initialization:    *
 **********************/
- initWithRows:(int)numberRows columns:(int)numberCols;
- free;
- (double *)matrixAllocate:(int)numberRows columns:(int)numberCols;
        
/**********************
 * Calculations:      *
 **********************/
- (double)determinant;
- (BOOL)luDecomp;
- (BOOL)linearSolve:(double *)b;
- (BOOL)invert;
- addMatrix:(double *)b;
- (BOOL)multMatrix:(double *)b rows:(int)bRows columns:(int)bCols;

/**********************
 * Set parameters:    *
 **********************/
- setHilbert;
- setMatrix: (double *)inputMatrix;

/**********************
 * Get parameters:    *
 **********************/
- (int)rows;
- (int)columns;
- (double *)matrix;
- (double *)luMatrix;
- (double *)xSolution;
- (double *)inverseMatrix;
- (double *)resultantMatrix;


/***********************
 * Reading, writing    *
 * matrices.           *
 ***********************/
- printOutput:(int)type;

@end