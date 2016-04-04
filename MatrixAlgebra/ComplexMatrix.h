/************************************************************************
 * This subclass of object implements a complex matrix.                 *
 *                                                                      *
 * File:ComplexMatrix.h                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/12/93  - Started                                              *
 ************************************************************************/
#import <objc/Object.h>
#import "c_headers.h"


@interface ComplexMatrix:Object

{
  int   rows, columns;          /* # of rows, cols in matrix    */
  int      *pivot;              /* pivot for LU decomp      */
  COMPLEX  *matrix;             /* matrix elements      */
  COMPLEX  *luMatrix;           /* LU decomposition     */
  COMPLEX  *xSolution;          /* Solution to Ax = b       */
}

/**********************
 * Initialization:    *
 **********************/
- initWithRows:(int)numberRows columns:(int)numberCols;
- free;
- (COMPLEX *)matrixAllocate:(int)numberRows columns:(int)numberCols;
- matrixFree:(COMPLEX *)inputMatrix;
        
/**********************
 * Calculations:      *
 **********************/
- (double)determinant;
- (BOOL)luDecomp;
- linearSolve:(COMPLEX *)b;


/**********************
 * Set parameters:    *
 **********************/
- setHilbert;
- setMatrix: (COMPLEX *)inputMatrix;

/**********************
 * Get parameters:    *
 **********************/
- (int)rows;
- (int)columns;
- (COMPLEX *)matrix;
- (COMPLEX *)luMatrix;
- (COMPLEX *)xSolution;


/***********************
 * Reading, writing    *
 * matrices.           *
 ***********************/
- printMatrix;
- printLUMatrix;


@end