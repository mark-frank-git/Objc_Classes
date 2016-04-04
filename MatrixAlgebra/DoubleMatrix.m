/************************************************************************
 * This subclass of object implements a double precision matrix.        *
 *                                                                      *
 * File:DoubleMatrix.m                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/93  - Started                                              *
 ************************************************************************/
 
#import "DoubleMatrix.h"
#import <stdlib.h>
#import <stdio.h>
#import <string.h>


/***
 *** in-line functions for use with 2D arrays:
 ***/
#define DOFOR(i,to) for(i=0;i<to;i++)
#define DFOR(i,from,to) for(i=from-1;i<to;i++)
#define DOBY(i,from,to,by) for(i=from-1;i<to;i+=by)
#define DOBYY(i,from,to,by) for(i=from;i<t;i+=by)
#define DOBYYY(i,from,to) for(i=from;i<to;i++)
#define DOV(i,to,by) for(i=0;i<to;i+=by)
#define INDEX(i,j) [j+(i)*number_columns]
#define VECTOR(i) [i-1]


/***
 *** BLAS functions
 ***/
int isamax(int n, double *sx, int incx);
void saxpy (int n, double sa, double *sx, int incx, double *sy, int incy);
double sdot(int n, double *sx, int incx, double *sy, int incy);
void sswap(int n, double *sx, int incx, double *sy, int incy);
void sscal(int n, double sa, double *sx, int incx);
double sasum(int n, double *sx, int incx);






@implementation DoubleMatrix


/*######################*
 * Initialization:  *
 *######################*/
- initWithRows:(int)numberRows columns:(int)numberCols
{
  [self init];

  rows      = numberRows;
  columns   = numberCols;
  resultantRows = resultantCols = 0;
  matrix    = [self matrixAllocate:rows columns:columns];
  luMatrix  = [self matrixAllocate:rows columns:columns];
  inverseMatrix = [self matrixAllocate:rows columns:columns];
  resultantMatrix = NULL;
  xSolution = (double *)malloc(rows*sizeof(double));
  pivot     = (int *)malloc(rows*sizeof(int));
  
  return self;
}

/*##############################*
 * Free up the matrix:      *
 *##############################*/
- free
{
  free(matrix);
  free(luMatrix);
  free(xSolution);
  free(inverseMatrix);
  free(resultantMatrix);
  [super free];
  return self;
}

/*##############################*
 * Allocate space for a double  *
 * matrix:          *
 *##############################*/
- (double *)matrixAllocate: (int)numberRows columns:(int)numberCols
{
  double *matrix_output;
  
  matrix_output = (double *)malloc(numberRows*numberCols*sizeof(double));
  return matrix_output;
}

/*##############################*
 * Calculate the determinant.   *
 * Assumes LU factorization has *
 * been performed.              *
 * Taken from: C Tools for      *
 * Engineers and Scientists.    *
 *##############################*/
- (double)determinant
{
  int number_columns, i, sign;
  double d;

  if(rows != columns)
  {
    printf("\ndeterminant ERROR: non-square, size = %d x %d\n",rows,columns);
    return 0.;
  }

  number_columns = columns;
  sign = 0;
  d    = 1.;
  for(i=0; i<number_columns; i++)
  {
    if(pivot[i] != i)
      sign++;
    d *= luMatrix INDEX(i,i);
  }
  sign = sign - ((sign>>1)<<1);
  if(sign)
    d = -d;
  return d;
}

/*##############################*
 * Find the LU decomposition    *
 * of the complex matrix.  Con- *
 * verted from "C Tools for     *
 * Engineers and Scientists."   *
 *##############################*/
- (BOOL)luDecomp
{
  int number_columns;
  int i,j,n,k,l,kp1,nm1;
  double *a, t;

  if(rows != columns)
  {
    printf("\nLU decomp ERROR: non-square, size = %d x %d\n",rows,columns);
    return NO;
  }

  a = luMatrix;
  for(i=0; i<rows*columns; i++)
    luMatrix[i] = matrix[i];
  n = number_columns = columns;
  nm1=n-1;
  if (nm1>=1)
  {/*nontrivial pblm*/
    DOFOR(k,nm1)
    {
      kp1=k+1;
/*partial pivoting ROW exchanges-search over column*/
/* in FORTRAN, the increment would be 1 not n in izamax call*/
      pivot [k]=l=isamax((n-k),&(a INDEX(k,k)),number_columns)+k;
      if(a INDEX(l,k) != 0.)
      {/*nonsingular pivot found*/
        if(l!=k)
        {/*interchange needed*/
          t=a INDEX(l,k);
          a INDEX(l,k)=a INDEX(k,k);
          a INDEX(k,k)=t;
        }
        t=-1./a INDEX(k,k);/*scale row*/
        sscal(nm1-k,t,&(a INDEX(k+1,k)),number_columns);
        DOBYYY(j,kp1,n)
        {
          t = a INDEX(l,j);
          if(l!=k)
          {
            a INDEX(l,j)=a INDEX(k,j);
            a INDEX(k,j)=t;
          }
          saxpy(nm1-k,t,&(a INDEX(k+1,k)), n, &(a INDEX(k+1,j)), n);
        }
      }
      else /*pivot singular*/
      { 
        printf("Singular matrix found in LU decomp\n");
        return NO;
      }
    }/*main loop over k*/
  }
  pivot [nm1]=nm1;
  if (a INDEX(nm1,nm1) ==0.0)
  {
    printf("Singular matrix found in LU decomp\n");
    return NO;
  }

  return YES;
}

/*##############################*
 * Find the solution to Ax=b,   *
 * using the LU decomposition.  *
 * This method calls the luDecomp*
 * method before proceeding.    *
 * Taken from "C Tools for  *
 * Engineers and Scientists."   *
 *##############################*/
- (BOOL)linearSolve: (double *)input
{
  int number_columns;
  int i,n,k,l,nm1;
  double *a, t, *b;

  if(![self luDecomp])
    return NO;
  a = luMatrix;
  n = number_columns = columns;
  nm1 = n-1;
/* solve ly=b first*/
  for(i=0; i<rows; i++)
    xSolution[i] = input[i];
  b = xSolution;
  DOFOR(k,nm1)
  {
    l=pivot[k];
    t=b[l];
    if(l!=k)
    {
      b [l] = b [k];
      b [k] = t;
    }
    saxpy( nm1-k,t, &(a INDEX(k+1,k)),n,&(b[k+1]),1);
  }

/* solve Ux=y*/
  DOFOR(l,n)
  {
    k = nm1-l;
    b [k] = b[k]/ a INDEX(k,k);
    t     = -b[k];
    saxpy(k,t,&(a INDEX(0,k)),n,b,1);
  }

  return YES;
}

/*##############################*
 * Finds the inverse matrix *
 * using the LU decomposition.  *
 * This method calls the luDecomp*
 * method before proceeding.    *
 * Taken from "C Tools for  *
 * Engineers and Scientists."   *
 *##############################*/
- (BOOL)invert
{
  int number_columns, n;
  int i,j,k,l,kb,kp1,nm1;
  double *work, *a, t;


  if(![self luDecomp])
    return NO;
  j = rows*columns;
  for(i=0; i<j; i++)
    inverseMatrix[i] = luMatrix[i];
  a = inverseMatrix;
  number_columns = n = columns;
  nm1=n-1;
  work = (double *)malloc(number_columns*sizeof(double));
  DOFOR(k,n)
  {
        a INDEX(k,k)=t=1./ a INDEX(k,k);
        t= -t;
        sscal(k,t,&(a INDEX(0,k)),number_columns);
        kp1=k+1;
        if (nm1>=kp1)
        {
            DOBYYY(j,kp1,n)
            {
                t=a INDEX(k,j);
                a INDEX(k,j)=0.0;
                saxpy(k+1,t,&(a INDEX(0,k)),number_columns,&(a INDEX(0,j)),
                      number_columns);
            }
        }

  }
/*inv(u)*inv(l)*/
  if (nm1>=1)
  {
        DOFOR(kb,nm1)
        {
            k=nm1-kb-1;
            kp1=k+1;
            DOBYYY(i,kp1,n)
            {
                work [i]=a INDEX(i,k);
                a INDEX(i,k)=0.0;
            }
            DOBYYY(j,kp1,n)
            {
                t=work [j];
                saxpy(n,t,&(a INDEX(0,j)),number_columns,&(a INDEX(0,k))
                           ,number_columns);
            }
            l=pivot [k];
            if(l!=k) sswap(n,&(a INDEX(0,k)),number_columns,&(a INDEX(0,l))
                          ,number_columns);
        }
  }
  free(work);
  return YES;
}

/*##############################*
 * Add the input matrix to      *
 * matrix -> resultantMatrix.   *
 *##############################*/
- addMatrix: (double *)b
{
  int i, size;

  if(resultantMatrix != NULL)
    free(resultantMatrix);
  resultantMatrix = [self matrixAllocate:rows columns:columns];
  resultantRows   = rows;
  resultantCols   = columns;
  size = rows*columns;
  for(i=0; i<size; i++)
    resultantMatrix[i] = matrix[i] + b[i];
  return self;
}


/*##############################*
 * Post multiply the input  *
 * matrix by matrix ->      *
 * resultantMatrix.     *
 *##############################*/
- (BOOL)multMatrix: (double *)b rows:(int)bRows columns:(int)bCols;
{
  int i, j, k, number_columns;
  double *a_row_ptr;

  if(bRows != columns)
  {
    printf("Bad multiplicand in multMatrix:\n");
    return NO;
  }
  if(resultantMatrix != NULL)
    free(resultantMatrix);
  resultantMatrix = [self matrixAllocate:rows columns:bCols];
  resultantRows   = rows;
  resultantCols   = bCols;
  number_columns = bCols;
  for(i=0; i<rows; i++)
  {
    a_row_ptr = &matrix[i*columns];
    for(j=0; j<bCols; j++)
    {
      resultantMatrix INDEX(i,j) = 0.;
      for(k=0; k<columns; k++)
        resultantMatrix INDEX(i,j) += a_row_ptr[k]*b INDEX(k,j);
    }
  }

  return YES;
}


/*##############################*
 * Set up a Hilbert matrix  *
 *##############################*/
- setHilbert
{
  int i, j, number_columns;

  if(rows != columns)
  {
    printf("\nsetHilbert ERROR: non-square, size = %d x %d\n",rows,columns);
    return 0;
  }

  number_columns = columns;
  for(i=0; i<rows; i++)
  {
    for(j=0; j<columns; j++)
      matrix INDEX(i,j) = 1.0/(i+j+1.);
  }
  return self;
}

/*############################*
 * Set matrix equal to input    *
 * matrix.          *
 *##############################*/
- setMatrix:(double *)inputMatrix
{
  int i;
  for(i=0; i<rows*columns; i++)
    matrix[i] = inputMatrix[i];
  return self;
}

/*##############################*
 * Return matrix parameters *
 *##############################*/
- (int)rows { return rows; }
- (int)columns { return columns; }
- (double *) matrix {return matrix;}
- (double *) luMatrix {return luMatrix;}
- (double *) xSolution {return xSolution;}
- (double *) inverseMatrix {return inverseMatrix;}
- (double *) resultantMatrix {return resultantMatrix;}

/*##############################*
 * Print out the matrix:    *
 *##############################*/
- printOutput:(int)type
{
  int i,j,cols_lin, matrix_rows, matrix_cols, number_columns;
  double *matrix_ptr;

  matrix_rows = rows;
  matrix_cols = columns;
  switch(type)
  {
    default:
    case INPUT_MATRIX:
      matrix_ptr = matrix;
      break;
    case LU_MATRIX:
      matrix_ptr = luMatrix;
      break;
    case INVERSE_MATRIX:
      matrix_ptr = inverseMatrix;
      break;
    case RESULT_MATRIX:
      matrix_ptr = resultantMatrix;
      matrix_rows = resultantRows;
      matrix_cols = resultantCols;
      break;
  }
  if(matrix_ptr == NULL)
    return self;

  number_columns = matrix_cols;
  cols_lin = 5;        /* 5 per line for double types */

  for(i = 0 ; i < matrix_rows ; i++)
  {
    for(j = 0 ; j<matrix_cols ; j++)
    {
      if(j%cols_lin == 0)
      {    /* newline every cols_lin */
        if(j == 0)           /* start of row */
          printf("\nRow%3d:",i);
        else
          printf("\n       ");
       }
       printf("%14.5g",matrix_ptr INDEX(i,j));
    }
  }
  printf("\n");
  return self;
}


@end