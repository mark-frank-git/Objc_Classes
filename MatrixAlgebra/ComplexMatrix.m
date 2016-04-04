/************************************************************************
 * This subclass of object implements a complex matrix.         *
 *                                  *
 * File:ComplexMatrix.m                         *
 *                                  *
 * Revision history:                            *
 *  1. 04/12/93  - Started                      *
 ************************************************************************/
 
#import "ComplexMatrix.h"
#import <stdlib.h>
#import <stdio.h>
#import <string.h>

   
@implementation ComplexMatrix


/*######################*
 * Initialization:  *
 *######################*/
- initWithRows:(int)numberRows columns:(int)numberCols
{
  [self init];

  rows      = numberRows;
  columns   = numberCols;
  matrix    = [self matrixAllocate:rows columns:columns];
  luMatrix  = [self matrixAllocate:rows columns:columns];
  xSolution = (COMPLEX *)malloc(rows*sizeof(COMPLEX));
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
  [super free];
  return self;
}

/*##############################*
 * Allocate space for a COMPLEX *
 * matrix:          *
 *##############################*/
- (COMPLEX *)matrixAllocate: (int)numberRows columns:(int)numberCols
{
  COMPLEX *matrix_output;
  
  matrix_output = (COMPLEX *)malloc(numberRows*numberCols*sizeof(COMPLEX));
  return matrix_output;
}

/*##############################*
 * Free up space for a double   *
 * matrix:          *
 *##############################*/
- matrixFree: (COMPLEX *)inputMatrix
{
  free(inputMatrix);
  return self;
}

/*##############################*
 * Calculate the determinant:   *
 * Taken from: C Language Al-   *
 * gorithms for Signal Proces-  *
 * sing.            *
 *##############################*/
- (double)determinant
{
  return(0.);
}

/*##############################*
 * Find the LU decomposition    *
 * of the complex matrix.  Con- *
 * verted from FORTRAN subrou-  *
 * tine, zgefa in LINPACK, and  *
 * LU decomp in "C Tools for    *
 * Engineers and Scientists."   *
 *##############################*/
- (BOOL)luDecomp
{
  int number_columns;
  int i,j,n,k,l,kp1,nm1;
  COMPLEX *a, t, ctemp;
  double dtemp;

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
      pivot [k]=l=izamax((n-k),&(a INDEX(k,k)),number_columns)+k;
      ctemp = a INDEX(l,k);
      dtemp = CABS2(ctemp);
      if(dtemp != 0.)
      {/*nonsingular pivot found*/
        if(l!=k)
        {/*interchange needed*/
          t=a INDEX(l,k);
          a INDEX(l,k)=a INDEX(k,k);
          a INDEX(k,k)=t;
        }
        CMPLX(ctemp, -1., 0.);
        CDIV(t, ctemp, a INDEX(k,k));/*scale row*/
        zscal(nm1-k,t,&(a INDEX(k+1,k)),number_columns);
        DOBYYY(j,kp1,n)
        {
          t = a INDEX(l,j);
          if(l!=k)
          {
            a INDEX(l,j)=a INDEX(k,j);
            a INDEX(k,j)=t;
          }
          zaxpy(nm1-k,t,&(a INDEX(k+1,k)), n, &(a INDEX(k+1,j)), n);
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
  ctemp = a INDEX(nm1,nm1);
  if (CABS2(ctemp) ==0.0)
  {
    printf("Singular matrix found in LU decomp\n");
    return NO;
  }

  return YES;
}

/*##############################*
 * Find the solution to Ax=b,   *
 * using the LU decomposition.  *
 * Taken from "C Tools for  *
 * Engineers and Scientists."   *
 *##############################*/
- linearSolve: (COMPLEX *)input
{
  int number_columns;
  int i,n,k,l,nm1;
  COMPLEX *a, t, ctemp, *b;

  if(![self luDecomp])
    return self;
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
    zaxpy( nm1-k,t, &(a INDEX(k+1,k)),n,&(b[k+1]),1);
  }

/* solve Ux=y*/
  DOFOR(l,n)
  {
    k = nm1-l;
    CDIV(ctemp, b[k], a INDEX(k,k));
    b [k] = ctemp;
    t.x = -b[k].x;
    t.y = -b[k].y;
    zaxpy(k,t,&(a INDEX(0,k)),n,b,1);
  }

  return self;
}


/*##############################*
 * Set up a Hilbert matrix  *
 *##############################*/
- setHilbert
{
  int i, j;
  int number_columns;


  if(rows != columns)
  {
    printf("\nsetHilbert ERROR: non-square, size = %d x %d\n",rows,columns);
    return 0;
  }
  number_columns = columns;
  for(i=0; i<rows; i++)
  {
    for(j=0; j<columns; j++)
    {
      matrix INDEX(i,j).x = 1.0/(i+j+1.);
      matrix INDEX(i,j).y = 0.;
    }

  }
  return self;
}

/*############################*
 * Set matrix equal to input    *
 * matrix.          *
 *##############################*/
- setMatrix:(COMPLEX *)inputMatrix
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
- (COMPLEX *) matrix {return matrix;}
- (COMPLEX *) luMatrix {return luMatrix;}
- (COMPLEX *) xSolution {return xSolution;}

/*##############################*
 * Print out the matrix:    *
 *##############################*/
- printMatrix
{
  int i,j,cols_lin;
  int number_columns;

  number_columns = columns;
  cols_lin = 5;        /* 5 per line for float types */

  for(i = 0 ; i < rows ; i++)
  {
    for(j = 0 ; j<columns ; j++)
    {
      if(j%cols_lin == 0)
      {    /* newline every cols_lin */
        if(j == 0)           /* start of row */
          printf("\nRow%3d:",i);
        else
          printf("\n       ");
       }
       printf("(%7.3g, %7.3g)",matrix INDEX(i,j).x,matrix INDEX(i,j).y);
    }
  }
  printf("\n");
  return self;
}


/*##############################*
 * Print out the LU matrix: *
 *##############################*/
- printLUMatrix
{
  int i,j,cols_lin;
  int number_columns;

  number_columns = columns;
  cols_lin = 5;        /* 5 per line for float types */

  for(i = 0 ; i < rows ; i++)
  {
    for(j = 0 ; j<columns ; j++)
    {
      if(j%cols_lin == 0)
      {    /* newline every cols_lin */
        if(j == 0)           /* start of row */
          printf("\nRow%3d:",i);
        else
          printf("\n       ");
       }
       printf("(%7.3g, %7.3g)",luMatrix INDEX(i,j).x,luMatrix INDEX(i,j).y);
    }
  }
  printf("\n");
  return self;
}

@end