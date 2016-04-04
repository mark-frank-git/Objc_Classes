/*********************
 * Include files:    *
 *********************/
#import "DesStat.h"
#import <stdio.h>
#import <math.h>
#import <stdlib.h>

#define MAXIMUM 10000
#define DIMENSION 3     /* if changed then modify a[], b[] */
/*********************
 * Constants:        *
 *********************/


void main()
{
   int      i,j,k,max=MAXIMUM,dim=DIMENSION;
   float    a[] = {0,0,0,0},
            b[] = {2,4,8,16};
   long double      data[DIMENSION+1],*ld;
   float            *f;
   id     desStat;              /* this is really a pointer! */
 
 
/*********************
 * Init statistics   *
 * object:           *
 *********************/
   desStat = [[DesStat alloc] initWithDimension:dim];
 
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

/*********************
 * Add data          *
 *********************/
   printf("\nNow we add %d data vectors.",max);
   for( i=0; i<max; i++ )
   {  for( j=0; j<dim; j++ )
         data[j] = (b[j]-a[j])*rand()/RAND_MAX + a[j];
      [desStat addData: data];
   }
   printf("\n     New data added.");
   
/*********************
 * print parameters  *
 *********************/
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

   ld = [desStat getVecSum];
   printf("\nThe unnormalized sums are:\n   ");
   for( i=0; i<dim; i++ )
      printf("%10.3Le",ld[i]);

   ld = [desStat getArraySum];
   printf("\nThe unnormalized sums of products are:");
   k = 0;
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for(j=0; j<i; j++)
         printf("%10.3Le",0.0);
      for( j=i; j<dim; j++ )
         printf("%10.3Le",ld[k++]);
   }

   f = [desStat getMean];
   printf("\nThe mean is:\n   ");
   for( i=0; i<dim; i++ )
      printf("%8.3f",f[i]);
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

   f = [desStat getStdDev];
   printf("\nThe standard deviation is:\n   ");
   for( i=0; i<dim; i++ )
      printf("%8.3f",f[i]);
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

   f = [desStat getCovar];
   printf("\nThe covariance matrix is:");
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for( j=0; j<dim; j++ )
         printf("%8.3f",f[i*dim+j]);
   }
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);
   printf("\n");

   f = [desStat getNormCovar];
   printf("\nThe normalized covariance matrix is:");
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for( j=0; j<dim; j++ )
         printf("%8.3f",f[i*dim+j]);
   }
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);
   printf("\n");
   
/*********************
 * reset dimension   *
 *********************/
   printf("\n*********** redo with dim = DIMENSION + 1 ***************\n");
   dim++;
   [desStat resetWithDimension:dim];
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

/*********************
 * Add data          *
 *********************/
   printf("\nNow we add %d data vectors.",max);
   for( i=0; i<max; i++ )
   {  for( j=0; j<dim; j++ )
         data[j] = (b[j]-a[j])*rand()/RAND_MAX + a[j];
      [desStat addData: data];
   }
   printf("\n     New data added.");
   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);
   
/*********************
 * print parameters  *
 *********************/
   ld = [desStat getVecSum];
   printf("\nThe unnormalized sums are:\n   ");
   for( i=0; i<dim; i++ )
      printf("%10.3Le",ld[i]);

   ld = [desStat getArraySum];
   printf("\nThe unnormalized sums of products are:");
   k = 0;
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for(j=0; j<i; j++)
         printf("%10.3Le",0.0);
      for( j=i; j<dim; j++ )
         printf("%10.3Le",ld[k++]);
   }

   f = [desStat getMean];
   printf("\nThe mean is:\n   ");
   for( i=0; i<dim; i++ )
      printf("%8.3f",f[i]);

   f = [desStat getStdDev];
   printf("\nThe standard deviation is:\n   ");
   for( i=0; i<dim; i++ )
      printf("%8.3f",f[i]);

   f = [desStat getCovar];
   printf("\nThe covariance matrix is:");
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for( j=0; j<dim; j++ )
         printf("%8.3f",f[i*dim+j]);
   }
   printf("\n");

   f = [desStat getNormCovar];
   printf("\nThe normalized covariance matrix is:");
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for( j=0; j<dim; j++ )
         printf("%8.3f",f[i*dim+j]);
   }
   
/*********************
 * reset             *
 *********************/
   printf("\n\n*********** reset ***************\n");
   [desStat reset];

   printf("\ndimension = %d, count = %ld, count at last mean = %ld",
              [desStat getDim],[desStat getCount],
              [desStat getCountAtMean]);
   printf("\n     count at last std dev = %ld, count at last covar = %ld",
              [desStat getCountAtStdDev], [desStat getCountAtCovar]);
   printf("\n     count at last norm covar = %ld.",
              [desStat getCountAtNormCovar]);

   ld = [desStat getVecSum];
   printf("\nThe unnormalized sums are:\n   ");
   for( i=0; i<dim; i++ )
      printf("%10.3Le",ld[i]);

   ld = [desStat getArraySum];
   printf("\nThe unnormalized sums of products are:");
   k = 0;
   for( i=0; i<dim; i++ )
   {  printf("\n   ");
      for(j=0; j<i; j++)
         printf("%10.3Le",0.0);
      for( j=i; j<dim; j++ )
         printf("%10.3Le",ld[k++]);
   }

   printf("\n");

}
