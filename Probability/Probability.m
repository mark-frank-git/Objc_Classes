/****************************************************************************
 * This subclass of object implements a general probability distribution    *
 * function.                                                                *
 *                                                                          *
 * File: /User/frank/Objc_Classes/Probability/Probability.m                 *
 *                                                                          *
 * Revision History:                                                        *
 *  1. 05/04/92 - Started                                                   *
 *  2. 06/10/92 - Free arrays in here                                       *
 *  3. 07/09/92 - Add gamma, log_gamma, log_normal, and correlation         *
 *                coefficient per Bill Kushner.                             *
 *  4. 11/19/92 - Add GAUSS_MARKOV, distributionFunction -> densityFunction *
 *  5. 02/10/94 - Add Rician, Rayleigh, noncentral chi squared              *
 ****************************************************************************/

#import "Probability.h"
#import "Random.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#define MIN_RHO         -0.99999999
#define MAX_RHO         0.99999999
#define MIN_VARIANCE    1.e-40
#define SEED1           2723
#define SEED2           177
#define SEED3           1139

@implementation Probability

/*###############################*
 * Initialization:               *
 *###############################*/
- init
{
  
  [super init];

// 
// Initialize instance variables:
//
  randomObject = [[Random alloc] initSeeds: SEED1: SEED2: SEED3];
  dimensions    = 1;
  mean          = variance = stdDev = NULL;
  distribution  = NORMAL;
  probability   = 1.;
  rho           = 0.;
  noncentralParameter   = 100.;
  degrees       = 100;
  doubleOutputs = NULL;
  [self initMeanVariance];
  
  return self;
}

/*################################*
 * Initialize means and variances *
 *################################*/
- initMeanVariance
{
  int i;
  if(mean != NULL)
  {
    free(mean);
    free(variance);
    free(stdDev);
  }
  mean          = (double *)malloc(dimensions*sizeof(double));
  variance      = (double *)malloc(dimensions*sizeof(double));
  stdDev    = (double *)malloc(dimensions*sizeof(double));
  for(i=0; i<dimensions; i++)
  {
    mean[i]     = 0.;
    variance[i] = stdDev[i] = 1.;
  }
  return self;
}

/*##############################*
 * Re-initialize seeds to       *
 * default values:              *
 *##############################*/
- initializeSeeds
{
  [randomObject setSeeds: SEED1: SEED2: SEED3];
  return self;
}

/*##############################*
 * These methods set parameters *
 *##############################*/
- setMean:(double *)newMean
{
  int i;
  for(i=0; i<dimensions; i++)
    mean[i] = newMean[i];
  return self;
}

- setVariance: (double *)newVariance
{
  int i;
  for(i=0; i<dimensions; i++)
  {
    variance[i] = MAX(MIN_VARIANCE, newVariance[i]);
    stdDev[i]   = sqrt(variance[i]);
  }
  return self;
}

- setCorrelationCoefficient: (double)newRho
{
  rho = MAX(MIN_RHO,newRho);
  rho = MIN(MAX_RHO,rho);
  return self;
}

- setProbDistribution: (int)newDistribution
{
  distribution = newDistribution;
  return self;
}

- setProbability: (double)newProbability
{
  probability = newProbability;
  return self;
}

- setDimensions: (int)newDimensions
{
  dimensions = newDimensions;
  [self initMeanVariance];
  return self;
}

- setNoncentralParameter:(double)parameter
{
  noncentralParameter = parameter;
  return self;
}
- setDegreesFreedom:(int)degree
{
  degrees = degree;
  return self;
}
/*##############################*
 * These methods get parameters *
 *##############################*/
- (double *)getMean{        return mean;}
- (double *)getVariance{    return variance;}
- (int)getProbDistribution{ return distribution;}
- (double)getProbability{   return probability;}
- (int) getDimensions{      return dimensions;}
- (int) degrees{            return degrees;}

/*###############################*
 * Return a new sample from our  *
 * distribution:                 *
 *###############################*/
- (double *)newSample
{
  int    i, j, k;
  double *sample, range;
  double u1,s,ln_s;
  double v1,v1_sq,v2_sq, alpha, a, tr;
  double phi, mu, sigma_sq;
  double n1, n2;
  static double v2, sqrt_lns, temp;
  static int odd = 1;
  static double markov_old = 0.;
  
  if(doubleOutputs != NULL)
    free(doubleOutputs);
  sample = doubleOutputs = (double *)malloc(dimensions*sizeof(double));
  switch(distribution)
  {
    case UNIFORM:
      for(i=0; i<dimensions; i++)
      {
        sample[i]  = [randomObject percent];        /* sample in [0,1] */
        sample[i] -= 0.5;
        range      = SQRT3*stdDev[i];
        sample[i] *= range + range;
        sample[i] += mean[i];
      }
      break;
    case NORMAL:
    default:
      for(i=0; i<dimensions; i++)
      {
        sample[i] = [self normal01];
        if(i==1)
        {
          sample[i] = stdDev[i]*(rho*temp + sample[i]*sqrt(1.-rho*rho));
          sample[i] += mean[i];
        }
        else
        {
         temp       = sample[i];
         sample[i] *= stdDev[i];
         sample[i] += mean[i];
        }
      }
      break;
    case EXPONENTIAL_DIST:
      for(i=0; i<dimensions; i++)
      {
        alpha   = mean[i] - stdDev[i];
        u1 = [randomObject percent];
        if(u1>0.)
          sample[i] = -stdDev[i]*log(u1);
        else
          sample[i] = 0.;
        sample[i] += alpha;
      }
      break;
    case ERLANG_DIST:
      for(i=0; i<dimensions; i++)
      {
        a   = mean[i]/variance[i];
        k   = ROUND(mean[i]*a);
        tr  = 1.;
        for(j=0; j<k; j++)
          tr  *= [randomObject percent];
        if(tr>0.)
          sample[i] = -log(tr)/a;
        else
          sample[i] = 0.;
      }
      break;
    case LOG_NORMAL:
/***************** Kushner's old stuff *******************
      for(i=0; i<dimensions; i++)
      {
        std_devx = stdDev[i];
        term    = log(variance[i]/mean[i]/mean[i] + 1.);
        std_devy = sqrt(term);
        if(mean[i]>0.)
          mean_y   = log(mean[i]) - 0.5*term;
        else
          mean_y = 0.;
        sum = -6.0;
        for(j=0; j<12; j++)
          sum += [randomObject percent];
        sample[i] = exp(mean_y + std_devy*sum);
      }
**********************************************************/
      for(i=0; i<dimensions; i++)
      {
        sample[i] = [self normal01];
        if(mean[i]>0.)
        {
          sigma_sq = log(variance[i]/mean[i]*mean[i] + 1.);
          mu       = log(mean[i]) - sigma_sq/2;
        }
        else mu = sigma_sq = 0.;
        sample[i] = exp(mu + sqrt(sigma_sq)*sample[i]);
      }
      break;
    case GAUSS_MARKOV:
      for(i=0; i<dimensions; i++)
      {
        if(odd)
        {
          s = 2.;
          while(s > 1.)
          {
            u1 = [randomObject percent];
            v2 = u1+u1-1.;
            v2_sq = v2*v2;
            u1 = [randomObject percent];
            v1 = u1 + u1 - 1.;
            v1_sq = v1*v1;
            s  = v1_sq + v2_sq;
          }
          ln_s = log(s);
          sqrt_lns = sqrt(-(ln_s+ln_s)/s);
          odd = 0;
          sample[i]  = rho*markov_old + sqrt(1.-rho*rho)*v1*sqrt_lns;
          markov_old = sample[i];
        }
        else
        {
          odd = 1;
          sample[i] = rho*markov_old + sqrt(1.-rho*rho)*v2*sqrt_lns;
          markov_old = sample[i];
        }
      }
      break;
    case RAYLEIGH:                  /* The variance is really the variance of Gaussian noise */
                                    /* See Whalen, p. 102                                   */
      for(i=0; i<dimensions; i++)
      {
        u1 = [randomObject percent];
        phi = variance[i]*2.;
        if(u1>0.)
          temp = -phi*log(u1);
        else
          temp = 1.;
        if(temp>0.)
          sample[i] = sqrt(temp);
        else
          sample[i] = 0.;
      }
      break;
    case RICIAN:                    /* The variance is equal to 1, A = mean[i]  in Whalen */
      for(i=0; i<dimensions; i++)
      {
        n1 = [self normal01];
        n2 = [self normal01];
        sample[i] = sqrt(n1*n1 + ipow((n2+mean[i]), 2));
      }
      break;
  }
  return sample;
}

/*###############################*
 * Return the value of our den-  *
 * sity function evaluated at x  *
 *###############################*/
- (double)densityFunctionAtX: (double *)x
{
  int    i, j, n, error;
  BOOL   in_range;
  double p_of_x, range, volume, temp, mu, sigma_sq, sum;
  double alpha, c, arg;
  
  switch(distribution)
  {
    case UNIFORM:
      in_range = YES;
      volume   = 1.;
      for(i=0; i<dimensions; i++)
      {
        range   = SQRT3*stdDev[i];
        volume *= range;
        if( (x[i]<(mean[i]+range)) && (x[i]>=(mean[i]-range)) )
          ;
        else
        {
          in_range = NO;
          break;
        }
      }
      if(in_range)
        p_of_x = 1./volume;
      else
        p_of_x = 0.;
      break;
    case NORMAL:
    default:
      temp   = (1-rho*rho);
      if(dimensions == 2)
        p_of_x = exp(rho*(x[0]-mean[0])*(x[1]-mean[1])/stdDev[0]/
                 stdDev[1]/temp);
      else
        p_of_x = 1.;
      temp *= 2.;
      for(i=0; i<dimensions; i++)
      {
        p_of_x *= exp(-(x[i]-mean[i])*(x[i]-mean[i])/temp/variance[i]);
        p_of_x /= SQRT2PI*stdDev[i];
      }
      break;
    case EXPONENTIAL_DIST:
      in_range = YES;
      for(i=0; i<dimensions; i++)
      {
        alpha   = mean[i] - stdDev[i];
        if( x[i]>=alpha )
          ;
        else
        {
          in_range = NO;
          break;
        }
      }
      if(in_range)
      {
        p_of_x = 1.;
        for(i=0; i<dimensions; i++)
        {
          p_of_x *= exp(-(x[i]-mean[i])/stdDev[i]);
          p_of_x /= stdDev[i];
        }
      }
      else
        p_of_x = 0.;
      break;
    case ERLANG_DIST:       /* see Papoulis p. 77 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        c   = mean[i]/variance[i];
        n   = ROUND(mean[i]*c);
       sum = 0.;
       for(j=1; j<n; j++)
         sum += log((double)j);
       temp = n*log(c) + (n-1)*log(x[i]) - c*x[i] - sum;
       p_of_x *= exp(temp);
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case LOG_NORMAL:        /* see Whalen p. 23 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(mean[i]>0.)
        {
          sigma_sq = log(variance[i]/mean[i]*mean[i] + 1.);
          mu       = log(mean[i]) - sigma_sq/2;
        }
        else mu = sigma_sq = 0.;
        temp = log(x[i]) - mu;
        p_of_x *= exp(-temp*temp/2./sigma_sq)/SQRT2PI/sqrt(sigma_sq)/x[i];
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case RICIAN:        /* see Whalen p. 105 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(mean[i]>0. && variance[i]>0.)
        {
          alpha = mean[i];                  /* This isn't correct, but is easy */
          p_of_x *= modbes(alpha*x[i]/variance[i], I0_SELECT, 0);
          p_of_x *= x[i]*exp(-(x[i]*x[i] + alpha*alpha)/2./variance[i])/variance[i];
        }
        else p_of_x = 0.;
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case RAYLEIGH:      /* see Whalen p. 102 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(variance[i]>0.)
          p_of_x *= x[i]*exp(-x[i]*x[i]/2./variance[i])/variance[i];
        else
          p_of_x = 0.;
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case NON_CENTRAL_CHI:       /* see Whalen p. 114 (4-75)         */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(variance[i]>0.)
        {
          arg  = x[i]/variance[i];      /* normalize, (4-76) */
          temp = ModBesselI((double)(degrees/2-1), sqrt(arg*noncentralParameter), &error);
          temp *= exp( -(noncentralParameter+arg)/2.);
          temp *= ipow(arg/noncentralParameter, (degrees-2)/4)/2./variance[i];
          p_of_x *= temp;
        }
        else
          p_of_x = 0.;
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case LOG_RICIAN:        /* see Whalen p. 105, use p(g(x)) */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(mean[i]>0. && variance[i]>0.)
        {
          alpha = mean[i];                  /* This isn't correct, but is easy */
          p_of_x *= modbes(alpha*exp(x[i])/variance[i], I0_SELECT, 1);
          p_of_x *= exp( 2*x[i]-(exp(2*x[i])-2.*alpha*exp(x[i]) +
                                   alpha*alpha)/2./variance[i])/variance[i];
        }
        else p_of_x = 0.;
      }
      if(!in_range)
        p_of_x = 0.;
      break;
  }
  return p_of_x;
}

/*######################################*
 * Return a N(0,1) random sample        *
 *######################################*/
- (double) normal01
{
  double u1,s,ln_s;
  double v1,v1_sq,v2_sq;
  static double v2, sqrt_lns;
  static int odd = 1;

  if(odd)
  {
    s = 2.;
    while(s > 1.)
    {
      u1 = [randomObject percent];
      v2 = u1+u1-1.;
      v2_sq = v2*v2;
      u1 = [randomObject percent];
      v1 = u1 + u1 - 1.;
      v1_sq = v1*v1;
      s  = v1_sq + v2_sq;
    }
    ln_s = log(s);
    sqrt_lns = sqrt(-(ln_s+ln_s)/s);
    odd = 0;
    return v1*sqrt_lns;
  }
  else
  {
    odd = 1;
    return  v2*sqrt_lns;
  }
}

  

@end
