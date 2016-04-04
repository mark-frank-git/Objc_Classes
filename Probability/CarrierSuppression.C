/****************************************************************************
 * This subclass of object calculates QAM carrier suppression statistics.	*
 *																			*
 * File: /User/frank/Objc_Classes/Probability/CarrierSuppression.h			*
 *																			*
 * Revision History:														*
 *  1. 07/11/97 - Started													*
 ***************************************************************************/

#import "CarrierSuppression.h"
#import "Numerical.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>
#import <iostream.h>

#define NUMBER_STD_DEVIATIONS	4.
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

#define DEFAULT_SAMPLES		100
#define DEFAULT_NOISE		1.
#define DEFAULT_OFFSET		1.
#define SIGNAL_POWER		10.
#define LARGE_BESSEL_ARG	5.

@implementation CarrierSuppression

- init
{
  double tol1;

  [super init];
  integrator		= [[Numerical alloc] init];
  numberSamples		= DEFAULT_SAMPLES;
  noiseVariance		= DEFAULT_NOISE;
  iqOffset			= DEFAULT_OFFSET;
  signalPower		= SIGNAL_POWER;

/*******************
 * Find machine eps*
 *******************/
  machineEps = 1.;
  tol1 = 1.1;
  while(tol1 > 1.)
  {
    machineEps /= 2.;
    tol1 = 1. + machineEps;
  }
  
  return self;
}

/*##############################*
 * Local methods				*
 *##############################*/
- findNoncentralVariance
{
  noncentralVariance	= 4.*noiseVariance*(noiseVariance/numberSamples + iqOffset)/numberSamples;
  return self;
}

/*##############################*
 * These methods set parameters *
 *##############################*/
- setNumberSamples:(int) samples
{
  numberSamples	= MAX(1, samples);
  return self;
}

- setNoiseVariance: (double) variance
{
  noiseVariance	= variance;
  return self;
}

- setIQOffset: (double) sSquared
{
  iqOffset	= sSquared;
  return self;
}

- setSignalPower: (double) power
{
  signalPower	= power;
  return self;
}

/*#################################*
 * Return the non-central var.		*
 *###################################*/
- (double) getNoncentralVariance
{
  [self findNoncentralVariance];
  return noncentralVariance;
}
  

#define NUMBER_STDEVS	6.					/* Number of std deviations to integrate over */
/*##################################*
 * Return the output variance by	*
 * integrating x^2 times the pdf.	*
 *##################################*/
- (double) outputVariance
{
  double non_central_variance, std_dev;
  double output_variance;

  [self findNoncentralVariance];
  std_dev	= sqrt(noncentralVariance);

  if(std_dev > 0.)
  {
    xMin	= MAX(machineEps, (iqOffset - NUMBER_STDEVS*std_dev) );
    xMax	= iqOffset + NUMBER_STDEVS*std_dev;
  }
  else
  {
    xMax	= 2.*iqOffset;
    xMin	= machineEps;
  }

  integrateVariance	= YES;
  output_variance	= [integrator integrateFrom:xMin to:xMax from:self];

  return output_variance;
}

/*##################################*
 * Return the output variance by	*
 * integrating x^2 times the pdf.	*
 *##################################*/
- (double) outputVarianceFrom:(double)xLo to:(double)xHi
{
  double non_central_variance, std_dev;
  double output_variance;

  xMin				= MAX(machineEps, xLo );
  xMax				= xHi;
  integrateVariance	= YES;
  output_variance	= [integrator integrateFrom:xMin to:xMax from:self];

  return output_variance;
}

/*##################################*
 * Return the output expected value	*
 * integrating x times the pdf.		*
 *##################################*/
- (double) outputExpectedValue
{
  double non_central_variance, std_dev;
  double output_expected;

  [self findNoncentralVariance];
  std_dev	= sqrt(noncentralVariance);

  if(std_dev > 0.)
  {
    xMin	= MAX(machineEps, (iqOffset - NUMBER_STDEVS*std_dev) );
    xMax	= iqOffset + NUMBER_STDEVS*std_dev;
  }
  else
  {
    xMax	= 2.*iqOffset;
    xMin	= machineEps;
  }

  integrateVariance	= NO;
  output_expected	= [integrator integrateFrom:xMin to:xMax from:self];

  return output_expected;
}

/*##################################*
 * Return the output expected value	*
 * integrating x times the pdf.		*
 *##################################*/
- (double) outputExpectedValueFrom:(double)xLo to:(double)xHi
{
  double non_central_variance, std_dev;
  double output_expected;

  xMin				= MAX(machineEps, xLo );
  xMax				= xHi;

  integrateVariance	= NO;
  integratePDF		= NO;
  integrateMeanPDF	= NO;
  output_expected	= [integrator integrateFrom:xMin to:xMax from:self];

  return output_expected;
}

/*##################################*
 * Return the output expected value	*
 * integrating x times the pdf.		*
 *##################################*/
- (double) outputMeanOfPDFFrom:(double)xLo to:(double)xHi
{
  double non_central_variance, std_dev;
  double output_expected;

  xMin				= MAX(machineEps, xLo );
  xMax				= xHi;

  integrateVariance	= NO;
  integratePDF		= NO;
  integrateMeanPDF	= YES;
  output_expected	= [integrator integrateFrom:xMin to:xMax from:self];

  return output_expected;
}
/*##################################*
 * Return the output variance by	*
 * integrating x^2 times the pdf.	*
 *##################################*/
- (double) pdfFrom:(double)xLo to:(double)xH
{
  double pdf;

  xMin				= MAX(machineEps, xLo );
  xMax				= xH;

  integrateVariance	= NO;
  integratePDF		= YES;
  pdf	= [integrator integrateFrom:xMin to:xMax from:self];

  return pdf;
}

/*##################################*
 * Return the pdf of the non-		*
 * central chi squared				*
 * evaluated at the point, x.		*
 *###############################*/
- (double)pdf:(double)x
{
  double	n_over_sigma2, bessel_fn;
  double	bessel_arg, exp_fn, exp_arg;
  double	pdf, sqrt_x_iq;

  sqrt_x_iq		= sqrt(x*iqOffset);
  n_over_sigma2	= numberSamples/noiseVariance;
  bessel_arg	= sqrt_x_iq*n_over_sigma2;
  exp_arg		= - 0.5*(x + iqOffset)*n_over_sigma2;
  if(bessel_arg > LARGE_BESSEL_ARG)
  {
    bessel_fn	= modbes(bessel_arg, 0, 1);							// exp(-arg) Io(arg) is returned
    exp_arg		+= bessel_arg;
  }
  else
    bessel_fn	= modbes(bessel_arg, 0, 0);
  exp_fn		= exp(exp_arg);

/////  cout << "bessel_arg = " << bessel_arg << ", exp_arg = " << exp_arg << "\n";
  pdf			= 0.5*n_over_sigma2*exp_fn*bessel_fn;

  return pdf;
}

/*##################################*
 * Return the pdf of the carrier	*
 * suppression estimator			*
 * evaluated at the point, x.		*
 *###############################*/
- (double)pdfCS:(double)x
{
  double	n_over_sigma2, bessel_fn;
  double	bessel_arg, exp_fn, exp_arg;
  double	pdf, sqrt_x1_iq, x1, g_prime;

  x1			= signalPower*pow(10., -x/10.);
  sqrt_x1_iq	= sqrt(x1*iqOffset);
  n_over_sigma2	= numberSamples/noiseVariance;
  bessel_arg	= sqrt_x1_iq*n_over_sigma2;
  exp_arg		= - 0.5*(x1 + iqOffset)*n_over_sigma2;
  if(bessel_arg > LARGE_BESSEL_ARG)
  {
    bessel_fn	= modbes(bessel_arg, 0, 1);							// exp(-arg) Io(arg) is returned
    exp_arg		+= bessel_arg;
  }
  else
    bessel_fn	= modbes(bessel_arg, 0, 0);
  exp_fn		= exp(exp_arg);

  pdf			= 0.5*n_over_sigma2*exp_fn*bessel_fn;
//
//  Now, add 1/|g'(x)|
//
  g_prime		= 10.*log10(exp(1.))/x1;
  pdf			/= g_prime;

  return pdf;
}

/*##################################*
 * Return the  integrand			*
 * evaluated at the point, x.		*
 *###############################*/
- (double)functionToIntegrate:(double)x
{
  double	pdf, output;

  pdf			= [self pdfCS:x];

  if(integrateVariance)
    output		= x*x*pdf;										// Integral( g(x) p(x) )
  else if(integratePDF)
    output		= pdf;
  else if(integrateMeanPDF)
    output		= x*pdf;
  else
    output		= x*pdf;
  return output;
}

/*##################################*
 * Return the  integrand			*
 * evaluated at the point, x.		*
 *###############################*/
- (double)OldfunctionToIntegrate:(double)x
{
  double	pdf, output;

  pdf			= [self pdf:x];

  output		= log10(x/signalPower);								// add g(x)
  if(integrateVariance)
    output		= 100.*output*output*pdf;							// Integral( g(x) p(x) )
  else if(integratePDF)
    output		= pdf;
  else if(integrateMeanPDF)
    output		= x*pdf;
  else
    output		= -10.*output*pdf;
  return output;
}


@end
