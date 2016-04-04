/************************************************************************
 * This subclass of object implements an object for computing fast      *
 * Fourier transforms, as well as DFTs and DCTs                         *
 *                                                                      *
 * File: FFT.h                                                          *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 09/03/93  - Started                                              *
 *  2. 05/06/95  - DigitalFilter.h -> DigitalFilterOjbc.h.              *
 ************************************************************************/

/**************************
 * Include files:         *
 **************************/
#import "FFT.h"
#import "c_headers.h"
#import "DigitalFilter.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

@implementation FFT

- init
{
  outputSpectrum    = autocorrelation   = NULL;
  frequencyPoints   = twoSidedFrequency = NULL;
  oldSpectrumPoints = oldAutoPoints = numberPoints = 0;
  samplingFreq      = 1.;
  meanIndexType     = MEAN_FROM_SCM;
  twoSidedBandwidth = 0.;
  digitalFilter     = [[DigitalFilter alloc] init];
  return [super init];
}

- initWithSamplingFreq:(float)frequency;
{
  [self init];
  samplingFreq = frequency;
  return self;
}

/*##############################*
 * Allocate memory needed:      *
 *##############################*/
- allocateMemory
{
  return self;
}

/*##############################*
 * Free up memory               *
 *##############################*/
- freeMemory
{
  return self;
}

/*##############################*
 * Free myself:                 *
 *##############################*/
- free
{
  [self freeMemory];
  [super free];
  return self;
}

/*##############################*
 * Set the new sampling freq.   *
 *##############################*/
- setSamplingFrequency:(float)frequency
{
  samplingFreq = frequency;
  return self;
}

/*##############################*
 * Set the new sampling freq.   *
 *##############################*/
- setWindowType:(int)type
{
  [digitalFilter setWindowingType:type];
  return self;
}

/*##############################*
 * Set the way to find the center*
 * of the PSD.                  *
 *##############################*/
- setMeanIndexType:(int)type
{
  meanIndexType = type;
  return self;
}

/*##############################*
 * Return the frequency points  *
 *##############################*/
- (float *)frequencyPoints
{
  int    i;
  static int points_store = 0;
  double frequency, frequency_step;

  if(points_store < numberPoints)
  {
    if(frequencyPoints != NULL)
      free(frequencyPoints);
    frequencyPoints = (float *)malloc(numberPoints*sizeof(float));
    points_store = numberPoints;
  }
  frequency_step = 0.;
  if(numberPoints > 0)
    frequency_step = samplingFreq/numberPoints;
  frequency = 0.;
  for(i=0; i<numberPoints; i++)
  {
    frequencyPoints[i] = frequency;
    frequency += frequency_step;
  }
  return frequencyPoints;  
} 

/*##############################*
 * Return the frequency points  *
 *##############################*/
- (float )deltaFrequency
{
  float frequency_step;
  frequency_step = 0.;
  if(numberPoints > 0)
    frequency_step = samplingFreq/numberPoints;
  return frequency_step;  
}

/*##############################*
 * Return the two sided PSD     *
 * bandwidth, first call        *
 * meanIndexOf                  *
 *##############################*/
- (float )twoSidedBandwidth{ return twoSidedBandwidth; }

/*##############################*
 * Return the frequency points  *
 *##############################*/
- (float *)twoSidedFrequency
{
  int    i;
  static int points_store = 0;
  double frequency, frequency_step;

  if(points_store < numberPoints)
  {
    if(twoSidedFrequency != NULL)
      free(twoSidedFrequency);
    twoSidedFrequency = (float *)malloc(numberPoints*sizeof(float));
    points_store = numberPoints;
  }
  frequency_step = 0.;
  if(numberPoints > 0)
    frequency_step = samplingFreq/numberPoints;
  frequency = -numberPoints*frequency_step/2.;
  for(i=0; i<numberPoints; i++)
  {
    twoSidedFrequency[i] = frequency;
    frequency += frequency_step;
  }
  return twoSidedFrequency;  
} 


/*##########################################################*
 * Calculate the power spectrum                             *
 * of real data.                                            *
 * Note: there are at least two ways to scale the power     *
 * spectrum, one way gets noise PSDs right, the other way   *
 * gets sinusoid magnitudes right.  I think this one gets   *
 * noise magnitudes right.                                  *
 *##########################################################*/
- (float *)powerSpectrum:(float *)input numberPoints:(int)number
{
  int i;
  float temp;
  float resolution_bw;

//
// Allocate space
//
  numberPoints = number;
  if(oldSpectrumPoints < numberPoints)
  {
    oldSpectrumPoints = numberPoints;
    if(outputSpectrum != NULL)
      free(outputSpectrum);
    outputSpectrum = (float *)malloc(numberPoints*sizeof(float));
  }
//
// FInd the FFT, input will contain the real part of FFT, outputSpectrum will
// contain the imaginary part
//
  fft(input, outputSpectrum, numberPoints, REAL_FFT);
  temp = 1.;
//
// Calculate scale factor
//
  if(numberPoints)
  {
    resolution_bw = samplingFreq/(float)numberPoints;
    temp          = resolution_bw*numberPoints*numberPoints;
  }
//
// Find power spectrum from real and imaginary plus scaling
//
  for(i=0; i<numberPoints; i++)
  {
    outputSpectrum[i] = input[i]*input[i] + outputSpectrum[i]*outputSpectrum[i];
    outputSpectrum[i] /= temp;
  }
  return outputSpectrum;
}

/*##############################*
 * Calculate the FFT magnitude  *
 * of real data.                *
 *##############################*/
- (float *)fftMagnitude:(float *)input numberPoints:(int)number
{
  int i;
  float temp;

  numberPoints = number;
  if(oldSpectrumPoints < numberPoints)
  {
    oldSpectrumPoints = numberPoints;
    if(outputSpectrum != NULL)
      free(outputSpectrum);
    outputSpectrum = (float *)malloc(numberPoints*sizeof(float));
  }
  fft(input, outputSpectrum, numberPoints, REAL_FFT);
  temp = (float)numberPoints;
  for(i=0; i<numberPoints; i++)
  {
    outputSpectrum[i] = input[i]*input[i] + outputSpectrum[i]*outputSpectrum[i];
    outputSpectrum[i] = sqrt(outputSpectrum[i])/temp;
  }
  return outputSpectrum;
}

/*##############################*
 * Calculate the cosine xform   *
 * of real data.                *
 *##############################*/
- (float *)cosineTransform:(float *)input numberPoints:(int)number
{

  numberPoints = number;
  if(oldSpectrumPoints < numberPoints)
  {
    oldSpectrumPoints = numberPoints;
    if(outputSpectrum != NULL)
      free(outputSpectrum);
    outputSpectrum = (float *)malloc(numberPoints*sizeof(float));
  }
  dct(input, outputSpectrum, numberPoints, numberPoints, 0);

  return outputSpectrum;
}

/*##############################*
 * Calculate the PSD from one   *
 * sided auto correlation.      *
 *##############################*/
- (float *)psdFromAuto:(float *)input numberPoints:(int)number
{
  int n, k;
  float  *window;
  static float *coefficients;
  static int old_coeff_points;
  double arg;

  numberPoints = 2*number-1;
  if(oldSpectrumPoints < numberPoints)  /* Should check if windowType has changed */
  {
    oldSpectrumPoints = numberPoints;
    if(outputSpectrum != NULL)
      free(outputSpectrum);
    outputSpectrum = (float *)malloc(numberPoints*sizeof(float));
  }
  if(old_coeff_points != numberPoints)
  {
    old_coeff_points = numberPoints;
    if(coefficients != NULL)
      free(coefficients);
    coefficients   = (float *)malloc(numberPoints*sizeof(float));
    arg = TWOPI/numberPoints;
    for(n=0; n<numberPoints; n++)
       coefficients[n] = cos(arg*n);
  }
  window = [digitalFilter windowFunction:numberPoints];
  for(k=0; k<numberPoints; k++)
  {
    outputSpectrum[k] = window[number-1]*input[0];
    for(n=1; n<number; n++)
      outputSpectrum[k] += window[number-1+n]*input[n]*
        (coefficients[ (k*n)%numberPoints] + 
         coefficients[ (k*(numberPoints-1))%numberPoints]);
  }
  return outputSpectrum;
}

/*##############################*
 * Calculate the power spectrum *
 * of complex data. Window the  *
 * data prior to FFT.  Return   *
 * power in watts.              *
 * Note: Input arrays are over- *
 * written.                     *
 *##############################*/
- (float *)avgPowerSpectrum:(float *)realInput imag:(float *)imagInput
               numberPoints:(int)number numberAvgs:(int)avgs
{
  int i, n;
  float temp, spectrum;
  float resolution_bw;
  float *real_ptr, *imag_ptr;
  float *window;

  numberPoints = number;
  if(oldSpectrumPoints < numberPoints)
  {
    oldSpectrumPoints = numberPoints;
    if(outputSpectrum != NULL)
      free(outputSpectrum);
    outputSpectrum = (float *)malloc(numberPoints*sizeof(float));
  }
  for(i=0; i<numberPoints; i++)
    outputSpectrum[i] = 0.;
  window = [digitalFilter windowFunction:numberPoints];
  real_ptr = realInput;
  imag_ptr = imagInput;
  temp     = 1.;
  if(numberPoints)
  {
    resolution_bw = samplingFreq/(float)numberPoints;
    temp          = avgs*resolution_bw*numberPoints*numberPoints;
  }
  for(n=0; n<avgs; n++)
  {
    for(i=0; i<numberPoints; i++)                   /* Window the input arrays  */
    {
      real_ptr[i] = window[i]*real_ptr[i];
      imag_ptr[i] = window[i]*imag_ptr[i];
    }
    fft(real_ptr, imag_ptr, numberPoints, COMPLEX_FFT);
    for(i=0; i<numberPoints; i++)
    {
      spectrum = real_ptr[i]*real_ptr[i] + imag_ptr[i]*imag_ptr[i];
      outputSpectrum[i] += spectrum/temp;
    }
    real_ptr += numberPoints;
    imag_ptr += numberPoints;
  }
  return outputSpectrum;
}

/*##############################*
 * Calculate the autocovariance *
 * function by first finding    *
 * the averaged power spectrum  *
 * Note: Input arrays are over- *
 * written.                     *
 *##############################*/
- (float *)autoFromPSD:(float *)inputSpectrum numberPoints:(int)number type:(int)type
          scale:(float *)scale
{
  int i;
  if(oldAutoPoints < number)
  {
    oldAutoPoints = number;
    if(autocorrelation != NULL)
      free(autocorrelation);
    autocorrelation = (float *)malloc(number*sizeof(float));
  }
/* Now take inverse FFT to find autocorrelation */
  for(i=0; i<number; i++)                   /* First zero out imaginary part    */
    autocorrelation[i] = 0.;
  
  fft(inputSpectrum, autocorrelation, number, INVERSE_FFT);
  *scale = sqrt(autocorrelation[0]*autocorrelation[0] + inputSpectrum[0]*inputSpectrum[0]);
  switch(type)
  {
    case COMPLEX_MAGNITUDE:
    default:
      for(i=0; i<number; i++)
        autocorrelation[i] = sqrt(autocorrelation[i]*autocorrelation[i] +
                                  inputSpectrum[i]*inputSpectrum[i]);
      break;
    case COMPLEX_PHASE:
      for(i=0; i<number; i++)
        autocorrelation[i] = atan2(autocorrelation[i], inputSpectrum[i]);
      break;
    case COMPLEX_REAL_PART:
      for(i=0; i<number; i++)
        autocorrelation[i] = inputSpectrum[i];
      break;
    case COMPLEX_IMAG_PART:
      for(i=0; i<number; i++)
        autocorrelation[i] = autocorrelation[i];
      break;
  }
  return autocorrelation;
}

/*##################################*
 * This  method finds the modified  * 
 * short-time autocorrelation func- *
 * tion from (4.30), page 143,      *
 * Rabiner and Schafer. A pointer   *
 * to the output autocorrelation is *
 * returned.                        *
 * Note:Result is divided by N      *
 *##################################*/
- (float *)autocorrelationOf:(float *)realInput imag:(float *)imagInput
          numberPoints:(int)number windowLength:(int)length lagLength:(int)lag type:(int)type
          scale:(float *)scale
{
  int    m,k;   
  float real_part, imag_part;
  float *window;
        
// Check # of points: 
  if( (number<length) || (length<=lag) )
    return NULL;                                             

  if(oldAutoPoints < lag)
  {
    oldAutoPoints = lag;
    if(autocorrelation != NULL)
      free(autocorrelation);
    autocorrelation = (float *)malloc(lag*sizeof(float));
  }
  window = [digitalFilter windowFunction:length];
/****************************
 * Now, calculate it:       *  
 * choose n = 0  in (4.36)  *
 ****************************/
  switch(type)
  {
    case COMPLEX_MAGNITUDE:
    default:
      for(k=0; k<lag; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<length-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        autocorrelation[k] = sqrt(real_part*real_part + imag_part*imag_part);
        autocorrelation[k] /= (double) (length-k);
      }
      *scale = autocorrelation[0];
      break;
    case COMPLEX_PHASE:
      for(k=0; k<lag; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<length-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        autocorrelation[k] = atan2(real_part, imag_part);
      }
      *scale = 1;
      break;
    case COMPLEX_REAL_PART:
      for(k=0; k<lag; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<length-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        autocorrelation[k] = real_part/(double) (length-k);
        if(k==0)
          *scale = sqrt(real_part*real_part + imag_part*imag_part)/(double) (length-k);
      }
      break;
    case COMPLEX_IMAG_PART:
      for(k=0; k<lag; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<length-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        autocorrelation[k] = imag_part/ (length-k);
        if(k==0)
          *scale = sqrt(real_part*real_part + imag_part*imag_part)/(double) (length-k);
      }
      break;
  }

  return autocorrelation;
}

/*##############################*
 * Find the center of mass of   *
 * the spectrum.                *
 *##############################*/
- (float  )meanIndexOf:(float *)inputSpectrum size:(int)points
{
  int i, min_index, max_index, flo, fhi;
  float min_spect, max_spect, bandwidth_threshold;
  double mean_sum, spectrum_sum;

// First find the two sided bandwidth by checking when spectrum falls below threshold:
  float_min_max(inputSpectrum, points, &min_spect, &max_spect, &min_index, &max_index);
  bandwidth_threshold = max_spect/pow(10., DB_THRESHOLD/10.);
  flo = fhi = 0;
  for(i=max_index; i>0; i--)
    if(inputSpectrum[i] < bandwidth_threshold)
    {
      flo = i;
      break;
    }
  for(i=max_index; i<points; i++)
    if(inputSpectrum[i] < bandwidth_threshold)
    {
      fhi = i;
      break;
    }
  twoSidedBandwidth = (fhi - flo)*[self deltaFrequency];
  if(meanIndexType == MEAN_FROM_THRESHOLD)
    mean_sum = (float)(fhi+flo)/2.;
// Else, find spectral center of mass
  else
  {
    mean_sum = spectrum_sum = 0.;
    for(i=1; i<points; i++)
    {
      mean_sum     += (double)i*inputSpectrum[i];
      spectrum_sum += inputSpectrum[i];
    }
    if(spectrum_sum>0.)
     mean_sum /= spectrum_sum;
  }
  return (float)mean_sum;
}

/*##############################*
 * Shift spectrum by half its   *
 * length for display.          *
 *##############################*/
- shiftRight:(float *)inputArray by:(int)shift size:(int)points
{
  int i, j;
  float *temp;

  shift %= points;
  temp = (float *)malloc(points*sizeof(float));
  j = points - shift;
  j %= points;
  for(i=0; i<points; i++)
  {
    temp[i] = inputArray[j++];
    j %= points;
  }
  for(i=0; i<points; i++)
    inputArray[i] = temp[i];
  free(temp);
  return self;
}

/*##############################*
 * Convert the power spectrum   *
 * to dB:                       *
 *##############################*/
- dBSpectrum:(float *)inputSpectrum size:(int)points
{
  int i;

  for(i=0; i<points; i++)
  {
    if(inputSpectrum[i]>0.)
      inputSpectrum[i] = 10.*log10(inputSpectrum[i]);
  }
  return self;
}

/*##############################*
 * Normalize the autocorrelation*
 *##############################*/
- normalizeAuto:(float *)autoc size:(int)points by:(float)scale
{
  int i;

  if(scale > 0.)
    for(i=0; i<points; i++)
      autoc[i] /= scale;
  return self;
}

@end