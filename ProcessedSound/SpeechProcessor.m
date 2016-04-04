/************************************************************************
 * This object is a subclass of Object.  It provides processing of      *
 * speech signals.  Some of the code originated from Bill Kushner of    *
 * CSRL.                                                                *
 *                                                                      *
 * File: /User/frank/Objc_Classes/ProcessedSound/SpeechProcessor.m      *
 *                                                                      *
 * The speech is processed on a frame-by-frame basis.  Each processed   *
 * frame of speech consists of the following data:                      *
 *   1. spectral data: e.g., cepstral, filter bank outputs, etc.        *
 *   2. zero crossing data                                              *
 *   3. peak to peak amplitude                                          *
 *   4. frame energy                                                    *
 *   5. speech/noise classification                                     *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/20/93 - Converted from perceptual.c                           *
 *  2. 05/20/93 - Calculate energyData[] in calculateSpectralData.      *
 *  3. 06/17/93 - Added lpcMagnitudeResponse                            *
 ************************************************************************/

/*************************************
 * Constants:                        *
 *************************************/
#define MAX_POWER       15                              /* 2**15 = max frame length     */
#define STEVENS_POWER   0.333333333333                  /* 1/3 power law                */
#define NOISE_FRAMES    10                              /* # of beginning noise frames  */
#define NUMBER_STDEV    4.5                             /* For speech/noise classify    */
#define SPECTRAL_THRESHOLD 0.38                         /* For speech/noise classify    */
#define SPECTRAL_FLOOR  0.2                             /* noise floor scale factor      */
#define NOISE_FACTOR    0.9                             /* noise subtraction scale factor */
#define TIME_CONSTANT   0.9                             /* noise update time constant     */
#define TIMEOUT         100                             /* # of frames of speech for timeout */
#define NNSUF           20

/* Two-dimensional indices */
#define INDEX_FB(i,j) [(j)+(i)*(numberFilterBanks+1)]
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

#import "SpeechProcessor.h"
#import "speechParams.h"
#import "c_headers.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

/*************************************
 * These strings must match the      *
 * #defines in plotParams.h      *
 *************************************/
char *spectral_names[] = {"LPC_Cepstrum", "Cepstrum", "LPC_Spectrum",
                               "Sonogram", "Energy"};
char *window_names[]   = {"Hamming", "Hanning", "Rectangular", "Triangular",
                             "Blackman", "Blackman_Harris"};
char *lifter_names[]  = {"No_liftering", "Raised_Sine", "Linear",
                                           "Gaussian", "Exponential"};
char *frequency_analysis_names[] = {"Perceptual", "Linear"};


@implementation SpeechProcessor

/*###################################*
 * Initialize instance variables     *
 *###################################*/
- initWithOrder:(int)order cepstralLength:(int)clength
               frameLength:(int)flength overlapLength:(int)olength
               numberBanks:(int)numberBanks samplingFrequency:(float)frequency
{
  int i;

//
// Default some instance variables:
//
  spectralType          = LPC_CEPSTRUM;
  windowType            = HAMMING;
  lifterType            = EXPONENTIAL;
  perceptualProcessing  = YES;
  spectralSubtraction   = NO;
  spectralFloor         = SPECTRAL_FLOOR;
  noiseFactor           = NOISE_FACTOR;
  numberFrames          = 0;
  window = lifter       = NULL;
  filterBanks           = loudnessWeighting = centerFrequencies = NULL;
//
// Set other instance variables:
//
  lpcOrder              = order;
  cepstralLength        = spectralLength = clength;
  frameLength           = flength;
  overlapLength         = olength;
  numberFilterBanks     = numberBanks;
  samplingFrequency     = frequency;
  foldingFrequency      = frequency/2.;
  powerOf2Length        = 2;
  oldCepstralLength     = oldFrameLength = 0;
  i = 0;
  while(powerOf2Length<frameLength && i<MAX_POWER)
    powerOf2Length *= 2;
  if(i==MAX_POWER)
    printf("Bad frame length in SpeechProcessor\n");
  [self initializeWindows];
  [self initializeFilterBanks];
  return self;
}

/*######################################*
 * Set the spectral data length based   *
 * on the spectralType, etc.        *
 *######################################*/
- initializeSpectralLength
{
  switch(spectralType)
  {
    case LPC_CEPSTRUM:
    case CEPSTRUM:
    default:
      spectralLength = cepstralLength;
      break;
    case LPC_SPECTRUM:
    case LPC_SONOGRAM:
      spectralLength = lpcOrder+1;
      break;
    case SONOGRAM:
      if(perceptualProcessing)
        spectralLength = numberFilterBanks;
      else
        spectralLength = powerOf2Length/2+1;
      break;
  }
  return self;
}

#define TWOSIGMASQ  72.     /* Used in Gaussian lifter */
#define EXPONENT    0.6     /* Used in Exponential lifter */
/*##################################################################*
 * Initialize the speech window and the liftering function      *
 * References:                                                      *
 * 1.  P. Meyer, J. Schroeter, M.M. Sondhi, "Design and evaluation  *
 * of optimal cepstral lifters for accessing articulatory codebooks,*
 * IEEE Trans. Signal Processing, vol. 39, no. 7, pp. 1493-1502,    *
 * July 1991.                                                       *
 * 2. B.-H. Juang, L.R. Rabiner, and J.G. Wilpon, "On the use of    *
 * bandpass liftering in speech recognition," IEEE Trans. ASSP, vol.*
 * ASSP-35, no. 7, July 1987.                                       *
 ********************************************************************/
- initializeWindows
{
  int i, k, q;
  double temp1, temp2, arg;

// First do the lifter:
  if(oldCepstralLength < cepstralLength)
  {
     if(lifter != NULL)
       free(lifter);
     lifter = (float *)malloc((cepstralLength+2)*sizeof(float));
     oldCepstralLength = cepstralLength;
  }
  q = cepstralLength + 1;
  switch(lifterType)
  {
    case RAISED_SINE:
      temp1 = PI/(double)cepstralLength;
      temp2 = 0.5*(double)cepstralLength;
      for(k=0; k<=q; k++)
        lifter[k] = 1. + temp2*sin(k*temp1);
      break;
    case LINEAR:
      for(k=0; k<=q; k++)
        lifter[k] = (float) k+1;
      break;
    case GAUSSIAN:
      for(k=0; k<=q; k++)
        lifter[k] = (double)k*exp(-k*k/TWOSIGMASQ);
      break;
    case EXPONENTIAL:
      for(k=0; k<=q; k++)
        lifter[k] = pow((double)(k+1),EXPONENT);
      break;
    default:
      for(k=0; k<=q; k++)
        lifter[k] = 1.;
      break;
  }

// Now do the window:
  if(oldFrameLength < frameLength)
  {
     if(window != NULL)
       free(window);
     window = (float *)malloc(frameLength*sizeof(float));
     oldFrameLength = frameLength;
  }
  switch(windowType)
  {
    case RECTANGULAR:
      for(i=0; i<frameLength; i++)
        window[i] = 1.;
      break;
    case HAMMING:
      arg = TWOPI/(frameLength-1.);
      for(i=0; i<frameLength; i++)
        window[i] = 0.54 - 0.46*cos(arg*i);
      break;
    case HANNING:
      arg = TWOPI/(frameLength-1.);
      for(i=0; i<frameLength; i++)
        window[i] = 0.50 - 0.50*cos(arg*i);
      break;
    case TRIANGULAR:
      arg = 2.0/(frameLength-1.);
      for(i=0; i<=(frameLength-1)/2; i++)
        window[i] = i*arg;
      for(   ; i<frameLength; i++)
        window[i] = 2. - i*arg;
      break;
    case BLACKMAN:
      arg = TWOPI/(frameLength-1.);
      for(i=0; i<frameLength; i++)
        window[i] = 0.42 - 0.50*cos(arg*i) + 0.08*cos((arg+arg)*i);
      break;
    case BLACKMAN_HARRIS:
      arg = TWOPI/(frameLength-1.);
      for(i=0; i<frameLength; i++)
        window[i] = 0.35875 - 0.48829*cos(arg*i) +
                      0.14128*cos((arg+arg)*i) - 0.01168*cos(3.*arg*i);
      break;
    default:      
      for(i=0; i<frameLength; i++)
        window[i] = 1.;
      break;
  }
  return self;
}

#define KLUGE   128.            /* I don't know what this should be? */
#define ZERO    1200.           /* zero at 1200 Hz */
#define POLE_1  400.            /* double pole at 400 Hz */
#define POLE_2  3100.           /* next pole at 3100 Hz */
#define POLE_3  5000.           /* triple pole at 5000 Hz */
#define GAIN    1.75901e27      /* to get unity gain @4000 */
/*######################################################*
 * Initialize the perceptual  filter banks.             *
 * purpose:                                             *
 * to calculate the auditory filter response curves and *
 * loudness compensation tables for use in the anrfv    *
 * filter bank routine.                                 *
 *                                                      *
 * This model is based on:                              *
 * 1) H. Hermansky, "Perceptual linear predictive (PLP) *
 * analysis of speech," J. Acoust. Soc. Am. 87, 4,      *
 * pp. 1738-1752, April '90                             *
 * 2) ghitza, o. "auditory nerve representation criteria*
 * for speech analysis/synthesis", ieee trans. acous.,  *
 * speech, and sig. proc., vol. assp-35, no.6, june,    *
 * 1987, pp.736-740.                                    *
 *######################################################*/
- initializeFilterBanks
{
  int   i,k,l;
  float scf,bfact,bfact2,ffact,f2;
  double temp, zero, pole1, pole2, pole3, max_filter;
  double fo, omega;

  int   half_length = powerOf2Length/2;
  float ee = exp(1.0);

// Allocate space for arrays:
  if(centerFrequencies != NULL)
  {
     free(centerFrequencies);
     free(loudnessWeighting);
     free(filterBanks);
  }
  centerFrequencies = (float *)malloc((numberFilterBanks+1)*sizeof(float));
  loudnessWeighting = (float *)malloc((numberFilterBanks+1)*sizeof(float));
  filterBanks   = (float *)malloc((numberFilterBanks+1)*
                                           (half_length+1)*sizeof(float));
  if( (foldingFrequency >= 9000.0) && 
                 (foldingFrequency <= 11000.0) )    /* for 10,000 */
   scf = 0.1555;
  else if (foldingFrequency >= 7000.0)          /* for 8,000  */
   scf = 0.1454;
  else if (foldingFrequency >= 4500.0)          /* for 5,000  */
   scf = 0.1241;
  else if (foldingFrequency >= 3500.0)          /* for 4,000  */
   scf = 0.114;
  else
  {
    printf("Undefined sampling frequency in anrfv\n");
    scf = 0.;
  }
  scf   *= KLUGE/half_length;
  bfact  = scf*half_length/(float)numberFilterBanks;
  bfact2 = bfact/2.0;
  ffact  = foldingFrequency/(float)half_length;

/*******************************************
 * calculate filter center frequencies (in *
 * bark scale) and equal loudness weighting*
 *******************************************/
  zero = ZERO*TWOPI;
  zero *= zero;
  pole1 = POLE_1*TWOPI;
  pole1 *= pole1;
  pole2 = POLE_2*TWOPI;
  pole2 *= pole2;
  pole3 = POLE_3*TWOPI;
  pole3 *= pole3;
  pole3 *= pole3*pole3;
  max_filter = 0.;
  for (k=numberFilterBanks; k>=0; k--)
  {
    centerFrequencies[k] = bfact*k;
    fo = 300.0*(pow(ee,(centerFrequencies[k]/6.0)) -
                pow(ee,(-centerFrequencies[k]/6.0)));
    f2 = fo*fo*TWOPI*TWOPI;
    loudnessWeighting[k] = f2*f2/(f2+pole1)/(f2+pole1);
    loudnessWeighting[k] = loudnessWeighting[k]*(f2+zero);
    loudnessWeighting[k] = loudnessWeighting[k]/(f2+pole2);
    loudnessWeighting[k] = GAIN*loudnessWeighting[k]/(f2*f2*f2+pole3);
    if(loudnessWeighting[k]>max_filter)
      max_filter = loudnessWeighting[k];
  }
/* Normalize max filter gain = 1 */
  for(k=0; k<=numberFilterBanks; k++)
   loudnessWeighting[k] /= max_filter;

/*
 Find the perceptual filter banks shaping factors:
*/
  for (l=0; l<=numberFilterBanks; l++)
  {
    for (i=0; i<=half_length; i++)
    {
      temp  = i*ffact/600.0;
      omega = 6.0*log(temp + sqrt(temp*temp + 1.0));

      if (omega <= (centerFrequencies[l]-bfact2))
        filterBanks INDEX_FB(i,l) =
                      pow(10.0,(omega-centerFrequencies[l]+bfact2));
      if ((omega > (centerFrequencies[l]-bfact2)) && 
                      (omega < (centerFrequencies[l]+bfact2)))
        filterBanks INDEX_FB(i,l) = 1.0;
      if (omega >= (centerFrequencies[l]+bfact2))
        filterBanks INDEX_FB(i,l) =
                       pow(10.0,(-2.5*(omega-centerFrequencies[l]-bfact2)));
      if(filterBanks INDEX_FB(i,l) < 1.5259e-5)
        filterBanks INDEX_FB(i,l) = 0.0;
    }
  }
  return self;
}


/*******************************
 * These methods return     *
 * instance variables:          *
 *******************************/
- (int)spectralType{ return spectralType;}
- (int)windowType{ return windowType; }
- (int)lifterType{ return lifterType;}
- (int)lpcOrder{ return lpcOrder;}
- (int)cepstralLength{ return cepstralLength;}
- (int)frameLength{ return frameLength;}
- (int)overlapLength{ return overlapLength;}
- (int)powerOf2Length{ return powerOf2Length;}
- (int)numberFilterBanks{ return numberFilterBanks;}
- (int)numberFrames{ return numberFrames;}
- (int)spectralLength{ return spectralLength;}
- (int *)speechOrNoise{ return speechOrNoiseFromSpectral;}
- (BOOL)perceptualProcessing{ return perceptualProcessing;}
- (BOOL)spectralSubtraction{ return spectralSubtraction;}
- (float *)windowData{ return window;}
- (float *)lifterData{ return lifter;}
- (float *)spectralData{ return spectralData;}
- (float *)ptpData{ return ptpData;}
- (float *)zeroCrossData{ return zeroCrossData;}
- (float *)energyData{ return energyData;}
- (float)samplingFrequency{ return samplingFrequency;}
- (float)spectralFloor{ return spectralFloor;}
- (float)noiseFactor{ return noiseFactor;}
- (char *)spectralName{ return spectral_names[spectralType];}
- (char *)windowName{  return window_names[windowType];}
- (char *)lifterName{  return lifter_names[lifterType];}
- (char *) frequencyAnalysisName
{
  if(perceptualProcessing)
    return frequency_analysis_names[0];
  else
    return frequency_analysis_names[1];
}

/*##############################*
 * These methods set instance   *
 * variables:           *
 *##############################*/
- setSpectralType: (int)type
{
  spectralType = type;
  [self initializeSpectralLength];
  return self;
}

- setWindow: (int)type
{
  windowType = type;
  [self initializeWindows];
  return self;
}

- setLifter: (int)type
{
  lifterType = type;
  [self initializeWindows];
  return self;
}

- setLPCFilterOrder: (int)order
{
  lpcOrder = order;
  [self initializeSpectralLength];
  return self;
}

- setCepLength: (int)length
{
  cepstralLength = length;
  [self initializeSpectralLength];
  [self initializeWindows];
  return self;
}

- setSpeechFrameLength: (int)length
{
  int i;
  frameLength = length;
  powerOf2Length = 2;
  i = 0;
  while(powerOf2Length<frameLength && i<MAX_POWER)
    powerOf2Length *= 2;
  if(i==MAX_POWER)
    printf("Bad frame length in SpeechProcessor\n");
  [self initializeSpectralLength];
  [self initializeWindows];
  [self initializeFilterBanks];
  return self;
}

- setOverlapLength: (int)length
{
  overlapLength = length;
  return self;
}

- setNumberFilterBanks: (int)number
{
  numberFilterBanks = number;
  [self initializeSpectralLength];
  [self initializeFilterBanks];
  return self;
}

- setPerceptualProcessing: (BOOL)flag
{
  perceptualProcessing = flag;
  [self initializeSpectralLength];
  return self;
}

- setSpectralSubtraction: (BOOL)flag
{
  spectralSubtraction = flag;
  return self;
}

- setSpectralFloor: (float)floor
{
  spectralFloor = floor;
  return self;
}

- setNoiseFactor: (float)factor
{
  noiseFactor = factor;
  return self;
}

- setSamplingFrequency: (float)frequency
{
  samplingFrequency = frequency;
  foldingFrequency  = frequency/2.;
  [self initializeFilterBanks];
  return self;
}


/*********************************
 * Find the spectral data which  *
 * depends on spectralType.      *
 * Do perceptual weighting if    *
 * perceptualProcessing == YES.  *
 *********************************/
- (float *)calculateSpectralData:(float *)speechInput start:(int)start
           length:(int)length
{
  int      i, j, k, n, n_outputs;
  int      window_size, end;
  int      total_length;
  static   int old_length = 0;
  static   int old_window = 0;
  float    spectral_center;
  float    *auto_correlation, *filter_outputs;
  float    *power_spectrum;
  float    *lpc_frame, *cepstral_frame;
  static float *windowed_speech = NULL;
  static float *noise_spectrum;

/**************************
 * Initialize variables:  *
 **************************/
  numberFrames   = (length-frameLength)/(frameLength-overlapLength);
  total_length    = numberFrames*spectralLength;
  window_size     = powerOf2Length;
  end             = start+length;
  filter_outputs  = NULL;

/********************************
 * Allocate arrays:             *
 ********************************/
  if(old_length < total_length)
  {
     if(spectralData != NULL)
     {
       free(spectralData);
       free(speechOrNoiseFromSpectral);
       free(energyData);
     }
     spectralData  = (float *)malloc(total_length*sizeof(float));
     energyData    = (float *)malloc(numberFrames*sizeof(float));
     old_length    = total_length;
     speechOrNoiseFromSpectral = (int *)malloc(numberFrames*sizeof(int));
  }
  if(old_window < window_size)
  {
    if(windowed_speech != NULL)
    {
       free(windowed_speech);
       free(noise_spectrum);
    }
    windowed_speech  = get_flt_spc(window_size);
    noise_spectrum   = get_flt_spc(window_size/2 + 1);
    for(i=frameLength; i<window_size; i++)  /* zero pad  */
      windowed_speech[i] = 0.;
    old_window = window_size;
    for(i=0; i<window_size/2+1; i++)
      noise_spectrum[i] = 0.;
  }
  n_outputs = 0;
  for(n=0; n<numberFrames; n++)
  {
    j = start;
    for(i=0; i<frameLength; i++)
    {
      if(j<end)
        windowed_speech[i] = window[i]*speechInput[j++];
      else
        windowed_speech[i] = 0.;
    }
    start += frameLength - overlapLength;
/***************************************
 * Find the power spectrum of the       *
 * current frame of speech              *
 ***************************************/
    power_spectrum = [self powerSpectrumOf:windowed_speech length:window_size
                             frameEnergy:&energyData[n] spectralCenter:&spectral_center];
    speechOrNoiseFromSpectral[n] = [self speechOrNoise:energyData[n] frame:n];
    if(spectralSubtraction)
      [self subtractNoiseFrom:power_spectrum noise:noise_spectrum 
                   speechFlag:speechOrNoiseFromSpectral[n] length:window_size
                  frameEnergy:&energyData[n]];
    if(perceptualProcessing)
    {
      filter_outputs = [self filterBankOutputs:power_spectrum
                        length:window_size];
      auto_correlation = [self autoFromFilterBanks:filter_outputs
                        length:window_size];
    }
    else
      auto_correlation = [self autoFromSpeech:windowed_speech
                         length:frameLength lagLength:lpcOrder+1];

/***************************************
 * Load data into output array:        *
 ***************************************/
    switch(spectralType)
    {
      case SONOGRAM:
      default:
        if(perceptualProcessing)
        {
          k = MIN(spectralLength, numberFilterBanks);
          for(i=0; i<k; i++)
            spectralData[n_outputs++] = filter_outputs[i+1];
        }
        else
        {
          for(i=0; i<spectralLength; i++)
            spectralData[n_outputs++] = power_spectrum[i];
        }
        break;
      case LPC_CEPSTRUM:
        lpc_frame      = [self lpcFromAuto:auto_correlation order:lpcOrder];
        cepstral_frame = [self cepstralFromAuto:auto_correlation
                          lpcOrder:lpcOrder length:spectralLength];
        for(i=0; i<spectralLength; i++)
          spectralData[n_outputs++] = lifter[i]*cepstral_frame[i+1];
        break;
      case CEPSTRUM:
        if(perceptualProcessing)
          power_spectrum = [self powerSpectrumFromAuto:auto_correlation
                                                length:window_size];
        cepstral_frame = [self cepstralFromSpectrum:power_spectrum
                                             length:cepstralLength];
        for(i=0; i<spectralLength; i++)
          spectralData[n_outputs++] = lifter[i]*cepstral_frame[i+1];
        break;
      case LPC_SPECTRUM:
      case LPC_SONOGRAM:
        lpc_frame      = [self lpcFromAuto:auto_correlation order:lpcOrder];
        k = MIN(spectralLength, (lpcOrder+1));
        for(i=0; i<k; i++)
          spectralData[n_outputs++] = lpc_frame[i];
        break;
    }
  }
  return spectralData;
}

/***********************************
 * This method finds the peak-     *
 * to-peak value of an input array *
 * on a frame by frame basis.      *
 ***********************************/
- (float *)calculatePTPData:(float *)speechInput start:(int)start
           length:(int)length
{
  int    i, j, end;
  int    n;
  static int old_size = 0;
  float  min_peak, max_peak;
    
  end    = start+length;
  numberFrames = (length-frameLength)/(frameLength-overlapLength);
  if(old_size < numberFrames)
  {
     if(ptpData != NULL)
       free(ptpData);
     ptpData  = (float *)malloc(numberFrames*sizeof(float));
     old_size = numberFrames;
  }
  for(n=0; n<numberFrames; n++)
  {
    min_peak = 1.e10;
    max_peak = -1.e10;
    j = start;
    for(i=0; i<frameLength; i++)
    {
      min_peak = min_peak<speechInput[j] ? min_peak : speechInput[j];
      max_peak = max_peak>speechInput[j] ? max_peak : speechInput[j];
      if(++j >= end)
        break;
    }
    start += frameLength - overlapLength;
    ptpData[n] = max_peak - min_peak;
  }
  return ptpData;
}

/***********************************
 * This method finds the zero      *
 * crossing data of an input array *
 * on a frame by frame basis.      *
 ***********************************/
- (float *)calculateZeroCrossData:(float *)speechInput start:(int)start
           length:(int)length
{
  int    j, k, end;
  int    n, crossings;
  static int old_size = 0;
  float  old_sound;
    
  end    = start+length;
  numberFrames = (length-frameLength)/(frameLength-overlapLength);
  if(old_size < numberFrames)
  {
     if(zeroCrossData != NULL)
       free(zeroCrossData);
     zeroCrossData  = (float *)malloc(numberFrames*sizeof(float));
     old_size = numberFrames;
  }
  old_sound = 0.;
  for(n=0; n<numberFrames; n++)
  {
    crossings = 0;
    k = start;
    for(j=0; j<frameLength; j++)
    {
      if(SGN(speechInput[k]) != SGN(old_sound))
        crossings++;
      old_sound = speechInput[k++];
      if(k >= end)
        break;
    }
    start += frameLength - overlapLength;
    zeroCrossData[n]  = (float)crossings/(float)frameLength;
  }
  return zeroCrossData;
}

/***********************************
 * This method finds the energy     *
 * data of an input array           *
 * on a frame by frame basis.       *
 * NOTE: the frame energies are     *
 * also calculated in calculate-    *
 * SpectralData in the frequency    *
 * domain.                          *
 ***********************************/
- (float *)calculateEnergyData:(float *)speechInput start:(int)start
           length:(int)length
{
  int    i, j, end;
  int    n;
  static int old_size = 0;
  float  sum, speech_sample;
    
  end    = start+length;
  numberFrames = (length-frameLength)/(frameLength-overlapLength);
  if(old_size < numberFrames)
  {
     if(energyData != NULL)
       free(energyData);
     energyData  = (float *)malloc(numberFrames*sizeof(float));
     old_size = numberFrames;
  }
  for(n=0; n<numberFrames; n++)
  {
    sum = 0.;
    j = start;
    for(i=0; i<frameLength; i++)
    {
      speech_sample = speechInput[j];
      sum += speech_sample*speech_sample;
      if(++j >= end)
        break;
    }
    start += frameLength - overlapLength;
    energyData[n] = sum;
  }
  return energyData;
}

/************************************
 * This method classifies a set     *
 * of frames as speech or noise     *
 * based on energy data.            *
 ************************************/
- (int *)calculateSpeechOrNoiseFromEnergy:(float *)energyInput numberFrames:(int)frames
{
  int    frame;
  static int ncnt;
  static int frame_count = TIMEOUT+1;
  static int old_size = 0;
  static float noise_average = 0.;
  static float noise_variance = 0.;
  static float reference = 1.e20;
  float  diff;
    
  if(old_size < frames)
  {
     if(speechOrNoiseFromWave != NULL)
       free(speechOrNoiseFromWave);
     speechOrNoiseFromWave  = (int *)malloc(frames*sizeof(float));
     old_size = frames;
  }
  for(frame=0; frame<frames; frame++)
  {

/********************************
 check noise timeout flag; reset*
 ********************************/
    if (frame_count > TIMEOUT)
      frame_count = ncnt = 0;

/********************************
 classify the frame             *
 ********************************/
    if ((ncnt >= NNSUF) && energyInput[frame] > reference)
    {
      speechOrNoiseFromWave[frame] = 1;
      frame_count++;
    }
    else
      speechOrNoiseFromWave[frame] = frame_count = 0;

/********************************
 calculate the mean and variance*
 ********************************/
    if (ncnt < NNSUF)
      ncnt++;

    if (!speechOrNoiseFromWave[frame])
    {
/********************************
 calculate the average frame energy
 ********************************/
      noise_average = TIME_CONSTANT*noise_average + (1.0 - TIME_CONSTANT)*energyInput[frame];

/********************************
 calculate the noise variance   *
 ********************************/
      diff = noise_average - energyInput[frame];
      noise_variance = TIME_CONSTANT*noise_variance + (1.0 - TIME_CONSTANT)*diff*diff;

/********************************
 determine the energy reference *
 ********************************/
      reference = noise_average + NUMBER_STDEV*sqrt(noise_variance);
    }
  }
  return speechOrNoiseFromWave;
}

/*######################################*
 * Window the input speech using the    *
 * previously selected window.  Return  *
 * a pointer to the windowed speech.    *
 *######################################*/
- (float *)windowSpeechData:(float *)speech length:(int)length
{
  int i, window_length;
  static int old_length         = 0;
  static float *windowed_speech = NULL;

  if(old_length < frameLength)
  {
     if(windowed_speech != NULL)
       free(windowed_speech);
     windowed_speech  = (float *)malloc(frameLength*sizeof(float));
     old_length = frameLength;
  }
  window_length = MIN(frameLength, length);
  for (i=0; i<window_length; i++)
    windowed_speech[i] = speech[i]*window[i];
  return windowed_speech;
}

/*######################################*
 * Find the power spectrum of a section *
 * of speech.                           *
 * This is also called the periodogram  *
 * in Oppenheim and Schafer and is      *
 * actually a biased estimate of the    *
 * power spectrum.                      *
 * Also return the frame energy and the *
 * spectral center of mass.             *
 *######################################*/
- (float *)powerSpectrumOf:(float *)speech length:(int)length
               frameEnergy:(float *)frameEnergy spectralCenter:(float *)scm
{
  int i, half_length;
  static int old_length = 0;
  static float *xi = NULL;
  static float *power_spectrum = NULL;
  float  sum, temp;

  if(old_length < length)
  {
     if(xi != NULL)
     {
       free(xi);
       free(power_spectrum);
     }
     xi              = (float *)malloc(length*sizeof(float));
     power_spectrum  = (float *)malloc(length*sizeof(float));
     old_length = length;
  }
  
/*
 calculate the periodogram:
 --------------------------
*/
  for (i=0; i<length; i++)
    power_spectrum[i] = speech[i];
  fft( power_spectrum, xi, length, REAL_FFT);
  half_length = length/2;
  temp = (float)length;
  sum = *scm = 0.0;
  for (i=0; i<=half_length; i++)
  {
    power_spectrum[i] = power_spectrum[i]*power_spectrum[i] + xi[i]*xi[i];
    power_spectrum[i] /= temp;
    sum  += power_spectrum[i];
    *scm += i*power_spectrum[i];
  }
  *frameEnergy = 2.*sum - power_spectrum[0] - power_spectrum[half_length];
  temp = sum - power_spectrum[0];
  if(temp>0.)
    *scm /= half_length*temp;

  return power_spectrum;
}

/*#########################################*
 * This  program finds the power spectrum  *
 * from the autocorrelation function.      *
 * NOTE: length must be a power of 2.      *
 *#########################################*/
- (float *)powerSpectrumFromAuto:(float *)autocorr length:(int)length;
{
  int     i, half_length;
  static  int old_length = 0;
  static  float *spectrum = NULL;
  static  float *imag = NULL;

  if(old_length < length)
  {
     if(spectrum != NULL)
     {
       free(spectrum);
       free(imag);
     }
     spectrum = (float *)malloc(length*sizeof(float));
     imag     = (float *)malloc(length*sizeof(float));
  }
  half_length = length/2;
  spectrum[0] = autocorr[0];
  for(i=1; i<=half_length; i++)
    spectrum[length-i] = spectrum[i] = autocorr[i];
  fft(spectrum, imag, length, REAL_FFT);
  return spectrum;
}

/*######################################*
 * This method finds the energy out of  *
 * the perceptual filters given the *
 * input power spectrum.        *
 *######################################*/
- (float *)filterBankOutputs:(float *)powerSpectrum length:(int)length
{
  int i, k, half_length;
  float  sum, scale, temp;
  static int old_banks = 0;
  static float *filter_outputs = NULL;

  if(old_banks < numberFilterBanks)
  {
     if(filter_outputs != NULL)
       free(filter_outputs);
     filter_outputs  = (float *)malloc((numberFilterBanks+1)*sizeof(float));
     old_banks = numberFilterBanks;
  }
/* 
 Scale the power spectrum per Kushner:
 */
  half_length = length/2;
  temp = length;
  sum = 0.;
  for(i=0; i<=half_length; i++)
     sum += powerSpectrum[i]*temp;
  if(sum>0.)
    scale = temp*temp/sum;
  else
    scale = 1.;
/*
 calculate the bark domain filter energies; weight each filter by the
 number of synched filters; store enhanced spectrum filter energies in q[ ]
 --------------------------------------------------------------------------
*/
  for(i=0; i<=numberFilterBanks; i++)
    filter_outputs[i] = 0.;
  for (i=0; i<=half_length; i++)
    for(k=0; k<=numberFilterBanks; k++)
      filter_outputs[k] +=  filterBanks INDEX_FB(i,k)*scale*powerSpectrum[i];

/*
 perform equal loudness weighting and perform intensity to loudness 
 domain conversion using steven's power law {l = i**(1/3)}.  
 set filter_outputs[0] = filter_outputs[1] per Hermansky.
 ---------------------------------------------------------------------
*/
  for (k=0; k<=numberFilterBanks; k++)
  {
    filter_outputs[k] *= loudnessWeighting[k];
    filter_outputs[k] = pow(filter_outputs[k],STEVENS_POWER);
  }
  filter_outputs[0] = filter_outputs[1];
  
  return filter_outputs;
}

/*######################################*
 * This method finds the autocorrelation*
 * function from the filter bank        *
 * energies.                            *
 *######################################*/
- (float *)autoFromFilterBanks:(float *)filterEnergies length:(int)length
{
  int i, j, k, half_length;
  static int old_length = 0;
  static float *auto_correlation = NULL;
  static float *xi = NULL;
  float  domega, omega;

  if(old_length < length)
  {
     if(auto_correlation != NULL)
     {
       free(auto_correlation);
       free(xi);
     }
     auto_correlation = (float *)malloc(length*sizeof(float));
     xi               = (float *)malloc(length*sizeof(float));
     old_length = length;
  }
/*
 interpolate warped spectrum to get a new uniform bark-based energy spectrum
 load real values symetrically for an inverse fft; load hz pwr spec into
 the imaginary array (NOTE: this is no longer done);
 this trick gets you 2 auto-corrs with 1 fft
 ---------------------------------------------------------------------------
*/
  half_length = length/2;
  domega = centerFrequencies[numberFilterBanks]/(float)half_length;
  k                   = 0;
  omega               = 0.0;
  auto_correlation[0] = filterEnergies[0];
  xi[0]               = 0.;

  for (j=1; j<=half_length; j++)
  {
    if (omega >= centerFrequencies[k+1])
      k++;
    i = length-j;
    auto_correlation[i] = auto_correlation[j] = filterEnergies[k] +
             ((filterEnergies[k+1]-filterEnergies[k])/
             (centerFrequencies[k+1]-centerFrequencies[k]))*
             (omega-centerFrequencies[k]);
    xi[i] = xi[j] = 0.;
    omega += domega;
  }

/*
 take the inverse fft of power spectra to get the
 autocorrelation function;
*/
  fft( auto_correlation, xi, length, -1);

/*
 normalize the pfb auto-corr
 ---------------------------
*/
  if(auto_correlation[0] > 0.)
    for (i=1; i<=half_length; i++)
      auto_correlation[i] /= auto_correlation[0];
  auto_correlation[0] = 1.;
  return auto_correlation;
}

/*######################################*
 * This  method finds the  short-time   *
 * autocorrelation function from (4.30) *
 * page 143, Rabiner and Schafer. A     *
 * pointer to the output autocorrelation*
 * is returned.                         *
 *######################################*/
- (float *)autoFromSpeech:(float *)speech length:(int)length
                                    lagLength:(int)lagLength
{
  int     k, m;
  static int old_length;
  static float *rn;
  double sum;

  if(old_length < lagLength)
  {
     if(rn != NULL)
       free(rn);
     rn = (float *)malloc(lagLength*sizeof(float));
     old_length = lagLength;
  }
        
/****************************
 * Now, calculate it:       *  
 * choose n = 0  in (4.30)  *
 ****************************/
  if(length > lagLength)
    for(k=0; k<lagLength; k++)
    { 
      sum = 0.;
      for(m=0; m<(length-k); m++)
        sum +=  speech[m] * speech[m+k];
      rn[k] = sum/length;
    }

  return rn;
}

/*#########################################*
 * This  method finds the LPC coefficients *
 * from the short-time autocorrelation     *
 * function using Durbin's algorithm,      *
 * eqs. (8.67)->(8.71), page 411, Rabiner  *
 * and Schafer. The gain is returned in the*
 * zeroth array position, the rest of the  *
 * coefficients are returned in positions  *
 * 1->lpcOrder                             *
 ****NOTE: The alphas found here have the  *
 * opposite sign of those in Rabiner and   *
 * Schafer, i.e.:                          *
 *                                         *
 *              gain                       *
 * H(z) =  ----------------                *
 *         1 + alpha(1)*z^-1 + ...         *
 *                                         *
 *#########################################*/
- (float *)lpcFromAuto:(float *)autoCorrelation order:(int)order
{
  int    i, ii, j, lpc_length;
  static int old_order;
  static float *lpc, *k, *alpha_old, *rn_normalized;
  float  gain, error;

  lpc_length = order+1;
  if(old_order < order)
  {
     if(lpc != NULL)
     {
       free(lpc);
       free(k);
       free(alpha_old);
       free(rn_normalized);
     }
     lpc           = (float *)malloc(lpc_length*sizeof(float));
     k             = (float *)malloc(lpc_length*sizeof(float));
     alpha_old     = (float *)malloc(lpc_length*sizeof(float));
     rn_normalized = (float *)malloc(lpc_length*sizeof(float));
     old_order = order;
  }
        
/****************************
 * Initialize:              *
 ****************************/
  for(ii=1; ii<=order; ii++)
    rn_normalized[ii] = (double)autoCorrelation[ii]/(double)autoCorrelation[0];
  error = 1.;

/****************************
 * Loop over filter order:  *
 ****************************/
  for(i=1; i<=order; i++)
  {
    k[i] = 0.;
    for(j=1; j<i; j++)
      k[i] += alpha_old[j]*rn_normalized[i-j];
    lpc[i]  = k[i] = -(rn_normalized[i]+k[i])/error;
    error *= (1. - k[i]*k[i]);
    for(j=1; j<i; j++)
      lpc[j] = alpha_old[j] + k[i]*alpha_old[i-j];
    for(j=1; j<=i; j++)
      alpha_old[j] = lpc[j];

  }     

/****************************
 * Calculate the gain from  *
 * Eqn 8.38:                *
 ****************************/
  gain = 1.;
  for(i=1; i<=order; i++)
    gain -= lpc[i]*autoCorrelation[i];
  gain = sqrt(gain);

/****************************
 * Save gain in alpha[0]    *
 ****************************/
  lpc[0]   = gain;

  return lpc;
}

/*#########################################*
 * This  program finds the cepstral        *
 * coefficients from the LPC coefficients. *
 *                                         *
 * References:                             *
 * M. Basseville, "Distance measures for   *
 * signal proc. and pattern recognition,"  *
 * Signal Processing, 18, pp. 349-369, 1989*
 *                                         *
 * NOTE: on output: cepstral[0] = 0,       *
 * cepstral[1] ->cepstral[cepstral_length] *
 * are valid                               *
 *#########################################*/
- (float *)cepstralFromLPC:(float *)lpc lpcOrder:(int)order 
                    length:(int)length;
{
  int    n, k, lpc_length, len;
  static int old_length;
  static float *cepstral;
  float  sum;

  lpc_length = order+1;
  if(old_length < length)
  {
     if(cepstral != NULL)
       free(cepstral);
     cepstral = (float *)malloc((length+1)*sizeof(float));
     old_length = length;
  }
  
/**************************************
 * Calculate Cepstral Coefficients,   *
 * Use (77) in Basseville:            *
 **************************************/
  cepstral[0] = 0.;
  cepstral[1] = -lpc[1];
  len = MIN(order, length);
  for (n=2; n<=len; n++)
  {
    sum = 0.0;
    for(k=1; k<n; k++)
      sum += k*cepstral[k]*lpc[n-k];
    cepstral[n] = -(lpc[n] + sum/(double)n);
  }
  for(n=lpc_length; n<=length; n++)
  {
    sum = 0.;
    for(k=1; k<=order; k++)
      sum -= (n-k)*cepstral[n-k]*lpc[k];
    cepstral[n] = sum/(double)n;
  }
  
  return cepstral;
}

/*#########################################*
 * This method finds the cepstral          *
 * coefficients from the LPC coefficients. *
 * The LPC coefficients are obtained from  *
 * the autocorrelation function.           *
 *                                         *
 * References:                             *
 * M. Basseville, "Distance measures for   *
 * signal proc. and pattern recognition,"  *
 * Signal Processing, 18, pp. 349-369, 1989*
 *                                         *
 * NOTE: on output: cepstral[0] = 0,       *
 * cepstral[1] ->cepstral[cepstral_length] *
 * are valid                               *
 *#########################################*/
- (float *)cepstralFromAuto:(float *)autoCorrelation lpcOrder:(int)order
                    length:(int)outputLength;
{
  float *lpc;

  lpc = [self lpcFromAuto:autoCorrelation order:order];
  return [self cepstralFromLPC:lpc lpcOrder:order length:outputLength];
}


/*###############################################*
 * This  program finds the cepstral coefficients *
 * from the inverse transform of the logarithm   *
 * the magnitude of the Fourier transform of the *
 * of the speech, see Rabiner and Schafer,       *
 * Fig. 7.7.                                     *
 *                                               *
 * References:                                   *
 *  1. A.M. Noll, "Cepstrum pitch determination,"*
 * J. Acoust. Soc. Am., pp. 293-309, Feb. 1967.  *
 * Reprinted in IEEE Speech Analysis.            *
 *###############################################*/
- (float *)cepstralFromSpectrum:(float *)spectrum length:(int)outputLength
{
  int     i, len, window_length, cepstral_length, half_length;
  static  int old_clength = 0;
  static  int old_wlength = 0;
  float   temp;
  static  float *cepstral = NULL;
  static  float *real = NULL;

  if(old_clength < outputLength)
  {
     if(cepstral != NULL)
       free(cepstral);
     cepstral = (float *)malloc((outputLength+1)*sizeof(float));
     old_clength = outputLength;
  }
  if(old_wlength < powerOf2Length)
  {
     if(real != NULL)
       free(real);
     real = (float *)malloc(powerOf2Length*sizeof(float));
     old_wlength = powerOf2Length;
  }
  window_length = powerOf2Length;
  cepstral_length = outputLength;
  
/**************************************
 * Now inverse transform the log      *
 * magnitude function:                *
 **************************************/
  temp = (float)window_length;
  half_length = window_length/2;
  real[0] = log(temp*spectrum[0]);
  for(i=1; i<=half_length; i++)
    real[window_length-i] = real[i] = log(temp*spectrum[i]);
  len = MIN(cepstral_length+1, window_length);
  idct(real, cepstral, window_length, len);  
  cepstral[0] = 0.;
  for(i=1; i<len; i++)
    cepstral[i] /= 2.;
  
  return cepstral;
}


/*#########################################*
 * This method finds the magnitude res-    *
 * ponse of the input LPC filter.          *
 *#########################################*/
- (float *)lpcMagnitudeResponse:(float *)alpha lpcOrder:(int)order
                    digitalFrequencies:(float *)omega numberPoints:(int)length;
{
  int    i, n, lpc_length;
  static int old_length = 0;
  static float  *response = NULL;
  double  temp;
  COMPLEX den,exponent,ctemp,carg, cprod, csum;

  if(old_length < length)
  {
     if(response != NULL)
       free(response);
     response = (float *)malloc((length)*sizeof(float));
     old_length = length;
  }
  
  lpc_length = order + 1;
  for(n=0; n<length; n++)
  {
    CMPLX(exponent,0.,-omega[n])
    CEXP(carg,exponent)

/****************************
 * Find denominator value   *
 * use Horner's method:     *
 ****************************/ 
    CMPLX(ctemp,alpha[lpc_length-1],0.);
    CMULT(cprod,ctemp,carg);
    for(i=(lpc_length-2); i>0; i--)
    {
     CMPLX(ctemp,alpha[i],0.);
     CADD(csum,ctemp,cprod);
     CMULT(cprod,csum,carg); 
    }
    CMPLX(ctemp,1.,0.);
    CADD(den,cprod,ctemp);

    temp = sqrt(den.x*den.x + den.y*den.y);
    if(temp > 0.)
      response[n] = 1./temp;
  }
  return(response);
}
/********************************************************************************
*
* program:  snclass_FIL.c
*
* programmer:   william m. kushner
*
* purpose:  to generate mean and variance codebook vectors for noise
*               "on-the-fly". This code produces a filtered running average
*               based on the contents of the sample buffer.
*
********************************************************************************
* rev    date    by   reason
* ---  --------  ---  ----------------------------------------------------------
* org  02/21/93  wmk  original C version
*
********************************************************************************/
- (int)speechOrNoise:(float)frameEnergy frame:(int)frame;
{
  int           speech_or_noise;
  static int    ncnt;
  static int   frame_count = TIMEOUT+1;
  static float noise_average = 0.;
  static float noise_variance = 0.;
  static float reference = 1.e20;
  float        diff;

/********************************
 check noise timeout flag; reset*
 ********************************/
  if (frame_count > TIMEOUT)
    frame_count = ncnt = 0;

/********************************
 classify the frame             *
 ********************************/
  if ((ncnt >= NNSUF) && frameEnergy > reference)
  {
    speech_or_noise = 1;
    frame_count++;
  }
  else
    speech_or_noise = frame_count = 0;

/********************************
 calculate the mean and variance*
 ********************************/
  if (ncnt < NNSUF)
   ncnt++;

  if (!speech_or_noise)
  {
/********************************
 calculate the average frame energy
 ********************************/
    noise_average = TIME_CONSTANT*noise_average + (1.0 - TIME_CONSTANT)*frameEnergy;

/********************************
 calculate the noise variance   *
 ********************************/
    diff = noise_average - frameEnergy;
    noise_variance = TIME_CONSTANT*noise_variance + (1.0 - TIME_CONSTANT)*diff*diff;

/********************************
 determine the energy reference *
 ********************************/
    reference = noise_average + NUMBER_STDEV*sqrt(noise_variance);
  }

  return speech_or_noise;
}

/*#########################################*
 * This method subtracts an estimate of the*
 * noise spectrum from the speech power    *
 * spectrum.                               *
 *#########################################*/
- subtractNoiseFrom:(float *)spectrum noise:(float *)noiseSpectrum
          speechFlag:(int)speech length:(int)length
          frameEnergy:(float *)frameEnergy;
{
  int i, half_length;
  float ne_temp, diff, noise_floor;
  float sum;

  half_length = length/2;
  if(!speech)
  {
    for(i=0; i<=half_length; i++)
      if(spectrum[i]>=0.)
        noiseSpectrum[i] = TIME_CONSTANT*noiseSpectrum[i] + 
                (1.0 - TIME_CONSTANT)*sqrt(sqrt(spectrum[i]))*noiseFactor;
  }
/************************
 * perform subtraction  *
 * with whitening.  *
 ************************/
  ne_temp = sum = 0.;
  for(i=0; i<=half_length; i++)
  {
    if(spectrum[i]>=0.)
      ne_temp     = sqrt(sqrt(spectrum[i]));
    diff        = ne_temp - noiseSpectrum[i];
    noise_floor = spectralFloor*noiseSpectrum[i];
    spectrum[i] = (diff > noise_floor) ? diff : noise_floor;
    spectrum[i] *= spectrum[i];
    spectrum[i] *= spectrum[i];
    sum += spectrum[i];
  }
  *frameEnergy = sum + sum - spectrum[0] - spectrum[half_length];
  return self;
}
 

@end