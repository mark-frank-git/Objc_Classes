/************************************************************************
 * This object is a subclass of the Sound object.  In addition to the   *
 * to the normal sound methods, it provides for processing of the sound.*
 *                                                                      *
 * File: /User/frank/Objc_Classes/ProcessedSound/ProcessedSound.m       *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 06/17/92 - Started                                               *
 *  2. 07/01/92 - Added NIST capability                                 *
 *  3. 09/22/92 - Added DigitalFilter object                            *
 *  4. 11/17/92 - Add numberBytes for output data                       *
 *  5. 12/07/92 - if findEndpoints = null, speech = total sound?        *
 *  6. 02/04/93 - Add extendEndpointsBy:                                *
 *  7. 04/16/93 - Move savePrcFile to NoisySound                        *
 *  8. 04/26/93 - Add SpeechProcessor object                            *
 *                Delete strobe stuff.                                  *
 *  9. 05/28/93 - Combined maxs and mins into arrays.                   *
 * 10. 06/12/93 - Add SNDSwapSoundToHost()                              *
 ************************************************************************/

#import "ProcessedSound.h"
#import "SpeechProcessor.h"
#import "DigitalFilter.h"
#import "PlotSoundView.h"
#import "speechParams.h"
#import "c_headers.h"
#import <appkit/Control.h>
#import <soundkit/soundkit.h>
#import <dpsclient/wraps.h>
#import <string.h>
#import <math.h>
#import <time.h>
#import <ctype.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

/***************************************
 * The following define the processing *
 * of the sub vectors in the feature   *
 * vector.  The order of these must    *
 * be consistent.              *
 * NOTE: If these change, make the     *
 * same changes in NoisySound.m        *
 ***************************************/
char *sub_names[SUB_VECTORS] = {"cepstral", "zero cross", "ptp",
                                "frame energy", "speech or noise"};

/***************************************
 * If forced_range[i] = YES then force *
 * the range of the subvector, e.g., if*
 * forced_range[0] = YES and           *
 * minVector[0] = -0.8, maxVector[0]   *
 * = 0.8, then any cepstral coefficient*
 * >= 0.8 gets mapped to 1.0, and any  *
 * cepstral coefficient <= -0.8 gets   *
 * mapped to -1.0                      *
 ***************************************/
BOOL forced_range[SUB_VECTORS] = {YES, YES, NO, NO, YES};

/***************************************
 * If linear_map[i] = YES, then map    *
 * using a linear mapping to [-1, 1],  *
 * else, map using an affine mapping.  *
 ***************************************/
BOOL linear_mapping[SUB_VECTORS] = {YES, YES, YES, NO, YES};

/*************************************
 * These strings must match radio    *
 * button definitions in plotParams.h*
 * Also note that the cepstral       *
 * coefficients are assumed to be    *
 * first in the vector array by many *
 * methods.                          *
 *************************************/
char *frequency_type_names[] = {"Perceptual", "Linear"};

char *window_type_names[] = {"Hamming", "Hanning", "Rectangular", "Triangular",
                             "Blackman", "Blackman_Harris"};
char *lifter_type_names[] = {"No_liftering", "Raised_Sine", "Linear",
                                           "Gaussian", "Exponential"};

/*************************************
 * Constants:                        *
 *************************************/
#define MAX_FILE    120     /* Maximum file name length */

@implementation ProcessedSound

/*###################################*
 * Initialize the instance variables *
 *###################################*/
- initWithSoundWindow: (id)newSoundWindow
{
  int i;
  
  [super init];
  soundWindow = newSoundWindow;
  digitalFilter     = [[DigitalFilter alloc] init];
  speechProcessor   = [[SpeechProcessor alloc] initWithOrder:LPC_DEFAULT
                           cepstralLength:LPC_DEFAULT frameLength:FRAME_LENGTH
                           overlapLength:OVERLAP_LENGTH 
                           numberBanks:FILTER_BANKS_DEFAULT
                           samplingFrequency:SND_RATE_CODEC];
// Initialize instance variables
  numberBytes               = ELEMENT_BYTES;
  intScaleFactor            = ONE_BYTE_SCALE;   /* must match above */
  inputSoundType            = -1;
  numberSoundSamples        = 0;
  preemphasisCoefficient    = PREEMPHASIS_COEFFICIENT;
  minVector[0]              = MIN_CEPSTRAL;
  maxVector[0]              = MAX_CEPSTRAL;
  minVector[1]              = MIN_ZERO_CROSS;
  maxVector[1]              = MAX_ZERO_CROSS;
  minVector[2]              = MIN_PTP;
  maxVector[2]              = MAX_PTP;
  minVector[3]              = MIN_FRAME_ENERGY;
  maxVector[3]              = MAX_FRAME_ENERGY;
  floatVectorsNeeded        = YES;
  selectionStart            = selectionLength = featureVectorLength = 0;
  rawSoundData              = NULL;
  soundData                 = NULL;
  nistHeader                = NULL;
  charFeatureVectors        = NULL;
  intFeatureVectors         = NULL;
  floatFeatureVectors       = NULL;
  outputSamplingRate        = 0.;
  forcedRange               = forced_range;
  linearMapping             = linear_mapping;
  
  for(i=0; i<SUB_VECTORS; i++)
  {
    subvectorCount[i]  = 1;
    subvectorMin[i]    = -100;
    subvectorMax[i]    = 100;
    subvectorScale[i]  = 1.;
    subvectorOffset[i] = 0.;
  }
  subvectorCount[0] = [speechProcessor cepstralLength];
   
  return self;
}

/*#############################*
 * These methods respond to    *
 * other objects request for   *
 * information:                *
 *#############################*/
- (id)speechProcessor { return speechProcessor; }
- (int)getInputSoundType { return inputSoundType; }
- (int)getVectorLength { return featureVectorLength; }
- (int)getNumberBytes { return numberBytes; }
- (int)getNumberFrames
{
  if(inputSoundType==SPEECH_FROM_PRC)
    return numberPRCFrames;
  else
    return [speechProcessor numberFrames];
}
- (int)getSelectionLength { return selectionLength; }
- (int)getSelectionStart { return selectionStart; }
- (int)getFrameLength { return [speechProcessor frameLength]; }
- (int)getOverlapLength { return [speechProcessor overlapLength]; }
- (int)getCepstralLength { return [speechProcessor cepstralLength]; }
- (int)getLPCOrder { return [speechProcessor lpcOrder]; }
- (int)getFilterBanks { return [speechProcessor numberFilterBanks]; }
- (int)getIntScaleFactor { return intScaleFactor; }
- (int *)getIntFeatureVector { return intFeatureVectors; }
- (float) getMinVector:(int)index { return minVector[MIN(index,SUB_VECTORS)]; }
- (float) getMaxVector:(int)index { return maxVector[MIN(index,SUB_VECTORS)]; }
- (float) getSamplingFrequency { return [self samplingRate]; }
- (float) getPreemphasisCoefficient { return preemphasisCoefficient; }
- (float *)getSoundData { return soundData; }
- (float *)getFloatFeatureVectors { return floatFeatureVectors; }
- (char *)getUtteranceId { return utteranceId; }
- (char *) windowName{ return [speechProcessor windowName];}
- (char *) spectralName{  return [speechProcessor spectralName];}
- (char *) lifterName{  return [speechProcessor lifterName];}
- (char *) frequencyAnalysisName
           { return [speechProcessor frequencyAnalysisName];}


/*#####################################*
 * The following methods set instance   *
 * variables.               *
 *#####################################*/
- setMinVector:(int)index to:(float)value
{
  minVector[index] = value;
  return self;
}

- setMaxVector:(int)index to:(float)value
{
  maxVector[index] = value;
  return self;
}

- setNumberBytes:(int)number
{
  numberBytes = number;
  if(numberBytes == 1)
    intScaleFactor = ONE_BYTE_SCALE;
  else
    intScaleFactor = TWO_BYTE_SCALE;
  return self;
}

- setPreemphasis:(float)coefficient
{
    int sound_format;
    
    sound_format           = [self dataFormat];
    preemphasisCoefficient = coefficient;
    if(rawSoundData != NULL)
      soundData = [self rawSoundToFloat];
    return self;
}

- setFrequencyType:(int)type
{
  [speechProcessor setPerceptualProcessing:type];
  return self;
}

- setLPCOrder:(int)order
{
  [speechProcessor setLPCFilterOrder:order];
  return self;
}

- setLiftering:(int)type
{
  [speechProcessor setLifter:type];
  return self;
}

- setWindowType:(int)type
{
  [speechProcessor setWindow:type];
  return self;
}

- setProcessType:(int)type
{
  [speechProcessor setSpectralType:type];
  return self;
}

- setFilterBanks:(int)number
{
  [speechProcessor setNumberFilterBanks:number];
    return self;
}

- setFrameLength:(int)length
{
  [speechProcessor setSpeechFrameLength:length];
   return self;
}

- setOverlapLength:(int)length
{
  [speechProcessor setOverlapLength:length];
  return self;
}

- setCepstralLength:(int)length
{
  [speechProcessor setCepLength:length];
  return self;
}

- setInputSoundType:(int)new_type
{
  inputSoundType = new_type;
  return self;
}

- setFloatVectorsNeeded: (BOOL)flag
{
  floatVectorsNeeded = flag;
  return self;
}

/*###############################*
 *  Read the sound file:         *
 *###############################*/
-(BOOL) readSndFile:(const char *)fileName
{
   
// read in new sound:
   if( [self readSoundfile:(char *)fileName] != SND_ERR_NONE)
    return NO;
   strncpy(inputSpeechFile, fileName, STRING_LEN);
   [self setDataFromSoundStructure];
   inputSoundType = SPEECH_FROM_SND;
   [speechProcessor setSamplingFrequency:[self samplingRate]];
// set sound in the plotting window:
   if(soundWindow != nil)
     [soundWindow setSound:self];
// Set selection to total sound,
   [self setSelectionRange: 0: numberSoundSamples];
       
   return YES;
}

/*###############################*
 *  Read the sound file:         *
 *  This is the same as the      *
 *  method above, except we set  *
 *  inputSoundType = SPEECH_FROM *
 *  _MIKE.                       *
 *###############################*/
-(BOOL) readMikeFile:(const char *)fileName
{
   
// read in new sound:
   if( [self readSoundfile:(char *)fileName] != SND_ERR_NONE)
    return NO;

   strncpy(inputSpeechFile, fileName, STRING_LEN);
   [self setDataFromSoundStructure];
   inputSoundType = SPEECH_FROM_MIKE;
   [speechProcessor setSamplingFrequency:[self samplingRate]];
// set sound in the plotting window:
   if(soundWindow != nil)
     [soundWindow setSound:self];
// Set selection to total sound,
   [self setSelectionRange: 0: numberSoundSamples];
       
   return YES;
}

/*###############################*
 *  Read in the processed        *
 *  sound file in NIST format    *
 *  Return YES if read successful*
 *###############################*/
-(BOOL) readPrcFile:(const char *)fileName
{
  char *err_msg;
  FILE *file_ptr;
  int  data_size, i, j, k;
  char  *char_input_data;
  short *short_input_data;
  
// read in new processed data:
  if( (file_ptr = fopen(fileName, "r")) == NULL)
    return NO;
  if(nistHeader != NULL)
    sp_close_header(nistHeader);

  if( (nistHeader = sp_open_header(file_ptr, 1, &err_msg))== HDRNULL)
    return NO;

  strncpy(inputSpeechFile, fileName, STRING_LEN);
  inputSoundType = SPEECH_FROM_PRC;
// Get data from NIST header:
  [self getNISTHeaderData];

// Read in feature vector data:
  if(charFeatureVectors != NULL)
    free(charFeatureVectors);
  if(intFeatureVectors != NULL)
    free(intFeatureVectors);
  if(floatFeatureVectors != NULL)
    free(floatFeatureVectors);

  data_size     = numberPRCFrames*featureVectorLength;
  charFeatureVectors    = get_char_spc(numberBytes*data_size);
  intFeatureVectors = get_int_spc(data_size);
  if(floatVectorsNeeded)
    floatFeatureVectors = get_flt_spc(data_size);
  if(fread(charFeatureVectors, sizeof(char), numberBytes*data_size, file_ptr) )
  {
    j = 0;
    switch(numberBytes)
    {
      case 1:
        char_input_data = charFeatureVectors;
        for(i=0; i<numberPRCFrames; i++)
        {
          for(k=0; k<featureVectorLength; k++)
          {
            intFeatureVectors[j] = char_input_data[j];
            if(floatVectorsNeeded)
              floatFeatureVectors[j] = char_input_data[j]/ONE_BYTE_SCALE;
            j++;
          }
        }
        break;
      case 2:
        short_input_data = (short *)charFeatureVectors;
        for(i=0; i<numberPRCFrames; i++)
        {
          for(k=0; k<featureVectorLength; k++)
          {
            intFeatureVectors[j] = short_input_data[j];
            if(floatVectorsNeeded)
              floatFeatureVectors[j] = short_input_data[j]/TWO_BYTE_SCALE;
            j++;
          }
        }
        break;
      default:
        printf("Illegal # of bytes in readPrcFile\n");
        break;
    }
    fclose(file_ptr);
    return YES;
  }
  return NO;
}

/*###############################*
 * Set the instance variables   *
 * from the sound structure:    *
 *##############################*/
- setDataFromSoundStructure
{
   int    sound_format;
   double sampling_frequency;

// Init data pointer, sample size, sampling frequency
   sound_format         = [self dataFormat];
   sampling_frequency   = [self samplingRate];
   if( (sound_format==SND_FORMAT_MULAW_8) || 
       (sound_format==SND_FORMAT_MULAW_SQUELCH) )
   {
     sound_format = SND_FORMAT_LINEAR_16;
     [self convertToFormat:sound_format
                       samplingRate:(double)sampling_frequency channelCount:1];
   }
   numberSoundSamples= [self sampleCount];
   rawSoundData      = [self data];
// Swap bytes if necessary:
   SNDSwapSoundToHost(rawSoundData, rawSoundData, numberSoundSamples, [self channelCount],
                      sound_format);
   
// Convert to floating point sound, result is stored in soundData:
   soundData = [self rawSoundToFloat];
  
   return self;
}

/************************************
 * Get the NIST header information  *
 * into local variables:            *
 ************************************/
- getNISTHeaderData
{
  int i;
  
  for(i=0; i<SUB_VECTORS; i++)
    [self getSubVectorMinMax: 0: &minPRC[i]: &maxPRC[i]];
  
  [self getNISTInteger: "element_n_bytes":      &numberBytes];
  [self getNISTInteger: "vector_element_count":     &featureVectorLength];
  [self getNISTInteger: "number_frames":        &numberPRCFrames];
  [self getNISTFloat  : "preemphasis_coefficient":  
                                                   &preemphasisCoefficient];
  [self getNISTChar   : "utterance_id":         utteranceId];

  if(numberBytes == 1)
    intScaleFactor = ONE_BYTE_SCALE;
  else
    intScaleFactor = TWO_BYTE_SCALE;
  
  return self;
}

/*********************************
 * Get the min and max values    *
 * for a subvector:              *
 *********************************/
- getSubVectorMinMax: (int) vectorNumber: (int *)min: (int *)max
{
  int length;
  char subvector_min_max[40];
  
  sprintf(subvector_min_max,"subvector_%1d_min",vectorNumber);
  length = sp_get_size(nistHeader, subvector_min_max);
  sp_get_data(nistHeader, subvector_min_max, min, &length);
  sprintf(subvector_min_max,"subvector_%1d_max",vectorNumber);
  sp_get_data(nistHeader, subvector_min_max, max, &length);
  
  return self;
}

/*********************************
 * Get an integer variable from  *
 * the NIST header:              *
 *********************************/
- getNISTInteger: (char *)nistName: (int *)data
{
  int length;
  length = sizeof(int);
  sp_get_data(nistHeader, nistName, data, &length);
  return self;
}

/*********************************
 * Get a float    variable from  *
 * the NIST header:              *
 *********************************/
- getNISTFloat: (char *)nistName: (float *)data
{
  int length;
  double d_data;
  
  length = sp_get_size(nistHeader, nistName);
  sp_get_data(nistHeader, nistName, &d_data, &length);
  *data = (float)d_data;
  return self;
}

/*********************************
 * Get a char     variable from  *
 * the NIST header:              *
 *********************************/
- getNISTChar: (char *)nistName: (char *)data
{
  int length;
  length = strlen(data);
  while(length--)
    data[length] = '\0';
  length = sp_get_size(nistHeader, nistName);
  sp_get_data(nistHeader, nistName, data, &length);
  return self;
}

/*###############################################*
 * This method finds the endpoints of the speech *
 * waveform, and sets the selection range ac-    *
 * cordingly.                                    *
 *###############################################*/
- findEndpoints
{
    int    start, length, end;
    int    speech_start, speech_end;
    double sampling_frequency;
    
// Get current selection
    [self getSelectionRange];
    sampling_frequency = [self samplingRate];
    if( (numberSoundSamples>0) && (sampling_frequency>0.) && 
        (selectionLength>0) )
    {
     start  = selectionStart;
     length = selectionLength;  
     end    = start+length; 
     endpoint_detection( &soundData[start], length, sampling_frequency,
                              &speech_start, &speech_end);
     speech_start += start;
     speech_end   += start;
     if(speech_end > speech_start)
       [self setSelectionRange: speech_start: (speech_end-speech_start)];
     [self getSelectionRange];  
    }
       
    return self;
}

/*###############################################*
 * This method sets the current sound selection  *
 * to the entire range:                          *
 *###############################################*/
- resetEndpoints
{

// Set selection to total sound,
   [self setSelectionRange: 0: numberSoundSamples];
       
   return self;
}

/*###############################################*
 * This method sets the current sound selection  *
 * to the entire range:                          *
 *###############################################*/
- extendEndpointsBy:(int)samples
{
  int start, end, length;

// Get current selection
  [self getSelectionRange];
  start  = selectionStart;
  length = selectionLength;
  end    = start + length + samples;
  end    = MIN(numberSoundSamples-1, end);
  start  = MAX(0, start-samples);
  if(end > start)
       [self setSelectionRange: start: (end-start)];
  [self getSelectionRange];  
       
  return self;
}
  
/*#####################################*
 * This method gets the selection      *
 * range from the sound subview:       *
 * If selection range has changed      *
 * since last call, return YES         *
 *#####################################*/
- (BOOL) getSelectionRange
{
  static int old_start = 0;
  static int old_length = 0;
  double input_rate;
  
  switch(inputSoundType)
  {
    case SPEECH_FROM_PRC:
      return NO;
    case SPEECH_FROM_WAV:
      if(outputSamplingRate<=0.)
        return NO;
      input_rate = [self samplingRate];
      [soundWindow getSelection: &selectionStart: &selectionLength];
      selectionStart  = ROUND(selectionStart*input_rate/outputSamplingRate);
      selectionLength = ROUND(selectionLength*input_rate/outputSamplingRate);
      break;
    default:
      [soundWindow getSelection: &selectionStart: &selectionLength];
      break;
  }
  selectionStart = (selectionStart>=0) ? selectionStart : 0;
  if( (selectionStart+selectionLength) > numberSoundSamples)
     selectionLength = numberSoundSamples - selectionStart;
  if( (old_start!=selectionStart) || (old_length!=selectionLength) )
  {
    old_start = selectionStart;
    old_length = selectionLength;
    return YES;
  }
  return NO;
}

/*#####################################*
 * This method sets the selection      *
 * range for  the sound subview:       *
 *#####################################*/
- setSelectionRange: (int)start: (int)length
{
  double input_rate;
  
  if(inputSoundType == SPEECH_FROM_WAV)
  {
    input_rate = [self samplingRate];
    if(input_rate>0.)
    {
      start  = ROUND(start*outputSamplingRate/input_rate);
      length = ROUND(length*outputSamplingRate/input_rate);
    }
  }
  selectionStart = start;
  selectionLength = length;
  [soundWindow setSelection: start: length];

  return self;
}

/*######################################*
 * Find feature vector on a frame-by-   *
 * frame basis:                         *
 * NOTE: speechProcessor does the work  *
 * The following sub vectors are        *
 * currently calculated: spectral, zero-*
 * cross, PTP, log(frame_energy), speech*
 * or noise.                            *
 *######################################*/
- (BOOL)calculateFloatFeatureVectors
{
  int frame_length, overlap_length;
  int   n, i, j, k;
  int   number_frames, spectral_length, length;
  int   *speech_or_noise;
  static int old_length = 0;
  float *sub_vector[SUB_VECTORS-1];
  float temp, temp1;
         
  if( (numberSoundSamples<=0)  || (inputSoundType == SPEECH_FROM_PRC) )
    return NO;
  frame_length   = [speechProcessor frameLength];
  overlap_length = [speechProcessor overlapLength];
  if( (selectionLength<frame_length) || (frame_length<=overlap_length) )
  {
    printf("Frame length>data length or frame<overlap in processSound\n");
    return NO;
  }

  sub_vector[0]     = [speechProcessor calculateSpectralData:soundData
                             start:selectionStart length:selectionLength];
  sub_vector[1]     = [speechProcessor calculateZeroCrossData:soundData
                             start:selectionStart length:selectionLength];
  sub_vector[2]     = [speechProcessor calculatePTPData:soundData
                             start:selectionStart length:selectionLength];
  sub_vector[3]     = [speechProcessor energyData];
  speech_or_noise   = [speechProcessor speechOrNoise];
  number_frames     = [speechProcessor numberFrames];
  spectral_length   = [speechProcessor spectralLength];
  featureVectorLength = spectral_length + SUB_VECTORS - 1;
     
// Allocate space if needed:
  length = featureVectorLength*number_frames;
  if(old_length < length)
  {
    if(floatFeatureVectors != NULL)
      free(floatFeatureVectors);
    floatFeatureVectors = get_flt_spc(length);
    old_length = length;
  }

// Take log of frame energies
  for(n=0; n<number_frames; n++)
    if(sub_vector[3][n] > 0.)
       sub_vector[3][n] = log10(sub_vector[3][n]);
               
// Find the max and min values if range is not forced:
  for(i=0; i<(SUB_VECTORS-1); i++)
  {
    if(!forced_range[i])
    {
      if(i==0)
        length = spectral_length*number_frames;
      else
        length = number_frames;
      float_min_max(sub_vector[i], length, &minVector[i], &maxVector[i],
                    &k, &j);
    }
  }
/**************************************
 * Find a scale factor to map feature *
 * vectors to [-1,1] or [0,1].  For   *
 * frame_energies use affine map.     *
 **************************************/
  for(i=0; i<(SUB_VECTORS-1); i++)
  {
    if(i==0)
      length = spectral_length*number_frames;
    else
      length = number_frames;
    if(linear_mapping[i])
    {
      linear_map(sub_vector[i], length, minVector[i], maxVector[i], 1.,
                 &temp);
      subvectorScale[i]  = temp;
      subvectorOffset[i] = 0.;
    }
    else
    {
      affine_map(sub_vector[i], length, minVector[i], maxVector[i], 0., 1.,
                 &temp, &temp1);
      subvectorScale[i]  = temp;
      subvectorOffset[i] = temp1;
    }
  }

  j = k = 0;
  for(n=0; n<number_frames; n++)
  {
    for(i=0; i<spectral_length; i++)
      floatFeatureVectors[k++] = sub_vector[0][j++];
    floatFeatureVectors[k++]   = sub_vector[1][n];
    floatFeatureVectors[k++]   = sub_vector[2][n];
    floatFeatureVectors[k++]   = sub_vector[3][n];
    floatFeatureVectors[k++]   = speech_or_noise[n];
  }
   
  return YES;   
}

/*##########################################*
 * The following routine converts the raw   *
 * sound data to floating point, and pre-   *
 * emphasizes if requested.                 *
 *##########################################*/
- (float *)rawSoundToFloat
{
   int     i;
   static  float *emphasized_data = NULL;
   double  b[2], a[2];
   char    *cdata;
   short   *sdata;
   int     *idata;
   
   
   if(numberSoundSamples)
   {
     if(emphasized_data != NULL)
       free(emphasized_data);
     if( (emphasized_data=(float*)malloc(numberSoundSamples*sizeof(float)))
                    != NULL)
     {
       switch([self dataFormat])
       {
         case SND_FORMAT_LINEAR_8:
       cdata = (char *)rawSoundData;
           for(i=0; i<numberSoundSamples; i++)
             emphasized_data[i] = cdata[i];
       break;
     case SND_FORMAT_LINEAR_16:
           sdata = (short *)rawSoundData;
           for(i=0; i<numberSoundSamples; i++)
             emphasized_data[i] = sdata[i];
       break;
     case SND_FORMAT_LINEAR_32:
       idata = (int *)rawSoundData;
           for(i=0; i<numberSoundSamples; i++)
             emphasized_data[i] = idata[i];
       break;
     default:
       printf("Not an implemented sound format in rawSoundToFloat\n");
       return emphasized_data;
       break;
       }
       if(preemphasisCoefficient < EPS)
         return emphasized_data;
       [digitalFilter setFilterOrder:1];
       b[1] = 1.;
       b[0] = -preemphasisCoefficient;
       a[1] = 1.;
       a[0] = 0.;
       [digitalFilter setFilterACoeff:a bCoeff:b];
       [digitalFilter filterFloatArray:emphasized_data
                        numberPts:numberSoundSamples];
       
       return emphasized_data;
     }
   }
   
   return NULL;
}

/*######################################*
 * The following method converts    *
 * floatFeatureVectors -> intFeature    *
 * Vectors.                             *
 *######################################*/
- floatVectorsToInt
{
  int    i, j, k, number_frames;
   
// Find the strobes and distance functions:
  if(featureVectorLength == 0)
     return self;
    
  if( intFeatureVectors != NULL)
     free(intFeatureVectors);
  number_frames = [speechProcessor numberFrames];
  intFeatureVectors = get_int_spc(number_frames*featureVectorLength);
  k = 0;
  for(i=0; i<number_frames; i++)
  {
     for(j=0; j<featureVectorLength; j++)
     {
       intFeatureVectors[k] = intScaleFactor*floatFeatureVectors[k];
       k++;
     }
  }
  return self;
}

/*######################################*
 * This method converts the feature     *
 * vector into float format for         *
 * saving into a .prc file.             *
 *######################################*/
- floatVectorsToChar
{
  int   i, j, k, n, number_frames, spectral_length;
  char  *char_data;
  short *short_data;

  number_frames   = [speechProcessor numberFrames];
  spectral_length = [speechProcessor spectralLength];

  if(charFeatureVectors != NULL)
    free(charFeatureVectors);
  charFeatureVectors =
               get_char_spc(numberBytes*featureVectorLength*number_frames);
  k = n = 0;
  switch(numberBytes)
  {
    case 1:
      char_data = charFeatureVectors;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<featureVectorLength; j++)
          char_data[k++] = (char) intScaleFactor*floatFeatureVectors[n++];
      }

// Set min and max values:
      subvectorCount[0] = spectral_length;
      for(i=0; i<SUB_VECTORS; i++)
      {
        subvectorMin[i] = ONE_BYTE_SCALE;
        subvectorMax[i] = -(ONE_BYTE_SCALE);
      }
      k = 0;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<spectral_length; j++)
        {
          subvectorMin[0] = MIN(subvectorMin[0], char_data[k]);
          subvectorMax[0] = MAX(subvectorMax[0], char_data[k]);
          k++;
        }
        for(n=1; n<SUB_VECTORS; n++)
        {
          subvectorMin[n] = MIN(subvectorMin[n], char_data[k]);
          subvectorMax[n] = MAX(subvectorMax[n], char_data[k]);
          k++;
        }
      }
      break;
    case 2:
      short_data = (short *)charFeatureVectors;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<featureVectorLength; j++)
          short_data[k++] = (short) intScaleFactor*floatFeatureVectors[n++];
      }

// Set min and max values:
      subvectorCount[0] = spectral_length;
      for(i=0; i<SUB_VECTORS; i++)
      {
        subvectorMin[i] = TWO_BYTE_SCALE;
        subvectorMax[i] = -(TWO_BYTE_SCALE);
      }
      k = 0;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<spectral_length; j++)
        {
          subvectorMin[0] = MIN(subvectorMin[0], short_data[k]);
          subvectorMax[0] = MAX(subvectorMax[0], short_data[k]);
          k++;
        }
        for(n=1; n<SUB_VECTORS; n++)
        {
          subvectorMin[n] = MIN(subvectorMin[n], short_data[k]);
          subvectorMax[n] = MAX(subvectorMax[n], short_data[k]);
          k++;
        }
      }
      break;
    default:
      printf("Illegal # of bytes in intVectorsToChar\n");
      break;
  }
  
  return self;
}

/*######################################*
 * This method converts the feature     *
 * vector into integer format for       *
 * compatibility with MANNsim input.    *
 * Also set the min and max array vals. *
 *######################################*/
- intVectorsToChar
{
  int   i, j, k, n, number_frames, spectral_length;
  char  *char_data;
  short *short_data;

  number_frames   = [speechProcessor numberFrames];
  spectral_length = [speechProcessor spectralLength];

  if(charFeatureVectors != NULL)
    free(charFeatureVectors);
  charFeatureVectors =
               get_char_spc(numberBytes*featureVectorLength*number_frames);
  k = n = 0;
  switch(numberBytes)
  {
    case 1:
      char_data = charFeatureVectors;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<featureVectorLength; j++)
          char_data[k++] = (char) intFeatureVectors[n++];
      }

// Set min and max values:
      subvectorCount[0] = spectral_length;
      for(i=0; i<SUB_VECTORS; i++)
      {
        subvectorMin[i] = ONE_BYTE_SCALE;
        subvectorMax[i] = -(ONE_BYTE_SCALE);
      }
      k = 0;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<spectral_length; j++)
        {
          subvectorMin[0] = MIN(subvectorMin[0], char_data[k]);
          subvectorMax[0] = MAX(subvectorMax[0], char_data[k]);
          k++;
        }
        for(n=1; n<SUB_VECTORS; n++)
        {
          subvectorMin[n] = MIN(subvectorMin[n], char_data[k]);
          subvectorMax[n] = MAX(subvectorMax[n], char_data[k]);
          k++;
        }
      }
      break;
    case 2:
      short_data = (short *)charFeatureVectors;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<featureVectorLength; j++)
          short_data[k++] = (short) intFeatureVectors[n++];
      }

// Set min and max values:
      subvectorCount[0] = spectral_length;
      for(i=0; i<SUB_VECTORS; i++)
      {
        subvectorMin[i] = TWO_BYTE_SCALE;
        subvectorMax[i] = -(TWO_BYTE_SCALE);
      }
      k = 0;
      for(i=0; i<number_frames; i++)
      {
        for(j=0; j<spectral_length; j++)
        {
          subvectorMin[0] = MIN(subvectorMin[0], short_data[k]);
          subvectorMax[0] = MAX(subvectorMax[0], short_data[k]);
          k++;
        }
        for(n=1; n<SUB_VECTORS; n++)
        {
          subvectorMin[n] = MIN(subvectorMin[n], short_data[k]);
          subvectorMax[n] = MAX(subvectorMax[n], short_data[k]);
          k++;
        }
      }
      break;
    default:
      printf("Illegal # of bytes in intVectorsToChar\n");
      break;
  }
  
  return self;
}


@end
