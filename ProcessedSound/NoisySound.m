/************************************************************************
 * This object is a subclass of the NISTSound object. It provides for   *
 *  adding noise to a sound object.  The data from noise.snd is in:     *
 *  *floatNoise.  The data from signal.snd is in: *soundData.           *
 *  (Integer data for both is elsewhere).                               *
 *                                                                      *
 * File path: /User/frank/Objc_Classes/ProcessedSound                   *
 * File name: NoisySound.m                                              *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 04/13/93  - Started with skeleton by MF                          *
 *  2. 04/16/93  - First version by SG                                  *
 *  3. 07/12/93  - Added spectral subtraction headers.                  *
 ***********************************************************************/

#import "NoisySound.h"
#import "SpeechProcessor.h"
#import "c_headers.h"
#import "speechParams.h"
#import <math.h>
#import <string.h>

#define DEFAULT_THRESHOLD 12.5  /* = 2000/nfrsz, in BK's code nfrsz=160? */

/****************************************
 * Names for putting in .prc header *
 * NOTE: If these change, make the same *
 * changes in ProcessedSound.m      *
 ****************************************/
char *feature_subvectors[SUB_VECTORS] = 
               {" ", "Zero_Cross", "PTP", "Frame_Energy", "Speech_Or_Noise"};


@implementation NoisySound

/********************************
 * Initialization & free:       *
 * - initialize instance vars   *
 * and call super's init.       *
 ********************************/
- initWithSoundWindow: (id)newSoundWindow
{
  [super initWithSoundWindow: (id)newSoundWindow];
// Initialize instance variables
  noise             = [Sound new];
  playSound         = nil;
  floatNoise        = NULL;
  cleanClassifier   = NULL;
  noiseSamples      = noiseOffset = 0;
  snr               = 100.;
  noiseEnergy       = 0;
  signalEnergy      = -1;
  threshold         = DEFAULT_THRESHOLD;
  noiseFile[0]      = '\0';
  addNoise          = NO;
  hiEnergyFrameCt   = 0;
  speechNoiseClassifierType = NOISELESS_CLASSIFIER;
  
  return self;
}

- free
{
  free(floatNoise);
  return [super free];
}

/*******************************
 * These methods set            *
 * instance variables           *
 *******************************/
- setAddNoise: (BOOL)flag
{
  addNoise = flag;
  return self;
}

- setSpeechNoiseClassifier: (int)type
{
  speechNoiseClassifierType = type;
  return self;
}

- setThreshold: (float)thold
{
  threshold = thold;
  return self;              /* what's standard here? */
}

- setSignalToNoiseRatio: (float)signalToNoise
{
  snr = signalToNoise;
  return self;
}

- setNoiseOffset:(int)newOffset
{
  noiseOffset = newOffset;
  return self;
}

- setSpectralSubtraction:(BOOL)flag
{
  [speechProcessor setSpectralSubtraction:flag];
  return self;
}

- setSpectralFloor:(float)floor
{
  [speechProcessor setSpectralFloor:floor];
  return self;
}

- setNoiseFactor:(float)factor
{
  [speechProcessor setNoiseFactor:factor];
  return self;
}

/*##################################*
 * Find energies in the signal and  *
 * in the noise                     *
 *##################################*/
- (int)calculateEnergies
{
  int   i, j, j_limit, frame_length;
  float frame_energy_sum, hi_energy_sum=0.0; 
  float frame_threshold;

  frame_length = [speechProcessor frameLength];

/* Find noise energy */
  if(floatNoise != NULL)
  {
    i = noiseOffset%noiseSamples;
    noiseEnergy = 0.;
    for(j=0; j<numberSoundSamples; j++)
    {
     noiseEnergy += (floatNoise[i]*floatNoise[i]);
     if(++i == noiseSamples)
       i = 0;
    } 
    noiseEnergy /= ((double) numberSoundSamples);
  }

/* get ave energy of signals (over data from frames with energy > th) */
  hiEnergyFrameCt=0;
  frame_threshold = threshold*frame_length;
  j = 0;
  j_limit = numberSoundSamples - frame_length;
  while( j <= j_limit ) {
    i = j;
    frame_energy_sum = 0.;
    j += frame_length;
    while( i < j ) {
      frame_energy_sum += (soundData[i])*(soundData[i]);
      i++;
    }
    if( frame_energy_sum > frame_threshold ) {
      hi_energy_sum += frame_energy_sum;
      hiEnergyFrameCt++;
    }
  }
  if( hiEnergyFrameCt == 0 ) {
    printf("In NoisySound's addNoiseToSoundWithSNR method, hiEnergyCt == 0\n");
    return 1;
  }
  signalEnergy = hi_energy_sum/(frame_length*hiEnergyFrameCt);
  return 0;
}

/*******************************
 * These methods return     *
 * instance variables       *
 *******************************/
- (int)noiseSamples     {return noiseSamples;}
- (int)noiseOffset      {return noiseOffset;}
- (float)snr            {return snr;}
- (char *)noiseFile     {return noiseFile;}
- (float)threshold      {return threshold;}
- (float)noiseEnergy    {return noiseEnergy;}
- (float)signalEnergy   {return signalEnergy;}
- (int)hiEnergyFrameCt  {return hiEnergyFrameCt;}
- (float)spectralFloor  {return [speechProcessor spectralFloor];}
- (float)noiseFactor    {return [speechProcessor noiseFactor];}

/****************************************************************
 *    This method adds noise to the sound file.  We probably    *
 * want to overwrite the sound file, but should we save a       *
 * copy in case this method gets called more than once?         *
 ****************************************************************/
- (BOOL) addNoiseToSound
{
  int   i, j=0;
  float scale_factor;

  if( [self calculateEnergies] != 0 )
  {
    printf("In NoisySound's addNoiseToSoundWithSNR: method, no signal energy.\n");
    return NO;
  }
/* calculate scale factor */
  if( noiseEnergy <= 0 )
  {
    printf("In NoisySound's addNoiseToSoundWithSNR: method, bad energy.\n");
    return NO;
  }
  scale_factor = sqrt( signalEnergy/(noiseEnergy*pow(10.0,snr*0.1)));
  
/* scale noise and add to each sample */
  j = 0;
  i = noiseOffset%noiseSamples;
  while( j < numberSoundSamples )
  {
    soundData[j] += scale_factor*floatNoise[i];
    j++;
    if( ++i == noiseSamples ) i=0;  /* reuses noise if not enough */
  }
  return YES;
}

/*##############################*
 * Play the noisy sound:    *
 *##############################*/
- playNoisySound
{
  unsigned char *sound_data;
  int format, channel_count, size;
  double sampling_rate;

  if(playSound != nil)
    [playSound free];
  playSound = [Sound new];
  format        = [self dataFormat];
  sampling_rate = [self samplingRate];
  channel_count = [self channelCount];
  size          = format*numberSoundSamples;
  [playSound setDataSize:size dataFormat:format
             samplingRate:sampling_rate channelCount:channel_count
             infoSize:4];
  sound_data = [playSound data];
  [self convert:sound_data length:numberSoundSamples format:[self dataFormat]
                fromFloat:soundData];
  [playSound play];
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
 * NOTE: This overrides super's method. *
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
  if( (speechNoiseClassifierType==NOISELESS_CLASSIFIER) && (cleanClassifier!=NULL) )
    speech_or_noise = cleanClassifier;
  else
    speech_or_noise = [speechProcessor speechOrNoise];
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
    if(!forcedRange[i])
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
    if(linearMapping[i])
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

/*##############################*
 * Override super's method in   *
 * order to add in noise.       *
 *##############################*/
- (BOOL) readSndFile:(const char *)fileName
{
  float *energy;
  int   number_frames;

  if([super readSndFile:fileName])
  {
    if(speechNoiseClassifierType == NOISELESS_CLASSIFIER)
    {
      energy          = [speechProcessor calculateEnergyData:soundData
                                       start:selectionStart length:selectionLength];
      number_frames   = [speechProcessor numberFrames];
      cleanClassifier = [speechProcessor calculateSpeechOrNoiseFromEnergy:energy
                                                        numberFrames:number_frames];
    }
    if(addNoise)
      return([self addNoiseToSound]);
    return YES;
  }
  return NO;
}

/*##############################*
 * Override super's method in   *
 * order to add in noise.       *
 *##############################*/
- (BOOL) readMikeFile:(const char *)fileName
{
  float *energy;
  int   number_frames;

  if([super readMikeFile:fileName])
  {
    if(speechNoiseClassifierType == NOISELESS_CLASSIFIER)
    {
      energy          = [speechProcessor calculateEnergyData:soundData
                                       start:selectionStart length:selectionLength];
      number_frames   = [speechProcessor numberFrames];
      cleanClassifier = [speechProcessor calculateSpeechOrNoiseFromEnergy:energy
                                                        numberFrames:number_frames];
    }
    if(addNoise)
      return([self addNoiseToSound]);
    return YES;
  }
  return NO;
}

/*##############################*
 * Override super's method in   *
 * order to add in noise.       *
 *##############################*/
- (BOOL) readWavFile:(const char *)fileName
{
  float *energy;
  int   number_frames;

  if([super readWavFile:fileName])
  {
    if(speechNoiseClassifierType == NOISELESS_CLASSIFIER)
    {
      energy          = [speechProcessor calculateEnergyData:soundData
                                       start:selectionStart length:selectionLength];
      number_frames   = [speechProcessor numberFrames];
      cleanClassifier = [speechProcessor calculateSpeechOrNoiseFromEnergy:energy
                                                        numberFrames:number_frames];
    }
    if(addNoise)
      return([self addNoiseToSound]);
    return YES;
  }
  return NO;
}

/*********************************
 * This method reads in         *
 * the noise file:              *
 *********************************/
- (int) readNoiseFile:(const char *)nFile
{           /* changed parameter name to avoid msg */
  float *pFloat;
  short *pShort;

  if ( [noise readSoundfile: nFile] != 0 ) {
    printf("In NoisySound's readNoiseFile method, cannot read.\n");
    return 1;
  }
  strncpy(noiseFile, nFile, 255);
  if ( (noiseSamples = [noise sampleCount]) == 0 ) {
    printf("In NoisySound's readNoiseFile method, no data.\n");
    return 2;
  }
  if( [noise dataFormat] != SND_FORMAT_LINEAR_16 ) {
    printf("In NoisySound's readNoiseFile method, bad format\n");
    return 3;
  }
  if( [noise channelCount] != 1 ) {
    printf("In NoisySound's readNoiseFile method, bad channel\n");
    return 4;
  }
  if( fabs([noise samplingRate] - SND_RATE_CODEC)> 100. ) {
    printf("In NoisySound's readNoiseFile method, bad rate\n");
    return 5;
  }
  if( floatNoise != NULL ) free(floatNoise);
  if ( (floatNoise = malloc(sizeof(float)*noiseSamples)) == NULL ) {
    printf("In NoisySound's readNoiseFile method, cannot alloc\n");
    return 6;
  }
  if ( [noise compactSamples] != 0 ) {
    printf("In NoisySound's readNoiseFile method, cannot compact\n");
    return 7;
  }
// Swap bytes if necessary:
  SNDSwapSoundToHost([noise data], [noise data], noiseSamples, [noise channelCount],
                      [noise dataFormat]);
  pShort= (short*) [noise data];
  for(  pFloat=floatNoise; pFloat<floatNoise+noiseSamples; pFloat++)
  {
    *pFloat = (float) *pShort;
    pShort++;           /* assumes 16-bit data, shorts! */
  }
  return 0;
}

/*###############################*
 * Save the processed sound file *
 *###############################*/
- (BOOL) savePrcFile: (char *)processedFile
{
  int   vector_size, output_count, number_frames;
  FILE  *file_write;

  vector_size   = featureVectorLength;
  number_frames = [speechProcessor numberFrames];
  output_count  = numberBytes*vector_size*number_frames;
  if( (file_write = fopen(processedFile,"w")) == NULL)
    return NO;
  if( featureVectorLength )
  {
    if(charFeatureVectors != NULL)
      if([self writeProcessedHeader: file_write: number_frames: vector_size])
      {
       [self writeProcessedData: file_write: output_count: charFeatureVectors];
       fclose(file_write);
       return YES;
      }
  }
  fclose(file_write);
  
  return NO;
}


/*###############################*
 * Write the processed sound     *
 * header:                       *
 *###############################*/
- writeProcessedData: (FILE *)filePtr: (int)outputCount:
                              (char *)outputData
{
  int i;
  for(i=0; i<outputCount; i++)
    putc(outputData[i], filePtr);
  return self;
}
    
/*###############################*
 * Write the processed sound     *
 * header:                       *
 *###############################*/
- (BOOL) writeProcessedHeader: (FILE *)filePtr: (int)speechFrames:
                              (int) vectorSize
{
  char   database_id[STRING_LEN], database_version[STRING_LEN];
  char   utterance_id[STRING_LEN];
  char   name_subvector[STRING_LEN], subvector_count[STRING_LEN];
  char   subvector_min[STRING_LEN], subvector_max[STRING_LEN];
  char   subvector_scale[STRING_LEN], subvector_offset[STRING_LEN];
  char   low_pass_filter[STRING_LEN];
  char   *info_ptr, *window_name, *frequency_type, *lifter_name;
  char   *time_and_date;
  int    low_pass_order, channel_count, nist_sample_count;
  int    n, itemp;
  int    output_sample_count, sub_count;
  int    byte_format, element_bits;
  int   frame_length, overlap_length;
  long   h_bytes, data_bytes;
  float  temp1, temp2, temp3;
  double nist_sample_rate, output_sample_rate, low_pass_offset, pc, dtemp;
  struct header_t *output_header;
  time_t clock;

  frame_length   = [speechProcessor frameLength];
  overlap_length = [speechProcessor overlapLength];
  
// First, get header data from input sound info
// Note: this must be compatible with NISTSound.m

  info_ptr = [self info];
  if(inputSoundType!=SPEECH_FROM_MIKE)
    sscanf(info_ptr, "%s %s %s %f %f %s %f %d", database_id, database_version,
                   utterance_id, &temp1, &temp2, low_pass_filter,
           &temp3, &low_pass_order);
  else
  {
    strcpy(database_id, "input_speech");
    strcpy(database_version, "v0");
    strcpy(utterance_id, inputSpeechFile);
    strcpy(low_pass_filter, "none");
    low_pass_order = 0;
    temp1 = temp3 = 0.;
    temp2 = [self samplingRate];
  }
  nist_sample_rate = temp1;
  output_sample_rate = temp2;
  low_pass_offset  = temp3;
           
// Get other info from sound object:
  channel_count = [self channelCount];
  output_sample_count = [self sampleCount];
  if(output_sample_rate > 0.)
    nist_sample_count = output_sample_count *
                        (nist_sample_rate/output_sample_rate);
  else
    nist_sample_count = 0;
  speechFrames  = (selectionLength-frame_length)/(frame_length-overlap_length);
  window_name    = [self windowName];
  frequency_type = [self frequencyAnalysisName];
  lifter_name    = [self lifterName];
  pc             = (double) preemphasisCoefficient;
  sub_count  = SUB_VECTORS;
  byte_format    = ELEMENT_BYTE_FORMAT;
  element_bits   = ELEMENT_BITS;
  time(&clock);
  time_and_date  = ctime(&clock);
  info_ptr       = strchr(time_and_date,'\n');
  *info_ptr      = '\0';
  
// Initialize the header and add the fields to it:
  output_header = sp_create_header();
  sp_add_field(output_header, "database_id", T_STRING, database_id);
  sp_add_field(output_header, "database_version", T_STRING, database_version);
  sp_add_field(output_header, "utterance_id", T_STRING, utterance_id);
  sp_add_field(output_header, "channel_count", T_INTEGER, &channel_count);
  sp_add_field(output_header, "sample_count", T_INTEGER, &nist_sample_count);
  sp_add_field(output_header, "sample_rate", T_REAL, &nist_sample_rate);
  sp_add_field(output_header, "output_sample_rate", T_REAL,
                                                &output_sample_rate);
  sp_add_field(output_header, "low_pass_filter", T_STRING, low_pass_filter);
  sp_add_field(output_header, "low_pass_offset", T_REAL, &low_pass_offset);
  sp_add_field(output_header, "low_pass_order", T_INTEGER, &low_pass_order);
  sp_add_field(output_header, "output_sample_count", T_INTEGER,
                                                &output_sample_count);
  sp_add_field(output_header, "frame_sample_count", T_INTEGER, &frame_length);
  sp_add_field(output_header, "overlap_sample_count", T_INTEGER,
                                                &overlap_length);
  sp_add_field(output_header, "number_frames", T_INTEGER, &speechFrames);
  sp_add_field(output_header, "speech_start_sample", T_INTEGER,
                                                &selectionStart);
  sp_add_field(output_header, "speech_length", T_INTEGER, &selectionLength);
  sp_add_field(output_header, "preemphasis_coefficient", T_REAL, &pc);
  sp_add_field(output_header, "window_type", T_STRING, 
                                        [speechProcessor windowName]);
    itemp = [speechProcessor numberFilterBanks];
  sp_add_field(output_header, "number_filter_banks", T_INTEGER,
                                                &itemp);
    itemp = [speechProcessor lpcOrder];
  sp_add_field(output_header, "lpc_order", T_INTEGER, &itemp);
  sp_add_field(output_header, "frequency_analysis", T_STRING, frequency_type);
  sp_add_field(output_header, "lifter_type", T_STRING, 
                                     [speechProcessor lifterName]);
  sp_add_field(output_header, "noise_file", T_STRING, noiseFile);
  sp_add_field(output_header, "noise_offset", T_INTEGER, &noiseOffset);
  sp_add_field(output_header, "signal_threshold", T_REAL, &threshold);
  sp_add_field(output_header, "snr", T_REAL, &snr);
    itemp = [speechProcessor spectralSubtraction];
  sp_add_field(output_header, "spectral_subtraction", T_INTEGER, &itemp);
    dtemp = [speechProcessor spectralFloor];
  sp_add_field(output_header, "spectral_floor", T_REAL, &dtemp);
    dtemp = [speechProcessor noiseFactor];
  sp_add_field(output_header, "noise_factor", T_REAL, &dtemp);
  sp_add_field(output_header, "analysis_date", T_STRING, time_and_date);
  sp_add_field(output_header, "vector_element_count", T_INTEGER,
                                                &vectorSize);
  itemp = SUB_VECTORS;
  sp_add_field(output_header, "subvector_count", T_INTEGER, &itemp);
  sp_add_field(output_header, "element_n_bytes", T_INTEGER,
                                                &numberBytes);
  sp_add_field(output_header, "element_byte_format", T_INTEGER,
                                                &byte_format);
  sp_add_field(output_header, "element_sig_bits", T_INTEGER,
                                                &element_bits);
  for(n=0; n<SUB_VECTORS; n++)
  {
    sprintf(name_subvector, "name_subvector_%d", n);
    sprintf(subvector_count, "subvector_%d_count",n);
    sprintf(subvector_min, "subvector_%d_min",n);
    sprintf(subvector_max, "subvector_%d_max",n);
    sprintf(subvector_scale, "subvector_%d_scale",n);
    sprintf(subvector_offset, "subvector_%d_offset",n);
    if(n == 0)
      sp_add_field(output_header, name_subvector, T_STRING, 
                                  [speechProcessor spectralName]);
    else
      sp_add_field(output_header, name_subvector, T_STRING, 
                                  feature_subvectors[n]);
    sp_add_field(output_header, subvector_count, T_INTEGER, 
                                  &subvectorCount[n]);
    sp_add_field(output_header, subvector_min, T_INTEGER, 
                                  &subvectorMin[n]);
    sp_add_field(output_header, subvector_max, T_INTEGER, 
                                  &subvectorMax[n]);
    sp_add_field(output_header, subvector_scale, T_REAL, 
                                  &subvectorScale[n]);
    sp_add_field(output_header, subvector_offset, T_REAL, 
                                  &subvectorOffset[n]);
  }
/*****************************************
 * To write header to stdout:            *
  sp_print_lines(output_header, stdout); *
 *****************************************/
  n = sp_write_header(filePtr, output_header, &h_bytes, &data_bytes);
  sp_close_header(output_header);
  if(n == 0 )
    return YES;
  else
    return NO;
}




@end
