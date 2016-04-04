/************************************************************************
 * This object is a subclass of the processedSound object. It provides  *
 * additional processing of NIST .wav files.                            *
 *                                                                      *
 * File: /User/frank/Objc_Classes/ProcessedSound/NISTSound.m            *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 07/02/92  - Started                                              *
 *  2. 12/14/92  - Corrected error in convertSound to indicate that     *
 *                 output sample rate = NIST sample rate.               *
 *  3. 01/27/93  - Add playSound and saveSndFile: methods.              *
 *               - Fix -convertSound to store info in convertedSound    *
 *  4. 04/26/93  - Delete -findStrobes with addition of speechProcessor *
 *                                                                      *
 ************************************************************************/

#define INFO_PAD_SIZE   103     /* Extro info size */

#import "NISTSound.h"
#import "NISTPhonemeMap.h"
#import "DigitalFilter.h"
#import "c_headers.h"
#import <string.h>
#import <math.h>


@implementation NISTSound

/*****************************************
 * Create the NISTSound and initialize   *
 * instance variables:                   *
 *****************************************/
- initWithSoundWindow: (id)newSoundWindow
{
  [super initWithSoundWindow: (id)newSoundWindow];
// Initialize instance variables
  phonemeMap            = [[NISTPhonemeMap alloc] init];
  convertedSound        = nil;
  filterType            = HAMMING_SINC;
  numberPoles           = DEFAULT_POLES;
  outputSamplingType    = SAMPLE_CODEC;
  numberPhonemes        = 0;
  offset                = 0.01;
  phonemeLabelsAvailable = NO;
  
  return self;
}

/***********************************************
 * The following methods get information from  *
 * the NIST Header:                            *
 ***********************************************/
#define MAX_FIELD 80

/*****************************************
 * Get the format of the NIST data       *
 *****************************************/
- (int) NISTdataFormat
{
  int   field_size;
  int   output;
  
  field_size = sp_get_size(nistHeader, "sample_n_bytes");
  if(field_size == sizeof(int))
  {
    sp_get_data(nistHeader, "sample_n_bytes", &output, &field_size);
    switch(output)
    {
      case 1:
        return SND_FORMAT_LINEAR_8;
    break;
      case 2:
        return SND_FORMAT_LINEAR_16;
    break;
      case 3:
        return SND_FORMAT_LINEAR_24;
    break;
      case 4:
        return SND_FORMAT_LINEAR_32;
    break;
      default:
        break;
    }
  }
  return SND_FORMAT_UNSPECIFIED;
}

/*****************************************
 * Get the number of bytes in each       *
 * sample                                *
 *****************************************/
- (int) NISTdataSize
{
  int   data_size;
  
  [super getNISTInteger:"sample_n_bytes" :&data_size];
  return data_size;
}

/*****************************************
 * Get the sampling rate:                *
 *****************************************/
- (int) NISTsamplingRate
{
  int   sample_rate;
  
  [super getNISTInteger:"sample_rate" :&sample_rate];
  return sample_rate;
}

/*****************************************
 * Get the channel count:                *
 *****************************************/
- (int) NISTchannelCount
{
  int   channel_count;
  
  [super getNISTInteger:"channel_count" :&channel_count];
  return channel_count;
}

/*****************************************
 * Get the info size:                    *
 * Info consists of database_id, database*
 * version and utterance size.           *
 *****************************************/
- (int) NISTinfoSize
{
  int   field_size;
  
  field_size  = sp_get_size(nistHeader, "database_id");
  field_size += sp_get_size(nistHeader, "database_version");
  field_size += sp_get_size(nistHeader, "utterance_id");
  field_size += 4;      /* For spaces */
  return field_size;
}


/*****************************************
 * Get the number of samples:            *
 *****************************************/
- (int) NISTnumberSamples
{
  int   sample_count;
  
  [super getNISTInteger:"sample_count" :&sample_count];
  return sample_count;
}

/*****************************************
 * Get the sound info:                   *
 * Get info from three fields, pack into *
 * output:                               *
 *****************************************/
- (char *) NISTsoundInfo
{
  int   i, j, field_size, output_size;
  char  *field, *output;
  
  output_size = [self NISTinfoSize];
  output = (char  *)malloc(output_size*sizeof(char));
  if(output == NULL)
    return output;
  for(j=0; j<output_size; j++)
    output[j]= '\0';
  
  i = 0;
// field 1:
  field_size  = sp_get_size(nistHeader, "database_id");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "database_id", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
// field 2:
  field_size  = sp_get_size(nistHeader, "database_version");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "database_version", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
// field 3:
  field_size  = sp_get_size(nistHeader, "utterance_id");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "utterance_id", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
  return output;
}

/*****************************************
 * Get the minimum sound value:          *
 *****************************************/
- (int) NISTsampleMin
{
  int   sample_min;
  
  [super getNISTInteger:"sample_min" :&sample_min];
  return sample_min;
}

/*****************************************
 * Get the maximum sound value:          *
 *****************************************/
- (int) NISTsampleMax
{
  int   sample_max;
  
  [super getNISTInteger:"sample_max" :&sample_max];
  return sample_max;
}

/*****************************************
 * Get the info size:                    *
 * Info consists of database_id, database*
 * version and utterance size.           *
 *****************************************/
- (int) convertInfoSize
{
  int   field_size;
  
  field_size  = sp_get_size(nistHeader, "database_id");
  field_size += sp_get_size(nistHeader, "database_version");
  field_size += sp_get_size(nistHeader, "utterance_id");
  field_size += INFO_PAD_SIZE;      /* For extra info */
  return field_size;
}


/*****************************************
 * Get info from three fields, pack into *
 * output:                               *
 * Also input_sample_rate, output_sample_*
 * rate, filterType, offset, filter order*
 *****************************************/
- (char *) convertSoundInfo:(float) sampleFreq
{
  int   i, j, field_size, output_size;
  char  *field, *output, *filter;
  
  output_size = [self convertInfoSize];
  output = (char  *)malloc(output_size*sizeof(char));
  if(output == NULL)
    return output;
  
  for(j=0; j<output_size; j++)
    output[j] = '\0';
    
  i = 0;
// field 1:
  field_size  = sp_get_size(nistHeader, "database_id");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "database_id", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
// field 2:
  field_size  = sp_get_size(nistHeader, "database_version");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "database_version", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
// field 3:
  field_size  = sp_get_size(nistHeader, "utterance_id");
  if(field_size)
  {
    field = (char *)malloc(field_size*sizeof(char));
    if(field != NULL)
      sp_get_data(nistHeader, "utterance_id", field, &field_size);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    output[i++] = ' ';
    free(field);
  }
  
// field 4:
  field_size  = (INFO_PAD_SIZE-3)/4;
  field = (char *)malloc(field_size*sizeof(char));
  if(field != NULL)
  {
    sprintf(field,"%d ", [self NISTsamplingRate]);
    field_size = strlen(field);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
  }
  
// field 5:
  if(field != NULL)
  {
    sprintf(field,"%-7.1f ", sampleFreq);
    field_size = strlen(field);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
  }

// field 6:
  if(filterType == BUTTERWORTH)
  {
    field_size = 12;
    filter = "butterworth ";
  }
  else
  {
    field_size = 13;
    filter      = "hamming_sinc ";
  }
  for(j=0; j<field_size; j++)
      output[i++] = filter[j];
  
  
// field 7:
  if(field != NULL)
  {
    sprintf(field,"%-7.3f ", offset);
    field_size = strlen(field);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
  }
  
// field 8:
  if(field != NULL)
  {
    sprintf(field,"%-1d ", numberPoles);
    field_size = strlen(field);
    for(j=0; j<field_size; j++)
      output[i++] = field[j];
    free(field);
  }

  return output;
}

/********************************************
 * Get the min sample value in the sound    *
 * data                                     *
 ********************************************/
- (int) soundSampleMin;
{
  int      min;
  int      i, number_samples;
  unsigned char *sound_data;
  short   *sdata;
  int     *idata;
  
  min        = 10000;
  sound_data = [self data];
  number_samples = [self sampleCount];
  switch([self dataFormat])
  {
    case SND_FORMAT_LINEAR_8:
      for(i=0; i<number_samples; i++)
      {
         if(sound_data[i] < min)
           min = sound_data[i];
      }
      break;
    case SND_FORMAT_LINEAR_16:
      sdata = (short *)sound_data;
      for(i=0; i<number_samples; i++)
      {
         if(sdata[i] < min)
           min = sdata[i];
      }
      break;
    case SND_FORMAT_LINEAR_32:
      idata = (int *)sound_data;
      for(i=0; i<number_samples; i++)
      {
         if(idata[i] < min)
           min = idata[i];
      }
      break;
    default:
      printf("Not an implemented sound format.\n");
      break;
  }
  return min;
}

/********************************************
 * Get the max sample value in the sound    *
 * data                                     *
 ********************************************/
- (int) soundSampleMax;
{
  int      max;
  int      i, number_samples;
  unsigned char *sound_data;
  short   *sdata;
  int     *idata;
  
  max        = -10000;
  sound_data = [self data];
  number_samples = [self sampleCount];
  switch([self dataFormat])
  {
    case SND_FORMAT_LINEAR_8:
      for(i=0; i<number_samples; i++)
      {
         if(sound_data[i] > max)
           max = sound_data[i];
      }
      break;
    case SND_FORMAT_LINEAR_16:
      sdata = (short *)sound_data;
      for(i=0; i<number_samples; i++)
      {
         if(sdata[i] > max)
           max = sdata[i];
      }
      break;
    case SND_FORMAT_LINEAR_32:
      idata = (int *)sound_data;
      for(i=0; i<number_samples; i++)
      {
         if(idata[i] > max)
           max = idata[i];
      }
      break;
    default:
      printf("Not an implemented sound format.\n");
      break;
  }
  return max;
}

/*######################################*
 * The following methods set instance   *
 * variables:               *
 *######################################*/
- setFilterType: (int)type
{
  filterType = type;
  return self;
}

- setNumberPoles: (int)poles
{
  numberPoles = poles;
  return self;
}

- setOutputSamplingType: (int)newType
{
  outputSamplingType = newType;
  return self;
}

- setOffset: (float)newOffset
{
  offset = newOffset;
  return self;
}
 
/*######################################*
 * Convert the sound.  If input sound   *
 * was read from a .wav file, convert   *
 * to NeXT format.          *
 *######################################*/
- convertSound
{
  int    i;
  int    interpolate, decimate;
  int    output_samples, sound_format;
  int    channel_count, info_size, output_size;
  unsigned char    *sound_data;
  char   *info, *info_ptr;
  float  *interp_sound, *dec_sound;
  double input_frequency, output_rate;
  double ratio;
  double fc_norm, fs_norm;
  
  if(convertedSound != nil)
    [convertedSound free];

// Convert the sound to NeXT format:
  convertedSound  = [Sound new];
  input_frequency = (double) [self samplingRate];
  switch(outputSamplingType)
  {
    case SAMPLE_CODEC:
      output_rate = SND_RATE_CODEC;
      break;
    case SAMPLE_LOW:
      output_rate = SND_RATE_LOW;
      break;
    case SAMPLE_HIGH:
      output_rate = SND_RATE_HIGH;
      break;
    default:
      outputSamplingRate = output_rate = input_frequency;
      break;
  }
// Find interpolation/decimation for sampling conversion:
  if(input_frequency>0.)
    ratio = output_rate/input_frequency;
  else
    ratio = 1.;
  if(fabs(ratio-1.)<EPS)        /* no warping necessary */
  {
    dec_sound      = soundData;
    output_samples = numberSoundSamples;
  }
  else
  {
    [digitalFilter findWarpFactorsFor:ratio];
    interpolate = [digitalFilter getInterpolateFactor];
    decimate    = [digitalFilter getDecimateFactor];
    outputSamplingRate = input_frequency*interpolate/decimate;

// Find butterworth filter coefficients
    fs_norm = TWOPI;
    i = interpolate>decimate ? interpolate : decimate;
    fc_norm = PI*(1./(double)i-offset);
    [digitalFilter setFilterFo:0. fc:fc_norm fs:fs_norm];
    [digitalFilter setFilterOrder:numberPoles];
    if(filterType == BUTTERWORTH)
    {
      [digitalFilter findPolesZeros];
      [digitalFilter findTransferFunction];
    }
    else
      [digitalFilter transferForFIRWindow];

// Interpolate sound:
    interp_sound = [digitalFilter interpolateInput:soundData 
                                numberPts:numberSoundSamples];
  
// LPF:
    [digitalFilter filterFloatArray:interp_sound
                        numberPts:numberSoundSamples*interpolate];

// Decimate:
    dec_sound = [digitalFilter decimateInput:interp_sound
                       numberPts:interpolate*numberSoundSamples]; 
    output_samples = numberSoundSamples*interpolate/decimate;
  }

//  Store in converted sound object:
  channel_count = [self channelCount];
  info_size     = [self convertInfoSize];
  sound_format  = [self dataFormat];
  info          = [self convertSoundInfo:outputSamplingRate];

  output_size = [self dataWidth: sound_format]*output_samples;
  [convertedSound setDataSize:output_size dataFormat: sound_format 
                     samplingRate: output_rate
                     channelCount: channel_count infoSize: info_size];
  sound_data = [convertedSound data];
  [self convert: sound_data length:output_samples format:sound_format
                         fromFloat:dec_sound];
  info_ptr = [convertedSound info];
  for(i=0; i<info_size; i++)
    info_ptr[i] = info[i];
  free(info);
  
  return self;
}

/********************************
 * Play the converted sound:    *
 ********************************/
- playConvertedSound
{
  if(convertedSound != nil)
    [convertedSound play];
  return self;
}

/********************************
 * Read the .wav file           *
 ********************************/
- (BOOL) readWavFile: (const char *)fileName
{
  int  i;
  int  sound_format, sound_size, sampling_rate, channel_count;
  int  info_size, sound_error, data_size, data_width;
  int  input_samples, c, extension_location;
  float temp;
  char *errmsg;
  char *sound_data;
  char *sound_info, *info_ptr;
  static char phone_file[STRING_LEN];
  FILE *fp;
  
  fp = fopen(fileName,"r");
  if (fp == FPNULL)
  {
    printf("Can't open .wav file\n");
    return NO;
  }

  if(nistHeader != NULL)
    sp_close_header(nistHeader);
  nistHeader = sp_open_header(fp,TRUE,&errmsg);
  if (nistHeader == HDRNULL) 
    return NO;
  strncpy(inputSpeechFile, fileName, STRING_LEN);
  sound_format  = [self NISTdataFormat];
  data_size     = [self NISTdataSize];
  sampling_rate = [self NISTsamplingRate];
  channel_count = [self NISTchannelCount];
  info_size     = [self NISTinfoSize];
  input_samples = [self NISTnumberSamples];
  sound_info    = [self NISTsoundInfo];
  sound_size    = input_samples*data_size;
  
  if(soundStruct != NULL)
    SNDFree(soundStruct);
  sound_error = SNDAlloc(&soundStruct, sound_size,
                sound_format,sampling_rate, channel_count, info_size);
  if(sound_error != SND_ERR_NONE)
    return NO;
  info_ptr = soundStruct->info;
  for(i=0; i<info_size; i++)
    info_ptr[i] = sound_info[i];
  free(sound_info);
  
// read data into data buffer
  SNDGetDataPointer(soundStruct, &sound_data, &data_size, &data_width);
  for(i=0; i<sound_size; i++)
  {
    if ((c = getc(fp)) != EOF)
      sound_data[i] = c;
    else
      break;
  }
  fclose(fp);
// swap bytes, get instance variables, convert to NeXT format
  [self swapBytes: sound_data: input_samples: sound_format];
  [self setDataFromSoundStructure];
  [self convertSound];
   inputSoundType = SPEECH_FROM_WAV;
// set sound in the plotting window:
   if(soundWindow != nil)
   {
     [soundWindow setSound:convertedSound];
     temp = (float)sampling_rate;
     [soundWindow setSamplingFrequency:temp];
   }
// Set selection to total sound,
   [self setSelectionRange: 0: numberSoundSamples];

// Now, open phoneme file
  strcpy(phone_file, fileName);
  extension_location = sindex(phone_file,WAVE_FILE);
  strcpy(&phone_file[extension_location], PHONE_FILE);
  fp = fopen(phone_file,"r");
  if (fp == FPNULL)
  {
    numberPhonemes = 0;
    phonemeStartSamples[0] = 100000;
    phonemeLabelsAvailable = NO;
    return YES;
  }
  phonemeLabelsAvailable = YES;
  numberPhonemes = 0;
  while(fscanf(fp, "%d %d %s", &phonemeStartSamples[numberPhonemes], 
                               &phonemeEndSamples[numberPhonemes],
                                phonemeList[numberPhonemes]) != EOF)
  {
     phonemeMiddleSamples[numberPhonemes] = (phonemeEndSamples[numberPhonemes]
                             + phonemeStartSamples[numberPhonemes])/2;
     if(numberPhonemes++ == MAX_PHONEMES)
       break;
  }
  fclose(fp);

  return YES;
}

/********************************
 * Save the .snd file           *
 ********************************/
- (BOOL) saveSndFile: (const char *)fileName
{
  if(convertedSound != nil)
    if( [convertedSound writeSoundfile:fileName]==SND_ERR_NONE )
      return YES;
  return NO;
}

/*###############################################*
 * This method finds the endpoints of the speech *
 * waveform, and sets the selection range ac-    *
 * cordingly.                                    *
 *###############################################*/
- findEndpoints
{
  int    i, start;
  int    speech_start, speech_end;
    
  if( (inputSoundType!=SPEECH_FROM_WAV) || (!phonemeLabelsAvailable) )
    return([super findEndpoints]);

  speech_start = speech_end = start = 0;
  for(i=0; i<numberPhonemes-1; i++)
  {
    if([phonemeMap isNoisePhoneme:phonemeList[i]])
    {
      start        = i;
      speech_start = phonemeEndSamples[i];
      break;
    }
  }
  for(i=numberPhonemes-1; i>start; i--)
  {
    if([phonemeMap isNoisePhoneme:phonemeList[i]])
    {
      speech_end = phonemeStartSamples[i];
      break;
    }
  }
  if(speech_end > speech_start)
     [self setSelectionRange: speech_start: (speech_end-speech_start)];
  [self getSelectionRange];
       
  return self;
}

/********************************************
 * The following routine swaps the bytes in *
 * the input sound data:                    *
 ********************************************/
- swapBytes: (char *)data: (int)no_points:
                            (int)data_format
{
   int     i;
   char    *cdata;
   
   if(no_points)
   {
     switch(data_format)
     {
       case SND_FORMAT_LINEAR_8:
     break;
       case SND_FORMAT_LINEAR_16:
         cdata = (char *)malloc(no_points*sizeof(short));
     if(cdata == NULL)
       break;
     swapbytes( data, cdata, 2*no_points);
     for(i=0; i<2*no_points; i++)
       data[i] = cdata[i];
     free(cdata);
     break;
       default:
         printf("Not an implemented sound format.\n");
     break;
     }
   }
   
   return self;
}


/********************************************
 * The following method converts the sound  *
 * floating point data to integer:          *
 ********************************************/
- convert: (UCHAR *)data length:(int)no_points format:(int)data_format
                         fromFloat:(float *)fdata
{
   int     i;
   char    *cdata;
   short   *sdata;
   int     *idata;
   
   if(no_points)
   {
     switch(data_format)
     {
         case SND_FORMAT_LINEAR_8:
       cdata = (char *)data;
           for(i=0; i<no_points; i++)
             cdata[i] = fdata[i];
       break;
     case SND_FORMAT_LINEAR_16:
       sdata = (short *)data;
       for(i=0; i<no_points; i++)
         sdata[i] = fdata[i];
       break;
     case SND_FORMAT_LINEAR_32:
       idata = (int *)data;
           for(i=0; i<no_points; i++)
             idata[i] = fdata[i];
       break;
     default:
           printf("Not an implemented sound format.\n");
       break;
     }
       
   }
   
   return self;
}

/*************************************
 * Return the width of the sound data*
 * That is, # of bytes/word          *
 *************************************/
- (int) dataWidth: (int)soundFormat
{
  switch(soundFormat)
  {
    case SND_FORMAT_LINEAR_8:
      return 1;
      break;
    case SND_FORMAT_LINEAR_16:
      return 2;
      break;
    case SND_FORMAT_LINEAR_32:
      return 4;
      break;
    default:
      printf("Not an implemented sound format.\n");
       break;
  }
  return 0;
}

/********************************
 * Dummy prototype              *
 ********************************/
- setSamplingFrequency: (float)frequency { return self; }

@end
