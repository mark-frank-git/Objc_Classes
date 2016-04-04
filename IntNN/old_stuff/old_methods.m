
#define ETA_INTEGER 10
#define ETA_SCALE   100
/*##############################*
 * Update the weights of    *
 * a single neuron      *
 *##############################*/
- intupdateWeightsOfNeuron: (int) neuronNumber
{
  int    i;
  static int old_number = 0;
  int    *old_weights;
  static int *delta_weights = NULL;
  int   delta;
  
  if( old_number<numberWeightsPerNeuron )
  {
    if(delta_weights!=NULL)
      NXZoneFree(myZone,delta_weights);
    delta_weights =(int *)NXZoneMalloc(myZone,
                                 (numberWeightsPerNeuron*sizeof(int)));
    old_number = numberWeightsPerNeuron;
  }
  old_weights = [mapNeurons[neuronNumber] getAllWeights];
  for(i=0; i<numberInputWeights; i++)
  {
    delta        = ETA_INTEGER*(mapInputs[i] - old_weights[i]);
    if( (delta_weights[i] = delta/ETA_SCALE) == 0 )
      delta_weights[i] = ISGN(delta);
  }
  for(; i<numberWeightsPerNeuron; i++)
  {
    delta        = ETA_INTEGER*(mapInputs[i] - old_weights[i]);
    if( (delta_weights[i] = delta/ETA_SCALE) == 0 )
      delta_weights[i] = ISGN(delta);
  }
  [mapNeurons[neuronNumber] changeAllWeightsBy: delta_weights];

  return self;
}

