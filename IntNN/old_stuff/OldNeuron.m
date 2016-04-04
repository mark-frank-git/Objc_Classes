/************************************************************************
 *  Neural Network Classes for the NeXT Computer            *
 *  Written by: Ralph Zazula                    *
 *  University of Arizona                       *
 *  zazula@bonehead.tucson.az.us (NeXT Mail)            *
 *                                  *
 * Revisions:                               *
 *   1.3  92/01/02  14:04:31  zazula                    *
 *       Faster linked-list for connections             *
 *       No more Storage object                     *
 *   1.31 04/01/92  added Sigmoid type = None   (frank)         *
 *                  add neuronType                  *
 *        04/06/92  added numberInputs to speed up weight calculations. *
 *        04/23/92  added stepStart:length:             *
 *   1.4  04/24/92  Changed to new class: IntNeuron (Integer calcs.)    *
 *                  deleted things I wasn't using (temperature, etc.)   *
 *    06/01/92  - Free output arrays on next call           *
 *                                  *
 ************************************************************************/
#import "IntNeuron.h"
#import <appkit/nextstd.h>
#import "math.h"

#define INT_SCALE   100 // Divide by this prior to activation
#define RANDOM_WEIGHT  5    // Random weights [-5,5]


//----------------------------------------------------------

@implementation IntNeuron

- inputs { return inputs; }
- setType:(int)type { nodeType = type; return self; }
- (int)getType { return nodeType; }
- setNeuronType:(int)type {neuronType = type; return self; }
- (int)getNeuronType { return neuronType; }
- setRandom:theRandom { random = theRandom; return self; }
- setRandomMin:(int)min max:(int)max 
         {randomMin=min; randomMax=max; return self;}

//----------------------------------------------------------

- (int)activation:(int)net
{
   int    output;
   double temp;

   switch (nodeType) {
   case Binary :
      temp = (double)net/INT_SCALE; 
      output = (temp > 0.5) ? INT_SCALE : 0;
      break;
   case Sigmoid :
      temp = (double)net/INT_SCALE; 
      output = INT_SCALE/(1.0+exp(-temp));
      break;
   case Sign : 
      temp = (double)net/INT_SCALE; 
      output = (temp > 0.0) ? INT_SCALE : -INT_SCALE;
      break;
   case Tanh :
      temp = (double)net/INT_SCALE; 
      output = INT_SCALE*tanh(temp);
      break;
   case None:
   default:
      output = net;
      break;
   }
   
   return output;
}

//----------------------------------------------------------

- init
{
   [super init];
   lastOutput = 0;
   nodeType = None;           // default node type
   numberInputs = 0;
   neuronType= EUCLIDEAN_DISTANCE;  // default step calculation
   head = tail = NULL;        // initialize the linked-list of connections
   randomMin = -RANDOM_WEIGHT;
   randomMax = RANDOM_WEIGHT;
   intOutputs = NULL;
    
   return self;
}

//-----------------------------------------------------------

- initializeWeights
{
   int        i;
   connection *C;
          
   C = head;
   i = 0;
   while((C != NULL) && (i<numberInputs))
   {
     C->weight = [random randMin:randomMin max:randomMax];
     C = (connection *)C->next; 
   }
   return self;
}

- freeStorage
{
  connection *C, *next;
  
  if(intOutputs!=NULL)
    free(intOutputs);
  C = head;
  while(C != NULL)
  {
    next = (connection *)C->next;
    free(C);
    C = next;
  }
  return self;
}

//-----------------------------------------------------------

- step
// update the output value based on our inputs
{
   connection *C;
   int    output;
   int    temp = 0;
   double norm_input, norm_weight;
   
   norm_input = norm_weight = 0.;
   C = head;
   switch(neuronType)
   {
     case DOT_PRODUCT:
       while(C != NULL)
       {
         temp += C->weight*[C->source lastOutput];
         C = (connection *)C->next;
       }
       break;
     case NORMALIZED_DOT_PRODUCT:
       while(C != NULL)
       {
         output       = [C->source lastOutput];
         temp        += C->weight*output;
         norm_input  += output*output;
         norm_weight += (C->weight)*(C->weight);
         C = (connection *)C->next;
       }
       temp /= sqrt(norm_input*norm_weight);
       break;
     case EUCLIDEAN_DISTANCE:
       while(C != NULL)
       {
         output       = [C->source lastOutput];
     output      -= C->weight;
         temp        += output*output;
         C = (connection *)C->next;
       }
       break;
     default:
       printf("Unknown neuron type\n");
       break;
   }
    
   lastOutput = [self activation:temp];
   
   return self;
}

//-----------------------------------------------------------

- stepStart:(int)start length:(int)length
// update the output value based on our inputs
{
   int  i;
   connection *C;
   int    output;
   int    temp=0;
   double norm_input, norm_weight;
   
   norm_input = norm_weight = 0.;
   C = head;
   i = 0;
   while(C != NULL)
   {
     if(++i>start)
       break;
     C = (connection *)C->next;
   }
   switch(neuronType)
   {
     case DOT_PRODUCT:
       for(i=0; i<length; i++)
       {
         temp += C->weight*[C->source lastOutput];
         C = (connection *)C->next;
     if(C == NULL)
       break;
       }
       break;
     case NORMALIZED_DOT_PRODUCT:
       for(i=0; i<length; i++)
       {
         output       = [C->source lastOutput];
         temp        += C->weight*output;
         norm_input  += output*output;
         norm_weight += (C->weight)*(C->weight);
         C = (connection *)C->next;
     if(C == NULL)
       break;
       }
       temp /= sqrt(norm_input*norm_weight);
       break;
     case EUCLIDEAN_DISTANCE:
       for(i=0; i<length; i++)
       {
         output       = [C->source lastOutput];
     output      -= C->weight;
         temp        += output*output;
         C = (connection *)C->next;
     if(C == NULL)
       break;
       }
       break;
     default:
       printf("Unknown neuron type\n");
       break;
   }
    
   lastOutput = [self activation:temp];
   
   return self;
}

//-----------------------------------------------------------

- (int)lastOutput
{
   return lastOutput;
}

//-----------------------------------------------------------

- connect:sender
{
   if(random == nil) random = [[Random alloc] init];
   return [self connect:sender withWeight:
          [random randMin:randomMin max:randomMax]];
}

//-----------------------------------------------------------

- connect:sender withWeight:(int)weight
//
// adds sender to the list of inputs
// we should check to make sure sender is a Neruon
// also need to check if it is already in the list
//
{
   connection *C;
   
   C = (connection *)malloc(sizeof(connection));
   if(head == NULL)
     head = C;
   else 
     tail->next = C;
   tail = C;
   C->source = sender;
   C->weight = weight;
   C->next   = NULL;
   numberInputs++;
   return self;
}

//-----------------------------------------------------------

- (int)getWeightFor:source
{
   connection *C;

   C = head;
   while((C != NULL) && (C->source != source))
     C = (connection *)C->next;
        
   if(C != NULL) {            // if C==NULL, source isn't an input
      return C->weight;
   }
   else {
      fprintf(stderr,"connection not found in getWeightFor:\n");
      return NAN;
   }
    
}

//-----------------------------------------------------------

- (int *)getAllWeights
{
   int        i;
   connection *C;

   if(intOutputs!=NULL)
     free(intOutputs);
   intOutputs = (int *)malloc(numberInputs*sizeof(int));
   C = head;
   i = 0;
   while((C != NULL) && (i<numberInputs))
   {
     intOutputs[i++] = C->weight;
     C = (connection *)C->next;
        
   }
   return intOutputs;
    
}

//-----------------------------------------------------------
- setAllWeightsTo: (int *)newWeights
{
   int        i;
   connection *C;
   
   C = head;
   i = 0;
   while((C != NULL) && (i<numberInputs))
   {
     C->weight = newWeights[i++];
     C = (connection *)C->next;
   }
   return self;
}

//-----------------------------------------------------------
- setWeightFor:source to:(int)weight
{
   connection *C;
   
   C = head;
   while((C != NULL) && (C->source != source))
     C = (connection *)C->next;
        
   if(C != NULL) {            // if C==NULL, source isn't an input
      C->weight = weight;
      return self;
   }
   else {
      fprintf(stderr,"connection not found in setWeightFor:to:\n");
      return nil;
   }
}

//-----------------------------------------------------------

- setOutput:(int)output
{
   lastOutput = output;
   
   return self;
}
//-----------------------------------------------------------

- changeWeightFor:source by:(int)delta
{
   connection *C;
   
   C = head;
   while((C != NULL) && (C->source != source))
     C = (connection *)C->next;
        
   if(C != NULL)    // if C==NULL, source isn't an input
   {            
      C->weight += delta;
      return self;
   }
   else 
   {
      fprintf(stderr,"connection not found in changeWeightfor:by:\n");
      return nil;
   }
}

//-----------------------------------------------------------

- changeAllWeightsBy: (int *)delta
{
   int        i;
   connection *C;
   
   C = head;
   i = 0;
   while((C != NULL) && (i<numberInputs))
   {
     C->weight += delta[i++];
     C = (connection *)C->next;
   }
   return self;
}

@end
