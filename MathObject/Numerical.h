/************************************************************************
 * This subclass of object implements certain numerical math methods    *
 *                                                                      *
 * File:Numerical.h                                                     *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/26/92  - Started                                              *
 *  2. 11/23/93  - Add numerical integration                            *
 ************************************************************************/
#import <objc/Object.h>


@interface Numerical:Object

{
  double  machineEps;           /* Machine epsilon */
}

/**********************
 * Initialization:    *
 **********************/
- init;
        
/**********************
 * Running:           *
 **********************/
- (double) zeroIn:(double)ax   to:(double)bx   withTol:(double)tol
             from:(id)sender   error:(BOOL *)error;
- (double) integrateFrom:(double) ax to:(double)bx from:(id)sender;
- (double) twoDimIntegrateFrom:(double) ax to:(double)bx from:(id)sender;

/**********************
 * Set parameters:    *
 **********************/

/**********************
 * Get parameters:    *
 **********************/


@end



/*
 * We use the objective-C respondsTo: mechanism to see if we can send the
 * message to get the function to operate on.  This dummy interface
 * declaration declares those messages (so that even if they don't exists,
 * we can at least use them to check with respondsTo:).
 */

@interface PossibleFunction : Object

- (double)functionToFindRoot:(double)x;
- (double)functionToIntegrate:(double)x;
- (double)outsideFunctionToIntegrate:(double)x;

@end
