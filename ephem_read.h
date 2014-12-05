/**==========================================================================**/
/**                                                                          **/
/**  SOURCE FILE: ephem_read.h                                               **/
/**                                                                          **/
/**     This file contains C prototype declarations for the functions that   **/
/**     are defined in the source file ephem_read.c; the purpose of each of  **/
/**     them is explained there.                                             **/
/**                                                                          **/
/**  Programmer:  David Hoffman/EG5                                          **/
/**               NASA, Johnson Space Center                                 **/
/**               Houston, TX 77058                                          **/
/**               e-mail: david.a.hoffman1@jsc.nasa.gov                      **/
/**                                                                          **/
/**==========================================================================**/

#ifndef TYPES_DEFINED
#include "ephem_types.h"
#endif

/*----------------------------------------------------------------------------*/
/*  Initialize_Ephemeris                                                      */
/*----------------------------------------------------------------------------*/

    int Initialize_Ephemeris( char *fileName );

/*----------------------------------------------------------------------------*/
/*  Interpolate_Libration                                                     */
/*----------------------------------------------------------------------------*/

    void Interpolate_Libration( long double Time, int Target, long double Libration[3] );

/*----------------------------------------------------------------------------*/
/*  Interpolate_Nutation                                                      */
/*----------------------------------------------------------------------------*/

    void Interpolate_Nutation( long double Time , int Target , long double Nutation[2] );
   
/*----------------------------------------------------------------------------*/
/*  Interpolate_Position                                                      */
/*----------------------------------------------------------------------------*/

    void Interpolate_Position( long double Time , int Target , long double Position[3] );

/*----------------------------------------------------------------------------*/
/*  Interpolate_State                                                         */
/*----------------------------------------------------------------------------*/

    void Interpolate_State( long double Time , int Target , stateType *Planet );

long double Find_DataInHeader( char    *szTarget);