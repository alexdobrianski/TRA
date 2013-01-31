/**==========================================================================**/
/**                                                                          **/
/**  SOURCE FILE: rdeph.c                                                    **/
/**                                                                          **/
/**      Purpose: This program interpolates a state from data in a JPL       **/
/**               ephemeris file for a given Julian date. It is designed     **/
/**               to be called from a Tcl/Tk script, so it has a very        **/
/**               unfriendly user interface.                                 **/
/**                                                                          **/
/**               For convenience, here is the mapping between the command   **/
/**               line inputs and the global variables found in read.tcl:    **/
/**                                                                          **/
/**                        argv[1]  <-->  EphName                            **/
/**                        argv[2]  <-->  JD                                 **/
/**                        argv[3]  <-->  TN                                 **/
/**                        argv[4]  <-->  FmtOpt                             **/
/**                        argv[5]  <-->  Cnt                                **/
/**                                                                          **/
/**   Programmer: David Hoffman/EG5                                          **/
/**               NASA, Johnson Space Center                                 **/
/**               Houston, TX 77058                                          **/
/**               e-mail: david.a.hoffman1@jsc.nasa.gov                      **/
/**                                                                          **/
/**==========================================================================**/

#include "stdafx.h"
#define _USE_MATH_DEFINES 1
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ephem_read.h"
#ifndef TYPES_DEFINED
#include "ephem_types.h"
#endif
        


int _tmain(int argc, _TCHAR* argv[])

{
  stateType  State;
  stateType  StateMoon;
  stateType  StateEarthMoon;
  stateType  StateSun;
  stateType  StateJupiter;
  stateType  StateSaturn;
  char       ephemFileName[12] , *tgtName;
  double     Position[3] , Time;
  int        i , Target;
  double dGMB, dEMRAT, dAU, dGMS, dGMSaturn, dGMJupiter;

  if (argc != 6)
  {
      printf("\n Usage:");
      printf("\n rdeph.exe <data-file> <Julian Date> <Target-planet> <Format> <print case>"); 
      printf("\n  <data-file> - files from NASA with all coef");
      printf("\n <Julian date> - date and time injulian calendar");
      printf("\n <target planet> - 0:Mercury;  1:Venus;  2:Earth; 3:Mars");
      printf("\n                   4:Jupiter   5:Saturn  6:Uranus 7:Neptune");
      printf("\n                   8:Pluto     9:Moon   10:Sun    default:Mercury");
      printf("\n <Format> = PosVel - position and velocity");
      printf("\n <print case> to case for printing");
      return 0;
  }
  /* Convert command line arguments to numeric values.........................*/
  
  Time   = atof(argv[2]);
  Target = atoi(argv[3]);

  /* Initialize the ephemeris.................................................*/

  Initialize_Ephemeris(argv[1]);

  /* Compute the desired ephemeris data.......................................*/

  if ( !strcmp(argv[4],"PosVel") ) 
  {
      Interpolate_State( Time , Target , &State );
      Interpolate_State( Time , 10 , &StateSun );
      Interpolate_State( Time , 2 , &StateEarthMoon );
      Interpolate_State( Time , 9 , &StateMoon );
      Interpolate_State( Time , 4 , &StateJupiter );
      Interpolate_State( Time , 5 , &StateSaturn );

  }
  else
      Interpolate_Position( Time , Target , Position );

  /* Print the answer.........................................................*/

  if ( !strcmp(argv[4],"PosVel") ) 
     {
       /* Print seperator.....................................................*/
       
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       printf("  Case %2s  ",argv[5]);
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       
       /* Print planet name...................................................*/
       
       switch (Target) {
          case  0: tgtName = "Mercury";
                   break;
          case  1: tgtName = "Venus";
                   break;
          case  2: tgtName = "Earth";
                   break;
          case  3: tgtName = "Mars";
                   break;
          case  4: tgtName = "Jupiter";
                   break;
          case  5: tgtName = "Saturn";
                   break;
          case  6: tgtName = "Uranus";
                   break;
          case  7: tgtName = "Neptune";
                   break;
          case  8: tgtName = "Pluto";
                   break;
          case  9: tgtName = "Moon";
                   break;
          case 10: tgtName = "Sun";
                   break;
          default: tgtName = "Mercury";
                   break;
       }
       
       printf("\n\n  Target:  %s",tgtName);
       printf("\n      JD:  %8.2f",Time);

       /* Print position......................................................*/
       
       printf("\n\n  Position (km):     [1] =  % 22.15e",State.Position[0]);
       printf("\n                     [2] =  % 22.15e",State.Position[1]);
       printf("\n                     [3] =  % 22.15e",State.Position[2]);

       /* Print velocity......................................................*/

       printf("\n\n  Velocity (km/sec): [1] =  % 22.15e",State.Velocity[0]);
       printf("\n                     [2] =  % 22.15e",State.Velocity[1]);
       printf("\n                     [3] =  % 22.15e\n\n\n",State.Velocity[2]);

       dGMB = Find_DataInHeader("GMB   ");
       dGMS = Find_DataInHeader("GMS   ");
       dEMRAT = Find_DataInHeader("EMRAT ");
       dGMJupiter = Find_DataInHeader("GM5   ");
       dGMSaturn = Find_DataInHeader("GM6   ");
       dAU = Find_DataInHeader("AU    ")*1000.0;
       printf("\n  <CT:setting name=\"AU\" value=\"%23.15e\" />",dAU);
       printf("\n  <CT:setting name=\"GMSun\"  value=\"%23.15e\" />",dGMS*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0));
       printf("\n  <CT:setting name=\"GMEarthMoon\"  value=\"%23.15e\" />",dGMB*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0));
       printf("\n  <CT:setting name=\"GMSaturn\"  value=\"%23.15e\" />",dGMSaturn*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0));
       printf("\n  <CT:setting name=\"GMJupiter\"  value=\"%23.15e\" />",dGMJupiter*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0));
       printf("\n  <CT:setting name=\"GMEarth\"  value=\"%23.15e\" />",dGMB*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0)*(dEMRAT)/(dEMRAT+1.0));
       printf("\n  <CT:setting name=\"GMMoon\"  value=\"%23.15e\" />",dGMB*dAU*dAU*dAU/((24.0*60.0*60.0)*24.0*60.0*60.0)/(dEMRAT+1.0));

       printf("\n  <CT:setting name=\"MassRatioEarthToMoon\"  value=\"%23.15e\" />",dEMRAT);

       printf("\n  <CT:setting name=\"SunNASA_X\" value=\"%23.15e\" />",StateSun.Position[0]*1000.0);
       printf("\n  <CT:setting name=\"SunNASA_Y\" value=\"%23.15e\" />",StateSun.Position[1]*1000.0);
       printf("\n  <CT:setting name=\"SunNASA_Z\" value=\"%23.15e\" />",StateSun.Position[2]*1000.0);
       printf("\n  <CT:setting name=\"SunNASA_VX\" value=\"%23.15e\" />",StateSun.Velocity[0]*1000.0);
       printf("\n  <CT:setting name=\"SunNASA_VY\" value=\"%23.15e\" />",StateSun.Velocity[1]*1000.0);
       printf("\n  <CT:setting name=\"SunNASA_VZ\" value=\"%23.15e\" />",StateSun.Velocity[2]*1000.0);

       printf("\n  <CT:setting name=\"EarthNASA_X\" value=\"%23.15e\" />",StateEarthMoon.Position[0]*1000.0);
       printf("\n  <CT:setting name=\"EarthNASA_Y\" value=\"%23.15e\" />",StateEarthMoon.Position[1]*1000.0);
       printf("\n  <CT:setting name=\"EarthNASA_Z\" value=\"%23.15e\" />",StateEarthMoon.Position[2]*1000.0);
       printf("\n  <CT:setting name=\"EarthNASA_VX\" value=\"%23.15e\" />",StateEarthMoon.Velocity[0]*1000.0);
       printf("\n  <CT:setting name=\"EarthNASA_VY\" value=\"%23.15e\" />",StateEarthMoon.Velocity[1]*1000.0);
       printf("\n  <CT:setting name=\"EarthNASA_VZ\" value=\"%23.15e\" />",StateEarthMoon.Velocity[2]*1000.0);

       printf("\n  <CT:setting name=\"MoonNASA_X\" value=\"%23.15e\" />",StateMoon.Position[0]*1000.0);
       printf("\n  <CT:setting name=\"MoonNASA_Y\" value=\"%23.15e\" />",StateMoon.Position[1]*1000.0);
       printf("\n  <CT:setting name=\"MoonNASA_Z\" value=\"%23.15e\" />",StateMoon.Position[2]*1000.0);
       printf("\n  <CT:setting name=\"MoonNASA_VX\" value=\"%23.15e\" />",StateMoon.Velocity[0]*1000.0);
       printf("\n  <CT:setting name=\"MoonNASA_VY\" value=\"%23.15e\" />",StateMoon.Velocity[1]*1000.0);
       printf("\n  <CT:setting name=\"MoonNASA_VZ\" value=\"%23.15e\" />",StateMoon.Velocity[2]*1000.0);

       printf("\n  <CT:setting name=\"JupiterNASA_X\" value=\"%23.15e\" />",StateJupiter.Position[0]*1000.0);
       printf("\n  <CT:setting name=\"JupiterNASA_Y\" value=\"%23.15e\" />",StateJupiter.Position[1]*1000.0);
       printf("\n  <CT:setting name=\"JupiterNASA_Z\" value=\"%23.15e\" />",StateJupiter.Position[2]*1000.0);
       printf("\n  <CT:setting name=\"JupiterNASA_VX\" value=\"%23.15e\" />",StateJupiter.Velocity[0]*1000.0);
       printf("\n  <CT:setting name=\"JupiterNASA_VY\" value=\"%23.15e\" />",StateJupiter.Velocity[1]*1000.0);
       printf("\n  <CT:setting name=\"JupiterNASA_VZ\" value=\"%23.15e\" />",StateJupiter.Velocity[2]*1000.0);

       printf("\n  <CT:setting name=\"SaturnNASA_X\" value=\"%23.15e\" />",StateSaturn.Position[0]*1000.0);
       printf("\n  <CT:setting name=\"SaturnNASA_Y\" value=\"%23.15e\" />",StateSaturn.Position[1]*1000.0);
       printf("\n  <CT:setting name=\"SaturnNASA_Z\" value=\"%23.15e\" />",StateSaturn.Position[2]*1000.0);
       printf("\n  <CT:setting name=\"SaturnNASA_VX\" value=\"%23.15e\" />",StateSaturn.Velocity[0]*1000.0);
       printf("\n  <CT:setting name=\"SaturnNASA_VY\" value=\"%23.15e\" />",StateSaturn.Velocity[1]*1000.0);
       printf("\n  <CT:setting name=\"SaturnNASA_VZ\" value=\"%23.15e\" />",StateSaturn.Velocity[2]*1000.0);
     }
  else
     {
       /* Print seperator.....................................................*/
       
       for ( i=0 ; i<21 ; i++ ) putchar('-');
       printf("  Case %2s  ",argv[5]);
       for ( i=0 ; i<21 ; i++ ) putchar('-');

       /* Print planet name...................................................*/
       
       switch (Target) {
          case  1: tgtName = "Mercury";
                   break;
          case  2: tgtName = "Venus";
                   break;
          case  3: tgtName = "Earth";
                   break;
          case  4: tgtName = "Mars";
                   break;
          case  5: tgtName = "Jupiter";
                   break;
          case  6: tgtName = "Saturn";
                   break;
          case  7: tgtName = "Uranus";
                   break;
          case  8: tgtName = "Neptune";
                   break;
          case  9: tgtName = "Pluto";
                   break;
          case 10: tgtName = "Moon";
                   break;
          case 11: tgtName = "Sun";
                   break;
          default: tgtName = "Mercury";
                   break;
       }
       
       printf("\n\n  Target:  %s",tgtName);
       printf("\n      JD:  %8.2f",Time);

       /* Print position......................................................*/
       
       printf("\n\n  Position (km):     [1] =  % 22.15e",Position[0]);
       printf("\n                     [2] =  % 22.15e",Position[1]);
       printf("\n                     [3] =  % 22.15e\n\n\n",Position[2]);
     }

  /* Exit normally............................................................*/
  
  exit(0);
}
