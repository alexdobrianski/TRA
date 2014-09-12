/***********************************************************************
     
    2009-2014 (C) Alex Dobrianski tra.cpp oprbit calculation app

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
************************************************************************/


#include "stdafx.h"
#include <afxinet.h>
#include <afxsock.h>

#define _USE_MATH_DEFINES 1
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include "ephem_read.h"

//////////////////////////////////////////////////////////////////////////////
//   predefine vaiable to build different flavor
//#define FIND_IMPULSE_TIME 1
//////////////////////////////////////////////////////////////////////////////

#ifdef _DO_VISUALIZATION
#include "JPEGLIB.H"
void
write_JPEG_file (char * filename, 
				 int quality, 
				 int SizeW, 
				 int SizeH, 
				 int SizeB, 
				 unsigned char *bArray, 
				 J_COLOR_SPACE ColorCode
				 );

#endif

#define MERCURY 0
#define VENUS 1
#define EARTH 2
#define MARS 3
#define JUPITER 4
#define SATURN 5
#define URANUS 6
#define NEPTUNE 7
#define PLUTO 8
#define MOON  9
#define SUN 10

#define PLANET_COUNT 11
    // Earth gravitation model
    // earth GEM-T3 3.98600436E5 6378.137
    // earth 1082.6260745913E-6 -2.5325160653E-6 -1.6185636000E-6 -0.2266690830E-6 0.5390785906E-6
    // C21 -0.0002194691E-6 S21 +0:0015362834E-6 C22 +1.5744102040E-6 S22 0.9037571782E-6
    // earth JGM3 3.986004415E5 6378.1363
	// earth 1082.6360229830E-6 -2.5324353458E-6 -1.6193312052E-6 -0.2277161017E-6 +0.5396484905E-6
    // C21 -0.0002414000E-6 S21 +0.0015431000E-6 C22 +1.5745360428E-6 S22 -0.9038680730E-6
    // see http://www.csr.utexas.edu/publications/statod/TabD.3.new.txt
	// New table containing conventional gravity coefficients generated from
	// Table D.1 (normalized gravity coefficients)
	//Table D.3 JGM-3 Earth Gravity Field (Conventional Coefficients)
	//l	m	            C                                S
	//2	0	  -0.10826360229840D-02	                       0.0
	//3	0	   0.25324353457544D-05	                       0.0
	//4	0	   0.16193312050719D-05	                       0.0
	//5	0	   0.22771610163688D-06	                       0.0
	//6	0	  -0.53964849049834D-06	                       0.0
	//7	0	   0.35136844210318D-06	                       0.0
	//8	0	   0.20251871520885D-06 	               0.0

	//2	1	  -0.24140000522221D-09	       0.15430999737844D-08
	//3	1	   0.21927988018965D-05	       0.26801189379726D-06
	//4	1	  -0.50872530365024D-06	      -0.44945993508117D-06
	//5	1	  -0.53716510187662D-07	      -0.80663463828530D-07
	//6	1	  -0.59877976856303D-07	       0.21164664354382D-07
	//7	1	   0.20514872797672D-06	   	   0.69369893525908D-07
	//8	1	   0.16034587141379D-07	   	   0.40199781599510D-07 

	//2	2	   0.15745360427672D-05	  	  -0.90386807301869D-06
	//3	2	   0.30901604455583D-06	 	  -0.21140239785975D-06
	//4	2	   0.78412230752366D-07	       0.14815545694714D-06
	//5	2	   0.10559053538674D-06	   	  -0.52326723987632D-07
	//6	2	   0.60120988437373D-08		  -0.46503948132217D-07
	//7	2	   0.32844904836492D-07	       0.92823143885084D-08
	//8	2	   0.65765423316743D-08		   0.53813164055056D-08

	//3	3 	   0.10055885741455D-06	       0.19720132389889D-06
	//4	3	   0.59215743214072D-07	  	  -0.12011291831397D-07
	//5	3	  -0.14926153867389D-07	  	  -0.71008771406986D-08
	//6	3	   0.11822664115915D-08	       0.18431336880625D-09
	//7	3	   0.35285405191512D-08	  	  -0.30611502382788D-08
	//8	3	  -0.19463581555399D-09	  	  -0.87235195047605D-09

	//4	4	  -0.39823957404129D-08	   	   0.65256058113396D-08
	//5	4	  -0.22979123502681D-08	   	   0.38730050770804D-09
	//6	4	  -0.32641389117891D-09	   	  -0.17844913348882D-08
	//7	4	  -0.58511949148624D-09	   	  -0.26361822157867D-09
	//8	4	  -0.31893580211856D-09	       0.91177355887255D-10

	//5	5	   0.43047675045029D-09	      -0.16482039468636D-08
	//6	5	  -0.21557711513900D-09	      -0.43291816989540D-09
	//7	5	   0.58184856030873D-12	       0.63972526639235D-11
	//8	5	  -0.46151734306628D-11	       0.16125208346784D-10

	//6	6	   0.22136925556741D-11	       -0.55277122205966D-10
	//7	6	  -0.24907176820596D-10	        0.10534878629266D-10
	//8	6	  -0.18393642697634D-11	        0.86277431674150D-11

	//7	7	   0.25590780149873D-13 	   0.44759834144751D-12 
	//8	7	   0.34297618184624D-12 	   0.38147656686685D-12 

	//8	8	  -0.15803322891725D-12	           0.15353381397148D-12
    	
    // see http://vadimchazov.narod.ru/lepa_zov/lesat.pdf
    // P0 == P[0](x) = 1
    // P1 == P[1](x) = x
    // PN == P[N+1](x) = (x*P[N](x)*(2*n+1) - n*P[N-1](x))/(n+1)
    // sample:
    // P2(x) == P[2](x) = 1/2 *(-1+3*x*x)
    // P3(x) == P[3](x) = 1/2 *(-3*x+5*x*x*x)
    // PNK(x) = pow((1-x*x),k/2)* dK(PN(x)/dxK
    // = pow((1-x*x),k/2)* dK((x*P[N](x)*(2*n+1) - n*P[N-1](x))/(n+1))/dxK
    // P21(x) = 3*x*sqrt(1-x*x)
    // P22(x) = 3*(1-x*x)
    // P31(x) = 3/2*(-1+5*x*x)*sqrt(1-x*x)
    // P32(x) = 15*x*(1-x*x)
    // P33(x) = 15*x*pow((1-x*x),3/2)
    // P41(x) = 5/2*(-3*x+7*x*x*x*)*sqrt(1-x*x)
    // P42(x) = 15/2*(-1+7*x*x*)*sqrt(1-x*x)
	// JN = -sqrt(2*n+1) Cl0 from Table D.3
    // CNK = sqrt(2*(2*n+1)) * sqrt(((n-k)!/(n+k)!) * Clm
    // SNK = sqrt(2*(2*n+1) * sqrt((n-k)!/(n+k)!)) * Slm
    // but coeffs from http://www.csr.utexas.edu/publications/statod/TabD.3.new.txt
    // just needs to copy it
	// FM/r * 
	//( 1 
    //	- J2 * (r0/r)**2 * P2(sinPHI)  
    //	- J3 * (r0/r)**3 * P3(sinPHI) 
    //  - J4 * (r0/r)**4 * P4(sinPHI)
    //  - J5 * (r0/r)**5 * P5(sinPHI)
    //  - J6 * (r0/r)**6 * P6(sinPHI)
    //  - J7 * (r0/r)**7 * P7(sinPHI)
    //  - J8 * (r0/r)**8 * P8(sinPHI)
    //  + SUM2=( (r0/r)**2 * P21(sinPHI) * (C21 * cos(1*Lambda) + S21*sin(1*Lanbda)) +
    //           (r0/r)**2 * P22(sinPhi) * (C22 *cos(2*Lambda) + S22*sin(2*Lambda))     )
    //    SUM3=( (r0/r)**3 *P31(sinPHI) * (C31 * cos(1*Lambda) + S31*sin(1*Lambda)) +
    //           (r0/r)**3 *P32(sinPHI) * (C32 * cos(2*Lambda) + S32*sin(2*Lambda)) +
    //           (r0/r)**3 *P33(sinPHI) * (C33 * cos(3*Lambda) + S33*sin(3*Lambda))      )
    //  SUM4= ( (r0/r)**4 *P41(sinPHI) * (C41 * cos(1*Lambda) + S41*sin(1*Lambda)) +
    //          (r0/r)**4 *P42(sinPHI) * (C42 * cos(2*Lambda) + S42*sin(2*Lambda)) +
    //          (r0/r)**4 *P43(sinPHI) * (C43 * cos(3*Lambda) + S43*sin(3*Lambda)) +
    //          (r0/r)**4 *P44(sinPHI) * (C44 * cos(4*Lambda) + S44*sin(3*Lambda))     )
    //  SUM5= ( (r0/r)**5 *P51(sinPHI) * (C51 * cos(1*Lambda) + S51*sin(1*Lambda)) +
    //          (r0/r)**5 *P52(sinPHI) * (C52 * cos(2*Lambda) + S52*sin(2*Lambda)) +
    //          (r0/r)**5 *P53(sinPHI) * (C53 * cos(3*Lambda) + S53*sin(3*Lambda)) +
    //          (r0/r)**5 *P54(sinPHI) * (C54 * cos(4*Lambda) + S54*sin(4*Lambda)) +
    //          (r0/r)**5 *P55(sinPHI) * (C55 * cos(5*Lambda) + S55*sin(5*Lambda))         )
    //  SUM6= ( (r0/r)**6 *P61(sinPHI) * (C61 * cos(1*Lambda) + S61*sin(1*Lambda)) +
    //          (r0/r)**6 *P62(sinPHI) * (C62 * cos(2*Lambda) + S62*sin(2*Lambda)) +
    //          (r0/r)**6 *P63(sinPHI) * (C63 * cos(3*Lambda) + S63*sin(3*Lambda)) +
    //          (r0/r)**6 *P64(sinPHI) * (C64 * cos(4*Lambda) + S64*sin(4*Lambda)) +
    //          (r0/r)**6 *P65(sinPHI) * (C65 * cos(5*Lambda) + S65*sin(5*Lambda)) +
    //          (r0/r)**6 *P66(sinPHI) * (C66 * cos(6*Lambda) + S66*sin(6*Lambda))         )
    //  SUM7= ( (r0/r)**7 *P71(sinPHI) * (C71 * cos(1*Lambda) + S71*sin(1*Lambda)) +
    //          (r0/r)**7 *P72(sinPHI) * (C72 * cos(2*Lambda) + S72*sin(2*Lambda)) +
    //          (r0/r)**7 *P73(sinPHI) * (C73 * cos(3*Lambda) + S73*sin(3*Lambda)) +
    //          (r0/r)**7 *P74(sinPHI) * (C74 * cos(4*Lambda) + S74*sin(4*Lambda)) +
    //          (r0/r)**7 *P75(sinPHI) * (C75 * cos(5*Lambda) + S75*sin(5*Lambda)) +
    //          (r0/r)**7 *P76(sinPHI) * (C76 * cos(6*Lambda) + S76*sin(6*Lambda)) +
    //          (r0/r)**7 *P76(sinPHI) * (C77 * cos(7*Lambda) + S77*sin(7*Lambda))        )
    //  SUM8= ( (r0/r)**8 *P81(sinPHI) * (C81 * cos(1*Lambda) + S81*sin(1*Lambda)) +
    //          (r0/r)**8 *P82(sinPHI) * (C82 * cos(2*Lambda) + S82*sin(2*Lambda)) +
    //          (r0/r)**8 *P83(sinPHI) * (C83 * cos(3*Lambda) + S83*sin(3*Lambda)) +
    //          (r0/r)**8 *P84(sinPHI) * (C84 * cos(4*Lambda) + S84*sin(4*Lambda)) +
    //          (r0/r)**8 *P85(sinPHI) * (C85 * cos(5*Lambda) + S85*sin(5*Lambda)) +
    //          (r0/r)**8 *P86(sinPHI) * (C86 * cos(6*Lambda) + S86*sin(6*Lambda)) +
    //          (r0/r)**8 *P87(sinPHI) * (C87 * cos(7*Lambda) + S86*sin(7*Lambda)) +
    //          (r0/r)**8 *P88(sinPHI) * (C88 * cos(8*Lambda) + S86*sin(8*Lambda))         )

    // moon 0.0002027
#define MAX_COEF_J 9
double GM_MODEL = 3.986004415E5;
//double R0_MODEL = 6378136.30; // in m
#define R0_MODEL 6378136.30
//#define R0_MODEL 6378000.00
//#define R0_MODEL 6379137.0
double Clm[MAX_COEF_J][MAX_COEF_J] = {
    0.0,                  0.0,-0.10826360229840e-02, 0.25324353457544E-05, 0.16193312050719e-05, 0.22771610163688E-06,-0.53964849049834e-06, 0.35136844210318e-06, 0.20251871520885e-06,
	0.0,                  0.0,-0.24140000522221e-09, 0.21927988018965e-05,-0.50872530365024e-06,-0.53716510187662e-07,-0.59877976856303e-07, 0.20514872797672e-06, 0.16034587141379e-07,
	0.0,                  0.0, 0.15745360427672e-05, 0.30901604455583e-06, 0.78412230752366e-07, 0.10559053538674e-06, 0.60120988437373e-08, 0.32844904836492e-07, 0.65765423316743e-08,
	0.0,                  0.0,                  0.0, 0.10055885741455e-06, 0.59215743214072e-07,-0.14926153867389e-07, 0.11822664115915e-08, 0.35285405191512e-08,-0.19463581555399e-09,
    0.0,                  0.0,                  0.0,                  0.0,-0.39823957404129e-08,-0.22979123502681e-08,-0.32641389117891e-09,-0.58511949148624e-09,-0.31893580211856e-09,
    0.0,                  0.0,                  0.0,                  0.0,                  0.0, 0.43047675045029e-09,-0.21557711513900e-09, 0.58184856030873e-12,-0.46151734306628e-11, 
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0, 0.22136925556741e-11,-0.24907176820596e-10,-0.18393642697634e-11, 
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0, 0.25590780149873e-13, 0.34297618184624e-12,
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,-0.15803322891725e-12};
double Slm[MAX_COEF_J][MAX_COEF_J] = {
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,
    0.0,                  0.0, 0.15430999737844e-08, 0.26801189379726e-06,-0.44945993508117e-06,-0.80663463828530e-07, 0.21164664354382e-07, 0.69369893525908e-07, 0.40199781599510e-07,
    0.0,                  0.0,-0.90386807301869e-06,-0.21140239785975e-06, 0.14815545694714e-06,-0.52326723987632e-07,-0.46503948132217e-07, 0.92823143885084e-08, 0.53813164055056e-08,
    0.0,                  0.0,                  0.0, 0.19720132389889e-06,-0.12011291831397e-07,-0.71008771406986e-08, 0.18431336880625e-09,-0.30611502382788e-08,-0.87235195047605e-09, 
    0.0,                  0.0,                  0.0,                  0.0, 0.65256058113396e-08, 0.38730050770804e-09,-0.17844913348882e-08,-0.26361822157867e-09, 0.91177355887255e-10,
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,-0.16482039468636e-08,-0.43291816989540e-09, 0.63972526639235e-11, 0.16125208346784e-10,
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,-0.55277122205966e-10, 0.10534878629266e-10, 0.86277431674150e-11, 
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0, 0.44759834144751e-12, 0.38147656686685e-12, 
    0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0,                  0.0, 0.15353381397148e-12};
    
typedef struct TraObj
{
    int Elem;
    int flInUse[PLANET_COUNT];
    long double X[PLANET_COUNT];
    long double Y[PLANET_COUNT];
    long double Z[PLANET_COUNT];
    long double VX[PLANET_COUNT];
    long double VY[PLANET_COUNT];
    long double VZ[PLANET_COUNT];

    long double X_[PLANET_COUNT];
    long double VX_[PLANET_COUNT];
    long double FX[PLANET_COUNT];
    long double X__[PLANET_COUNT];
    long double Y_[PLANET_COUNT];
    long double VY_[PLANET_COUNT];
    long double FY[PLANET_COUNT];
    long double Y__[PLANET_COUNT];
    long double Z_[PLANET_COUNT];
    long double VZ_[PLANET_COUNT];
    long double FZ[PLANET_COUNT];
    long double Z__[PLANET_COUNT];

    long double GM[PLANET_COUNT];
    long double M[PLANET_COUNT];
    long double GMxM[PLANET_COUNT][PLANET_COUNT];
    long double Distance[PLANET_COUNT][PLANET_COUNT];
    long double Distance2[PLANET_COUNT][PLANET_COUNT];
    long double ForceDD[PLANET_COUNT][PLANET_COUNT];
    char Kepler1[PLANET_COUNT][100];
    char Kepler2[PLANET_COUNT][100];
    char Kepler3[PLANET_COUNT][100];
    // 3 punch card calculation helper vars
    //
    long double ProbEpoch[PLANET_COUNT];
    long double ProbEpochS[PLANET_COUNT];
    long double ProbMeanMotion[PLANET_COUNT];
    long double ProbFirstDervMeanMotion[PLANET_COUNT];
    long double ProbSecondDervmeanMotion[PLANET_COUNT];
    long double ProbDragterm[PLANET_COUNT];
    unsigned char ProbElementSetType[PLANET_COUNT];
    long double ProbIncl[PLANET_COUNT];
    long double ProbAscNode[PLANET_COUNT];
    long double ProbEcc[PLANET_COUNT];
    long double ProbArgPer[PLANET_COUNT];
    long double ProbMeanAnom[PLANET_COUNT];
    long double ProbTPeriod[PLANET_COUNT];
    long double ProbTDays[PLANET_COUNT];
    long double ProbTSec[PLANET_COUNT];
    long double ProbRevAtEpoch[PLANET_COUNT];
    // satelitte close to body
    // it can be only one body
    int iLeg;
    int iLeg_longit;
    int LegBody;
    long double J[MAX_COEF_J];
    long double CNK[MAX_COEF_J][MAX_COEF_J];
    long double SNK[MAX_COEF_J][MAX_COEF_J];
    long double CosTetta[MAX_COEF_J];
    long double SinTetta[MAX_COEF_J];
    long double P[MAX_COEF_J];
    long double Rn1divR[MAX_COEF_J];
    long double Ptilda[MAX_COEF_J];
    long double Pnk_tilda[MAX_COEF_J][MAX_COEF_J];
    //double Pnk[MAX_COEF_J][MAX_COEF_J];
    long double Qnk[MAX_COEF_J][MAX_COEF_J];
    long double R0divR[MAX_COEF_J];
    long double ForceDD_;
    long double Lambda;
    long double DeltaVX[PLANET_COUNT];
    long double DeltaVY[PLANET_COUNT];
    long double DeltaVZ[PLANET_COUNT];
    long double Xk[MAX_COEF_J];
    long double Yk[MAX_COEF_J];
    long double XkDxr[MAX_COEF_J];
    long double YkDxr[MAX_COEF_J];
    long double XkDyr[MAX_COEF_J];
    long double YkDyr[MAX_COEF_J];
    long double tempValX;
    long double tempValY;
    long double tempValZ;
    long double tempValR;
    void CalcCosK(void)
    {
        
    }
    void CalcP(double sinTetta, double XdivR, double YdivR, double ValX, double ValY, double ValZ, double ValR)
    {
        int n,k;
        tempValX = ValX; tempValY = ValY; tempValZ = ValZ; tempValR = ValR;
        // sinTetta is a Z/R - tetta is Latitude
        // XdivR, YdivR - must be rotated from original X,Y based on a time (and angle Lambda) 
        //
        //
        //Lambda +=1.0471975511965977461542144610932;//1.075; //0.183;//0.36;//1.72944494;
        //Lambda += 1.5707963267948966192313216916398;//1.2337005501361698273543113749845;
        //Lambda = -Lambda;
        //if (Lambda != -2)
        {
            XdivR = cos(Lambda) * XdivR - sin(Lambda) * YdivR;
            YdivR = sin(Lambda) * XdivR + cos(Lambda) * YdivR;
        }

        SinTetta[1] = sinTetta;
        if (ValX > 0.0)
            CosTetta[1] = cos(asin(sinTetta));
        else 
            CosTetta[1] = -cos(asin(sinTetta));
        // power of a cos
        for (k = 2; k <= iLeg; k++)
        {
            CosTetta[k] = CosTetta[k-1]*CosTetta[1];
            SinTetta[k] = SinTetta[k-1]*SinTetta[1];
        }
        // legandr functions from sinTetta
        P[0] = 1.0;
        P[1] = sinTetta;
        Pnk_tilda[0][0] = P[0];
        Pnk_tilda[1][0] = P[1];
        for (n = 2; n <=iLeg; n++)
        {
            // page 90- formula 8
            // (n + 1)Pn+1(z) - (2n + 1)zPn(z) + nPn-1(z) = 0 =>
            // (n + 1)Pn+1(z) = (2n + 1)zPn(z) - nPn-1(z) =>
            // Pn+1(z) = ((2n + 1)zPn(z) - nPn-1(z))/(n + 1) =>
            // or Pn(z) = ((2(n-1) + 1)zPn-1(z) - (n-1)Pn-2(z))/(n-1 + 1)
            // or Pn(z) = ((2n-1)  zPn-1(z) - (n-1)Pn-2(z))/n
            

            //P[n] = ((2.0* (n-1) +1) *sinTetta * P[n-1] - (n-1)*P[n-2])/((n-1)+1);
            P[n] = ((2.0* n-1.0) *sinTetta * P[n-1] - (n-1)*P[n-2])/n;
            Pnk_tilda[n][0] = P[n];
            //break;
            //}
        }
        // derivetiveas from legandr fucntions
        Ptilda[0] = 0; // P0 was constant == P0' == 0
        Ptilda[1] = 1; // P1 was a X == P1' == 1

        Pnk_tilda[0][1] = Ptilda[0];
        Pnk_tilda[1][1] = Ptilda[1];
        for (n=2;n<=iLeg;n++)
        {   // all derivatives
            // for P2' == 2* P1 + sin(tetta) * 1
            Ptilda[n] = n * P[n-1] + sinTetta * Ptilda[n-1];
            Pnk_tilda[n][1] = Ptilda[n];
        }
        // derivatives for dK(Pn(sinTetta))/d(sinTetta)**K
        for (n= 2;n <=iLeg;n++)
        {
            for (k = 1; k<=iLeg;k++)
            {
                Pnk_tilda[n][k] = (2*n-1) * Pnk_tilda[n-1][k-1] - Pnk_tilda[n-2][k];
            }
        }
        //for (n= 2;n <=iLeg;n++)
        //{
        //    for (k = 1; k<=iLeg;k++)
        //    {
        //        //Pnk[n][k] = pow((1 - X*X),k/2)*Pnk_tilda[n][k];
        //        Pnk[n][k] = CosTetta[k]*Pnk_tilda[n][k];
        //    }
        //}
        //if (iLeg_longit) // this is case when need to acount Longitude
        // last formula on page 92
        Xk[0] = 1; Yk[0] =0; 
        Xk[1] = XdivR; Yk[1] = YdivR;
        for (k = 2; k<=iLeg; k++)
        {
            Xk[k] = Xk[k-1]*XdivR - Yk[k-1] * YdivR;
            Yk[k] = Yk[k-1]*XdivR + Xk[k-1] * YdivR;
        }
        XkDxr[0] = 0; XkDyr[0] = 0; YkDxr[0] = 0; YkDyr[0] = 0;
        for (k = 1; k <=iLeg; k++)
        {
            XkDxr[k] = XkDxr[k-1]*XdivR + Xk[k-1] - YkDxr[k-1]*YdivR;
            XkDyr[k] = XkDyr[k-1]*XdivR - YkDyr[k-1]*YdivR - Yk[k-1];
            YkDxr[k] = YkDxr[k-1]*XdivR + Yk[k-1] + XkDxr[k-1]*YdivR;
            YkDyr[k] = YkDyr[k-1]*XdivR + XkDyr[k-1]*YdivR+Xk[k-1];
        }
        {
            // formula 8 on page 92
             for (n = 2; n <=iLeg; n++)
             {
                 for (k = 0; k <=iLeg; k++)
                 {
                     Qnk[n][k] = CNK[n][k] * Xk[k] + SNK[n][k] * Yk[k];
                 }
             }
        }
    }
    //void CalcPNK(double X)
    //{
    //    LastPNK[2][1] = 3*X*sqrt(1-X*X);
    //    LastPNK[2][2] = 3*(1-X*X);
    //    // rest skiped
    //}
    void SummXYZ(double &X, double &Y, double &Z)
    {
        int n,k;
        X=0; Y=0; Z=0;
        for (n= 2;n <=iLeg;n++)
        {
            // implementation of derivetive (page 91 on Aksenov lectures)
            // Potential Unk:
            // Unk = fm/r (r0/r)**n Pnk(sinTetta) [Cnk*cos(k*lambda) + Snk*sin(k*lambda)
            // then Unk separated into 3 parts:
            // (1) Rn(1/r) = fm/r *(r0/r)**n         => Rn(1/r)' = (n+1) * fm * r0**n * (1/r)**n = (n+1) * fm * (r0/r)**n
            //
            //  on a page 12 is (kaind hard to write formulas in C code)
            // Pnk(sinTetta) = (1-sinTetta**2) ** k/2 * dk (Pn(sinTetta))/ d(sinTetta)**k ==
            //   cosTetta **k * dK(Pn(sinTetta))/ d(sinTtta)**k
            // now cosTetta**k separated from Pnk(sinTetta) and a second part is:
            // (2) Znk(z/r) = dK(Pn(sinTetta)/ d(sinTetta)**k == Pnk_tetta[n][k]  => Znk(z/r)' = Pnk_tilda[n][k+1]
            // also (2) for K == 0 is:
            // Zn0 = Pn(sinTetta)
            // and cosTetta goes to a last part:
            // (3) Qnk(x/r,y/r) = cosTetta**k * [Cnk*cos(k*lambda) + Snk * sin(k*lambda)] ==
            // [Cnk*cosTetta**k * cos(k*lambda) + Snk * cosTetta**k * sin(k*lambda)] =
            // [Cnk* Xk + Snk*Yk] and 
            // Xk = (cosTetta)**k * cos(k*Lambda) 
            // Yk = (costetta)**k * sin(k*Lambda)
            // or recurcevly:
            // X0 =1, Y0 = 0
            // X1 = cosTeatta * Cos(Lambda) = x/r; Y1 = cosTetta * sin(Lambda) = y/r
            // X[k+1] = Xk * x/r - Yk* y/r  Y[k+1] = Yk* x/r  + Xk * y/r
            //  
            // X2 = X1 * x/r - Y1 * y/r       = (x/r)**2 - (y/r)**2;
            // Y2 = Y1 *x/r + X1 *y/r         = y/r * x/r + x/r* y/r = 2 (y/r) * (x/r)
            //
            // Q20 = C20 (constant) => D(Q20)/D(x/r) = 0 ; D(Q20)/ D (y/r) = 0
            // Q21 = C21 * (x/r) + S21 (y/r) => 
            //       D(Q21)/D(x/r) = C21; 
            //       D(Q21)/D(y/r) = S21
            // Q22 = C22 * [(x/r)*(x/r) - (y/r)*(y/r)] + S22 * [2*x/r*y/r] 
            //       D(Q22)/D(x/r) = 2 * C22 * x/r + 2 * S22 * y/r = 2*[C22*x/r + S22*y/r]
            //       D(Q22)/D(y/r) =-2 * C22 * y/r + 2 * S22 * x/r = 2*[C22*y/r + S22*x/r]
            // Q30 = C30 (constant) => D(Q30)/D(x/r) = 0; D(Q30)/D(y/r) = 0;
            // Q31 = C31 * (x/r) + S31 (y/r) => 
            //       D(Q31)/D(x/r) = C31; 
            //       D(Q31)/D(y/r) = S31
            // Q32 = C32 * [(x/r)*(x/r) - (y/r)*(y/r)] + S32 * [2*x/r*y/r] 
            //       D(Q32)/D(x/r) = 2 * C32 * x/r + 2 * S32 * y/r = 2*[C22*x/r + S22*y/r]
            //       D(Q32)/D(y/r) =-2 * C32 * y/r + 2 * S32 * x/r = 2*[C22*y/r + S22*x/r]
            // Q33 = C33 * [(x/r)*(x/r) - (y/r)*(y/r)] + S33 * [2*x/r*y/r] 
            //       D(Q32)/D(x/r) = 2 * C32 * x/r + 2 * S32 * y/r = 2*[C22*x/r + S22*y/r]
            //       D(Q32)/D(y/r) =-2 * C32 * y/r + 2 * S32 * x/r = 2*[C22*y/r + S22*x/r]
            // chastnie proizvodnie (simbol D) :
            // D(Unk)/D(x) = d(Rn) / d(1/r)    * D(1/r)/D(x)    * Znk    * Qnk +
            //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(x) * Qnk +
            //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(x) + D(Qnk)/D(y/r)*D(y/r)/D(x)]
            //            = (n+1)* fm * (r0/r)**n * (-x/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-x*z/r**3) * Qnk +
            //              Rn * Znk * [D(Qnk)/D(x/r) * (1/r - x**2/r**3) + D(Qnk)/D(y/r)*(-x*y/r**3)]
            //            = (n+1)* fm * (r0/r)**n * (-x/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-x*z/r**3) * Qnk +
            //              Rn * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (1/r - x**2/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-x*y/r**3)]
            // as a result for example n=2 and k = 0
            // D(U20)/D(X) = (2+1)* fm * (r0/r)**2 * (-x/r**3) * d(P2(sinTetta))/d(sinTetta) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(0+1)(P2(sinTetta)/d(sintetta)**(0+1) *(-x*z/r**3) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(Pn(sinTetta))/d(sinTetta) * [D(Q20)/D(x/r) * (1/r - x**2/r**3) + D(Q20)/D(y/r)*(-x*y/r**3)]
            //            = - 3 * fm * (ro/r)**2 * (x/r**3) * Pnk_tilda[2][0] * J2 
            //              -     fm * (r0/r)**2 * (1/r) * Pnk_tilda[2][1] * x*z/r**3 *J2 
            //                    fm * (r0/r)**2 * (1/r) Pnk_tilda[2][0] * [0*(1/r-x**2/r**3) + 0*(-x*y/r**3)]
            //X += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//CNK[n][0] 
            //     - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            X += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]// * tempValX/(tempValR)
                 - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]// * tempValX* tempValZ/(tempValR*tempValR)
                 ;

            // D(Unk)/D(y) = d(Rn) / d(1/r)    * D(1/r)/D(y)    * Znk    * Qnk +
            //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(y) * Qnk +
            //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)*D(y/r)/D(y)]
            //            = (n+1)* fm * (r0/r)**n * (-y/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-y*z/r**3) * Qnk +
            //              Rn * Znk * [D(Qnk)/D(x/r) * (- x*y/r**3) + D(Qnk)/D(y/r)*(1/r-y**2/r**3)]
            //            = (n+1)* fm * (r0/r)**n * (-y/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-y*z/r**3) * Qnk +
            //              Rn * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x*y/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(1/r-y**2/r**3)]
            // as a result for example n=2 and k = 0
            // D(U20)/D(X) = (2+1)* fm * (r0/r)**2 * (-y/r**3) * d(P2(sinTetta))/d(sinTetta) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(0+1)(P2(sinTetta)/d(sintetta)**(0+1) * (-y*z/r**3) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(Pn(sinTetta))/d(sinTetta) * [D(Q20)/D(x/r) * ( - x*y/r**3) + D(Q20)/D(y/r)*(1/r-y**2/r**3)]
            //            = - 3 * fm * (ro/r)**2 * (y/r**3) * Pnk_tilda[2][0] * J2 
            //              -     fm * (r0/r)**2 * (1/r) * Pnk_tilda[2][1] * y*z/r**3 *J2 
            //                    fm * (r0/r)**2 * (1/r) Pnk_tilda[2][0] * [0*(-x*y/r**3) + 0*(1/r-y**2/r**3)]
            //Y += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//CNK[n][0]
            //     - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            Y += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//* tempValY/(tempValR)
                 - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]//* tempValY* tempValZ/(tempValR*tempValR)
            ;

            // D(Unk)/D(z) = d(Rn) / d(1/r)    * D(1/r)/D(z)    * Znk    * Qnk +
            //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(z) * Qnk +
            //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(z) + D(Qnk)/D(y/r)*D(y/r)/D(z)]
            //            = (n+1)* fm * (r0/r)**n * (-z/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (1/r-z**2/r**3) * Qnk +
            //              Rn * Znk * [D(Qnk)/D(x/r) * (- x*z/r**3) + D(Qnk)/D(y/r)*(-y*z/r**3)]
            //            = (n+1)* fm * (r0/r)**n * (-z/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (1/r-z**2/r**3) * Qnk +
            //              Rn * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x*z/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-y*z/r**3)]
            // as a result for example n=2 and k = 0
            // D(U20)/D(z) = (2+1)* fm * (r0/r)**2 * (-z/r**3) * d(P2(sinTetta))/d(sinTetta) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(0+1)(P2(sinTetta)/d(sintetta)**(0+1) * (1/r-z**2/r**3) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(Pn(sinTetta))/d(sinTetta) * [D(Q20)/D(x/r) * ( - x*z/r**3) + D(Q20)/D(y/r)*(-y*z/r**3)]
            //            = - 3 * fm * (ro/r)**2 * (z/r**3) * Pnk_tilda[2][0] * J2 
            //                   fm * (r0/r)**2 * (1/r) * Pnk_tilda[2][1] * (1/r - z**2/r**3 *J2 
            //                    fm * (r0/r)**2 * (1/r) Pnk_tilda[2][0] * [0*(-x*z/r**3) + 0*(-y*z/r**3)]
            //Z += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//CNK[n][0] 
            //+ R0divR[n] * CosTetta[2]/ SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            Z += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//* tempValZ/(tempValR)
            + R0divR[n] * CosTetta[2]/ SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]//* (1 - tempValZ* tempValZ/(tempValR*tempValR))
            ;
        }
        // second round
        
        for (n= 2;n <=iLeg;n++)
        {
            // k=0 already done
            for (k = 1; k<=n; k++)
            {
                //            = (n+1)* fm * (r0/r)**n * (-x/r**3) * Znk * Qnk +
                //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-x*z/r**3) * Qnk +
                //              Rn * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (1/r - x**2/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-x*y/r**3)]
                X += -(n+1) * R0divR[n] * Pnk_tilda[n][k] * Qnk[n][k]// * tempValX/(tempValR)
                    - R0divR[n] * SinTetta[1] * Pnk_tilda[n][k+1] * Qnk[n][k]// * tempValX* tempValZ/(tempValR*tempValR)
                    ; 
                Y += -(n+1) * R0divR[n] * Pnk_tilda[n][k] * Qnk[n][k]// * tempValY/(tempValR)
                    - R0divR[n] * SinTetta[1] * Pnk_tilda[n][k+1] * Qnk[n][k]//* tempValY* tempValZ/(tempValR*tempValR)
                ; 
                Z += -(n+1) * R0divR[n] * Pnk_tilda[n][k] * Qnk[n][k]// * tempValZ/(tempValR)
                    + R0divR[n] * CosTetta[2]/ SinTetta[1] * Pnk_tilda[n][k+1] * Qnk[n][k]//* (1 - tempValZ* tempValZ)/(tempValR*tempValR)
                ; 

            }
        }
        
    }
    double SummJ(void)
    {
        // this is a formula for gravitation potential (NOT a force)
        // to get a force needs to get analiticaly derivitive from Gravitation potencial
        // page 90 on Aksenov lectures : http://vadimchazov.narod.ru/lepa_zov/lesat.pdf
        //	- J2 * (r0/r)**2 * P2(sinPHI)  
        //	- J3 * (r0/r)**3 * P3(sinPHI) 
        //  - J4 * (r0/r)**4 * P4(sinPHI)
    
        double Summ = 1.0;
        for (int n = 2; n <=iLeg; n++)
        {
            Summ += J[n]*R0divR[n]*P[n];
        }
        // for now no earth rotation acounted
        //+ SUM2=( (r0/r)**2 * P21(sinPHI) * (C21 * cos(1*Lambda) + S21*sin(1*Lanbda)) +
        //           (r0/r)**2 * P22(sinPhi) * (C22 *cos(2*Lambda) + S22*sin(2*Lambda))     )
        //Summ += R0divR[2] * LastPNK[2][1] * (CNK[2][1] * cos(Lambda) + SNK[2][1]*sin(Lambda)) +
        //       R0divR[2]*LastPNK[2][2]*(CNK[2][2]*cos(2.0*Lambda) + SNK[2][2]*sin(2.0*Lambda));
        for (int n= 2; n <= iLeg; n++)
        {
            for (int k =1; k <=n; k++)
            {
                Summ += R0divR[n] * Ptilda[n] *(CNK[n][k]*cos((double)k*Lambda) +SNK[n][k]*sin((double)k*Lambda));
            }
        }

        return Summ;
    }

} TRAOBJ, *PTRAOBJ;

#define MAX_IMPULSE_LINES 100
typedef struct TraImplObj
{
    int iLine;
    int iEngineOnSatellite;
    int EngineOn;
    int EngineDone;
    int ImplsPointer;
    int iCalculate;
    double TotalImpulse;
    double Weight;
    double DeltaTime;
    int IteraPerSec;
    double FireTime;
    double Ang1;
    double Ang2;
    double XVec;
    double YVec;
    double ZVec;
    double ValImpl[MAX_IMPULSE_LINES];
    int NearBody;
    double TotalWeight;
    int AngleType;
    int AngleOnBody;
    double OptimizationInitialStep;
    double OptimizationDecCoef;
    double OptimizationInitialStepCopy;
    double OptimizationDecCoefCopy;
    double OptimizationStop;
    int OptimizationFirstDirectionSwitch;
    double SeartchForPeriod;
    int iCountApogPerig;
} TRAIMPLOBJ, *PTRAIMPLOBJ;

typedef struct TraOptimObj
{
    double FireTime; // firing time of an engine
    double Ang1; // first angle
    double Ang2; // second angle
    double XVec; // direction firing vector (X component)
    double YVec; // direction (Y component)
    double ZVec; // direction (Z component)
    int NearBody; 
    int AngleType;
    int AngleOnBody;
    //int OptimizationFirstDirectionSwitch;
    int EngineToOptimize;
    int TrajectoryOptimizationType;
    int LastEngine;
    int Calculate;
    double OptimizationInitialStep;
    double OptimizationDecCoef;
    double OptimizationStop;
    int iNumberOfTryValues;
#define _VAL_TRY 30*24
    double dValTry[_VAL_TRY]; 
    double dValTryMaxMin[_VAL_TRY]; 
    //double SeartchForPeriod;
    double Period;
} TRAOPTIMOBJ, *PTRAOPTIMOBJ;

int EngineToOptimize = 4; // 0 == first engine firing - search for apogee 
                       // 1 == second engine firing - search for perigee
                       // 2 == third engine firing 
                       // 3 == 4th impulse
                       // 4 == 5th impulse
    
    
int TrajectoryOptimizationType = 0; 
          //1 - search for a minimum by adjusting time of firing
          //2 - search for a maximum by adjusting time of firing
          //3 - search for minimum by adjusting time of firing
          //4 - search fo minimum by adjusting angle of firing 

#define MINIMUM_BY_TIME 1
#define MAXIMUM_BY_TIME 2
#define MINIMUM_BY_ANGLE 4
#define MINIMUM_BY_WEIGHT 6

#define CALC_PERIGEE 0
#define CALC_APOGEE 1
#define CALC_TARGET_PRACTICE 3
#define CALC_TARGET_POINT 5
#define CALC_PERIOD 6 
#define CALC_MAX_FROM_EARTH 7
#define CALC_INIT_PERIOD 8
#define CALC_FIRE_FIRST_ENGINE_TIME 9
#define CALC_FIRE_SECOND_ENGINE_TIME 10
#define CALC_FIRE_THIRD_ENGINE_TIME_TRY_ONE 11
#define CALC_AT_APOGEE_DIFF_TO_3_4_DIST 12

#define COPYKEPLER(a1,b1,c1) memset(szTempo, 0, sizeof(szTempo)); memcpy(szTempo, b1, c1);  if (szTempo[0] == ' ') {szTempo[0] = '0';if (szTempo[1] == ' ') {szTempo[1] ='0';if (szTempo[2] == ' '){szTempo[2] ='0';}}}a1 = atof(szTempo);


long OldCurentTimeSec;
int OldCurentTimeTD;
int OldCurentIteraPerSec;

double dMinFromNow = 3.0;
double dStartJD = 0.0;//2451544.5; // if value dStartJD not set (==0.0) then use value from keplers elements of a satelite 0
// if it will be more then one satellite needs to set this value to a last epoch of all satellites

double SunX =.0;
double SunY =.0;
double SunZ = .0;
double SunM = 1.9891E30;
double GMSun;
double SunR;
double GMEarth;
double GMEarthMoon;
double GMMoon;

double AU;
double AUcalc;
double EarthSmAxAU;

TRAOBJ SolarSystem;
TRAOBJ Sat;
#define MAX_OPTIM 30
int MaxOptim = MAX_OPTIM;
int StartOptim = 0;
TRAOPTIMOBJ Opt[MAX_OPTIM];
int iOptPtr = 0;
int iOptimizationStep = 0;


double MinMaxX = .0;
double MinMaxY = .0;
double MinMaxZ = .0;

int iFirstMinMax = 0;
double FirstMinMaxX2 = .0;
double FirstMinMaxY2 = .0;
double FirstMinMaxZ2 = .0;

double DeltaMinMaxD = 1.0E60;


double EarthX = .0;
double EarthY = .0;
double EarthZ = .0;
double EarthX2 = .0;
double EarthY2 = .0;
double EarthZ2 = .0;
double EarthR = 6371000.0;
double EarthRP;
double EarthRE;
double EarthM = 5.9736E24;
double MassRatioSunToEarthPlusMoon;
double EarthTSolSec;
double EarthCurTime;
double EarthCurTimeS;
double EarthSmAx;
double EarthTDays;
double EarthTSec;
double EarthCalcKepler;

double EarthVX = .0;
double EarthVY = .0;
double EarthVZ = .0;

long double GST,SLONG,SRASN, SDEC;
long double GreenwichA = 0.0;



double MoonX = 363104000.0;
double MoonY = .0;
double MoonZ = .0;

double MoonR = 1737100.0;
double MoonRP;
double MoonRE;
double MoonM = 7.3477E22;
double MassRatioEarthToMoon;

double MoonVX = .0;
double MoonVY = 1022.0;
double MoonVZ = .0;
double MoonTSec;
double MoonCurTime;
double MoonCurTimeS;

double MoonTPeriod;
int MoonKeplerDone = 0;


double TimeSl = 0;//0.01;
double TimeSlOld = 0;//0.01;
double TimeSl_2 = 0;
double StartLandingIteraPerSec = 0.0;
int iStartLandingIteraPerSec = 0;

double Gbig = 0;//6.6725E-11;

BOOL OutLast = FALSE;
double TotalDays;
long iTotalSec;
double IterPerSec;

char szMoonKeplerLine1[1024];
char szMoonKeplerLine2[1024];
char szMoonKeplerLine3[1024];

int EnginesCount = 0;
#define MAX_ENGINES 6
TRAIMPLOBJ Engine[MAX_ENGINES];
int LastEngine = 0;
int iItaration = 0;


// ground stations
double GrLat[10];
double GrLong[10];
typedef struct tagPulsars
{
    int N;
    char Name[20];
    double ELONG;
    double ELAT;
    double P0;
    double S400mJy;

} PULSARS, *PPULSARS;

#define NPULSARS 150
int nPulsars= 0;
PULSARS Pulsars[150];

void IteraSolarSystem(int TimeDirection, TRAOBJ * SlS)
{
    int i;
    int j;
	//SUN_08 (1950,1,0,0,0,GST,SLONG,SRASN,SDEC);
	//GreenwichA = GreenwichAscension(50001.0);

#define PlanetsCount SlS->Elem
//#define PlanetsCount PLANET_COUNT
    // calculation of a forces
    for (i = 0; i < SlS->Elem; i++)
    {
        if (SlS->flInUse[i] ==0)
            continue;
        for (j = i; j < SlS->Elem; j++)
        {
            if (i == j) 
                continue;
            if (SlS->flInUse[j] ==0)
                continue;
            double tD_Obj1Obj2 = (SlS->X[i] - SlS->X[j])*(SlS->X[i] - SlS->X[j]) + 
					    (SlS->Y[i] - SlS->Y[j])*(SlS->Y[i] - SlS->Y[j]) + 
					    (SlS->Z[i] - SlS->Z[j])*(SlS->Z[i] - SlS->Z[j]);

            if (tD_Obj1Obj2 != SlS->Distance2[i][j])
            {
                SlS->Distance2[i][j] = tD_Obj1Obj2;
                SlS->Distance2[j][i] = tD_Obj1Obj2;
                double tD_ = sqrt(tD_Obj1Obj2);

                if (tD_ != SlS->Distance[i][j])
                {
                    SlS->Distance[i][j] = tD_;
                    SlS->Distance[j][i] = tD_;
#if 0
                    SlS->ForceDD[i][j] = SlS->GM[i] * SlS->M[j] / SlS->Distance2[i][j];
                    SlS->ForceDD[j][i] = SlS->GM[j] * SlS->M[i] / SlS->Distance2[j][i];
#else
                    SlS->ForceDD[i][j] = SlS->GMxM[i][j] / SlS->Distance2[i][j];
                    SlS->ForceDD[j][i] = SlS->GMxM[j][i] / SlS->Distance2[j][i];

#endif
                }
            }
        }
    }
    for (i = 0; i < SlS->Elem; i++)
    {
        SlS->FX[i] = 0.0;
        SlS->FY[i] = 0.0;
        SlS->FZ[i] = 0.0;
    }
    // calculation of a XYZ forces

    for (i = 0; i < SlS->Elem; i++)
    {
        if (SlS->flInUse[i] ==0)
            continue;
        for (j = 0; j < SlS->Elem; j++)
        {
            if (i == j) continue;
            if (SlS->flInUse[j] ==0)
            continue;
            SlS->FX[i] += -( SlS->X[i] - SlS->X[j]) * SlS->ForceDD[i][j]/SlS->Distance[i][j];
            SlS->FY[i] += -( SlS->Y[i] - SlS->Y[j]) * SlS->ForceDD[i][j]/SlS->Distance[i][j];
            SlS->FZ[i] += -( SlS->Z[i] - SlS->Z[j]) * SlS->ForceDD[i][j]/SlS->Distance[i][j];
        }
    }
    // calculation of velocities and positions 

    for (i = 0; i < SlS->Elem; i++)
    {
        if (SlS->flInUse[i] ==0)
            continue;
#if 1
        // this is original formula
        SlS->VX[i] += SlS->FX[i] * TimeSl / SlS->M[i];
        SlS->VY[i] += SlS->FY[i] * TimeSl / SlS->M[i];
        SlS->VZ[i] += SlS->FZ[i] * TimeSl / SlS->M[i];

        SlS->X[i] += SlS->VX[i]*TimeSl;
        SlS->Y[i] += SlS->VY[i]*TimeSl;
        SlS->Z[i] += SlS->VZ[i]*TimeSl;

#else
        SlS->VX_[i] += SlS->FX[i];
        SlS->VY_[i] += SlS->FY[i];
        SlS->VZ_[i] += SlS->FZ[i];

#ifdef PROP_VELOCITY
        // this can be skipped to make calculation faster
        // VX is different from VX_ by coef = TimeS1*Mass
        // 
        SlS->VX[i] = SlS->VX_[i]*TimeSl / SlS->M[i];
        SlS->VY[i] = SlS->VY_[i]*TimeSl / SlS->M[i];
        SlS->VZ[i] = SlS->VZ_[i]*TimeSl / SlS->M[i];
#endif
        SlS->X_[i] += SlS->VX_[i];
        SlS->Y_[i] += SlS->VY_[i];
        SlS->Z_[i] += SlS->VZ_[i];

        SlS->X[i] = SlS->X_[i]*TimeSl*TimeSl / SlS->M[i];
        SlS->Y[i] = SlS->Y_[i]*TimeSl*TimeSl / SlS->M[i];
        SlS->Z[i] = SlS->Z_[i]*TimeSl*TimeSl / SlS->M[i];
#endif
    }
}
long double GreenwichAscension(long double EP);
void IteraSat(int TimeDirection, TRAOBJ * SlS, TRAOBJ * Sat, long double TimeOfCalc)
{
    int i;
    int j;
    double Summ;
    double Temp;
    double DX,DY,DZ;

    Sat->Lambda = GreenwichAscension(TimeOfCalc);
    // calculation of a forces
    // loop for all calulated satellites
    for (i = 0; i < Sat->Elem; i++)
    {
       // loop for all selestial bodies
        for (j = 0; j < SlS->Elem; j++)
        {
            // do we account that? Sun == yes Moon == yes but Saturn == probably no
            if (SlS->flInUse[j] ==0)
                continue;
            // that is the distance*distance from the satellite to selestial body
            double tD_Obj1Obj2 = (Sat->X[i] - SlS->X[j])*(Sat->X[i] - SlS->X[j]) + 
					    (Sat->Y[i] - SlS->Y[j])*(Sat->Y[i] - SlS->Y[j]) + 
					    (Sat->Z[i] - SlS->Z[j])*(Sat->Z[i] - SlS->Z[j]);
            
            // if distance will changed = then force will change => otherwice force will be like in prevous iteration
            // that is not correct - enable this line only to speedup calculations
#if 0
            if (tD_Obj1Obj2 != Sat->Distance2[i][j])
#endif
            {
                Sat->Distance2[i][j] = tD_Obj1Obj2;
                // not distance in m
                double tD_ = sqrt(tD_Obj1Obj2);
                // enable this line only to speedup calculations
#if 0
                if (tD_ != Sat->Distance[i][j])
#endif
                {
                    Sat->Distance[i][j] = tD_;
                    Sat->ForceDD[i][j] = SlS->GM[j] /* Sat->M[i]*/ / Sat->Distance2[i][j]; // to get real force need to multiply on mass of the satellite
                    if (j == Sat->LegBody)        // check ? j == EARTH or MOON ?? or anothe selestial body
                    {
                        Sat->R0divR[1] = R0_MODEL/tD_;
                        for (int n = 2; n <= (Sat->iLeg); n++)
                        {
                            Sat->R0divR[n] = Sat->R0divR[1]*Sat->R0divR[n-1];
                        }
                        Sat->ForceDD_ = Sat->ForceDD[i][j];
                    }
                }
            }
            if (j == Sat->LegBody)
            {
                if (Sat->Distance[i][j] < 10*R0_MODEL)
                {
                    // is this a metter how to measure Tetta?
                    // sin Tetta(colatitude) = Z/R
                    // and tetta(Latitude) btw equator and vector
                    // this means that sinTetta(colatitude) equal cosTetta(latitude)
                    Sat->SinTetta[1] = ((Sat->Z[i] - SlS->Z[j]) /Sat->Distance[i][j]);
                    
                    //if (Sat->LastSinTetta <0)
                    //    Sat->LastCosTetta = cos(asin(Sat->LastSinTetta)-M_PI/2.0);
                    //else
                    //Sat->CosTetta[1] = cos(asin(Sat->LastSinTetta[1]));

                    //Sat->LastCosTetta = sin(acos(Sat->LastSinTetta));
                    //Sat->LastCosTetta = sqrt(Sat->Distance[i][j]*Sat->Distance[i][j] - (Sat->Z[i] - SlS->Z[j])*(Sat->Z[i] - SlS->Z[j]))/Sat->Distance[i][j];
                    // first just use only J
                    
                    Sat->CalcP(Sat->SinTetta[1],((Sat->X[i] - SlS->X[j]) /Sat->Distance[i][j]),((Sat->Y[i] - SlS->Y[j]) /Sat->Distance[i][j]),
                        (Sat->X[i] - SlS->X[j]),(Sat->Y[i] - SlS->Y[j]),(Sat->Z[i] - SlS->Z[j]),Sat->Distance[i][j]);//>LastCosTetta);
                    //Sat->CalcPNK(Sat->LastCosTetta);
                    //Summ = Sat->SummJ();
                    Sat->SummXYZ(DX,DY,DZ);
                    
                    //Summ *=1.0001;//0.999905;

                    //Sat->ForceDD[i][j] = Sat->ForceDD_;// * Summ;
                    Sat->DeltaVX[i] =(1-DX);
                    Sat->DeltaVY[i] =(1-DY);
                    Sat->DeltaVZ[i] =(1-DZ);

                    // is this a WGS84??:
                    //Temp = SlS->GM[j] / (R0_MODEL*R0_MODEL);
                    //Summ = (1- 0.00193185138639*Sat->LastSinTetta*Sat->LastSinTetta)/sqrt(1.0 - 0.00669437999013*Sat->LastSinTetta*Sat->LastSinTetta);
                    //Summ = (1- 0.00193185138639*Sat->LastCosTetta*Sat->LastCosTetta)/sqrt(1.0 - 0.00669437999013*Sat->LastCosTetta*Sat->LastCosTetta);
                    //Sat->ForceDD[i][j] = Sat->ForceDD_*Summ;
                    //GM_MODEL 3.986004415E5;
                }
                else
                {
                    Sat->DeltaVX[i] =1;
                    Sat->DeltaVY[i] =1;
                    Sat->DeltaVZ[i] =1;
                }
            }
            //else // switch off sun and moon == 2 m/s difference in 1 hour
            //{
            //    if (Sat->LegBody != -1)
            //    {
            //        if ((j == MOON) || (j == SUN))
            //            Sat->ForceDD[i][j] = 0.0;
            //
            //    }
            //}
        }
    }
    for (i = 0; i < Sat->Elem; i++)
    {
        Sat->FX[i] = 0.0;
        Sat->FY[i] = 0.0;
        Sat->FZ[i] = 0.0;
    }
    // calculation of a XYZ forces

    for (i = 0; i < Sat->Elem; i++)
    {
        for (j = 0; j < SlS->Elem; j++)
        {
            if (SlS->flInUse[j] ==0)
            continue;
            if (j == Sat->LegBody)
            {
                Sat->FX[i] += -( Sat->X[i] - SlS->X[j]) *Sat->DeltaVX[i] * Sat->ForceDD[i][j]/Sat->Distance[i][j];
                Sat->FY[i] += -( Sat->Y[i] - SlS->Y[j]) *Sat->DeltaVY[i]* Sat->ForceDD[i][j]/Sat->Distance[i][j];
                Sat->FZ[i] += -( Sat->Z[i] - SlS->Z[j]) *Sat->DeltaVZ[i] * Sat->ForceDD[i][j]/Sat->Distance[i][j];
            }
            else
            {
                Sat->FX[i] += -( Sat->X[i] - SlS->X[j]) * Sat->ForceDD[i][j]/Sat->Distance[i][j];
                Sat->FY[i] += -( Sat->Y[i] - SlS->Y[j]) * Sat->ForceDD[i][j]/Sat->Distance[i][j];
                Sat->FZ[i] += -( Sat->Z[i] - SlS->Z[j]) * Sat->ForceDD[i][j]/Sat->Distance[i][j];
            }
        }
    }
    // calculation of velocities and positions 

    for (i = 0; i < Sat->Elem; i++)
    {
#if 0
        // this is original formula
        Sat->VX[i] += Sat->FX[i] * TimeSl /* Sat->M[i]*/;
        Sat->VY[i] += Sat->FY[i] * TimeSl /* Sat->M[i]*/;
        Sat->VZ[i] += Sat->FZ[i] * TimeSl /* Sat->M[i]*/;

        Sat->X[i] += Sat->VX[i]*TimeSl;
        Sat->Y[i] += Sat->VY[i]*TimeSl;
        Sat->Z[i] += Sat->VZ[i]*TimeSl;


#else
        Sat->VX_[i] += Sat->FX[i];
        Sat->VY_[i] += Sat->FY[i];
        Sat->VZ_[i] += Sat->FZ[i];
//#ifdef PROP_VELOCITY
        // this can be skipped to make calculation faster
        // VX is different from VX_ by coef = TimeS1*Mass
        // 
        Sat->VX[i] = Sat->VX_[i]*TimeSl /* Sat->M[i]*/;
        Sat->VY[i] = Sat->VY_[i]*TimeSl /* Sat->M[i]*/;
        Sat->VZ[i] = Sat->VZ_[i]*TimeSl /* Sat->M[i]*/;
//#endif
        Sat->X_[i] += Sat->VX_[i];
        Sat->Y_[i] += Sat->VY_[i];
        Sat->Z_[i] += Sat->VZ_[i];

        Sat->X[i] = Sat->X_[i]*TimeSl*TimeSl /* Sat->M[i]*/;
        Sat->Y[i] = Sat->Y_[i]*TimeSl*TimeSl /* Sat->M[i]*/;
        Sat->Z[i] = Sat->Z_[i]*TimeSl*TimeSl /* Sat->M[i]*/;
#endif
    }
}

void FireEngine(int TimeDirection, int DeltaTimeFromStart, TRAIMPLOBJ * Engine, TRAOBJ *Sat, int iSat, double DirX, double DirY, double DirZ)
{
    double VectorValue = sqrt(DirX*DirX + DirY*DirY +DirZ*DirZ);
    double coeff = 1.0/VectorValue;
    double dirX = DirX*coeff;
    double dirY = DirY*coeff;
    double dirZ = DirZ*coeff;
    // direction of a force is oposit the direction of a firing engine
    double ForceX = -dirX * Engine->ValImpl[DeltaTimeFromStart];
    double ForceY = -dirY * Engine->ValImpl[DeltaTimeFromStart];
    double ForceZ = -dirZ * Engine->ValImpl[DeltaTimeFromStart];
    double VX = ForceX* Engine->DeltaTime / Sat->M[iSat];
    double VY = ForceY* Engine->DeltaTime / Sat->M[iSat];
    double VZ = ForceZ* Engine->DeltaTime / Sat->M[iSat];
    
    Sat->VX[iSat] += VX;
    Sat->VY[iSat] += VY;
    Sat->VZ[iSat] += VZ;

    Sat->X[iSat] += VX * Engine->DeltaTime;
    Sat->Y[iSat] += VY * Engine->DeltaTime;
    Sat->Z[iSat] += VZ * Engine->DeltaTime;
    Sat->M[iSat] -= Engine->Weight * Engine->ValImpl[DeltaTimeFromStart]*Engine->DeltaTime / Engine->TotalImpulse;
}
double GetRadius(TRAOBJ * SlS, int iBody, double LongOnMoon,double LatiOnMoon)
{
    double dR = 1000000.0; // just for init == bogus radius 1000km
    // TBD earth NOT round
    if (iBody == EARTH)
    {
        dR = EarthR;
    }
    // TBD with the moon the same 
    if (iBody == MOON)
    {
        dR = MoonR;
    }
    return dR;
}
double GetRadius(TRAOBJ * SlS, int iBody, TRAOBJ * Sat, int iSat)
{
    double dR = 1000000.0; // just for init == bogus radius 1000km
    // TBD earth NOT round
    if (iBody == EARTH)
    {
        dR = EarthR;
    }
    // TBD with the moon the same 
    if (iBody == MOON)
    {
        dR = MoonR;
    }
    return dR;
}
double Targetlongitude = -15.0; // dolgota
double Targetlatitude = -2.0; // shirota
// view from north pole to equator plane 
/// x to right
//  y to up
//  z perpendicular xy and to north
// x to west (negative axe to east) 
// y to grinvich meridian ( 0 longitude) west is a negative angle, east positive
// z to north
//   LAT = latitude * pi/180    // shirota
//   LON = longitude * pi/180   // dolgota
//   Y =  R * cos(LAT) * cos(LON)
//   Z =  R * sin(LAT) 
//   X = -R * cos(LAT) * sin(LON)
//
// (1) sin(LAT) = z/R  => LAT = arcsin(z/R)
//
// (2) sin(LON) = - X / (R * cos(LAT))
//     LON = arcsin(- X / (R * cos(LAT)))
//     cos(LON) = Y / (R * cos(LAT))
//                                             AY
//      0-> PI/2     sin(LON)<0 cos(LON)>0     | sin(LON) >0 cos(LON) >0      0 -> -PI/2
//      -------------------------------------------------------------------------------> X
//      PI/2->PI     sin(LON)<0 cos(LON)<0     | sin(LON)>0 cos(LON) <0    -PI/2 -> -PI

void getLongLatiMoon(double &LongOnMoon,double &LatiOnMoon, TRAOBJ *SlS, int iRefMOON,int iRefEARTH, TRAOBJ *Sat, int iSat)
{
    double dREM = sqrt((SlS->X[iRefMOON] - SlS->X[iRefEARTH])*(SlS->X[iRefMOON] - SlS->X[iRefEARTH])+
                       (SlS->Y[iRefMOON] - SlS->Y[iRefEARTH])*(SlS->Y[iRefMOON] - SlS->Y[iRefEARTH])+
                       (SlS->Z[iRefMOON] - SlS->Z[iRefEARTH])*(SlS->Z[iRefMOON] - SlS->Z[iRefEARTH]));
    double dRM = sqrt((SlS->X[iRefMOON] - Sat->X[iSat])*(SlS->X[iRefMOON] - Sat->X[iSat])+
                       (SlS->Y[iRefMOON] - Sat->Y[iSat])*(SlS->Y[iRefMOON] - Sat->Y[iSat])+
                       (SlS->Z[iRefMOON] - Sat->Z[iSat])*(SlS->Z[iRefMOON] - Sat->Z[iSat]));
    double SatZSin = (Sat->Z[0] -SlS->Z[iRefMOON])/ dRM;
    double EarthZSin = (SlS->Z[iRefEARTH] -SlS->Z[iRefMOON])/ dREM;
    double SatLatiAsin = asin(SatZSin);
    double EarthLatiAsin = asin(EarthZSin);
    LatiOnMoon = (SatLatiAsin-EarthLatiAsin)*180.0/M_PI;
	double SatZCos = cos(SatLatiAsin);
	double EarthZCos = cos(EarthLatiAsin);

    double SatXSin = -(Sat->X[0] -SlS->X[iRefMOON])/ dRM/ SatZCos;
    double EarthXSin = -(SlS->X[iRefEARTH] -SlS->X[iRefMOON])/ dREM/ EarthZCos;

    double SatYCos = (Sat->Y[0] -SlS->Y[MOON])/ dRM/ SatZCos ;
    double EarthYCos = (SlS->Y[iRefEARTH] -SlS->Y[iRefMOON])/ dREM/ EarthZCos;


    double SatLongAsin = asin(SatXSin);
    // check all quarters
    if ((SatXSin >= 0.0) && (SatYCos>=0.0))
        ;
    else if ((SatXSin < 0.0) && (SatYCos>=0.0))
        ;
    else if ((SatXSin < 0.0) && (SatYCos < 0.0))
        SatLongAsin += M_PI/2.0;
    else
        SatLongAsin -= M_PI/2.0;

    double EarthLongAsin = asin(EarthXSin);

    if ((EarthXSin >= 0.0) && (EarthYCos>=0.0))
        ;
    else if ((EarthXSin < 0.0) && (EarthYCos>=0.0))
        ;
    else if ((EarthXSin < 0.0) && (EarthYCos < 0.0))
        EarthLongAsin += M_PI/2.0;
    else
        EarthLongAsin -= M_PI/2.0;

    LongOnMoon = (SatLongAsin-EarthLongAsin)*180.0/M_PI;
}
// PosXMoon,PosYMoon,PosZMoon - is a point on a moon surface ; coordinates are common
void getXYZMoon(double LongOnMoon,double LatiOnMoon,double &PosXMoon,double &PosYMoon,double &PosZMoon, TRAOBJ *SlS,int iRefMOON,int iRefEARTH, double dRadius)
{
    double dREM = sqrt((SlS->X[iRefMOON] - SlS->X[iRefEARTH])*(SlS->X[iRefMOON] - SlS->X[iRefEARTH])+
                       (SlS->Y[iRefMOON] - SlS->Y[iRefEARTH])*(SlS->Y[iRefMOON] - SlS->Y[iRefEARTH])+
                       (SlS->Z[iRefMOON] - SlS->Z[iRefEARTH])*(SlS->Z[iRefMOON] - SlS->Z[iRefEARTH]));
    double EarthZSin = (SlS->Z[iRefEARTH] -SlS->Z[iRefMOON])/ dREM;
    double EarthLatiAsin = asin(EarthZSin);
	double EarthZCos = cos(EarthLatiAsin);
    double EarthLatiDeg = EarthLatiAsin * 180.0/M_PI;
    double SatLatideg = LatiOnMoon + EarthLatiDeg;
    double SatLatiAsin = SatLatideg * M_PI/180.0;
    double SatZSin = sin(SatLatiAsin);
	double SatZCos = cos(SatLatiAsin);

    double EarthXSin = - (SlS->X[iRefEARTH] -SlS->X[iRefMOON])/ dREM/EarthZCos;
    double EarthYCos = (SlS->Y[iRefEARTH] -SlS->Y[iRefMOON])/ dREM/EarthZCos;
    double EarthLongAsin = asin(EarthXSin);
    if ((EarthXSin >= 0.0) && (EarthYCos>=0.0))
        ;
    else if ((EarthXSin < 0.0) && (EarthYCos>=0.0))
        ;
    else if ((EarthXSin < 0.0) && (EarthYCos < 0.0))
        EarthLongAsin += M_PI/2.0;
    else
        EarthLongAsin -= M_PI/2.0;

    double SatLongAsin = LongOnMoon + EarthLongAsin* 180.0/M_PI;
    double SatXSin = sin(SatLongAsin* M_PI/180.0);
    double SatYCos = cos(SatLongAsin* M_PI/180.0);

    
    PosZMoon =   dRadius * SatZSin + SlS->Z[iRefMOON];
    PosXMoon = - dRadius * SatXSin * SatZCos + SlS->X[iRefMOON];
    PosYMoon =   dRadius * SatYCos * SatZCos + SlS->Y[iRefMOON];
}
void GetXYZfromLatLong(double Long,double Lat,double &PosX,double &PosY,double &PosZ, double dRadius)
{
//   LAT = latitude * pi/180    // shirota
//   LON = longitude * pi/180   // dolgota
//   Y =  R * cos(LAT) * cos(LON)
//   Z =  R * sin(LAT) 
//   X = -R * cos(LAT) * sin(LON)

    PosZ =   dRadius * sin(Lat* M_PI/180.0);
    PosX = - dRadius * cos(Lat* M_PI/180.0) * sin(Long* M_PI/180.0);
    PosY =   dRadius * cos(Lat* M_PI/180.0) * cos(Long* M_PI/180.0);
}

#ifdef _DO_VISUALIZATION
#define IMAGE_W 1280
#define IMAGE_H 720
int RGBReferenceBody = EARTH;
unsigned char bRGBImage[IMAGE_W*IMAGE_H*3];
int bRGBImageW = IMAGE_W;
int bRGBImageH = IMAGE_H;
double dRGBScale = 1000000000.0;
int iProfile = 0; // 0 == XY , 1 == YZ, 2 == XZ 3 == -YZ 4 == -XZ 5==-XY
// 0 or XY is a view from North to south, 5 (- XY) is a view from south to north
// 1 or YZ is a view to easter

int iMaxSeq = 128;
int iCounterToSkip = 0;
int iMaxCounter = 24*60*60;

                                  // mer      ven       ear      mar      jup         sat       urn      nep       plt       moon        sun
unsigned char PlanetColors[16*3] ={50,50,50, 0,255,0, 0,0,255, 127,0,0, 127,127,0, 127,0,127, 0,127,127, 0,64,64, 64,0,64, 127,127,127, 255,0,0};
unsigned char SatColors[16*3] ={0,0,0, 20,20,20, 30,30,30,};
char szXYZ[6][5] ={"XY", "YZ","XZ", "_YZ", "_XZ", "_XY"};
void putpixel(unsigned char *bRGB, int X, int Y)
{
    int iRow = (X + bRGBImageW/2) + (bRGBImageH - (Y + bRGBImageH/2))*bRGBImageW;
    if ((iRow*3 >=0) && (iRow*3 < (sizeof(bRGBImage) -3)))
    {
        memcpy(&bRGBImage[iRow*3],bRGB, 3);
    }
}
void DrawFinalBody(TRAOBJ *SlS, int iBodyReference, TRAOBJ *Sat, int iSec,char *szInitName, TRAOBJ *XYZReference, int iXYZReference, double Scale, int &StartSequence)
{
    double X1 = SlS->X[iBodyReference] - XYZReference->X[iXYZReference];
    double Y1 = SlS->Y[iBodyReference] - XYZReference->Y[iXYZReference];
    double Z1 = SlS->Z[iBodyReference] - XYZReference->Z[iXYZReference];
    double X2;
    double Y2;
    int i;
    double dRadius = GetRadius(SlS, iBodyReference, Sat, 0);
    switch(iProfile)
    {
    case 0: X2 = X1; Y2 = Y1;break; 
    case 1: X2 = Y1; Y2 = Z1;break; 
    case 2: X2 = X1; Y2 = Z1;break; 
    case 3: X2 = -Y1; Y2 = Z1;break;
    case 4: X2 = -X1; Y2 = Z1;break;
    case 5: X2 = -X1; Y2 = Y1;break; 
    default: break;
    }
    // draw circle
    for (i = 0; i < 365; i++)
    {
        double dDeg = ((double)i) * M_PI /180;
        double X2C = X2 + cos(dDeg)*dRadius;
        double Y2C = Y2 + sin(dDeg)*dRadius;

        if (-Scale < X2C  && X2C < Scale)
        {
            if (((-Scale * bRGBImageH/ bRGBImageW) < Y2C)  && Y2C < ((Scale* bRGBImageH/ bRGBImageW)))
            {
                putpixel(&PlanetColors[iBodyReference*3], (int)(X2C/Scale*bRGBImageH), (int)(Y2C/Scale*bRGBImageH));
            }
        }
    }
    // draw sun day & night lights
    double X1S = SlS->X[SUN] - XYZReference->X[iXYZReference];
    double Y1S = SlS->Y[SUN] - XYZReference->Y[iXYZReference];
    double Z1S = SlS->Z[SUN] - XYZReference->Z[iXYZReference];
    double DS = sqrt(X1S*X1S + Y1S*Y1S + Z1S*Z1S);
    switch(iProfile)
    {
    case 0: X1S = X1S; Y1S = Y1S;break; 
    case 1: X1S = Y1S; Y1S = Z1S;break; 
    case 2: X1S = X1S; Y1S = Z1S;break; 
    case 3: X1S = -Y1S; Y1S = Z1S;break;
    case 4: X1S = -X1S; Y1S = Z1S;break;
    case 5: X1S = -X1S; Y1S = Y1S;break; 
    default: break;
    }
    double dAngleSin = X1S/DS;
    double dAngleCos = Y1S/DS;
    for (i = -(int)dRadius; i < (int)dRadius; i+=1000)
    {
        double X2C = X2 + -dAngleCos*((double)i);
        double Y2C = Y2 + dAngleSin*((double)i);
        if (-Scale < X2C  && X2C < Scale)
        {
            if (((-Scale * bRGBImageH/ bRGBImageW) < Y2C)  && Y2C < ((Scale* bRGBImageH/ bRGBImageW)))
            {
                putpixel(&PlanetColors[iXYZReference*3], (int)(X2C/Scale*bRGBImageH), (int)(Y2C/Scale*bRGBImageH));
            }
        }

    }

    // draw direction to earth
    X1S = SlS->X[EARTH] - XYZReference->X[iXYZReference];
    Y1S = SlS->Y[EARTH] - XYZReference->Y[iXYZReference];
    Z1S = SlS->Z[EARTH] - XYZReference->Z[iXYZReference];
    DS = sqrt(X1S*X1S + Y1S*Y1S + Z1S*Z1S);
    switch(iProfile)
    {
    case 0: X1S = X1S; Y1S = Y1S;break; 
    case 1: X1S = Y1S; Y1S = Z1S;break; 
    case 2: X1S = X1S; Y1S = Z1S;break; 
    case 3: X1S = -Y1S; Y1S = Z1S;break;
    case 4: X1S = -X1S; Y1S = Z1S;break;
    case 5: X1S = -X1S; Y1S = Y1S;break; 
    default: break;
    }
    dAngleSin = X1S/DS;
    dAngleCos = Y1S/DS;
    for (i = 0; i < dRadius; i+=1000)
    {
        double X2C = X2 + dAngleSin*((double)i);
        double Y2C = Y2 + dAngleCos*((double)i);
        if (-Scale < X2C  && X2C < Scale)
        {
            if (((-Scale * bRGBImageH/ bRGBImageW) < Y2C)  && Y2C < ((Scale* bRGBImageH/ bRGBImageW)))
            {
                putpixel(&PlanetColors[iXYZReference*3], (int)(X2C/Scale*bRGBImageH), (int)(Y2C/Scale*bRGBImageH));
            }
        }

    }
    // draw last satelite point
    X1 = Sat->X[0] - XYZReference->X[iXYZReference];
    Y1 = Sat->Y[0] - XYZReference->Y[iXYZReference];
    Z1 = Sat->Z[0] - XYZReference->Z[iXYZReference];
    switch(iProfile)
    {
    case 0: X2 = X1; Y2 = Y1;break; 
    case 1: X2 = Y1; Y2 = Z1;break; 
    case 2: X2 = X1; Y2 = Z1;break; 
    case 3: X2 = -Y1; Y2 = Z1;break;
    case 4: X2 = -X1; Y2 = Z1;break;
    case 5: X2 = -X1; Y2 = Y1;break; 
    default: break;
    }

    if (-Scale < X2  && X2 < Scale)
    {
        if (((-Scale * bRGBImageH/ bRGBImageW) < Y2)  && Y2 < ((Scale* bRGBImageH/ bRGBImageW)))
        {
            putpixel(&SatColors[0], (int)(X2/Scale*bRGBImageH), (int)(Y2/Scale*bRGBImageH));
        }
    }
    if (StartSequence++ < iMaxSeq)
    {
        char szName[256];
        sprintf(szName, "%s%s%05d.jpg", szInitName, &szXYZ[iProfile],StartSequence-1);
        write_JPEG_file (szName, 80, bRGBImageW, bRGBImageH, 3,  bRGBImage, JCS_RGB);
    }
    // draw North pole on a moon
    if (iBodyReference == MOON)
    {

    }
}
void DrawAnimationSequence(TRAOBJ *SlS, TRAOBJ *Sat, int iSec,char *szInitName, TRAOBJ *XYZReference, int iXYZReference, double Scale, int &StartSequence, int storeNow)
{
    int i,j;

    if ((StartSequence == 0) && (storeNow ==0))
    {
        memset(bRGBImage, 255, sizeof(bRGBImage));
        iCounterToSkip = iMaxCounter;
        StartSequence++;
    }
    for (i = 0; i < SlS->Elem; i++)
    {
        double X1 = SlS->X[i] - XYZReference->X[iXYZReference];
        double Y1 = SlS->Y[i] - XYZReference->Y[iXYZReference];
        double Z1 = SlS->Z[i] - XYZReference->Z[iXYZReference];
        double X2;
        double Y2;
        switch(iProfile)
        {
        case 0: X2 = X1; Y2 = Y1;break; 
        case 1: X2 = Y1; Y2 = Z1;break; 
        case 2: X2 = X1; Y2 = Z1;break; 
        case 3: X2 = -Y1; Y2 = Z1;break;
        case 4: X2 = -X1; Y2 = Z1;break;
        case 5: X2 = -X1; Y2 = Y1;break; 
        default: break;
        }
        if (-Scale < X2  && X2 < Scale)
        {
            if (((-Scale * bRGBImageH/ bRGBImageW) < Y2)  && Y2 < ((Scale* bRGBImageH/ bRGBImageW)))
            {
                putpixel(&PlanetColors[i*3], (int)(X2/Scale*bRGBImageH), (int)(Y2/Scale*bRGBImageH));
            }
        }
    }
    for (i = 0; i < Sat->Elem; i++)
    {
        double X1 = Sat->X[i] - XYZReference->X[iXYZReference];
        double Y1 = Sat->Y[i] - XYZReference->Y[iXYZReference];
        double Z1 = Sat->Z[i] - XYZReference->Z[iXYZReference];
        double X2;
        double Y2;
        switch(iProfile)
        {
        case 0: X2 = X1; Y2 = Y1;break; 
        case 1: X2 = Y1; Y2 = Z1;break; 
        case 2: X2 = X1; Y2 = Z1;break; 
        case 3: X2 = -Y1; Y2 = Z1;break;
        case 4: X2 = -X1; Y2 = Z1;break;
        case 5: X2 = -X1; Y2 = Y1;break;
        default: break;
        }
        if (-Scale < Y2  && Y2 < Scale)
        {
            if ((-Scale * bRGBImageW/ bRGBImageH) < X2  && X2 < (Scale* bRGBImageW/ bRGBImageH))
            {
                putpixel(&SatColors[i*3], (int)(X2/Scale*bRGBImageH), (int)(Y2/Scale*bRGBImageH));
            }
        }
    }
    // store last image
    if ((--iCounterToSkip <= 0) || (storeNow == 1))
    {
        if (storeNow == 1)
        {
            iCounterToSkip++;
        }
        else
            iCounterToSkip = iMaxCounter;
        if (StartSequence++ < iMaxSeq)
        {
            char szName[256];
            sprintf(szName, "%s%s%05d.jpg", szInitName, &szXYZ[iProfile],StartSequence-1);
            //DeleteFile(szName);
            write_JPEG_file (szName, 80, bRGBImageW, bRGBImageH, 3,  bRGBImage, JCS_RGB);

        }
    }

}
#endif

int RunOrVoidEngine(int TimeDirection, TRAIMPLOBJ * Engines, TRAOBJ * SlS, TRAOBJ *Sat, 
                     long &CurentTimeSec, int &CurentTimeTD, int &CurentIteraPerSec, double StartTimeSec)
{
    int iRet = 0;
    int i,j;
    double SecCalc;
    double DTempo;
    for(i = 0; i < Sat->Elem; i++)
    {
        if (Sat->flInUse)
        {
            for (j = 0; j < EnginesCount; j++)
            {
                if (Engines[j].iEngineOnSatellite == i)
                {
                    if (Engines[j].EngineDone == 0)
                    {
                        if (Engines[j].EngineOn)
                        {
                            // engine is ON
                            iRet = 1;
                            if ((Engines[j].ImplsPointer) < Engines[j].iLine)
                            {
                                FireEngine(1,Engines[j].ImplsPointer, &Engines[j], Sat, i, Engines[j].XVec, Engines[j].YVec, Engines[j].ZVec);
                                double VX = Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[EARTH];
                                double VY = Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[EARTH];
                                double VZ = Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[EARTH];
                                DTempo = sqrt(VX*VX+VY*VY+VZ*VX);
                                
                                VX = Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[MOON];
                                VY = Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[MOON];
                                VZ = Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[MOON];
                                DTempo = sqrt(VX*VX+VY*VY+VZ*VZ);
                                if (j ==4)
                                    printf("\n velocity = %f ", DTempo);
                                double X = Sat->X[Engines[j].iEngineOnSatellite] - SlS->X[MOON];
                                double Y = Sat->Y[Engines[j].iEngineOnSatellite] - SlS->Y[MOON];
                                double Z = Sat->Z[Engines[j].iEngineOnSatellite] - SlS->Z[MOON];
                                DTempo = sqrt(X*X+Y*Y+Z*Z);
                                if (j ==4)
                                    printf(" distance = %f ", DTempo- GetRadius(SlS, MOON, Sat, Engines[j].iEngineOnSatellite));

                            }
                            if ((++(Engines[j].ImplsPointer)) >= Engines[j].iLine)
                            {
                                if (Engines[j].ImplsPointer % Engines[j].IteraPerSec ==0)
                                {
                                    printf("\n cutoff engine=====> %d at= %f",j,Engines[j].FireTime+Engines[j].DeltaTime*Engines[j].iLine);
                                    // engine done
                                    iRet = 0;
                                    Engines[j].EngineDone = 1;
                                    Engines[j].iCountApogPerig = 0;
                                    TimeSl = TimeSlOld;
                                    CurentTimeSec = OldCurentTimeSec + Engines[j].ImplsPointer / Engines[j].IteraPerSec;
                                    CurentTimeTD = OldCurentTimeTD;
                                    CurentIteraPerSec = OldCurentIteraPerSec;
                                }
                            }
                        }
                        else
                        {
                            // skip firing engines which we are not interesting
                            if (j > LastEngine)
                                continue;
                            // engine is OFF
                            // SecCalc = CurentTimeSec * CurentIteraPerSec + CurentTimeTD;
                            SecCalc = ((double)CurentTimeSec) + ((double)CurentTimeTD) / ((double)CurentIteraPerSec);
                            if ((Engines[j].FireTime) <= SecCalc)
                            {
								unsigned char EngineColor[3]= {255,0,0};
								double X1 = Sat->X[j] - SolarSystem.X[EARTH];
								double Y1 = Sat->Y[j] - SolarSystem.Y[EARTH];
								double Z1 = Sat->Z[j] - SolarSystem.Z[EARTH];
								double X2;
								double Y2;
#ifdef _DO_VISUALIZATION
								switch(iProfile)
								{
									case 0: X2 = X1; Y2 = Y1;break; 
									case 1: X2 = Y1; Y2 = Z1;break; 
									case 2: X2 = X1; Y2 = Z1;break; 
									case 3: X2 = -Y1; Y2 = Z1;break;
									case 4: X2 = -X1; Y2 = Z1;break;
									case 5: X2 = -X1; Y2 = Y1;break; 
									default: break;
								}
								if (-dRGBScale < X2  && X2 < dRGBScale)
								{
									if (((-dRGBScale * bRGBImageH/ bRGBImageW) < Y2)  && Y2 < ((dRGBScale* bRGBImageH/ bRGBImageW)))
									{
										putpixel(&EngineColor[i*3], (int)(X2/dRGBScale*bRGBImageH), (int)(Y2/dRGBScale*bRGBImageH));
									}
								}
#endif
								printf("\n fire engine=====> %d at= %f",j,Engines[j].FireTime);
                                iRet = 1;
                                Engines[j].EngineOn = 1;
                                Engines[j].ImplsPointer = 0;
                                TimeSlOld = TimeSl;
                                TimeSl = Engines[j].DeltaTime;
                                OldCurentTimeSec = CurentTimeSec;
                                OldCurentTimeTD = CurentTimeTD;
                                OldCurentIteraPerSec = CurentIteraPerSec;
                                CurentTimeTD = 0;
                                CurentIteraPerSec = Engines[j].IteraPerSec;

                                // for now engine fire just in direction of a velocity around nearbody
                                double VX = Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[Engines[j].AngleOnBody];
                                double VY = Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[Engines[j].AngleOnBody];
                                double VZ = Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[Engines[j].AngleOnBody];
                                DTempo = sqrt(VX*VX+VY*VY+VZ*VX);
								// setting both vectors to 1.0 will fire engine in direction
								// ortogonal vector (not in direction of velocity)
								//if ((Engines[j].Ang1 == 1.0) && (Engines[j].Ang2 == 1.0))
                                if (Engines[j].AngleType == 4) // 4 - same direction as vector of velocity to nearbody // this is a brake impulse to land on the moon
                                {
                                    Engines[j].XVec = (Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[Engines[j].AngleOnBody]);
                                    Engines[j].YVec = (Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[Engines[j].AngleOnBody]);
                                    Engines[j].ZVec = (Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[Engines[j].AngleOnBody]);

                                }
                                if (Engines[j].AngleType == 3) // 3 - oposit vector of velocity
                                {
                                    Engines[j].XVec = -(Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[Engines[j].AngleOnBody]);
                                    Engines[j].YVec = -(Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[Engines[j].AngleOnBody]);
                                    Engines[j].ZVec = -(Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[Engines[j].AngleOnBody]);
                                }
                                if (Engines[j].AngleType == 0) // 0 - tangent line to orbit (elipse) oposit velocity  to nearbody
                                {
                                    //Engines[j].XVec = -VX;
                                    //Engines[j].YVec = -VY;
                                    //Engines[j].ZVec = -VZ;
                                    // equation of a plane perpendicular to a vector from center of earth to setelite
                                    // first calculates radius vector (R=(Xs,Ys,Zs)) of a satellite's point reletive to earth (P=(Xs,Yx,Zs))
                                    double Xs = Sat->X[Engines[j].iEngineOnSatellite] - SlS->X[Engines[j].AngleOnBody];
                                    double Ys = Sat->Y[Engines[j].iEngineOnSatellite] - SlS->Y[Engines[j].AngleOnBody];
                                    double Zs = Sat->Z[Engines[j].iEngineOnSatellite] - SlS->Z[Engines[j].AngleOnBody];
                                    // N = (A,B,C) is a normal to a plane (the needs to be in same direction of a radius vector)
                                    DTempo = sqrt(Xs*Xs+Ys*Ys+Zs*Zs);
                                    double Xn = Xs/DTempo;
                                    double Yn = Ys/DTempo;
                                    double Zn = Zs/DTempo;
                                    // PLane : A*(x-Xs) + B*(y-Ys) + C*(z-Zs) = 0
                                    // or Xn*(x-Xs) + Yn*(y-Ys) + Zn*(z-Zs) = 0
                                    // vector (Q=(Xq,Yq,Zq)) perpendicular to R and starting point is P, 4 equation:
                                    //(1)  (x - x1)(y2-y1) = (x2-x1)(y-y1) 
                                    //(2)  (y - y1)(z2-z1) = (y2-y1)(z-z1)
                                    //(3)  (z - z1)(x2-x1) = (z2-z1)(x-x1)
                                    //(4)  Xn*(x-Xs) + Yn*(y-Ys) + Zn*(z-Zs) = 0
                                    // for calculations:
                                    // first calculate 2 vector M1 and M2 ortogonal X and VX
                                    // (1) Xs*XM + Ys*YM + Zs*ZM = 0 => ZM = - (Xs*XM + Ys*YM)/ Zs
                                    // 
                                    // (2) VX*XM + VY*YM + VZ*ZM = 0 => VY*YM = - (VX*XM +VZ*ZM) => VY*YM = - (VX*XM -VZ*(Xs*XM + Ys*YM)/ Zs)
                                    // VY*YM = - ((VX*XM*Zs/Zs -VZ*(Xs*XM + Ys*YM)/ Zs)) => VY*YM = - (VX*XM*Zs -VZ*(Xs*XM + Ys*YM))/ Zs => VY*YM = - (VX*XM*Zs - VZ*Xs*XM - VZ*Ys*YM)/ Zs =>
                                    // VY*YM = - (VX*XM*Zs - VZ*Xs*XM)/ Zs + VZ*Ys*YM/ Zs => VY*YM - VZ*Ys*YM/ Zs = - (VX*XM*Zs - VZ*Xs*XM)/ Zs => YM*(VY*Zs - VZ*Ys)/ Zs = - (VX*XM*Zs - VZ*Xs*XM)/ Zs
                                    // YM*(VY*Zs - VZ*Ys) = - (VX*XM*Zs - VZ*Xs*XM) = > YM = - (VX*XM*Zs - VZ*Xs*XM)/(VY*Zs - VZ*Ys) => YM = - XM * (VX*Zs - VZ*Xs)/(VY*Zs - VZ*Ys)
                                    double XM =  1.0;
                                    double YM = - XM * (VX*Zs - VZ*Xs)/(VY*Zs - VZ*Ys);
                                    double ZM = - (Xs*XM + Ys*YM)/ Zs;
                                    // second calcualtes 2 vectors D1 D2 ortogonal to X and M1
                                    // (3) Xs*Xd + Ys*Yd + Zs*Zd = 0 => Zd = - (Xs*Xd + Ys*Yd)/Zs 
                                    //
                                    // (4) XM*Xd + YM*Yd + ZM*Zd = 0 => YM*Yd = -(XM*Xd + ZM*Zd) => YM*Yd = -(XM*Xd - ZM*(Xs*Xd + Ys*Yd)/Zs) => YM*Yd = -(XM*Xd*Zs - ZM*Xs*Xd - ZM*Ys*Yd)/Zs =>
                                    // YM*Yd = -(XM*Xd*Zs - ZM*Xs*Xd)/Zs + ZM*Ys*Yd/Zs => YM*Yd - ZM*Ys*Yd/Zs = -(XM*Xd*Zs - ZM*Xs*Xd)/Zs => Yd*(YM*Zs/Zs - ZM*Ys/Zs) = -(XM*Xd*Zs - ZM*Xs*Xd)/Zs
                                    // Yd*(YM*Zs - ZM*Ys) = -(XM*Xd*Zs - ZM*Xs*Xd) => Yd = -(XM*Xd*Zs - ZM*Xs*Xd)/(YM*Zs - ZM*Ys) => Yd = -Xd*(XM*Zs - ZM*Xs)/(YM*Zs - ZM*Ys)
                                    double Xd1 =  1.0;
                                    double Yd1 = -Xd1*(XM*Zs - ZM*Xs)/(YM*Zs - ZM*Ys);
                                    double Zd1 = - (Xs*Xd1 + Ys*Yd1)/Zs ;
                                    double Xd2 =  -1.0;
                                    double Yd2 = -Xd2*(XM*Zs - ZM*Xs)/(YM*Zs - ZM*Ys);
                                    double Zd2 = - (Xs*Xd2 + Ys*Yd2)/Zs ;
                                    // conformation: D and X are ortogonal and M and D is ortogonal
                                    DTempo = sqrt(Xs*Xd1+Ys*Yd1+Zs*Zd1);
                                    DTempo = sqrt(Xs*Xd2+Ys*Yd2+Zs*Zd2);
                                    DTempo = sqrt(Xs*XM+Ys*YM+Zs*ZM);
                                    // calculation - which of a vectors D1 or D2 close to vector of velocity
                                    // cos(alfa) = ( (x2 -x1)(x4-x3) + (y2-y1)(y4-y3) + (z2-z1)(z4-z3) )/ ( sqrt( (x2-x1)^2 + (y2-y1)^2 +(z2-z1)^2 )*sqrt( (x4-x3)^2 + (y4-y3)^2 +(z4-z3)^2 )
                                    double cosAlfa1 = (Xd1*VX + Yd1*VY + Zd1*VZ) / sqrt(Xd1*Xd1 + Yd1*Yd1+ Zd1*Zd1) / sqrt(VX*VX + VY*VY +VZ*VZ);
                                    double cosAlfa2 = (Xd2*VX + Yd2*VY + Zd2*VZ) / sqrt(Xd2*Xd2 + Yd2*Yd2+ Zd2*Zd2) / sqrt(VX*VX + VY*VY +VZ*VZ);

                                    if (cosAlfa1 < cosAlfa2)
                                    {
                                        Engines[j].XVec = -Xd2;
                                        Engines[j].YVec = -Yd2;
                                        Engines[j].ZVec = -Zd2;
                                    }
                                    else
                                    {
                                        Engines[j].XVec = -Xd1;
                                        Engines[j].YVec = -Yd1;
                                        Engines[j].ZVec = -Zd1;
                                    }
                                    if (j==4) // for a brake impulse into a vector of velociti around moon
                                    {
                                        //Engines[j].XVec = -(Sat->X[Engines[j].iEngineOnSatellite] - SlS->X[MOON]);
                                        //Engines[j].YVec = -(Sat->Y[Engines[j].iEngineOnSatellite] - SlS->Y[MOON]);
                                        //Engines[j].ZVec = -(Sat->Z[Engines[j].iEngineOnSatellite] - SlS->Z[MOON]);

                                        Engines[j].XVec = (Sat->VX[Engines[j].iEngineOnSatellite] - SlS->VX[MOON]);
                                        Engines[j].YVec = (Sat->VY[Engines[j].iEngineOnSatellite] - SlS->VY[MOON]);
                                        Engines[j].ZVec = (Sat->VZ[Engines[j].iEngineOnSatellite] - SlS->VZ[MOON]);
                                    }
                                }
								//else if ((Engines[j].Ang1 != 0.0) && (Engines[j].Ang2 != 0.0))
                                else if (Engines[j].AngleType == 1) // 1 - two angles set with reference to NearBody centre direction
								{
								}
                                else if (Engines[j].AngleType == 2) // 2 - three angles set vector fire (constant all fire time) 
								{
									//Engines[j].XVec = 8;
                                    //Engines[j].YVec = 0;
                                    //Engines[j].ZVec = -16;
                                    Engines[j].XVec = Engines[j].XVec;
                                    Engines[j].YVec = Engines[j].YVec;
                                    Engines[j].ZVec = Engines[j].ZVec;
								}

                            }
                        }
                    }
                }
            }
        }
    }
    return (iRet);
}


#define CALC_SOLAR_SYSTEM 1
void AssignFromNASAData(TRAOBJ * SlS, double JDSec)
{
    double aproxim;
    //double RBarisMoon;
    //double RBarisEarth;
    //double VEarthOrbit;
    //double VMoonOrbit;
    //double VDelta;
    double Temp0;
    double Temp1;
    double Temp2;
    double Temp3;
    double Temp4;

    stateType  State;
    double dEMRAT = Find_DataInHeader("EMRAT ");
    double dAU = Find_DataInHeader("AU    ")*1000.0;
	// store
	AU = dAU;
    SlS->Elem = 11;
    aproxim = 0.58782958984375;
    aproxim = 1.0;
    aproxim = 1.0;//1.0496393116944882060883002708885;
#define PlanetsCount SlS->Elem
//#define PlanetsCount PLANET_COUNT
    for (int i = 0; i < PlanetsCount; i++)
    {
        SlS->flInUse[i] = 1;        
        for (int j = 0; j < PlanetsCount; j++)
        {
            SlS->Distance2[i][j] = 0.0;
        }
        Interpolate_State( JDSec , i , &State );
        SlS->X[i] = State.Position[0]*1000.0;
        SlS->Y[i] = State.Position[1]*1000.0;
        SlS->Z[i] = State.Position[2]*1000.0;
        SlS->VX[i] = State.Velocity[0]*1000.0;
        SlS->VY[i] = State.Velocity[1]*1000.0;
        SlS->VZ[i] = State.Velocity[2]*1000.0;
        switch(i)
        {
        case 0: SlS->GM[i] = Find_DataInHeader("GM1   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 1: SlS->GM[i] = Find_DataInHeader("GM2   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 2: SlS->GM[i] = Find_DataInHeader("GMB   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0)*(dEMRAT)/(dEMRAT+1.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 3: SlS->GM[i] = Find_DataInHeader("GM4   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 4: SlS->GM[i] = Find_DataInHeader("GM5   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 5: SlS->GM[i] = Find_DataInHeader("GM6   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 6: SlS->GM[i] = Find_DataInHeader("GM7   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 7: SlS->GM[i] = Find_DataInHeader("GM8   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 8: SlS->GM[i] = Find_DataInHeader("GM9   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);SlS->M[i] = SlS->GM[i]/Gbig;break;
        case 9: SlS->GM[i] = Find_DataInHeader("GMB   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0)/(dEMRAT+1.0);SlS->M[i] = SlS->GM[i]/Gbig;
            {
                double BSX = SlS->X[EARTH];
                double BSY = SlS->Y[EARTH];
                double BSZ = SlS->Z[EARTH];
                double BSVX = SlS->VX[EARTH];
                double BSVY = SlS->VY[EARTH];
                double BSVZ = SlS->VZ[EARTH];
#if 0
                SlS->X[i] += SlS->X[EARTH];
                SlS->Y[i] += SlS->Y[EARTH];
                SlS->Z[i] += SlS->Z[EARTH];
                SlS->VX[i] += SlS->VX[EARTH];
                SlS->VY[i] += SlS->VY[EARTH];
                SlS->VZ[i] += SlS->VZ[EARTH];
#else
                SlS->X[EARTH] = BSX - (SlS->X[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                SlS->Y[EARTH] = BSY - (SlS->Y[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                SlS->Z[EARTH] = BSZ - (SlS->Z[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                Temp4 = sqrt(
                                (SlS->X[9])*(SlS->X[9])+
                                (SlS->Y[9])*(SlS->Y[9])+
                                (SlS->Z[9])*(SlS->Z[9])
                                );
                printf("\n D Earth - Moon %f ", Temp4);
                Temp1 = sqrt(
                                ((SlS->X[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9])))*
                                ((SlS->X[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9]))) +
                                ((SlS->Y[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9])))*
                                ((SlS->Y[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9]))) +
                                ((SlS->Z[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9])))*
                                ((SlS->Z[9]*SlS->M[2]/(SlS->M[2]+SlS->M[9])))
                                );

                 printf("\n D BC  - Moon %f ", Temp1);
                 Temp1 = sqrt(
                                SlS->GM[2] / Temp1
                                );

                 Temp2 = sqrt(
                                (SlS->VX[9])*(SlS->VX[9])+
                                (SlS->VY[9])*(SlS->VY[9])+
                                (SlS->VZ[9])*(SlS->VZ[9])
                                );
                 Temp3 = Find_DataInHeader("GMB   ")*dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0);
                 SlS->X[MOON] = SlS->X[EARTH] + SlS->X[MOON];
                 SlS->Y[MOON] = SlS->Y[EARTH] + SlS->Y[MOON];
                 SlS->Z[MOON] = SlS->Z[EARTH] + SlS->Z[MOON];

                 SlS->VX[EARTH] = BSVX - /*aproxim**/(SlS->VX[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                 SlS->VY[EARTH] = BSVY - /*aproxim**/(SlS->VY[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                 SlS->VZ[EARTH] = BSVZ - /*aproxim**/(SlS->VZ[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]));
                 Temp0 = (/*aproxim**/(SlS->VX[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON])))*
                                (/*aproxim**/(SlS->VX[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]))) +
                                (/*aproxim**/(SlS->VY[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON])))*
                                (/*aproxim**/(SlS->VY[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON]))) +
                                (/*aproxim**/(SlS->VZ[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON])))*
                                (/*aproxim**/(SlS->VZ[MOON]*SlS->M[MOON]/(SlS->M[EARTH]+SlS->M[MOON])));
                 SlS->VX[MOON] = SlS->VX[EARTH] + SlS->VX[MOON];
                 SlS->VY[MOON] = SlS->VY[EARTH] + SlS->VY[MOON];
                 SlS->VZ[MOON] = SlS->VZ[EARTH] + SlS->VZ[MOON];
#endif
            }
            break;
        case 10: SlS->GM[i] = Find_DataInHeader("GMS   ") *dAU*dAU*dAU/(24.0*60.0*60.0*24.0*60.0*60.0); SlS->M[i] = SlS->GM[i]/Gbig;break;
        }
    }

    printf("\n Earth velocity =%f", sqrt(Temp0));
    printf("\n Moon sqrt =%f", Temp1);
    printf("\n Moon velocity N =%f", Temp2);

    // this has to be set to properly calc helper variables
    TimeSl = 1.0 / IterPerSec;
            
    // this is done to reduce errors and avid unnessary 5 mul/div operations
    // temporary X_, VX_ will just added (in paralel can be done actualy) 
    for (int i = 0; i <SlS->Elem; i++)
    {
        SlS->VX_[i] = SlS->VX[i]* SlS->M[i] /TimeSl ;
        SlS->VY_[i] = SlS->VY[i]* SlS->M[i] /TimeSl;
        SlS->VZ_[i] = SlS->VZ[i]* SlS->M[i] /TimeSl;

        SlS->X_[i] = SlS->X[i]* SlS->M[i] /TimeSl/TimeSl ;
        SlS->Y_[i] = SlS->Y[i]* SlS->M[i] /TimeSl/TimeSl ;
        SlS->Z_[i] = SlS->Z[i]* SlS->M[i] /TimeSl/TimeSl ;
        for (int j = 0; j < SlS->Elem; j++)
        {
            SlS->GMxM[i][j] = SlS->GM[i]*SlS->M[j];
        }
    }
    // skip that planets to avoid 
    //SlS->flInUse[0] = 0; //Mer
    //SlS->flInUse[1] = 0; //Ven
    //SlS->flInUse[3] = 0; //Mar
    //SlS->flInUse[6] = 0; //
    //SlS->flInUse[7] = 0; //
    //SlS->flInUse[8] = 0; //
    //SlS->flInUse[SUN] = 0;
    //SlS->flInUse[MOON] = 0;
    

}

// main engine to dump keplers elements as it is (no drag)
void DumpKeplers(long double &T,long double &Ecc, long double &Incl, long double &AssNode, long double &ArgPer, long double &MeanAnm, long double Mass, long double mass, 
                    long double &X, long double &Y, long double &Z, long double &VX, long double &VY, long double &VZ)
{
	// see http://www.projectpluto.com/source.htm or http://www.projectpluto.com/lunar.zip
	// code based on classel.cpp
	//const double *v = r + 3;
	//const double r_dot_v = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
	long double r_dot_v = X*VX + Y*VY + Z*VZ;
	//const double dist = sqrt( r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
	long double R = sqrt(X*X + Y*Y + Z*Z);
	//const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	long double Vsq = VX*VX + VY*VY + VZ*VZ;
	//const double inv_major_axis = 2. / dist - v2 / gm;
	long double GbigMass = Gbig * Mass;
	long double InvMAxis = 2.0 / R - Vsq / GbigMass;
	long double h0, n0;
	long double h[3], e[3], ecc2;
	long double ecc, perihelion_speed, gm_over_h0;
	int i;

	//h[0] = r[1] * v[2] - r[2] * v[1];
	h[0]   = Y * VZ   - Z * VY;
	//h[1] = r[2] * v[0] - r[0] * v[2];
	h[1] = Z * VX - X * VZ;
	// h[2] = r[0] * v[1] - r[1] * v[0];
	h[2] = X * VY - Y * VX;
	n0 = h[0] * h[0] + h[1] * h[1];
	h0 = n0 + h[2] * h[2];
	n0 = sqrt( n0);
	h0 = sqrt( h0);

	//elem->asc_node = atan2( h[0], -h[1]);
	AssNode= atan2( h[0], -h[1]);
	//if (AssNode < 0.0)
	//	AssNode = 2*M_PI + AssNode;
    // elem->incl = asine( n0 / h0);
	Incl = asin( n0 / h0);
    if( h[2] < 0.)                   /* retrograde orbit */
         Incl = M_PI - Incl;
	//e[0] = (v[1] * h[2] - v[2] * h[1]) / gm - r[0] / dist;
	e[0] = (VY * h[2] - VZ * h[1]) / GbigMass - X / R;
	//e[1] = (v[2] * h[0] - v[0] * h[2]) / gm - r[1] / dist;
	e[1] = (VZ * h[0] - VX * h[2]) / GbigMass - Y / R;
	//e[2] = (v[0] * h[1] - v[1] * h[0]) / gm - r[2] / dist;
	e[2] = (VX * h[1] - VY * h[0]) / GbigMass - Z / R;
	ecc2 = 0.;
	for( i = 0; i < 3; i++)
	{
		ecc2 += e[i] * e[i];
	}
	//elem->minor_to_major = sqrt( fabs( 1. - ecc2));
	long double BdivA = sqrt( abs( 1. - ecc2));
	//ecc = elem->ecc = sqrt( ecc2);
	ecc = Ecc = sqrt( ecc2);
	for( i = 0; i < 3; i++)
	{
		e[i] /= ecc;
	}
	gm_over_h0 = GbigMass / h0;
	//perihelion_speed = gm_over_h0 + sqrt( gm_over_h0 * gm_over_h0 - inv_major_axis * gm);
	perihelion_speed = gm_over_h0 + sqrt( gm_over_h0 * gm_over_h0 - InvMAxis * GbigMass);
	//elem->q = h0 / perihelion_speed;
	long double q = h0 / perihelion_speed;
	long double major_axis;
	long double t0;
	long double perih_time;
	long double w0;
	if( InvMAxis)
    {
		//elem->major_axis = 1. / inv_major_axis;
		major_axis = 1. / InvMAxis;
		//elem->t0 = elem->major_axis * sqrt( fabs( elem->major_axis) / gm);
		t0 = major_axis * sqrt( abs( major_axis) / GbigMass);
		T = sqrt(4.0*M_PI*M_PI*abs(major_axis)*abs(major_axis)*abs(major_axis)/GbigMass);
		// or T = 2*PI*t0
		// from formula T*T = sqrt(4*pi^2*a^3/(G*m))
    }
	const long double cos_arg_per = (h[0] * e[1] - h[1] * e[0]) / n0;

	if( cos_arg_per < .7 && cos_arg_per > -.7)
		ArgPer = acos( cos_arg_per);
	else
    {
		const long double sin_arg_per = (e[0] * h[0] * h[2] + e[1] * h[1] * h[2] - e[2] * n0 * n0) / (n0 * h0);
		if (n0 * h0)
		{
			ArgPer = abs( asin( sin_arg_per));
			if( cos_arg_per < 0.)
				ArgPer = M_PI - ArgPer;
		}
		else
			ArgPer = 0.0;

	}
	if( e[2] < 0.)
		ArgPer = M_PI + M_PI - ArgPer;

	if( InvMAxis > 0.)         /* elliptical case */
    {
		const long double e_cos_E = 1. - R * InvMAxis;
		const long double e_sin_E = r_dot_v / sqrt( GbigMass * major_axis);
		const long double ecc_anom = atan2( e_sin_E, e_cos_E);

		MeanAnm = ecc_anom - ecc * sin( ecc_anom);
/*    elem->t0 = elem->major_axis * sqrt( elem->major_axis / gm);   */
		//elem->perih_time = t - elem->mean_anomaly * elem->t0;
		//perih_time = t - MeanAnm * t0;
	}
	else if( InvMAxis < 0.)         /* hyperbolic case */
    {
		const long double z = (1. - R * InvMAxis) / ecc;
		long double f = log( z + sqrt( z * z - 1.));

		if( r_dot_v < 0.)
			f = -f;
		MeanAnm = ecc * sinh( f) - f;
		//elem->perih_time = t - elem->mean_anomaly * fabs( elem->t0);
		//perih_time = t - MeanAnm * abs( t0);
		h0 = -h0;
	}
	else              /* parabolic case */
    {
		long double tau;

		tau = sqrt( R / q - 1.);
		if( r_dot_v < 0.)
			tau = -tau;
		w0 = (3. / M_SQRT2) / (q * sqrt( q / GbigMass));
/*    elem->perih_time = t - tau * (tau * tau / 3. + 1) *                   */
/*                                      elem->q * sqrt( 2. * elem->q / gm); */
	    //elem->perih_time = t - tau * (tau * tau / 3. + 1) * 3. / elem->w0;
		//perih_time = t - tau * (tau * tau / 3. + 1) * 3. / w0;
    }

	//for( i = 0; i < 3; i++)
	//	elem->perih_vec[i] = e[i];
	//elem->sideways[0] = (e[2] * h[1] - e[1] * h[2]) / h0;
	//elem->sideways[1] = (e[0] * h[2] - e[2] * h[0]) / h0;
	//elem->sideways[2] = (e[1] * h[0] - e[0] * h[1]) / h0;
	//elem->angular_momentum = h0;
}

//A brief header is given below:
//
//Des'n     H     G   Epoch     M        Peri.      Node       Incl.       e            n           a        Reference #Obs #Opp    Arc    rms  Perts   Computer
//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//00001    3.34  0.12 K118R 199.31747   72.44742   80.37910   10.58595  0.0785634  0.21418897   2.7665466  0 MPO110568  6063  94 1802-2006 0.61 M-v 30h MPCW       0000      (1) Ceres              20061025
//00002    4.13  0.11 K118R 181.70142  310.06331  173.12644   34.84160  0.2312416  0.21365561   2.7711488  0 MPO135420  7253  93 1827-2008 0.59 M-v 28h MPCW       0000      (2) Pallas             20080101
//00003    5.33  0.32 K118R 122.32357  248.18733  169.91025   12.98086  0.2551627  0.22576192   2.6711744  0 MPO135420  6278  91 1824-2008 0.60 M-v 38h MPCW       0000      (3) Juno               20080308
//00004    3.20  0.32 K118R  56.40276  149.88186  103.90428    7.13461  0.0883124  0.27164416   2.3612165  0 MPO135420  6449  82 1827-2007 0.60 M-v 18h MPCW       0000      (4) Vesta              20070808
//00005    6.85  0.15 K118R 334.02337  358.78351  141.61007    5.36742  0.1903560  0.23868879   2.5738390  0 MPO 41498  2178 100 1845-2002 0.56 M-v 38h Goffin     0000      (5) Astraea            20021114
//00006    5.71  0.24 K118R  76.23230  239.28896  138.72650   14.75022  0.2022525  0.26104401   2.4247124  0 MPO135420  5011  81 1848-2008 0.55 M-v 38h MPCW       0007      (6) Hebe               20080222
//00007    5.51  0.15 K118R 111.53561  145.26277  259.65080    5.52306  0.2305428  0.26737788   2.3862672  0 MPO135420  4443  76 1848-2008 0.60 M-v 38h MPCW       0000      (7) Iris               20080308
//00008    6.49  0.28 K118R  70.05922  285.18420  110.93525    5.88786  0.1564761  0.30182794   2.2010492  0 MPO 25662  2438 102 1847-2000 0.54 M-v 38h Goffin     0000      (8) Flora              20001115
//00009    6.28  0.17 K118R 248.44234    6.24383   68.94919    5.57547  0.1227453  0.26749937   2.3855446  0 MPO 25662  2131  96 1822-2001 0.54 M-v 38h Goffin     0000      (9) Metis              20011024
//00010    5.43  0.15 K118R  15.29657  313.01188  283.41502    3.84094  0.1164623  0.17722252   3.1389913  0 MPC 35055  2261 104 1849-1999 0.53 M-v 38h Goffin     0000     (10) Hygiea             19990306
//00011    6.55  0.15 K118R 281.05531  195.04403  125.60104    4.62566  0.0990377  0.25651561   2.4531655  0 MPO135420  3974  82 1850-2007 0.63 M-v 38h MPCW       0000     (11) Parthenope         20070725
//00012    7.24  0.22 K118R  79.04846   69.82359  235.48673    8.36699  0.2207436  0.27628481   2.3347016  0 MPO 57621  2635  86 1850-2003 0.56 M-v 38h Goffin     0000     (12) Victoria           20030921
//00013    6.74  0.15 K118R 240.86929   80.14640   43.27015   16.54560  0.0853678  0.23830803   2.5765798  0 MPO  2601   946  54 1850-2000 0.77 M-v 38h Williams   0000     (13) Egeria             20000131
//00014    6.30  0.15 K118R 221.96593   97.35232   86.20615    9.10849  0.1656321  0.23651727   2.5895689  0 MPO   595  1765  84 1851-2000 0.54 M-v 38h Goffin     0000     (14) Irene              20000402
//00015    5.28  0.23 K118R   2.13972   97.67316  293.21965   11.73525  0.1886022  0.22924203   2.6440716  0 MPO 25662  2071  91 1851-2001 0.52 M-v 38h Goffin     0000     (15) Eunomia            20010712
//00016    5.90  0.20 K118R  98.10135  226.95821  150.29920    3.09929  0.1374216  0.19745426   2.9207320  0 MPC 35055  1831  96 1852-1998 0.53 M-v 38h Goffin     0000     (16) Psyche             19980614

//this function transformed from original code to cross check another functions
void setup_orbit_vectors( long double &e_epoch, long double &e_ecc, long double &e_incl, long double &e_asc_node, long double &e_arg_per,long double &e_mean_anomaly,
	long double &e_q, long double &e_major_axis,long double &e_t0,long double &e_w0,long double &e_angular_momentum,long double &e_perih_time,
	long double &e_minor_to_major, long double &e_lon_per,
	long double &e_sideways_x, long double &e_sideways_y, long double &e_sideways_z,
	long double &vec_x, long double &vec_y,long double &vec_z,long double gm)
{
   long double sin_incl = sin( e_incl), cos_incl = cos( e_incl);
   //double *vec;
   long double vec_len;
   //double up[3];
   long double up_x, up_y, up_z;
   int i;

   e_minor_to_major = sqrt( abs( (long double)1. - e_ecc * e_ecc));
   e_lon_per = e_asc_node + atan2( sin( e_arg_per) * cos_incl,
                                       cos( e_arg_per));
   //vec = e->perih_vec;

   vec_x = cos( e_lon_per);
   vec_y = sin( e_lon_per);
   vec_z = (sin_incl / cos_incl) * sin( e_lon_per - e_asc_node);
   //vec_len = sqrt( 1. + vec[2] * vec[2]);
   vec_len = sqrt( 1. + vec_z * vec_z);
   //for( i = 0; i < 3; i++)
   //   vec[i] /= vec_len;
   vec_x /= vec_len;
   vec_y /= vec_len;
   vec_z /= vec_len;
            /* 'up' is a vector perpendicular to the plane of the orbit */
   up_x =  sin( e_asc_node) * sin_incl;
   up_y = -cos( e_asc_node) * sin_incl;
   up_z = cos_incl;

   e_sideways_x = up_y * vec_z - up_z * vec_y;
   e_sideways_y = up_z * vec_x - up_x * vec_z;
   e_sideways_z = up_x * vec_y - up_y * vec_x;
}
// this function transformed from original code to cross check another functions
void derive_quantities( long double &e_epoch, long double &e_ecc, long double &e_incl, long double &e_asc_node, long double &e_arg_per,long double &e_mean_anomaly,
	long double &e_q, long double &e_major_axis,long double &e_t0,long double &e_w0,long double &e_angular_momentum,long double &e_perih_time,
	long double &e_minor_to_major, long double &e_lon_per,
	long double &e_sideways_x, long double &e_sideways_y, long double &e_sideways_z,
	long double &vec_x, long double &vec_y,long double &vec_z,long double gm)
{
   if( e_ecc != 1.)    /* for non-parabolic orbits: */
      {
      e_major_axis = e_q / abs(1. - e_ecc);
      e_t0 = e_major_axis * sqrt( e_major_axis / gm);
      }
   else
      {
      e_w0 = (3. / sqrt((long double)2.0)) / (e_q * sqrt( e_q / gm));
      e_major_axis = e_t0 = 0.;
      }
  //setup_orbit_vectors( e);
}
// this function transformed from original code to cross check another functions
void do_element_setup( long double &e_epoch, long double &e_ecc, long double &e_incl, long double &e_asc_node, long double &e_arg_per,long double &e_mean_anomaly,
	long double &e_q, long double &e_major_axis,long double &e_t0,long double &e_w0,long double &e_angular_momentum,long double &e_perih_time,
	long double &e_minor_to_major, long double &e_lon_per,
	long double &e_sideways_x, long double &e_sideways_y, long double &e_sideways_z,
	long double &vec_x, long double &vec_y,long double &vec_z,long double gm, long double T)
{
	//elem->mean_anomaly  *= PI / 180.;
	//elem->arg_per       *= PI / 180.;
	//elem->asc_node      *= PI / 180.;
	//elem->incl          *= PI / 180.;
	long double Temp = (((gm * T * T / 4.0) / M_PI) / M_PI);
    long double A = pow(Temp, (long double)1.0 / (long double)3.0);

    //A= EarthSmAx;
    long double B = A * sqrt(1.0 - e_ecc * e_ecc);
	e_major_axis =  A;
	e_q = e_major_axis * (1. - e_ecc);
	//derive_quantities( elem, gm);
	if( e_ecc != 1.)    /* for non-parabolic orbits: */
	{
		e_major_axis = e_q / abs(1. - e_ecc);
		e_t0 = e_major_axis * sqrt( e_major_axis / gm);
    }
	else
    {
		e_w0 = (3. / sqrt((long double)2.0)) / (e_q * sqrt( e_q / gm));
		e_major_axis = e_t0 = 0.;
    }
	//setup_orbit_vectors( e);
	{
		long double sin_incl = sin( e_incl), cos_incl = cos( e_incl);
		//double *vec;
		long double vec_len;
		//double up[3];
		long double up_x, up_y, up_z;
		int i;

		e_minor_to_major = sqrt( abs( (long double)1. - e_ecc * e_ecc));
		e_lon_per = e_asc_node + atan2( sin( e_arg_per) * cos_incl,
                                       cos( e_arg_per));
		//vec = e->perih_vec;

		vec_x = cos( e_lon_per);
		vec_y = sin( e_lon_per);
		vec_z = (sin_incl / cos_incl) * sin( e_lon_per - e_asc_node);
		//vec_len = sqrt( 1. + vec[2] * vec[2]);
		vec_len = sqrt( 1. + vec_z * vec_z);
		//for( i = 0; i < 3; i++)
		//   vec[i] /= vec_len;
		vec_x /= vec_len;
		vec_y /= vec_len;
		vec_z /= vec_len;
            /* 'up' is a vector perpendicular to the plane of the orbit */
		up_x =  sin( e_asc_node) * sin_incl;
		up_y = -cos( e_asc_node) * sin_incl;
		up_z = cos_incl;

		e_sideways_x = up_y * vec_z - up_z * vec_y;
		e_sideways_y = up_z * vec_x - up_x * vec_z;
		e_sideways_z = up_x * vec_y - up_y * vec_x;
	}   
	e_angular_momentum = sqrt( gm * e_q * (1. + e_ecc));
	e_perih_time = e_epoch - e_mean_anomaly * e_t0;
	//e_is_asteroid = 1;
	//elem->central_obj = 0;
}


#define MAX_ITERATIONS 7
long double asinh( long double z)
{
   return( log( z + sqrt( z * z + 1.)));
}
long double near_parabolic( long double ecc_anom, long double e)
{
   long double anom2 = (e > 1. ? ecc_anom * ecc_anom : -ecc_anom * ecc_anom);
   long double term = e * anom2 * ecc_anom / 6.;
   long double rval = (1. - e) * ecc_anom - term;
   int n = 4;

   while( fabs( term) > 1e-15)
      {
      term *= anom2 / (double)(n * (n + 1));
      rval -= term;
      n += 2;
      }
   return( rval);
}

#define THRESH 1.e-8
#define MIN_THRESH 1.e-15
// from the same source - good sample of a keplers solution
long double kepler(long double ecc, long double mean_anom)
{
   long double curr, err, thresh, offset = 0.;
   long double delta_curr = 1.;
   int is_negative = 0, n_iter = 0;

   if( !mean_anom)
      return( 0.);

   if( ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
      {
      curr = atan2( sin( mean_anom), cos( mean_anom) - ecc);
            /* two correction steps,  and we're done */
      for( n_iter = 2; n_iter; n_iter--)
         {
         err = curr - ecc * sin( curr) - mean_anom;
         curr -= err / (1. - ecc * cos( curr));
         }
      return( curr);
      }

   if( ecc < 1.)
      if( mean_anom < -M_PI || mean_anom > M_PI)
         {
         double tmod = fmod( mean_anom, (long double)M_PI * (long double)2.);

         if( tmod > M_PI)             /* bring mean anom within -pi to +pi */
            tmod -= 2. * M_PI;
         else if( tmod < -M_PI)
            tmod += 2. * M_PI;
         offset = mean_anom - tmod;
         mean_anom = tmod;
         }

   if( mean_anom < 0.)
      {
      mean_anom = -mean_anom;
      is_negative = 1;
      }

   curr = mean_anom;
   thresh = THRESH * fabs( 1. - ecc);
               /* Due to roundoff error,  there's no way we can hope to */
               /* get below a certain minimum threshhold anyway:        */
   if( thresh < MIN_THRESH)
      thresh = MIN_THRESH;
   if( (ecc > .8 && mean_anom < M_PI / 3.) || ecc > 1.)    /* up to 60 degrees */
      {
      double trial = mean_anom / fabs( 1. - ecc);

      if( trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
         {
         if( mean_anom < M_PI)
            trial = pow( 6. * mean_anom,(long double)1.0/(long double)3.0);
         else        /* hyperbolic w/ 5th & higher-order terms predominant */
            trial = asinh( mean_anom / ecc);
         }
      curr = trial;
      }
   if( ecc > 1. && mean_anom > 4.)    /* hyperbolic, large-mean-anomaly case */
      curr = log( mean_anom);

   if( ecc < 1.)
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = near_parabolic( curr, ecc) - mean_anom;
         else
            err = curr - ecc * sin( curr) - mean_anom;
         delta_curr = -err / (1. - ecc * cos( curr));
         curr += delta_curr;
         }
   else
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = -near_parabolic( curr, ecc) - mean_anom;
         else
            err = ecc * sinh( curr) - curr - mean_anom;
         delta_curr = -err / (ecc * cosh( curr) - 1.);
         curr += delta_curr;
         }
   return( is_negative ? offset - curr : offset + curr);
}

// this function transformed from original code to cross check another functions
int posn_and_vel( long double &e_epoch, long double &e_ecc, long double &e_incl, long double &e_asc_node, long double &e_arg_per,long double &e_mean_anomaly,
	long double &e_q, long double &e_major_axis,long double &e_t0,long double &e_w0,long double &e_angular_momentum,long double &e_perih_time,
	long double &e_minor_to_major, long double &e_lon_per,
	long double &e_sideways_x, long double &e_sideways_y, long double &e_sideways_z,
	long double &e_perih_vec_x,long double &e_perih_vec_y,long double &e_perih_vec_z,
	long double &loc_x, long double &loc_y, long double &loc_z, long double &loc_r,long double &vel_x,long double &vel_y,long double &vel_z,long double t, long double gm)
{
	t -= e_perih_time;
	if( e_ecc != 1.)    /* not parabolic */
	{
		t /= e_t0;
		if( e_ecc < 1.)     /* elliptical case;  throw out extra orbits */
        {                    /* to fit mean anom between -PI and PI */
			t = fmod( (long double)t, (long double)M_PI * 2.);
			if( t < -M_PI) t += 2. * M_PI;
			if( t >  M_PI) t -= 2. * M_PI;
		}
		e_mean_anomaly = t;
	}
	//posn_part_ii( elem, t, loc, vel);
	{
		long double true_anom, r, x, y, r0;

		if( e_ecc == 1.)    /* parabolic */
		{
			long double g = e_w0 * t * .5;

			y = pow( g + sqrt( g * g + 1.), (long double)1.0/ (long double)3.0);
			true_anom = 2. * atan( y - 1. / y);
		}
		else           /* got the mean anomaly;  compute eccentric,  then true */
		{
			long double ecc_anom;

			ecc_anom = kepler( e_ecc, e_mean_anomaly);
			if( e_ecc > 1.)     /* hyperbolic case */
			{
				x = (e_ecc - cosh( ecc_anom));
				y = sinh( ecc_anom);
			}
			else           /* elliptical case */
			{
				x = (cos( ecc_anom) - e_ecc);
				y =  sin( ecc_anom);
			}
			y *= e_minor_to_major;
			true_anom = atan2( y, x);
		}

		r0 = e_q * (1. + e_ecc);
		r = r0 / (1. + e_ecc * cos( true_anom));
		x = r * cos( true_anom);
		y = r * sin( true_anom);
		loc_x = e_perih_vec_x * x + e_sideways_x * y;
		loc_y = e_perih_vec_y * x + e_sideways_y * y;
		loc_z = e_perih_vec_z * x + e_sideways_z * y;
		loc_r = r;
		//if( vel && (elem->angular_momentum != 0.))
		{
			long double angular_component = e_angular_momentum / (r * r);
			long double radial_component = e_ecc * sin( true_anom) *
                                e_angular_momentum / (r * r0);
			long double x1 = x * radial_component - y * angular_component;
			long double y1 = y * radial_component + x * angular_component;
			int i;

			//for( i = 0; i < 3; i++)
			//   vel[i] = elem->perih_vec[i] * x1 + elem->sideways[i] * y1;
			vel_x = e_perih_vec_x * x1 + e_sideways_x * y1;
			vel_y = e_perih_vec_y * x1 + e_sideways_y * y1;
			vel_z = e_perih_vec_z * x1 + e_sideways_z * y1;
		}
	}

   return( 0);
}
// today value is in Gbig * Mass = 398600441499999.87 in m**3/s**2
//#define BIG_XKE 7.43669161331734132e-2

// was based in SGP4 on G=398600.8 km**3/s**2
#define BIG_XKE .743669161E-1


// see http://www.amsat.org/amsat/keps/kepmodel.html
// or http://en.wikipedia.org/wiki/Orbital_elements
//    SatEpoch  - satellite epoch
//    T         - orbit period in sec 
//                 for Mean Motion (revolutions/day) needs T = 60*60*24 / Mean_Motion
//    Ecc       - Eccentricity
//    Incl      - Inclination
//    AssNode   - Longitude of ascending node (where the orbit passes upward through the reference plane)
//    ArgPer    - Argument of perihelion angle to per (lower point) point
//    MeanAnm   - Mean Anomaly (degrees)
//    Mass      - Mass of first body
//    mass      - Mass of second body
//                assuming:
//     Z  - pointed up
//     Y  - to right
//     X  - to point of view
//
//                or:
//     Z  - point north from ecliptic
//     X  - to easter
//     Y  - to point of view
// XY plain is equatorial plain for earth's satellites, or for planets it is a ecliptic 
// plane (the plane of the Earth's orbit around the Sun). In both cases X points to easter/
// this function now does not care about drag - which is bad 
void KeplerPosition(long double SatEpoch, long double CurTime, long double T,long double Ecc, long double Incl, long double AssNode, long double ArgPer, long double MeanAnm, long double BSTAR, long double GM, int doCorrection, 
                    long double &Xm, long double &Ym, long double &Zm, long double &VX, long double &VY, long double &VZ, long double ProbMeanMotion)
{
    long double XMNPDA_XMNPDA = 24.0*60.0*60.0;//1440.0; // XMNPDA time units(minutes) /day 1440.0
	long double XMNPDA = 24.0*60.0;//1440.0; // XMNPDA time units(minutes) /day 1440.0
	long double TEMP_=2*M_PI/XMNPDA/XMNPDA; // 2*pi / (1440 **2)
	long double TEMP_TEMP_=2*M_PI/XMNPDA_XMNPDA/XMNPDA_XMNPDA; // 2*pi / (1440 **2)
	long double XNO=ProbMeanMotion*TEMP_*XMNPDA; // rotation per day * 2*pi /1440 == rotation per day on 1 unit (1 min/sec)
	long double XNO_XNO=ProbMeanMotion*TEMP_TEMP_*XMNPDA_XMNPDA; // rotation per day * 2*pi /1440 == rotation per day on 1 unit (1 min/sec)
	long double XKE = BIG_XKE;//.743669161E-1;
	long double XKMPER = 6378.1350;
	//IF (IFLAG .EQ. 0) GO TO 100
	//* RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
	//* FROM INPUT ELEMENTS
	//A1=(XKE/XNO)**TOTHRD;
	long double XJ2 = 1.082616E-3;
	long double XJ3 = -.253881E-5;
	long double XJ4 = -1.65597E-6;
	long double AE = 1.0;
	long double QO =120.0;
	long double SO = 78.0;

	long double CK2=.5*XJ2*AE*AE;
	//CK4=-.375*XJ4*AE**4
	long double CK4=-.375*XJ4*AE*AE*AE*AE;
	long double QOMS2T=pow(((QO-SO)*AE/XKMPER),(long double)4.0);
	long double S=AE*(1.+SO/XKMPER);
	//A1=(XKE/XNO)**TOTHRD;
	long double A1=pow((XKE/XNO),(long double)2.0/(long double)3.0);
	long double A1_A1=pow((XKE/XNO_XNO),(long double)2.0/(long double)3.0);
	long double COSIO=cos(Incl);
	long double THETA2=COSIO*COSIO;
	long double X3THM1=3.*THETA2-1.;
	long double EOSQ=Ecc*Ecc;
	long double BETAO2=1.-EOSQ;
	long double BETAO=sqrt(BETAO2);
	long double DEL1=1.5*CK2*X3THM1/(A1*A1*BETAO*BETAO2);
	//AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
	long double AO=A1*(1.-DEL1*(.5*(2.0/3.0)+DEL1*(1.+134./81.*DEL1)));
	long double DELO=1.5*CK2*X3THM1/(AO*AO*BETAO*BETAO2);
    // ProbMeanMotion = XNO / (2*pi) * 1440.0
	long double XNODP=XNO/(1.+DELO);
    if (doCorrection)
    {
        T = 1.0/ (XNODP / M_PI/2.0 * 1440.0) * 24.0 * 60.0 * 60.0;
    }
    // and this is a semimajor axis
	long double AODP=AO/(1.-DELO);
    if (doCorrection)
    {
    }
	//* INITIALIZATION
	//* FOR PERIGEE LESS THAN 220 KILOMETERS, THE ISIMP FLAG IS SET AND
	//* THE EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN SQRT A AND
	//* QUADRATIC VARIATION IN MEAN ANOMALY. ALSO, THE C3 TERM, THE
	//* DELTA OMEGA TERM, AND THE DELTA M TERM ARE DROPPED.
	int ISIMP=0;
	//IF((AODP*(1.-EO)/AE) .LT. (220./XKMPER+AE)) ISIMP=1
	if((AODP*(1.-Ecc)/AE) < (220./XKMPER+AE)) 
		ISIMP=1;
	//* FOR PERIGEE BELOW 156 KM, THE VALUES OF
	//* S AND QOMS2T ARE ALTERED
	long double S4=S;
	long double QOMS24=QOMS2T;
	long double PERIGE=(AODP*(1.-Ecc)-AE)*XKMPER;
	//IF(PERIGE .GE. 156.) GO TO 10
	if (PERIGE >= 156.) 
		goto M_10;
	S4=PERIGE-78.;
	//IF(PERIGE .GT. 98.) GO TO 9
	if (PERIGE >=98.) 
		goto M_9;
	S4=20.;
M_9:
	//QOMS24=((120.-S4)*AE/XKMPER)**4;
	QOMS24=pow(((120.-S4)*AE/XKMPER),4);
	S4=S4/XKMPER+AE;
M_16:
M_10:
	long double PINVSQ=1./(AODP*AODP*BETAO2*BETAO2);
	long double TSI=1./(AODP-S4);
	long double ETA=AODP*Ecc*TSI;
	long double ETASQ=ETA*ETA;
	long double EETA=Ecc*ETA;
	long double PSISQ=abs(1.-ETASQ);
	//COEF=QOMS24*TSI**4
	long double COEF=QOMS24*TSI*TSI*TSI*TSI;
	//COEF1=COEF/PSISQ**3.5;
	long double COEF1=COEF/pow(PSISQ,(long double)3.5);

	long double C2=COEF1*XNODP*(AODP*(1.+1.5*ETASQ+EETA*(4.+ETASQ))+.75* CK2*TSI/PSISQ*X3THM1*(8.+3.*ETASQ*(8.+ETASQ)));
	long double C1=BSTAR*C2;
	long double SINIO=sin(Incl);
	//A3OVK2=-XJ3/CK2*AE**3;
	long double A3OVK2=-XJ3/CK2*(AE*AE*AE);
	long double C3=COEF*TSI*A3OVK2*XNODP*AE*SINIO/Ecc;
	long double X1MTH2=1.-THETA2;
	long double C4=2.*XNODP*COEF1*AODP*BETAO2*(ETA*(2.+.5*ETASQ)+Ecc*(.5+2.*ETASQ)-2.*CK2*TSI/
			(AODP*PSISQ)*(-3.*X3THM1*(1.-2.*EETA+ETASQ* (1.5-.5*EETA))+.75*X1MTH2*(2.*ETASQ-EETA* (1.+ETASQ))*cos(2.*ArgPer)));
	long double C5=2.*COEF1*AODP*BETAO2*(1.+2.75*(ETASQ+EETA)+EETA*ETASQ);
	long double THETA4=THETA2*THETA2;
	long double TEMP1=3.*CK2*PINVSQ*XNODP;
	long double TEMP2=TEMP1*CK2*PINVSQ;
	long double TEMP3=1.25*CK4*PINVSQ*PINVSQ*XNODP;
	long double XMDOT=XNODP+.5*TEMP1*BETAO*X3THM1+.0625*TEMP2*BETAO* (13.-78.*THETA2+137.*THETA4);
	long double X1M5TH=1.-5.*THETA2;
	long double OMGDOT=-.5*TEMP1*X1M5TH+.0625*TEMP2*(7.-114.*THETA2+ 395.*THETA4)+TEMP3*(3.-36.*THETA2+49.*THETA4);
	long double XHDOT1=-TEMP1*COSIO;
	long double XNODOT=XHDOT1+(.5*TEMP2*(4.-19.*THETA2)+2.*TEMP3*(3.-7.*THETA2))*COSIO;
	long double OMGCOF=BSTAR*C3*cos(ArgPer);
	long double XMCOF=-(2.0/3.0)*COEF*BSTAR*AE/EETA;
	long double XNODCF=3.5*BETAO2*XHDOT1*C1;
	long double T2COF=1.5*C1;
	long double XLCOF=.125*A3OVK2*SINIO*(3.+5.*COSIO)/(1.+COSIO);
	long double AYCOF=.25*A3OVK2*SINIO;
	//DELMO=(1.+ETA*COS(XMO))**3
	long double DELMO=pow((1.+ETA*cos(MeanAnm)),3);
	long double SINMO=sin(MeanAnm);
	long double X7THM1=7.*THETA2-1.;
	// IF(ISIMP .EQ. 1) GO TO 90
	if (ISIMP == 1) 
		goto M_90;
	long double C1SQ=C1*C1;
	long double D2=4.*AODP*TSI*C1SQ;
	long double TEMP=D2*TSI*C1/3.;
	long double D3=(17.*AODP+S4)*TEMP;
	long double D4=.5*TEMP*AODP*TSI*(221.*AODP+31.*S4)*C1;
M_17:
	long double T3COF=D2+2.*C1SQ;
	long double T4COF=.25*(3.*D3+C1*(12.*D2+10.*C1SQ));
	long double T5COF=.2*(3.*D4+12.*C1*D3+6.*D2*D2+15.*C1SQ*(2.*D2+C1SQ));
M_90:
	//IFLAG=0;
	//* UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
M_100:
    long double TSINCE = (CurTime - SatEpoch)*(24.0 * 60.0 * 60.0);
	long double XMDF=MeanAnm+XMDOT*TSINCE;
	long double OMGADF=ArgPer+OMGDOT*TSINCE;
	long double XNODDF=AssNode+XNODOT*TSINCE;
	long double OMEGA=OMGADF;
	long double XMP=XMDF;
	long double TSQ=TSINCE*TSINCE;
	long double XNODE=XNODDF+XNODCF*TSQ;
	long double TEMPA=1.-C1*TSINCE;
	long double TEMPE=BSTAR*C4*TSINCE;
	long double TEMPL=T2COF*TSQ;
	// IF(ISIMP .EQ. 1) GO TO 110
	if (ISIMP == 1) 
		goto M_110;
	long double DELOMG=OMGCOF*TSINCE;
	// DELM=XMCOF*((1.+ETA*COS(XMDF))**3-DELMO)
	long double DELM=XMCOF*(pow((1.+ETA*cos(XMDF)),3)-DELMO);
	TEMP=DELOMG+DELM;
	XMP=XMDF+TEMP;
	OMEGA=OMGADF-TEMP;
	long double TCUBE=TSQ*TSINCE;
	long double TFOUR=TSINCE*TCUBE;
	TEMPA=TEMPA-D2*TSQ-D3*TCUBE-D4*TFOUR;
	TEMPE=TEMPE+BSTAR*C5*(sin(XMP)-SINMO);
	TEMPL=TEMPL+T3COF*TCUBE+TFOUR*(T4COF+TSINCE*T5COF);
M_110:
	//A=AODP*TEMPA**2;
	//long double A=AODP*TEMPA*TEMPA;
	long double Eee=Ecc-TEMPE;



    long double Omega0 = ArgPer;

	long double TettaDelta = 0.0;

    // Tetta = contrclockwise angle from perihelion
	// Mean Anomaly
    // [aka "M0" or "MA" or "Phase"]
	// 	Mean anomaly is simply an angle that marches uniformly in time from 0 to 360 degrees during one revolution. It is defined to be 0 degrees at perigee, and therefore is 180 degrees at apogee. 
	// The mean anomaly increases uniformly from 0 to 2*PI() radians during each orbit. However, it is not an angle. 
	// Due to Kepler's second law, the mean anomaly is proportional to the area swept by the focus-to-body line since the last periapsis.
    // The mean anomaly is usually denoted by the letter M, and is given by the formula
	// cos(TrueAnomaly) = (cos(E) - e) / (1 - e*cos(E))
	// and E - is a essentric anomaly formula is:
	// M = E - e*sin(E) 
	// first calculates how many orbits was done from a epoche time till current moment of time
	TettaDelta = CurTime - SatEpoch;
	TettaDelta  *= (24.0 * 60.0 * 60.0);
	TettaDelta /= T;
	int iTettaDelta = TettaDelta;
	TettaDelta -= iTettaDelta;
	TettaDelta *= 2.0 * M_PI;


	long double M = MeanAnm + TettaDelta;
    if (doCorrection)
    {
        M = XMP;
        Ecc = Eee;
        ArgPer = OMEGA;
        AssNode = XNODE;
    }
    //// find E
    long double E = M;
    long double ENext = E - Ecc * sin(E);
    long double DeltaE = abs(ENext - M);
	long double PresM = abs(M) * 1.0E-15;
	int iAdditionalAttempts = 3;
    while (iAdditionalAttempts > 0)
    {
        ENext = E - Ecc * sin(E);
        DeltaE = abs(ENext - M);
		if (PresM >= DeltaE)
			iAdditionalAttempts--;
	
		E = E - ( ENext - M) /2;
    }

	long double CosPhi = (cos(E) - Ecc)/(1.0 -Ecc*cos(E));



	long double Tetta = acos(CosPhi);
	// for conformation see: http://www.jgiesen.de/kepler/kepler.html
    
	// first calulate a = (G*M1*T*T/(4*Pi*Pi))in power 1/3
    long double Temp = ((((GM * T * T) / 4.0) / M_PI) / M_PI);
    long double A = pow(Temp, (long double)1.0 / (long double)3.0);

    //A= EarthSmAx;
    long double B = A * sqrt(1.0 - Ecc * Ecc);

    printf("\n Calc A      = %f", A);
	// calulates R
    // X - pointed right and to perihelion
    // Y - pointed down
    // Z is not in use X-Y is a orbit plain
    // see pic1.bmp
    long double R = (B * B / ( 1.0 + Ecc * cos(Tetta))) / A;
    long double X0 = R*cos(Tetta);
	long double Y0 = R*sin(Tetta);
	long double Z0 = 0.0;

    long double VR = sqrt(((GM))*A/B/B) * Ecc * sin(Tetta);
    long double VN = sqrt(((GM))*A/B/B) * (1.0 +Ecc *cos(Tetta));

    int DirectionContrClock = 1;

    if (AssNode >=0.0 && AssNode < M_PI)
        DirectionContrClock = 0;
    long double VXm0;
    long double VYm0;
    long double VZm0;
    long double M_TrAnom[4]={-1,1,1,1};
	long double M_argper[4] = {1,-1,1,1};
	long double M_incl[4] = {1,-1,1,1};
	long double M_AsNode[4] = {1,-1,1,1};
    //if (DirectionContrClock)
    //{
    //    VXm0 =  -sin(Tetta) * VN - cos(Tetta) * VR;//- cos(Tetta) * VN + sin(Tetta) * VR;// - cos(Tetta) * VN + sin(Tetta) * VR;
    //    VYm0 = cos(Tetta) * VN - sin(Tetta) * VR;//- sin(Tetta) * VN - cos(Tetta) * VR;//  sin(Tetta) * VN - cos(Tetta) * VR;
    //    VZm0 = 0.0;
    //
    //}
    //else
    {
        VXm0 = M_TrAnom[0]*sin(Tetta) * VN + M_TrAnom[1]*cos(Tetta) * VR;
        VYm0 = M_TrAnom[2]*cos(Tetta) * VN + M_TrAnom[3]*sin(Tetta) * VR;
        VZm0 = 0.0;
    }

    // perihelium direction adjust:
    // rotation vector (X0,Y0) by angle omega clockwise will get vector (X1,Y1)
    //
    
    long double X1 = X0 * cos(ArgPer) - Y0 * sin(ArgPer);
    long double Y1 = X0 * sin(ArgPer) + Y0 * cos(ArgPer);
    long double Z1 = Z0;

    long double VXm1 = M_argper[0]*VXm0 * cos(ArgPer) + M_argper[1]* VYm0 * sin(ArgPer);
    long double VYm1 = M_argper[2]*VXm0 * sin(ArgPer) + M_argper[3]*VYm0 * cos(ArgPer);
    long double VZm1 = VZm0;

    // inclanation direction adjust 
    // incl (satellite) = incl
    long double X2 = X1;
    long double Y2 = Y1 * cos(Incl) - Z1 * sin(Incl);
    long double Z2 = Y1 * sin(Incl) + Z1 * cos(Incl);;

    long double VXm2 = VXm1;
    long double VYm2 = M_incl[0]*VYm1 * cos(Incl) + M_incl[1]* VZm1 * sin(Incl);
    long double VZm2 = M_incl[2]*VYm1 * sin(Incl) + M_incl[3]*VZm1 * cos(Incl);

    // Assending Node last opartion 
    long double X3 = X2 * cos(AssNode) - Y2 * sin(AssNode);
    long double Y3 = X2 * sin(AssNode) + Y2 * cos(AssNode);
    long double Z3 = Z2;

    long double VXm3 = M_AsNode[0]*VXm2 * cos(AssNode) + M_AsNode[1]*VYm2 * sin(AssNode);
    long double VYm3 = M_AsNode[2]*VXm2 * sin(AssNode) + M_AsNode[3]*VYm2 * cos(AssNode);
    long double VZm3 = VZm2;

    Xm = X3;
    Ym = Y3;
    Zm = Z3;
	
    VX = VXm3;
    VY = VYm3;
    VZ = VZm3;
}
// gone for now
void AjustKeplerPosition(long double &T,long double &Ap, long double &Ph, long double &SmAx,long double &Ecc, long double &Incl, long double &AssNode, long double &ArgPer, long double EpochS, long double CurEpochS)
{
    //TBD
    // this is basicaly SGP4 formula - but it gives 1km at epoch
    // may be is there another formula to calculate? then just original model
    // just simulation sun-earth-mooon with 1/8 sec delta time gives error 1km per 1 month
}
void ConvertJulianDayToDateAndTime(double JulianDay, SYSTEMTIME *ThatTime)
{
    long daysfrom2000 = JulianDay - 2451544.5;
    double flInDay = (JulianDay - 2451544.5) - (double)daysfrom2000; 
    int iYear = 0;
    daysfrom2000 += 1; // 1Jan must be 1;
    while (daysfrom2000 > 366)
    {
        switch(iYear)
        {
        case 24:daysfrom2000-=366;break;
        case 23:daysfrom2000-=365;break;
        case 22:daysfrom2000-=365;break;
        case 21:daysfrom2000-=365;break;
        case 20:daysfrom2000-=366;break;
        case 19:daysfrom2000-=365;break;
        case 18:daysfrom2000-=365;break;
        case 17:daysfrom2000-=365;break;
        case 16:daysfrom2000-=366;break;
        case 15:daysfrom2000-=365;break;
        case 14:daysfrom2000-=365;break;
        case 13:daysfrom2000-=365;break;
        case 12:daysfrom2000-=366;break;
        case 11:daysfrom2000-=365;break;
        case 10:daysfrom2000-=365;break;
        case 9:daysfrom2000-=365;break;
        case 8:daysfrom2000-=366;break;
        case 7:daysfrom2000-=365;break;
        case 6:daysfrom2000-=365;break;
        case 5:daysfrom2000-=365;break;
        case 4:daysfrom2000-=366;break;
        case 3:daysfrom2000-=365;break;
        case 2:daysfrom2000-=365;break;
        case 1:daysfrom2000-=365;break;
        case 0:daysfrom2000-=366;break;
        }
        iYear++;
    }
    ThatTime->wYear = iYear+2000;
    int iMonth =1;
    int iComp, iDecr;
    while(1)
    {
        switch(iMonth)
        {
        case 1: iComp = 31; iDecr = 31; break; // jan
        case 2: if (iYear%4 ==0) //leap year    //feb
                {
                    iComp = 29; iDecr = 29;
                }
                else
                {
                    iComp = 28; iDecr = 28;
                }
                break;
        case 3: iComp = 31; iDecr = 31; break; // mar
        case 4: iComp = 30; iDecr = 30; break; // apr
        case 5: iComp = 31; iDecr = 31; break; // may
        case 6: iComp = 30; iDecr = 30; break; // jun
        case 7: iComp = 31; iDecr = 31; break; //jul
        case 8: iComp = 31; iDecr = 31; break; // aug
        case 9: iComp = 30; iDecr = 30; break; //sep
        case 10: iComp = 31; iDecr = 31; break; // oct
        case 11: iComp = 30; iDecr = 30; break; // nov
        }
        if (daysfrom2000 < iComp)
                break;
        daysfrom2000 -=iDecr;
        if (++iMonth == 12) // what ?? getout!!!
            break;
    }
    int iDay  = daysfrom2000;//+1;
    
    int iHour = flInDay*24;
    int iMinutes = ((flInDay - ((double)iHour)/24.0))*(24.0*60.0);
    int iSec =   ((flInDay - ((double)iHour)/24.0) - ((double)iMinutes)/(24.0*60.0))*(24.0*60.0*60.0);
    int iMils = (flInDay - ((double)iHour)/(24.0) - ((double)iMinutes)/(24.0*60.0) - ((double)iSec)/(24.0*60.0*60.0))*(24.0*60.0*60.0*1000.0);
    ThatTime->wMonth = iMonth;
    ThatTime->wDay = iDay;
    ThatTime->wHour = iHour;
    ThatTime->wMinute = iMinutes;
    ThatTime->wSecond = iSec;
    ThatTime->wMilliseconds = iMils;
    ThatTime->wDayOfWeek = 0;
}
long double ConverEpochDate2JulianDay(long double KeplerDate)
{
    // TLE elements is 1 day based - needs to minus at the end one day
    int iYear = KeplerDate /1000;
    // date as it is = 2000/01/01     2451544.5, 2451910.5, 2452275.5, 2452640.5, 2453005.5, 2453371.5, 2453736.5, 2013-2456293.5
    // 
    double t2000_01_01_01 = 2451544.5;
    switch(iYear)
    {
        // add years = if you still alive !!! or just put formula if ((iYear-1)%4 == 0) t2000_01_01_01+=366; else t2000_01_01_01+=365;
    case 24:t2000_01_01_01+=365;
    case 23:t2000_01_01_01+=365;
    case 22:t2000_01_01_01+=365;
    case 21:t2000_01_01_01+=366;
    case 20:t2000_01_01_01+=365;
    case 19:t2000_01_01_01+=365;
    case 18:t2000_01_01_01+=365;
    case 17:t2000_01_01_01+=366;
    case 16:t2000_01_01_01+=365;
    case 15:t2000_01_01_01+=365;
    case 14:t2000_01_01_01+=365;
    case 13:t2000_01_01_01+=366;
    case 12:t2000_01_01_01+=365;
    case 11:t2000_01_01_01+=365;
    case 10:t2000_01_01_01+=365;
    case  9:t2000_01_01_01+=366;
    case  8:t2000_01_01_01+=365;
    case  7:t2000_01_01_01+=365;
    case  6:t2000_01_01_01+=365;
    case  5:t2000_01_01_01+=366;
    case  4:t2000_01_01_01+=365;
    case  3:t2000_01_01_01+=365;
    case  2:t2000_01_01_01+=365;
    case  1:t2000_01_01_01+=366;
    case  0:;
        // minus years = add if you interesting in anything from last century or use formala!!
    }
    //long it2000_01_01_01 = t2000_01_01_01;
    //double RestOfTheDay = t2000_01_01_01 - (double)it2000_01_01_01;
    return t2000_01_01_01// - RestOfTheDay 
            + KeplerDate - ((long double)(iYear*1000))
            -1; // epoch date is 1== 1 Jan - needs to adjust date.
}


///////////////////////////////////////////////////////////////////////////
// quick XML parser
///
char szSection[1024] =  {0};
#define XML_READ(XML_PARAM) if (CallXMLPars(szString, #XML_PARAM)) XML_PARAM = atof(pszQuo);
#define IF_XML_READ(XML_PARAM) if (CallXMLPars(szString, #XML_PARAM))
#define XML_BEGIN char *pszQuo;
#define XML_END ;
#define XML_SECTION(XML_SEC_NAME)     if (strcmp(szSection, #XML_SEC_NAME)==0)\
    {\
        pszQuo = strstr(szString, "value=\"");\
        if (pszQuo != NULL)\
        {\
            pszQuo += sizeof("value=\"") -1;

#define XML_SECTION_END         }\
    }\


int CallXMLPars(char *szString, char *XML_Params)
{
    char szFullComapre[1024] = {"name=\""};
    strcat(szFullComapre, XML_Params);
    strcat(szFullComapre, "\"");
    if (strstr(szString, szFullComapre) != NULL)   
        return 1;
    else
        return 0;

}
int iDayOfTheYearZeroBase(int iDay, int iMonth, int iYear)
{
	int iDays = iDay-1;
	switch(iMonth-1)
	{
	case 11:// november
		iDays+=30;
	case 10:// october
		iDays+=31;
	case 9:// september
		iDays+=30;
	case 8://august
		iDays+=31;
	case 7:// july
		iDays+=31;
	case 6:// june
		iDays+=30;
	case 5:// may
		iDays+=31;
	case 4:// april
		iDays+=30;
	case 3:// march
		iDays+=31;
	case 2:// february
		if ((iYear %4) ==0) // leap year
			iDays+=29;
		else
			iDays+=28;
	case 1:iDays+=31; // january
	case 0:iDays+=0;
		break;
	}
	return iDays;
}
double ConvertDateTimeToTLEEpoch(int iDay, int iMonth, int iYear, int iHour, int iMin, int iSec, int iMills)
{
    // An epoch of 98001.00000000 corresponds to 0000 UT on 1998 January 01—in other words, 
    // midnight between 1997 December 31 and 1998 January 01. 
    // An epoch of 98000.00000000 would actually correspond to the beginning of 1997 December 31—strange as that might seem. 
    // Note that the epoch day starts at UT midnight (not noon) and that all times are measured mean solar rather than sidereal time units.
    int mYear = iYear-2000;
    int mDays = iDayOfTheYearZeroBase(iDay, iMonth, iYear)+1;
	long mCurSec = iHour * 60*60;
    mCurSec += iMin *60;
    mCurSec += iSec;
	double dEpoch = mYear *1000.0 + mDays;
	dEpoch += (((double)mCurSec)+ ((double)iMills/1000.))/ (24.0*60.0*60.0);
    return dEpoch;
}

char szURLTraVisualFileName[3*_MAX_PATH];
char szURLTraVisualServer[3*_MAX_PATH];
char szTraVisualFileName[_MAX_PATH*3]={"travisual.xml"};
int UrlTraVisualPort=80;
BOOL VisualFileSet = FALSE;
BOOL ParsURL(char * URLServer, int *port, char* URL,  char * szParsingName)
{
    char sztemp[3*_MAX_PATH];
    *port=80;
    strcpy(URL,szParsingName);
    //strcpy(szTraVisualFileName, "@travisual.xml");
    strcpy(URLServer,URL);
    char *iFirst = strstr(URLServer,"http://");
    if (iFirst)
    {
        iFirst += 7;
        iFirst = strstr(iFirst,"/");
        if (iFirst)
        {
            *iFirst++=0;
            strcpy(URL,iFirst);
            iFirst = strstr(&URLServer[7],":");
            if (iFirst)  // found :8080 port number
            {
                *port = atoi(iFirst+1);
                *iFirst =  0;
            }
            strcpy(sztemp,&URLServer[7]);
            strcpy(URLServer,sztemp);
        }
        else
        {
            printf(" URL:%s is wrong", szParsingName);
            exit(3);
        }
        return TRUE;
    }
    return FALSE;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   TRA.XMl processing => common data
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int iGr = 0;

void ParamCommon(char *szString)
{
    XML_BEGIN;
    XML_SECTION(TraInfo);
    
        IF_XML_READ(dMinFromNow) 
        {
            dMinFromNow = atof(pszQuo);
        }
        IF_XML_READ(GRSTNLat) 
        {
            GrLat[iGr] = atof(pszQuo);
        }
        IF_XML_READ(GRSTNLong) 
        {
            GrLong[iGr++] = atof(pszQuo);
        }
        IF_XML_READ(dStartJD) 
        {
			TotalDays = 0;
            dStartJD = atof(pszQuo);    // format: <TRA:setting name="dStartJD0" value="2455625.1696833" />
			if (dStartJD <=0) // negativge value set current date munis amount of the minutes (negative -3 mean = total time == 3 minutes starting from curent time)
			{                 // i.e. -60 = total 60 minutes starting from -57 min, and 3 minutes in a future
            	SYSTEMTIME MyTime;
                TIME_ZONE_INFORMATION tmzone;
                //SYSTEMTIME ThatTime;

                int iYear;
				int iDays;
				int iCurSec;
                int Iret = GetTimeZoneInformation(&tmzone); 
                GetSystemTime(&MyTime);
                //double dEpoch = ConvertDateTimeToTLEEpoch(1, 1, 2013, 0, 0, 0, 0);
                //dStartJD = ConverEpochDate2JulianDay(dEpoch);
                double dEpoch = ConvertDateTimeToTLEEpoch(MyTime.wDay, MyTime.wMonth, MyTime.wYear, MyTime.wHour, MyTime.wMinute, MyTime.wSecond, MyTime.wMilliseconds);
                dStartJD = ConverEpochDate2JulianDay(dEpoch);
				
				TotalDays = (60.0*atof(pszQuo))/(24.0*60.0*60);
                // negative value !!! + 3 minutes
                dStartJD +=TotalDays + (60.0*dMinFromNow)/(24.0*60.0*60);
			}
			else
			{
				int iYear = dStartJD/1000;
				if ((iYear > 0) && (iYear < 24)) // this is a <TRA:setting name="dStartJD" value="11291.79166666" />
				{
					dStartJD = ConverEpochDate2JulianDay(dStartJD);
				}
				else
				{                     // 012345678901234567890
					if (dStartJD<=31) // DD/MM/YY HH:MM:SS:MLS format
					{
						int iDD = atoi(&pszQuo[0]);
						int iMO = atoi(&pszQuo[3]);
						int iYY = atoi(&pszQuo[6]);
						int iHH = atoi(&pszQuo[9]);
						int iMM = atoi(&pszQuo[12]);
						int iSS = atoi(&pszQuo[15]);
						int iMLS = atoi(&pszQuo[18]);
                        double dEpoch = ConvertDateTimeToTLEEpoch(iDD, iMO, iYY+2000, iHH, iMM, iSS,iMLS);
						//int iDays = iDayOfTheYearZeroBase(iDD, iMO, iYY+2000);
						//int iCurSec = iHH * 60*60;
						//iCurSec += iMM *60;
				        //iCurSec += iSS;
				        //dStartJD = iYY *1000.0 + iDays;
				        //dStartJD += (((double)iCurSec)+ ((double)iMLS/1000.))/ (24.0*60.0*60.0);
				        //dStartJD = ConverEpochDate2JulianDay(dStartJD);
                        dStartJD = ConverEpochDate2JulianDay(dEpoch);
					}
                    else // in a normal format
                    {
                    }
				}
			}
        }
        XML_READ(TimeSl);
        XML_READ(Gbig);
        XML_READ(IterPerSec);
        XML_READ(StartLandingIteraPerSec);
        IF_XML_READ(TotalDays) 
        {
			if (TotalDays<0)
				TotalDays = -TotalDays;
			else
				TotalDays = atof(pszQuo);
			iTotalSec = (int)(TotalDays * 24.0 * 60.0 * 60.0);

		}
        IF_XML_READ(EarthCurTime)  
        {
            EarthCurTime  = atof(pszQuo);
            // TBD
            EarthCurTimeS = EarthCurTime;
        }
        IF_XML_READ(EarthSmAxAU)
        {
            EarthSmAxAU = atof(pszQuo);
            AUcalc = EarthSmAx / EarthSmAxAU;
        }
        IF_XML_READ(TRAVisual)
        {
            strcpy(szTraVisualFileName,pszQuo);
            char * iQuot = strstr(szTraVisualFileName,"\"");
            if (iQuot)
                *iQuot=0;
            VisualFileSet = TRUE;
            if (ParsURL(szURLTraVisualServer, &UrlTraVisualPort, szURLTraVisualFileName,  szTraVisualFileName))
            {
                strcpy(szTraVisualFileName, "@travisual.xml");
            }
        }
#ifdef _DO_VISUALIZATION
        IF_XML_READ(RGBImageW)
        {
            bRGBImageW = atoi(pszQuo);
        }
        IF_XML_READ(RGBImageH)
        {
            bRGBImageH = atoi(pszQuo);
        }
        IF_XML_READ(RGBView)
        {
            iProfile = atoi(pszQuo);
        }
        IF_XML_READ(RGBMaxPictures)
        {
            iMaxSeq = atoi(pszQuo);
        }
        IF_XML_READ(RGBSecPerPictures)
        {
            iMaxCounter = atoi(pszQuo);
        }
        IF_XML_READ(RGBScale)
        {
            dRGBScale = atof(pszQuo);
        }

        IF_XML_READ(RGBReferenceBody)
        {
            RGBReferenceBody = atoi(pszQuo);
        }
#endif
        IF_XML_READ(EngineToOptimize)
        {
            EngineToOptimize = atoi(pszQuo);
        }
    
        IF_XML_READ(TrajectoryOptimizationType)
        {
            TrajectoryOptimizationType = atoi(pszQuo);
        }
        IF_XML_READ(LastEngine)
        {
            LastEngine = atoi(pszQuo);
        }
        IF_XML_READ(MaxOptim)
        {
            MaxOptim = atoi(pszQuo);
        }
        IF_XML_READ(StartOptim)
        {
            StartOptim = atoi(pszQuo);
            iOptimizationStep = StartOptim;
        }

            
    XML_SECTION_END;
    // processing pulsar coordinates and parameters
    XML_SECTION(pulsars);
        IF_XML_READ(N)
        {
            Pulsars[nPulsars].N = atoi(pszQuo);
        }
        IF_XML_READ(Name)
        {
            strcpy(Pulsars[nPulsars].Name,pszQuo);
        }
        IF_XML_READ(ELONG)
        {
            Pulsars[nPulsars].ELONG = atoi(pszQuo);
        }
        IF_XML_READ(ELAT)
        {
            Pulsars[nPulsars].ELAT = atoi(pszQuo);
        }
        IF_XML_READ(P0)
        {
            Pulsars[nPulsars].P0 = atoi(pszQuo);
        }
        IF_XML_READ(S400mJy)
        {
            Pulsars[nPulsars].S400mJy = atoi(pszQuo);
            if (++nPulsars >= NPULSARS)
                nPulsars = NPULSARS-1;
        }
    XML_SECTION_END;
    XML_END;
}
// rotate any point around vector on a angle 
void RotationAngleVect(long double &X, long double &Y, long double &Z, long double Angle, long double VectX, long double VectY, long double VectZ)
{
	// see http://en.wikipedia.org/wiki/Rotation_matrix

	// first make sure that Vector around rotation is a unit vector
	long double dLength = sqrt(VectX*VectX + VectY*VectY + VectZ*VectZ);
	VectX /= dLength;
	VectY /= dLength;
	VectZ /= dLength;
	long double RotX = (cos(Angle) + VectX*VectX*(1-cos(Angle)))*X +
		          (VectX*VectY*(1-cos(Angle)) - VectZ*sin(Angle)) * Y +
				  (VectX*VectZ*(1-cos(Angle)) + VectY*sin(Angle)) *Z;

	long double RotY = (VectY*VectX*(1-cos(Angle)) + VectZ*sin(Angle))*X +
		          (cos(Angle) + VectY*VectY*(1-cos(Angle)))*Y +
				  (VectY*VectZ*(1-cos(Angle)) - VectX*sin(Angle))*Z;

	long double RotZ = (VectZ*VectX*(1-cos(Angle)) - VectY*sin(Angle))*X+
		          (VectZ*VectY*(1-cos(Angle)) + VectX*sin(Angle))*Y+
				  (cos(Angle) + VectZ*VectZ*(1-cos(Angle)))*Z;
	X = RotX;
	Y = RotY;
	Z = RotZ;
}
// calculates angle btw vectors
long double AngleBtw(long double X1,long double Y1,long double Z1,long double X2,long double Y2,long double Z2)
{
	long double v1Len = sqrt(X1*X1 + Y1*Y1 + Z1*Z1);
	long double v2Len = sqrt(X2*X2 + Y2*Y2 + Z2*Z2);
	long double Cosv1v2 = (X1*X2 + Y1*Y2 + Z1*Z2)/v1Len/v2Len;
	long double Angle = acos(Cosv1v2);
	return (Angle * 180 /M_PI);
}
// calculates ortogonal vector
void Ort(long double &Xpr, long double &Ypr, long double &Zpr, long double u1, long double u2, long double u3, long double v1, long double v2, long double v3)
{
    Xpr = u2*v3 - u3*v2;
    Ypr = u3*v1 - u1*v3;
    Zpr = u1*v2-u2*v1;
    double prMod = sqrt(Xpr*Xpr + Ypr*Ypr +Zpr*Zpr);
    Xpr/=prMod;Ypr/=prMod;Zpr/=prMod;
}
void SUN_08 (int IYEAR,int IDAY,int IHOUR,int MIN,int ISEC,
	long double &GST,long double &SLONG,long double &SRASN,long double &SDEC)
{
	//C*******************************************************************
	//c
	//      SUBROUTINE SUN_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
	//C
	//C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
	//C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
	//C
	//C-------  INPUT PARAMETERS:
	//C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
	//C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
	//C
	//C-------  OUTPUT PARAMETERS:
	//C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
	//C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
	//C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
	//C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
	//C
	//C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
	//C
	//C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
	//C
    long double DJ,FDAY;
    long double RAD = 180.0/M_PI;//57.295779513;
	//C
    //IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
	if(IYEAR < 1901 || IYEAR>2099) 
		return;
    FDAY=(IHOUR*3600.0+MIN*60.0+ISEC)/86400.0;
    DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5+FDAY;
    long double T=DJ/36525.0;
    long double VL=fmod((long double)279.696678+(long double)0.9856473354*(long double)DJ,(long double)360.0);
    GST=fmod((long double)279.690983+(long double)0.9856473354*(long double)DJ+(long double)360.0*FDAY+(long double)180.0,(long double)360.0)/RAD;
    long double G=fmod((long double)358.475845+(long double)0.985600267*(long double)DJ,(long double)360.0)/RAD;
    SLONG=(VL+(1.91946-0.004789*T)*sin(G)+0.020094*sin(2.0*G))/RAD;
    //IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
	if(SLONG > (M_PI*2))//6.2831853) 
		SLONG=SLONG-(M_PI*2);//6.2831853;
    //IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
	if (SLONG < 0.0) 
		SLONG=SLONG+(M_PI*2);//6.2831853;
    long double OBLIQ=(23.45229-0.0130125*T)/RAD;
    long double SOB=sin(OBLIQ);
    long double SLP=SLONG-9.924E-5;
	//C
	//C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
	//C   EARTH'S ORBITAL MOTION
	//C
    long double SIND=SOB*sin(SLP);
    long double COSD=sqrt((long double)1.-SIND*SIND);
    long double SC=SIND/COSD;
    SDEC=atan(SC);
    SRASN=M_PI-atan2(cos(OBLIQ)/SOB*SC,-cos(SLP)/COSD);
    //RETURN
    //  END
}
/*
FUNCTION ACTAN(SINX,COSX)
COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
ACTAN=0.
IF (COSX.EQ.0. ) GO TO 5
IF (COSX.GT.0. ) GO TO 1
ACTAN=PI
GO TO 7
1 IF (SINX.EQ.0. ) GO TO 8
IF (SINX.GT.0. ) GO TO 7
ACTAN=TWOPI
GO TO 7
5 IF (SINX.EQ.0. ) GO TO 8
IF (SINX.GT.0. ) GO TO 6
ACTAN=X3PIO2
GO TO 8
6 ACTAN=PIO2
GO TO 8
7 TEMP=SINX/COSX
ACTAN=ACTAN+ATAN(TEMP)
8 RETURN
END
*/
/*
FUNCTION FMOD2P(X)
COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
FMOD2P=X
I=FMOD2P/TWOPI
FMOD2P=FMOD2P-I*TWOPI
IF(FMOD2P.LT.0) FMOD2P=FMOD2P+TWOPI
RETURN
END
*/
// The function subroutine THETAG is passed the epoch time exactly as it appears on the input element cards.
// The routine converts this time to days since 1950 Jan 0.0 UTC, stores this in the COMMON E1,
// and returns the right ascension of Greenwich at epoch (in radians).
long double GreenwichAscension(long double EP)
{
	//FUNCTION THETAG(EP)
	//COMMON /E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,BSTAR,
	//1 X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
	//DOUBLE PRECISION EPOCH,D,THETA,TWOPI,YR,TEMP,EP,DS50
	long double EPOCH,D,THETA,TWOPI,YR,TEMP,DS50;
	TWOPI=2.0*M_PI;//6.28318530717959D0
	YR=(EP+2.e-7)*1.e-3;
	int JY=YR;
	YR=JY;
	D=EP-YR*1.e3;
	//IF(JY.LT.10) JY=JY+80
	if(JY < 20) 
		JY=JY+80;
	int N=(JY-69)/4;
	//IF(JY.LT.70) N=(JY-72)/4
	if (JY < 70) 
		N=(JY-72)/4;
	DS50=7305.0 + 365.0*(JY-70) +N + D;
	THETA=1.72944494 + 6.3003880987*DS50;
	TEMP=THETA/TWOPI;
	int I=TEMP;
	TEMP=I;
	long double THETAG=THETA-TEMP*TWOPI;
	//IF(THETAG.LT.0.D0) THETAG=THETAG+TWOPI
	if(THETAG < 0.0) 
		THETAG=THETAG+TWOPI;
	return THETAG;
	//RETURN
	//END
}
// in FORTARN variable with a name begining with I is an integer variable 
// The function subroutine FMOD2P is
// passed an angle in radians and returns the angle in radians within the range of 0 to 2 PI.
long double FMOD2P(long double X)
{
	// COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
	//double DE2RA,PI,PIO2,TWOPI,X3PIO2;
	long double FMOD2P;
	int I;
	FMOD2P=X;
	//I=FMOD2P/TWOPI;
	I=FMOD2P/(2.0*M_PI);
	//FMOD2P=FMOD2P-I*TWOPI;
	FMOD2P=FMOD2P-I*(2.0*M_PI);
	//IF(FMOD2P.LT.0) FMOD2P=FMOD2P+TWOPI
	if (FMOD2P < 0.0) 
		FMOD2P=FMOD2P+2.0*M_PI;
	return FMOD2P;
//RETURN
//END
}
// strange function :The function subroutine ACTAN is passed the values of sine and cosine in that order and
// it returns the angle in radians within the range of 0 to 2 pi
long double ACTAN(long double SINX,long double COSX)
{
	long double ACTAN,TEMP;
	//COMMON/C2/DE2RA,PI,PIO2,TWOPI,X3PIO2
	long double X3PIO2 = 3.0*M_PI/2.0;
	long double PIO2 = M_PI/2.0;
	ACTAN=0.0;
	if (COSX == 0.0) goto M_5;
	if (COSX > 0.0) goto M_1;
	ACTAN=M_PI;
	goto M_7;
M_1:
	if (SINX == 0.0) goto M_8;
	if (SINX > 0.0) goto M_7;
	ACTAN=2.0*M_PI;
	goto M_7;
M_5:
	if (SINX == 0.0) goto M_8;
	if (SINX > 0.0) goto M_6;
	ACTAN=X3PIO2;
	goto M_8;
M_6:
	ACTAN=PIO2;
	goto M_8;
M_7:
	TEMP=SINX/COSX;
	ACTAN=ACTAN+atan(TEMP);
M_8:
	return ACTAN;
	//RETURN
	//END
}
// original and less acurate
void SGP(long double TS, long double XNDT2O,long double XNDD6O,/*double IEXP,*/long double nuBSTAR,/*double IBEXP,*/long double XINCL, long double XNODEO,long double EO, long double OMEGAO, long double XMO, long double XNO, 
	long double &X,long double &Y,long double &Z,long double &XDOT,long double &YDOT,long double &ZDOT)
{
	//* SGP 31 OCT 80
    //SUBROUTINE SGP(IFLAG,TSINCE)
	long double TSINCE = TS;
	//COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,BSTAR,
	//1 X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
	//COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
	//1 XJ3,XKE,XKMPER,XMNPDA,AE
	//DATA DE2RA,E6A,PI,PIO2,QO,SO,TOTHRD,TWOPI,X3PIO2,XJ2,XJ3,
	//1 XJ4,XKE,XKMPER,XMNPDA,AE/.174532925E-1,1.E-6,
	//2 3.14159265,1.57079633,120.0,78.0,.66666667,
	//4 6.2831853,4.71238898,1.082616E-3,-.253881E-5,
	//5 -1.65597E-6,.743669161E-1,6378.135,1440.,1./
	//double DE2RA= .174532925E-1;
	// double E6A = 1.E-6;
	//double PI= 3.14159265;
	//double PIO2 = 1.57079633;
	//double QO= 120.0;
	//double SO= 78.0;
	//double TOTHRD= .66666667;
	//double TWOPI = 6.2831853;
	//double X3PIO2 = 4.71238898;
	long double XJ2 = 1.082616E-3;
	long double XJ3 = -.253881E-5;
	//double XJ4 = -1.65597E-6;
	long double XKE = BIG_XKE;//.743669161E-1;
	//double XKMPER = 6378.135;
	//double XMNPDA= 1440.;
	long double AE =1.;

	int IFLAG;
    //double XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,XNDD6O,BSTAR,
    ///*X,Y,Z,XDOT,YDOT,ZDOT,*/EPOCH,DS50;
    //double CK2,CK4,E6A,QOMS2T,S,TOTHRD,
    //XJ3,XKE,XKMPER,XMNPDA,AE;
    //double EPOCH, DS50;
	long double C1,C2,C3,C4;
	long double COSIO,SINIO;
	long double A1,D1,AO,PO,QO,XLO,D1O,D2O,D3O,D4O,PO2NO,OMGDT,XNODOT,C5,C6,E,A,P,XNODES,OMGAS,XLS,AXNSL,AYNSL,XL,U;
	int ITEM3;
	long double EO1,TEM5;
	long double SINEO1,COSEO1,TEM2,ECOSE,ESINE,EL2,PL,PL2,R,RDOT,RVDOT,TEMP,SINU,COSU,SU,SIN2U,COS2U,RK,UK,XNODEK,XINCK,SINUK,COSUK,SINNOK,COSNOK;
	long double SINIK,COSIK,XMX,XMY,UX,UY,UZ,VX,VY,VZ;
	//CK2=.5*XJ2*AE**2
	long double CK2=.5*XJ2*AE*AE;
    //if (IFLAG == 0) goto M_19;
    // * INITIALIZATION
    C1= CK2*1.5;
    C2= CK2/4.0;
    C3= CK2/2.0;
    C4= XJ3*(AE*AE*AE)/(4.0*CK2);
    COSIO=cos(XINCL);
    SINIO=sin(XINCL);
    //A1=(XKE/XNO)**TOTHRD;
	A1 = pow((XKE/XNO),(long double)2.0/(long double)3.0);
	//D1= C1/A1/A1*(3.*COSIO*COSIO-1.)/(1.-EO*EO)**1.5;
	D1= C1/A1/A1*(3.*COSIO*COSIO-1.)/pow((1.-EO*EO),(long double)1.5);
	AO=A1*(1.-1./3.*D1-D1*D1-134./81.*D1*D1*D1);
    PO=AO*(1.-EO*EO);
	QO=AO*(1.-EO);
	XLO=XMO+OMEGAO+XNODEO;
	D1O= C3 *SINIO*SINIO;
	D2O= C2 *(7.*COSIO*COSIO-1.);
	D3O=C1*COSIO;
	D4O=D3O*SINIO;
	PO2NO=XNO/(PO*PO);
	OMGDT=C1*PO2NO*(5.*COSIO*COSIO-1.);
	XNODOT=-2.*D3O*PO2NO;
	C5=.5*C4*SINIO*(3.+5.*COSIO)/(1.+COSIO);
	C6=C4*SINIO;
	IFLAG=0;
	//* UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
M_19: 
	A = XNO+(2.*XNDT2O+3.*XNDD6O*TSINCE)*TSINCE;
	//A=AO*(XNO/A)**TOTHRD;
	A=pow(AO*(XNO/A),(long double)2.0/(long double)3.0);
	E=1.0E-6;
	//IF(A.GT.QO) E=1.-QO/A
	if(A>QO) 
		E=1.-QO/A;
	P=A*(1.-E*E);
	XNODES= XNODEO+XNODOT*TSINCE;
	OMGAS= OMEGAO+OMGDT*TSINCE;
	XLS=FMOD2P(XLO+(XNO+OMGDT+XNODOT+(XNDT2O+XNDD6O*TSINCE)*TSINCE)*TSINCE);
   // LONG PERIOD PERIODICS
M_7:
	AXNSL=E*cos(OMGAS);
	AYNSL=E*sin(OMGAS)-C6/P;
	XL=FMOD2P(XLS-C5/P*AXNSL);
	//* SOLVE KEPLERS EQUATION
	U=FMOD2P(XL-XNODES);
	ITEM3=0;
	EO1=U;
	TEM5=1.0;
M_20:
	SINEO1=sin(EO1);
	COSEO1=cos(EO1);
	//IF(ABS(TEM5).LT.E6A) GO TO 30
	if(abs(TEM5)<1E-6) 
		goto M_30;
	//IF(ITEM3.GE.10) GO TO 30
	if(ITEM3 > 30) 
		goto M_30;
	ITEM3=ITEM3+1;
	TEM5=1.-COSEO1*AXNSL-SINEO1*AYNSL;
	TEM5=(U-AYNSL*COSEO1+AXNSL*SINEO1-EO1)/TEM5;
	TEM2=abs(TEM5);
	//IF(TEM2.GT.1.) TEM5=TEM2/TEM5
	if(TEM2>1.0) 
		TEM5=TEM2/TEM5;
	EO1=EO1+TEM5;

	goto M_20;
	//* SHORT PERIOD PRELIMINARY QUANTITIES
M_30:
	ECOSE=AXNSL*COSEO1+AYNSL*SINEO1;
	ESINE=AXNSL*SINEO1-AYNSL*COSEO1;
	EL2=AXNSL*AXNSL+AYNSL*AYNSL;
	PL=A*(1.-EL2);
	PL2=PL*PL;
	R=A*(1.-ECOSE);
	RDOT=XKE*sqrt(A)/R*ESINE;
	RVDOT=XKE*sqrt(PL)/R;
	TEMP=ESINE/(1.+sqrt(1.-EL2));
	SINU=A/R*(SINEO1-AYNSL-AXNSL*TEMP);
	COSU=A/R*(COSEO1-AXNSL+AYNSL*TEMP);
	//SU=ACTAN(SINU,COSU);
	SU=ACTAN(SINU,COSU);
	//* UPDATE FOR SHORT PERIODICS
	SIN2U=(COSU+COSU)*SINU;
	COS2U=1.-2.*SINU*SINU;
	RK=R+D1O/PL*COS2U;
	UK=SU-D2O/PL2*SIN2U;
	XNODEK=XNODES+D3O*SIN2U/PL2;
	XINCK =XINCL+D4O/PL2*COS2U;
	//* ORIENTATION VECTORS
M_8:
	SINUK=sin(UK);
	COSUK=cos(UK);
	SINNOK=sin(XNODEK);
	COSNOK=cos(XNODEK);
	SINIK=sin(XINCK);
	COSIK=cos(XINCK);
	XMX=-SINNOK*COSIK;
	XMY=COSNOK*COSIK;
	UX=XMX*SINUK+COSNOK*COSUK;
	UY=XMY*SINUK+SINNOK*COSUK;
	UZ=SINIK*SINUK;
	VX=XMX*COSUK-COSNOK*SINUK;
	VY=XMY*COSUK-SINNOK*SINUK;
	VZ=SINIK*COSUK;
	//* POSITION AND VELOCITY
	X=RK*UX;
	Y=RK*UY;
	Z=RK*UZ;
	XDOT=RDOT*UX;
	YDOT=RDOT*UY;
	ZDOT=RDOT*UZ;
	XDOT=RVDOT*VX+XDOT;
	YDOT=RVDOT*VY+YDOT;
	ZDOT=RVDOT*VZ+ZDOT;
//RETURN
//END
}
// original, more acurate, takes account about drag == according Aksenov needs to use exact formulas from SGP4 to correctly process
// data. Backwards formulas should applied before everything 
// Read comments: (do) RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP) FROM INPUT ELEMENTS
void SGP4(long double TSINCE,/*double EPOCH,*/ 
          long double nuXNDT2O, 
          long double nuXNDD6O,/*IEXP,*/
          long double BSTAR,/*IBEXP,*/
	      long double XINCL,  // the “mean” inclination at epoch
          long double XNODEO, //the “mean” longitude of ascending node at epoch
          long double EO,     // the “mean” eccentricity at epoch
          long double OMEGAO, // the “mean” argument of perigee at epoch
          long double XMO,   // (M0) the “mean” mean anomaly at epoch
          long double XNO,   // (k0) the SGP type “mean” mean motion at epoch
		  long double &X,long double &Y,long double &Z,long double &XDOT,long double &YDOT,long double &ZDOT)
{
//COMMON/E1/XMO,XNODEO,OMEGAO,EO,XINCL,XNO,XNDT2O,
//1 XNDD6O,BSTAR,X,Y,Z,XDOT,YDOT,ZDOT,EPOCH,DS50
//COMMON/C1/CK2,CK4,E6A,QOMS2T,S,TOTHRD,
//1 XJ3,XKE,XKMPER,XMNPDA,AE
	//long double /*EPOCH,*/ DS50;
	long double XKE = BIG_XKE;//.743669161E-1;
	long double XKMPER = 6378.1350;
	//IF (IFLAG .EQ. 0) GO TO 100
	//* RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
	//* FROM INPUT ELEMENTS
	//A1=(XKE/XNO)**TOTHRD;
	long double XJ2 = 1.082616E-3; //the second gravitational zonal harmonic of the Earth
	long double XJ3 = -.253881E-5; // the third gravitational zonal harmonic of the Earth
	long double XJ4 = -1.65597E-6; // the fourth gravitational zonal harmonic of the Earth
	long double AE = 1.0;          // the equatorial radius of the Earth - actualy it is not true it is one == everything measuared in that radiuses
	long double QO =120.0;         // parameter for the SGP4/SGP8 density function
	long double SO = 78.0;         // parameter for the SGP4/SGP8 density function

	long double CK2=.5*XJ2*AE*AE;
	//CK4=-.375*XJ4*AE**4
	long double CK4=-.375*XJ4*AE*AE*AE*AE;
	long double QOMS2T=pow(((QO-SO)*AE/XKMPER),(long double)4.0);
	long double S=AE*(1.+SO/XKMPER);

    //The original mean motion (n0") and semimajor axis (a0") are first recovered from the input elements by the equations:
    //
    // a1 = (k0 /n0) ** (2/3)
    //
    // b1 = 3/2 * k2/a1 * (3 * (cos(i0))**2 -1)/(1-e0**2)**(3/2)
    //
    // a0 = a1 * (1 - 1/3 b1 - b1**2  - 134/81 b1**3)
    //
    // b0 - 3/2 * k2/(a0**2) * (3 * (cos(i0))**2 - 1)/(1- e9**2)3/2
    //
    // n0" = n0/ (1+b0)
    // a0" = a0/(1-b0)
    //
	//A1=(XKE/XNO)**TOTHRD;
	long double A1=pow((XKE/XNO),(long double)2.0/(long double)3.0);
	long double COSIO=cos(XINCL);
	long double THETA2=COSIO*COSIO;
	long double X3THM1=3.*THETA2-1.;
	long double EOSQ=EO*EO;
	long double BETAO2=1.-EOSQ;
	long double BETAO=sqrt(BETAO2);
	long double DEL1=1.5*CK2*X3THM1/(A1*A1*BETAO*BETAO2);
	//AO=A1*(1.-DEL1*(.5*TOTHRD+DEL1*(1.+134./81.*DEL1)))
	long double AO=A1*(1.-DEL1*(.5*(2.0/3.0)+DEL1*(1.+134./81.*DEL1)));
	long double DELO=1.5*CK2*X3THM1/(AO*AO*BETAO*BETAO2);
    // ProbMeanMotion = XNO / (2*pi) * 1440.0
	long double XNODP=XNO/(1.+DELO);
    // and this is a semimajor axis
	long double AODP=AO/(1.-DELO);
	//* INITIALIZATION
	//* FOR PERIGEE LESS THAN 220 KILOMETERS, THE ISIMP FLAG IS SET AND
	//* THE EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN SQRT A AND
	//* QUADRATIC VARIATION IN MEAN ANOMALY. ALSO, THE C3 TERM, THE
	//* DELTA OMEGA TERM, AND THE DELTA M TERM ARE DROPPED.
	int ISIMP=0;
	//IF((AODP*(1.-EO)/AE) .LT. (220./XKMPER+AE)) ISIMP=1
	if((AODP*(1.-EO)/AE) < (220./XKMPER+AE)) 
		ISIMP=1;
	//* FOR PERIGEE BELOW 156 KM, THE VALUES OF
	//* S AND QOMS2T ARE ALTERED
	long double S4=S;
	long double QOMS24=QOMS2T;
	long double PERIGE=(AODP*(1.-EO)-AE)*XKMPER;
	//IF(PERIGE .GE. 156.) GO TO 10
	if (PERIGE >= 156.) 
		goto M_10;
	S4=PERIGE-78.;
	//IF(PERIGE .GT. 98.) GO TO 9
	if (PERIGE >=98.) 
		goto M_9;
	S4=20.;
M_9:
	//QOMS24=((120.-S4)*AE/XKMPER)**4;
	QOMS24=pow(((120.-S4)*AE/XKMPER),4);
	S4=S4/XKMPER+AE;
M_16:
M_10:
	long double PINVSQ=1./(AODP*AODP*BETAO2*BETAO2);
	long double TSI=1./(AODP-S4);
	long double ETA=AODP*EO*TSI;
	long double ETASQ=ETA*ETA;
	long double EETA=EO*ETA;
	long double PSISQ=abs(1.-ETASQ);
	//COEF=QOMS24*TSI**4
	long double COEF=QOMS24*TSI*TSI*TSI*TSI;
	//COEF1=COEF/PSISQ**3.5;
	long double COEF1=COEF/pow(PSISQ,(long double)3.5);

	long double C2=COEF1*XNODP*(AODP*(1.+1.5*ETASQ+EETA*(4.+ETASQ))+.75* CK2*TSI/PSISQ*X3THM1*(8.+3.*ETASQ*(8.+ETASQ)));
	long double C1=BSTAR*C2;
	long double SINIO=sin(XINCL);
	//A3OVK2=-XJ3/CK2*AE**3;
	long double A3OVK2=-XJ3/CK2*(AE*AE*AE);
	long double C3=COEF*TSI*A3OVK2*XNODP*AE*SINIO/EO;
	long double X1MTH2=1.-THETA2;
	long double C4=2.*XNODP*COEF1*AODP*BETAO2*(ETA*(2.+.5*ETASQ)+EO*(.5+2.*ETASQ)-2.*CK2*TSI/
			(AODP*PSISQ)*(-3.*X3THM1*(1.-2.*EETA+ETASQ* (1.5-.5*EETA))+.75*X1MTH2*(2.*ETASQ-EETA* (1.+ETASQ))*cos(2.*OMEGAO)));
	long double C5=2.*COEF1*AODP*BETAO2*(1.+2.75*(ETASQ+EETA)+EETA*ETASQ);
	long double THETA4=THETA2*THETA2;
	long double TEMP1=3.*CK2*PINVSQ*XNODP;
	long double TEMP2=TEMP1*CK2*PINVSQ;
	long double TEMP3=1.25*CK4*PINVSQ*PINVSQ*XNODP;
	long double XMDOT=XNODP+.5*TEMP1*BETAO*X3THM1+.0625*TEMP2*BETAO* (13.-78.*THETA2+137.*THETA4);
	long double X1M5TH=1.-5.*THETA2;
	long double OMGDOT=-.5*TEMP1*X1M5TH+.0625*TEMP2*(7.-114.*THETA2+ 395.*THETA4)+TEMP3*(3.-36.*THETA2+49.*THETA4);
	long double XHDOT1=-TEMP1*COSIO;
	long double XNODOT=XHDOT1+(.5*TEMP2*(4.-19.*THETA2)+2.*TEMP3*(3.-7.*THETA2))*COSIO;
	long double OMGCOF=BSTAR*C3*cos(OMEGAO);
	long double XMCOF=-(2.0/3.0)*COEF*BSTAR*AE/EETA;
	long double XNODCF=3.5*BETAO2*XHDOT1*C1;
	long double T2COF=1.5*C1;
	long double XLCOF=.125*A3OVK2*SINIO*(3.+5.*COSIO)/(1.+COSIO);
	long double AYCOF=.25*A3OVK2*SINIO;
	//DELMO=(1.+ETA*COS(XMO))**3
	long double DELMO=pow((1.+ETA*cos(XMO)),3);
	long double SINMO=sin(XMO);
	long double X7THM1=7.*THETA2-1.;
	// IF(ISIMP .EQ. 1) GO TO 90
	if (ISIMP == 1) 
		goto M_90;
	long double C1SQ=C1*C1;
	long double D2=4.*AODP*TSI*C1SQ;
	long double TEMP=D2*TSI*C1/3.;
	long double D3=(17.*AODP+S4)*TEMP;
	long double D4=.5*TEMP*AODP*TSI*(221.*AODP+31.*S4)*C1;
M_17:
	long double T3COF=D2+2.*C1SQ;
	long double T4COF=.25*(3.*D3+C1*(12.*D2+10.*C1SQ));
	long double T5COF=.2*(3.*D4+12.*C1*D3+6.*D2*D2+15.*C1SQ*(2.*D2+C1SQ));
M_90:
	//IFLAG=0;
	//* UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
M_100:
	long double XMDF=XMO+XMDOT*TSINCE;
	long double OMGADF=OMEGAO+OMGDOT*TSINCE;
	long double XNODDF=XNODEO+XNODOT*TSINCE;
	long double OMEGA=OMGADF;
	long double XMP=XMDF;
	long double TSQ=TSINCE*TSINCE;
	long double XNODE=XNODDF+XNODCF*TSQ;
	long double TEMPA=1.-C1*TSINCE;
	long double TEMPE=BSTAR*C4*TSINCE;
	long double TEMPL=T2COF*TSQ;
	// IF(ISIMP .EQ. 1) GO TO 110
	if (ISIMP == 1) 
		goto M_110;
	long double DELOMG=OMGCOF*TSINCE;
	// DELM=XMCOF*((1.+ETA*COS(XMDF))**3-DELMO)
	long double DELM=XMCOF*(pow((1.+ETA*cos(XMDF)),3)-DELMO);
	TEMP=DELOMG+DELM;
	XMP=XMDF+TEMP;
	OMEGA=OMGADF-TEMP;
	long double TCUBE=TSQ*TSINCE;
	long double TFOUR=TSINCE*TCUBE;
	TEMPA=TEMPA-D2*TSQ-D3*TCUBE-D4*TFOUR;
	TEMPE=TEMPE+BSTAR*C5*(sin(XMP)-SINMO);
	TEMPL=TEMPL+T3COF*TCUBE+TFOUR*(T4COF+TSINCE*T5COF);
M_110:
	//A=AODP*TEMPA**2;
	long double A=AODP*TEMPA*TEMPA;
	long double E=EO-TEMPE;
	long double XL=XMP+OMEGA+XNODE+XNODP*TEMPL;
	long double BETA=sqrt(1.-E*E);
	//XN=XKE/A**1.5
	long double XN=XKE/pow(A,(long double)1.5);
	//* LONG PERIOD PERIODICS
	long double AXN=E*cos(OMEGA);
	TEMP=1./(A*BETA*BETA);
	long double XLL=TEMP*XLCOF*AXN;
	long double AYNL=TEMP*AYCOF;
	long double XLT=XL+XLL;
	long double AYN=E*sin(OMEGA)+AYNL;
	//* SOLVE KEPLERS EQUATION
	long double CAPU=FMOD2P(XLT-XNODE);
M_18:
	TEMP2=CAPU;
	// DO 130 I=1,10
	long double SINEPW;
	long double COSEPW;
	long double TEMP4;
	long double TEMP5;
	long double TEMP6;
	long double EPW;
	for (int i =1; i <=30; i++)
	{
		SINEPW=sin(TEMP2);
		COSEPW=cos(TEMP2);
		TEMP3=AXN*SINEPW;
		TEMP4=AYN*COSEPW;
		TEMP5=AXN*COSEPW;
		TEMP6=AYN*SINEPW;
		EPW=(CAPU-TEMP4+TEMP3-TEMP2)/(1.-TEMP5-TEMP6)+TEMP2;
		//IF(ABS(EPW-TEMP2) .LE. E6A) GO TO 140
		if(abs(EPW-TEMP2) <= 1.e-15) 
			goto M_140;
M_130:	TEMP2=EPW;
	}

	//* SHORT PERIOD PRELIMINARY QUANTITIES
M_140:
	long double ECOSE=TEMP5+TEMP6;
	long double ESINE=TEMP3-TEMP4;
	long double ELSQ=AXN*AXN+AYN*AYN;
	TEMP=1.-ELSQ;
	long double PL=A*TEMP;
	long double R=A*(1.-ECOSE);
	TEMP1=1./R;
	long double RDOT=XKE*sqrt(A)*ESINE*TEMP1;
	long double RFDOT=XKE*sqrt(PL)*TEMP1;
	TEMP2=A*TEMP1;
	long double BETAL=sqrt(TEMP);
	TEMP3=1./(1.+BETAL);
	long double COSU=TEMP2*(COSEPW-AXN+AYN*ESINE*TEMP3);
	long double SINU=TEMP2*(SINEPW-AYN-AXN*ESINE*TEMP3);
	long double U=ACTAN(SINU,COSU);
	long double SIN2U=2.*SINU*COSU;
	long double COS2U=2.*COSU*COSU-1.;
	TEMP=1./PL;
	TEMP1=CK2*TEMP;
	TEMP2=TEMP1*TEMP;
	// * UPDATE FOR SHORT PERIODICS
	long double RK=R*(1.-1.5*TEMP2*BETAL*X3THM1)+.5*TEMP1*X1MTH2*COS2U;
	long double UK=U-.25*TEMP2*X7THM1*SIN2U;
	long double XNODEK=XNODE+1.5*TEMP2*COSIO*SIN2U;
	long double XINCK=XINCL+1.5*TEMP2*COSIO*SINIO*COS2U;
	long double RDOTK=RDOT-XN*TEMP1*X1MTH2*SIN2U;
	long double RFDOTK=RFDOT+XN*TEMP1*(X1MTH2*COS2U+1.5*X3THM1);
	//* ORIENTATION VECTORS
	long double SINUK=sin(UK);
	long double COSUK=cos(UK);
M_19:
	long double SINIK=sin(XINCK);
	long double COSIK=cos(XINCK);
	long double SINNOK=sin(XNODEK);
	long double COSNOK=cos(XNODEK);
	long double XMX=-SINNOK*COSIK;
	long double XMY=COSNOK*COSIK;
	long double UX=XMX*SINUK+COSNOK*COSUK;
	long double UY=XMY*SINUK+SINNOK*COSUK;
	long double UZ=SINIK*SINUK;
	long double VX=XMX*COSUK-COSNOK*SINUK;
	long double VY=XMY*COSUK-SINNOK*SINUK;
	long double VZ=SINIK*COSUK;
	//* POSITION AND VELOCITY
	X=RK*UX;
	Y=RK*UY;
	Z=RK*UZ;
	XDOT=RDOTK*UX+RFDOTK*VX;
	YDOT=RDOTK*UY+RFDOTK*VY;
	ZDOT=RDOTK*UZ+RFDOTK*VZ;
//RETURN
//END
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   TRA.XML processing => all data about the MOON
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParamMoon(char *szString)
{
    XML_BEGIN;
    XML_SECTION(TraInfo);
        XML_READ(MoonX);
        XML_READ(MoonY);
        XML_READ(MoonZ);
        XML_READ(MoonVX);
        XML_READ(MoonVY);
        XML_READ(MoonVZ);

        XML_READ(MoonR);
        XML_READ(MoonRP);
        XML_READ(MoonRE);
        XML_READ(MoonM);
        XML_READ(GMMoon);
        XML_READ(GMEarthMoon);


        IF_XML_READ(MassRatioEarthToMoon)
        {
            MassRatioEarthToMoon = atof(pszQuo);
            //if (MassRatioSunToEarthPlusMoon)
            {
                double Temp = SunM / MassRatioSunToEarthPlusMoon;
                if (GMEarth != 0.0 && GMMoon != 0.0 && GMSun != 0.0)
                {
                    printf("\n calculation mass based on GMSun GMMoon and GMEarth");
                    printf("\n Was  GMSun = %f", GMSun);
                    printf("\n Was  GMEarth = %f", GMEarth);
                    printf("\n Was  GMMoon = %f", GMMoon);
                    printf("\n Was  SunM = %f", SunM);
                    SunM = GMSun /Gbig;
                    printf("\n Calc SunM = %f", SunM);                    
                    printf("\n Was  EarthM = %f", EarthM);
                    EarthM = GMEarth / Gbig;
                    printf("\n Calc EarthM = %f", EarthM);
                    printf("\n Was  MoonM = %f", MoonM);
                    MoonM = GMMoon / Gbig;
                    printf("\n Calc MoonM = %f", MoonM);
                }
                else
                {
                    printf("\n calculation mass based on mass ratio earth to moon");
                    printf("\n Was  MoonM = %f", MoonM);
                    MoonM = Temp / (1.0 + MassRatioEarthToMoon);
                    printf("\n Calc MoonM = %f", MoonM);
                    printf("\n Was  EarthM = %f", EarthM);
                    EarthM = MoonM * MassRatioEarthToMoon;
                    printf("\n Calc EarthM = %f", EarthM);
                }
            }
        }
    XML_SECTION_END;
    XML_END;
}
void MoonXYZCalc(double &flX, double &flY, double &flZ, double tsec)
{
    flX = (383.0*sin( 8399.685*tsec + 5.381)+
            31.5*sin(   70.990*tsec + 6.169)+
            10.6*sin(16728.377*tsec + 1.453)+
             6.2*sin( 1185.622*tsec + 0.481)+
             3.2*sin( 7143.070*tsec + 5.017)+
             2.3*sin(15613.745*tsec + 0.857)+
             0.8*sin( 8467.263*tsec + 1.010))*1.0E6;
    
    flY = (351.0*sin( 8399.687*tsec + 3.811)+
            28.9*sin(   70.997*tsec + 4.596)+
            13.7*sin( 8433.466*tsec + 4.766)+
             9.7*sin(16728.380*tsec + 6.165)+
             5.7*sin( 1185.667*tsec + 5.164)+
             2.9*sin( 7143.058*tsec + 0.300)+
             2.1*sin(15613.755*tsec + 5.565))*1.0E6;

    flZ = (153.2*sin( 8399.672*tsec + 3.807)+
            31.5*sin( 8433.464*tsec + 1.629)+
            12.5*sin(   70.996*tsec + 4.595)+
             4.2*sin(16728.364*tsec + 6.162)+
             2.5*sin( 1185.645*tsec + 5.167)+
             3.0*sin(  104.881*tsec + 2.555)+
             1.8*sin( 8399.116*tsec + 6.248))*1.0E6;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   TRA.XMl processing => all data abou SUN (now it is not really used)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParamSun(char *szString)
{
    XML_BEGIN;
    XML_SECTION(TraInfo);
        
        XML_READ(SunX);
        XML_READ(SunY);
        XML_READ(SunZ);
		XML_READ(SunR);
        XML_READ(SunM);
        IF_XML_READ(GMSun) 
        {
            GMSun = atof(pszQuo);
            printf("\n Was  SunM= %f", SunM);
            SunM = GMSun / Gbig;
            printf("\n Calc SunM= %f", SunM);
        }
    XML_SECTION_END;
    XML_END;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   TRA.XML processing => all earth data
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParamEarth(char *szString)
{
    XML_BEGIN;
    XML_SECTION(TraInfo);
        XML_READ(EarthX);
        XML_READ(EarthY);
        XML_READ(EarthZ);
        XML_READ(EarthVX);
        XML_READ(EarthVY);
        XML_READ(EarthVZ);
        XML_READ(EarthR);
        XML_READ(EarthRP);
        XML_READ(EarthRE);
        XML_READ(EarthM);
        XML_READ(AU);
        XML_READ(EarthTSolSec);
        XML_READ(EarthSmAx);
        double Temp = 0.0;
        IF_XML_READ(EarthTDays)  
        {
            
            EarthTDays  = atof(pszQuo);
            printf("\n calc EarthTDays=%f ", EarthTDays);
            EarthTSec = EarthTDays * 24.0*60.0*60.0;
            printf("\n calc EarthTSec =%f ", EarthTSec);
            Temp = EarthTSec / EarthTSolSec;

        }

        XML_READ(GMEarth);
        XML_READ(MassRatioSunToEarthPlusMoon);
        IF_XML_READ(EarthCalcKepler)  
        {
            // not used - left for verification in debug mode only            
            EarthCalcKepler  = atof(pszQuo);
            if (EarthCalcKepler == 1.0)
            {
                double tempCenterEarthMoonX;
                double tempCenterEarthMoonY;
                double tempCenterEarthMoonZ;
                double tempMoonX;
                double tempMoonY;
                double tempMoonZ;
                double tempMoonVX;
                double tempMoonVY;
                double tempMoonVZ;
                double tempCenterEarthMoonVX;
                double tempCenterEarthMoonVY;
                double tempCenterEarthMoonVZ;


                printf("\n GMSun = %f", GMSun);
                printf("\n GMEarth = %f", GMEarth);
                GMEarth = Gbig * EarthM; 
                printf("\n GMEarth = %f", GMEarth);
                printf("\n GMMoon = %f", GMMoon);
                GMMoon = Gbig * MoonM;
                printf("\n GMMoon = %f", GMMoon);

                printf("\n Was  A Earth= %f", EarthSmAx);
                printf("\n Was AU      = %f", AU);

                // adjust Moon + earth position at a point of a centre of gravity
                tempCenterEarthMoonX = EarthX;
                tempCenterEarthMoonY = EarthY;
                tempCenterEarthMoonZ = EarthZ;
                EarthX = tempCenterEarthMoonX + tempMoonX * (MoonM/(EarthM+MoonM));
                EarthY = tempCenterEarthMoonY + tempMoonY * (MoonM/(EarthM+MoonM));
                EarthZ = tempCenterEarthMoonZ + tempMoonZ * (MoonM/(EarthM+MoonM));
                MoonX =  EarthX - tempMoonX;
                MoonY =  EarthY - tempMoonY;
                MoonZ =  EarthZ - tempMoonZ;

                // adjust Moon + earth speed at a point of a centre of gravity
                //VM = sqrt( tempMoonVX* tempMoonVX +  tempMoonVY* tempMoonVY +  tempMoonVZ* tempMoonVZ);
                tempCenterEarthMoonVX = EarthVX;
                tempCenterEarthMoonVY = EarthVY;
                tempCenterEarthMoonVZ = EarthVZ;
                EarthVX = tempCenterEarthMoonVX - tempMoonVX * (MoonM/(EarthM+MoonM));
                EarthVY = tempCenterEarthMoonVY - tempMoonVY * (MoonM/(EarthM+MoonM));
                EarthVZ = tempCenterEarthMoonVZ - tempMoonVZ * (MoonM/(EarthM+MoonM));
                MoonVX = tempCenterEarthMoonVX + tempMoonVX * (EarthM/(EarthM+MoonM));
                MoonVY = tempCenterEarthMoonVY + tempMoonVY * (EarthM/(EarthM+MoonM));
                MoonVZ = tempCenterEarthMoonVZ + tempMoonVZ * (EarthM/(EarthM+MoonM));
            }
        }
    XML_SECTION_END;
    XML_END;
}
// just for convinience to proper orbit clculation
double EngCoeff = 1.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       TRA.XML processing => all data about satellite
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParamProb(char *szString)
{
    
    char szTemp[128] = {"0."};;
    XML_BEGIN;
    XML_SECTION(TraInfo);
    // pisition of the prob (active by firing engines) can be set by X,Y,Z   VX,VY,VZ
        IF_XML_READ(ProbX)
        {
            Sat.X[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbY)
        {
            Sat.Y[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbZ)
        {
            Sat.Z[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbVX)
        {
            Sat.VX[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbVY)
        {
            Sat.VY[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbVZ)
        {
            Sat.VZ[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbM)
        {
            Sat.M[Sat.Elem] = atof(pszQuo);
        }
        // or by 3 punch card
        IF_XML_READ(ProbKeplerLine1)
        {
            strcpy(Sat.Kepler1[Sat.Elem], pszQuo);
        }
        IF_XML_READ(ProbKeplerLine2)
        {
            strcpy(Sat.Kepler2[Sat.Elem], pszQuo);
        }
        IF_XML_READ(ProbKeplerLine3)
        {
            char szTempo[1024];
            int iYear;
            int iDays;
            double dflTemp;

            strcpy(Sat.Kepler3[Sat.Elem], pszQuo);
            //strcpy(Sat.Kepler1[0], szProbKeplerLine1);
            //strcpy(Sat.Kepler2[0], szProbKeplerLine2);
            //strcpy(Sat.Kepler3[0], szProbKeplerLine3);
			//0         1         2         3         4         5         6         7
			//01234567890123456789012345678901234567890123456789012345678901234567890
			//1 25544U 98067A   04236.56031392  .00020137  00000-0  16538-3 0  5135\
			//2 25544  51.6335 341.7760 0007976 126.2523 325.9359 15.70406856328903"
			// format readings from FORTRAN:
			// 702 FORMAT(18X,D14.8,1X,F10.8,2(1X,F6.5,I2),/,7X,2(1X,F8.4),1X,F7.7,2(1X,F8.4),1X,F11.8)
			//
		    // "04236.56031392" Element Set Epoch (UTC) D14.8
			//  04                   - year
			//    236.56031392       - day
            Sat.ProbEpoch[Sat.Elem] = atof(&Sat.Kepler2[Sat.Elem][18]);
            Sat.ProbEpochS[Sat.Elem] = Sat.ProbEpoch[Sat.Elem] =ConverEpochDate2JulianDay(Sat.ProbEpoch[Sat.Elem]);
			Sat.ProbEpochS[Sat.Elem] *=  60*60*24;
            // now for initial step asign emulation starting time:
			// if nothing is set then starting point is a last satellite epoch
            if (dStartJD == 0.0)
            {
                dStartJD = ConverEpochDate2JulianDay(Sat.ProbEpoch[Sat.Elem]);
            }
			//  "_.00020137"      1st Derivative of the Mean Motion with respect to Time F10.8
			//double ProbFirstDervMeanMotion; XNDT2O
			COPYKEPLER(Sat.ProbFirstDervMeanMotion[Sat.Elem],&Sat.Kepler2[Sat.Elem][33],10);
			// "_00000-0"        2nd Derivative of the Mean Motion with respect to Time (decimal point assumed) 1X,F6.5,I2
			//double ProbSecondDervmeanMotion;,XNDD6O,IEXP
			memset(szTempo, 0, sizeof(szTempo)); 
            szTempo[0] = Sat.Kepler2[Sat.Elem][44]; if (szTempo[0]==' ') szTempo[0] = '0';
			szTempo[1] = '.'; 
            memcpy(&szTempo[2], &Sat.Kepler2[Sat.Elem][45], 5);
            Sat.ProbSecondDervmeanMotion[Sat.Elem] = atof(szTempo);
			if (Sat.Kepler2[Sat.Elem][50] == '-')
				Sat.ProbSecondDervmeanMotion[Sat.Elem] *= pow(10., - (Sat.Kepler2[Sat.Elem][51] - '0'));
			else //if (Sat.Kepler2[Sat.Elem][50] == '+') //is it corect ??
				Sat.ProbSecondDervmeanMotion[Sat.Elem] *= pow(10., (Sat.Kepler2[Sat.Elem][51] - '0'));
			//// "_16538-3"        B* Drag Term 1X,F6.5,I2
			//double ProbDragterm; ,BSTAR,IBEXP
			memset(szTempo, 0, sizeof(szTempo)); 
            szTempo[0] = Sat.Kepler2[Sat.Elem][53]; if (szTempo[0]==' ') szTempo[0] = '0';
			szTempo[1] = '.'; 
            memcpy(&szTempo[2], &Sat.Kepler2[Sat.Elem][54], 5);
            Sat.ProbDragterm[Sat.Elem] = atof(szTempo);
			if (Sat.Kepler2[Sat.Elem][59] == '-')
				Sat.ProbDragterm[Sat.Elem] *= pow(10., - (Sat.Kepler2[Sat.Elem][60] - '0'));
                        else
                            Sat.ProbDragterm[Sat.Elem] *= pow(10., (Sat.Kepler2[Sat.Elem][60] - '0'));
			// "0"              Element Set Type
			Sat.ProbElementSetType[Sat.Elem] = Sat.Kepler2[Sat.Elem][62];
			// "_513"            Element Number
			// "5"              Checksum
			//                        The checksum is the sum of all of the character in the data line, modulo 10. 
			//                        In this formula, the following non-numeric characters are assigned the indicated values: 
			//                        Blanks, periods, letters, '+' signs -> 0
			//                        '-' signs -> 1
			// "2"             Line Number
			// "25544"         Object Identification Number
			// "_51.6335"       Orbit Inclination (degrees) 1X,F8.4
			COPYKEPLER(Sat.ProbIncl[Sat.Elem],&Sat.Kepler3[Sat.Elem][8],8);
            Sat.ProbIncl[Sat.Elem] = M_PI * Sat.ProbIncl[Sat.Elem]/180.0;
			// "341.7760"      Right Ascension of Ascending Node (degrees) 1X,F8.4
			COPYKEPLER(Sat.ProbAscNode[Sat.Elem],&Sat.Kepler3[Sat.Elem][17],8);
			Sat.ProbAscNode[Sat.Elem] = M_PI * Sat.ProbAscNode[Sat.Elem]/180.0;
			// "0007976"       Eccentricity (decimal point assumed) F7.7
			memset(szTempo, 0, sizeof(szTempo)); 
            szTempo[0] = '.'; 
            memcpy(&szTempo[1], &Sat.Kepler3[Sat.Elem][26], 7);
            Sat.ProbEcc[Sat.Elem] = atof(szTempo);
			// "126.2523"      Argument of Perigee (degrees) 1X,F8.4
			COPYKEPLER(Sat.ProbArgPer[Sat.Elem],&Sat.Kepler3[Sat.Elem][34],8);
			Sat.ProbArgPer[Sat.Elem] = M_PI * Sat.ProbArgPer[Sat.Elem]/180.0;
			// "325.9359"      Mean Anomaly (degrees) 1X,F8.4
			COPYKEPLER(Sat.ProbMeanAnom[Sat.Elem],&Sat.Kepler3[Sat.Elem][43],8);
			Sat.ProbMeanAnom[Sat.Elem] = M_PI * Sat.ProbMeanAnom[Sat.Elem]/180.0;
			// "15.70406856"    Mean Motion (revolutions/day)F11.8
			COPYKEPLER(Sat.ProbTPeriod[Sat.Elem],&Sat.Kepler3[Sat.Elem][52],11);
			Sat.ProbMeanMotion[Sat.Elem] = Sat.ProbTPeriod[Sat.Elem];
			
            printf("\n calc ProbTPeriod=%f ", Sat.ProbTPeriod[Sat.Elem]);
            Sat.ProbTDays[Sat.Elem] = 1.0/Sat.ProbTPeriod[Sat.Elem];
            printf("\n calc ProbTDays=%f ", Sat.ProbTDays[Sat.Elem]);
            Sat.ProbTSec[Sat.Elem] = Sat.ProbTDays[Sat.Elem] * 24. * 60. * 60.;
            printf("\n calc ProbTSec=%f ", Sat.ProbTSec[Sat.Elem]);
			// "328903"        Revolution Number at Epoch
			COPYKEPLER(Sat.ProbRevAtEpoch[Sat.Elem],&Sat.Kepler3[Sat.Elem][63],6);
            Sat.Elem++;  // next satellite
        }
        // finding the tag "UseJPLxyz" will force to retrive all solar sistem parameters from the JPL data
        IF_XML_READ(UseJPLxyz)
        {
            long double tempProbX;
            long double tempProbY;
            long double tempProbZ;
            long double tempProbVX;
            long double tempProbVY;
            long double tempProbVZ;

			long double tProbX=0.;
            long double tProbY=0.;
            long double tProbZ=0.;
            long double tProbVX=0.;
            long double tProbVY=0.;
            long double tProbVZ=0.;
            // V of the earth calulated = sqrt(GM * (2/R - 1/a)))
            // where GM - graviattional constant of a EarthMoon baricenter
            //       R - distance from baricenter to center of earth
            //       a - half of big axe of a orbit of a Earth around baricenter
            // first assumption a==R in this case V = sqrt(GMEarthMoon * 1/R)
            // 
            //VMoonOrbit = sqrt(GMEarthMoon / RBarisMoon);
            //VMoonOrbit = sqrt(GMEarthMoon * (2.0/ RBarisMoon - 1.0/(405696000.0*(EarthM/(EarthM+MoonM))) ));
            //VMoonOrbit = sqrt(GMEarthMoon * (2.0/ RBarisMoon - 1.0/(405696000.0*(EarthM/(EarthM+MoonM))) ));
            //aproxim = 155.23657531046074952275554885572/aproxim;
            //aproxim = -2.4009724978895439513029247826218;//1.0496393116944882060883002708885;//2.4;//aproxim*aproxim;
            //aproxim = -1.0496393116944882060883002708885;
            //aproxim = - 0.5728;
            
            AssignFromNASAData(&SolarSystem, dStartJD);
       
            printf("\n uses JPL coordinates and velocities");
            
            //AjustKeplerPosition(Sat.ProbTSec[Sat.Elem],ProbAph,ProbPer,ProbSmAx,Sat.ProbEcc[Sat.Elem], 
            //       Sat.ProbIncl[Sat.Elem], Sat.ProbAscNode[Sat.Elem], Sat.ProbArgPer[Sat.Elem], Sat.ProbEpochS[Sat.Elem], ProbCurTimeS);
            //printf("\n Was  A Prob = %f", ProbSmAx);
			//0         1         2         3         4         5         6         7
			//01234567890123456789012345678901234567890123456789012345678901234567890
			//1 25544U 98067A   04236.56031392  .00020137  00000-0  16538-3 0  5135\
			//2 25544  51.6335 341.7760 0007976 126.2523 325.9359 15.70406856328903"
			// format readings from FORTRAN:
			// 702 FORMAT(18X,D14.8,1X,F10.8,2(1X,F6.5,I2),/,7X,2(1X,F8.4),1X,F7.7,2(1X,F8.4),1X,F11.8)
			//
		    // "04236.56031392" Element Set Epoch (UTC) D14.8
			//  04                   - year
			//    236.56031392       - day
			// double ProbEpochS;
			//  "_.00020137"      1st Derivative of the Mean Motion with respect to Time F10.8
			//double ProbFirstDervMeanMotion; XNDT2O
			// "_00000-0"        2nd Derivative of the Mean Motion with respect to Time (decimal point assumed) 1X,F6.5,I2
			//double ProbSecondDervmeanMotion;,XNDD6O,IEXP
			//// "_16538-3"        B* Drag Term 1X,F6.5,I2
			//double ProbDragterm; ,BSTAR,IBEXP
			// "0"              Element Set Type
			//unsigned char ProbElementSetType;
		    // "2"             Line Number
			// "25544"         Object Identification Number
			// "_51.6335"       Orbit Inclination (degrees) 1X,F8.4
			//double ProbIncl;
			// "341.7760"      Right Ascension of Ascending Node (degrees) 1X,F8.4
			//double ProbAscNode;
			// "0007976"       Eccentricity (decimal point assumed) F7.7
			//double ProbEcc;
			// "126.2523"      Argument of Perigee (degrees) 1X,F8.4
			//double ProbArgPer;
			//double ProbPer;
			// "325.9359"      Mean Anomaly (degrees) 1X,F8.4
			//double ProbMeanAnom;
			// "15.70406856"    Mean Motion (revolutions/day)F11.8
			//double ProbMeanMotion; XNO
			// "328903"        Revolution Number at Epoch
			//double ProbRevAtEpoch;
			
			
			SUN_08 (1950,1,0,0,0,GST,SLONG,SRASN,SDEC);
			GreenwichA = GreenwichAscension(50001.0);
            // GreenwichA == GST but SUN_08 works for 2011
			//SRASN = SRASN * 180/ M_PI;
            for (int nSat = 0; nSat <Sat.Elem; nSat++)
            {
			    
			    long double AE = 1.0;
			    long double XKMPER = 6378.1350; //XKMPER kilometers/Earth radii 6378.135
			    long double XKE = BIG_XKE;//.743669161E-1;
                //XKE = sqrt(Gbig * SolarSystem.M[EARTH])*pow(AE/*(long double)6378.1350*//(long double)1440.0, (long double)3.0/(long double)2.0);;
			    long double XJ2 = 1.082616E-3;
			    long double CK2=.5*XJ2*AE*AE;
			    //double XMNPDA = 1440.0;
			    //double TEMP=2.0*M_PI/XMNPDA/XMNPDA;
			    //double XNO=ProbMeanMotion*TEMP*XMNPDA;

    			long double XMNPDA = 1440.0; // XMNPDA time units(minutes) /day 1440.0
	    		long double TEMP=2*M_PI/XMNPDA/XMNPDA; // 2*pi / (1440 **2)
		    	long double XNO=Sat.ProbMeanMotion[nSat]*TEMP*XMNPDA; // rotation per day * 2*pi /1440 == rotation per day on 1 unit (1 min)
			    long double XNDT2O=Sat.ProbFirstDervMeanMotion[nSat]*TEMP;
			    long double XNDD6O=Sat.ProbSecondDervmeanMotion[nSat]*TEMP/XMNPDA;

			    /*
			    long double A1=pow((XKE/XNO),(long double)2.0/(long double)3.0);
			    TEMP=1.5*CK2*(3.*cos(Sat.ProbIncl[nSat])*cos(ProbIncl)-1.)/pow((1.-Sat.ProbEcc[nSat]*Sat.ProbEcc[nSat]),(long double)1.5);
			    long double DEL1=TEMP/(A1*A1);
			    long double AO=A1*(1.-DEL1*(.5*(2.0/3.0)+DEL1*(1.+134./81.*DEL1)));
			    long double DELO=TEMP/(AO*AO);
			    long double XNODP=XNO/(1.+DELO);
			    //IF((TWOPI/XNODP/XMNPDA) .GE. .15625) IDEEP=1
                */
			    long double BSTAR=Sat.ProbDragterm[nSat]/AE;
			    //TSINCE=TS
			    //IFLAG=1
			    // first one does not use BSTAR
			    //SGP((dStartJD - Sat.ProbEpoch[nSat])*XMNPDA, XNDT2O,XNDD6O,BSTAR,Sat.ProbIncl[nSat], Sat.ProbAscNode[nSat],Sat.ProbEcc[nSat], 
                //               Sat.ProbArgPer[nSat], Sat.ProbMeanAnom[nSat],XNO, 
			    //	tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
			    //tProbX=tProbX*XKMPER/AE*1000.0;
			    //tProbY=tProbY*XKMPER/AE*1000.0;
			    //tProbZ=tProbZ*XKMPER/AE*1000.0;
			    //tProbVX=tProbVX*XKMPER/AE*XMNPDA/86400.*1000.0;
			    //tProbVY=tProbVY*XKMPER/AE*XMNPDA/86400.*1000.0;
			    //tProbVZ=tProbVZ*XKMPER/AE*XMNPDA/86400.*1000.0;
			    // second one does not use XNDT2O,XNDD6O
                // in next call XN0 and ProbMeanMotion connected by a formula:
                // Sat.ProbMeanMotion[nSat] = XNO / (2*pi) * 1440.0
                BSTAR = 0.0;
                XNDD6O = 0.0;
                XNDT2O =0.0;
			    SGP4((dStartJD - Sat.ProbEpoch[nSat])*XMNPDA, XNDT2O,XNDD6O,BSTAR,Sat.ProbIncl[nSat], Sat.ProbAscNode[nSat],Sat.ProbEcc[nSat], 
                    Sat.ProbArgPer[nSat], Sat.ProbMeanAnom[nSat],XNO, 
				    tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
			    tProbX=tProbX*XKMPER/AE*1000.0;
			    tProbY=tProbY*XKMPER/AE*1000.0;
			    tProbZ=tProbZ*XKMPER/AE*1000.0;
			    tProbVX=tProbVX*XKMPER/AE*XMNPDA/86400.*1000.0;
			    tProbVY=tProbVY*XKMPER/AE*XMNPDA/86400.*1000.0;
			    tProbVZ=tProbVZ*XKMPER/AE*XMNPDA/86400.*1000.0;

                KeplerPosition(Sat.ProbEpoch[nSat],dStartJD,      // prob epoch, and curent time
	    			    Sat.ProbTSec[nSat], // - orbit period in sec
		    		    Sat.ProbEcc[nSat],             // - Eccentricity
                        Sat.ProbIncl[nSat],            // - Inclination
                        Sat.ProbAscNode[nSat],         // - Longitude of ascending node
                        Sat.ProbArgPer[nSat],          // - Argument of perihelion
                        Sat.ProbMeanAnom[nSat],        // - Mean Anomaly (degrees)
                        BSTAR,
                        Gbig *SolarSystem.M[EARTH],0,
                        tempProbX,tempProbY,tempProbZ,tempProbVX,tempProbVY,tempProbVZ, Sat.ProbMeanMotion[nSat]);

			    long double tProbTSec = 0;
			    long double tProbEcc = 0;
                long double tProbIncl = 0;
                long double tProbAscNode = 0;
                long double tProbArgPer = 0;
                long double tProbMeanAnom = 0;
            
    			/*
	    		//long double e_ecc = Sat.ProbEcc[nSat];
		    	//long double e_incl = Sat.ProbIncl[nSat];
			    //long double e_asc_node = Sat.ProbAscNode[nSat];
			    //long double e_arg_per = Sat.ProbArgPer[nSat];
			    //long double e_mean_anomaly = Sat.ProbMeanAnom[nSat];
			    long double e_q=0.0;
			    long double e_major_axis=0.0;
			    long double e_t0=0.0;
			    long double e_w0=0.0;
			    long double e_angular_momentum=0.0;
			    long double e_perih_time=0.0;
			    long double e_minor_to_major=0.0;
			    long double e_lon_per=0.0;
			    long double e_sideways_x=0.0;
			    long double e_sideways_y=0.0;
			    long double e_sideways_z=0.0;
			    long double vec_x=0.0;
			    long double vec_y=0.0;
			    long double vec_z=0.0;
			    long double e_perih_vec_x=0.0,e_perih_vec_y=0.0,e_perih_vec_z=0.0;
			    long double loc_x=0.0, loc_y=0.0, loc_z=0.0, loc_r=0.0,vel_x=0.0,vel_y=0.0,vel_z=0.0,t=0.0;
    
	    		// first (mass) from NASA and second is just convinient constant
		    	long double gm = SolarSystem.M[EARTH] * Gbig;

    			long double e_epoch = Sat.ProbEpoch[Sat.Elem]*24.0*60.0*60.0;
                    
			    do_element_setup( e_epoch, Sat.ProbEcc[nSat], Sat.ProbIncl[nSat], Sat.ProbAscNode[nSat], Sat.ProbArgPer[nSat],Sat.ProbMeanAnom[nSat],
								e_q, e_major_axis,e_t0,e_w0,e_angular_momentum,e_perih_time,
								e_minor_to_major, e_lon_per,
								e_sideways_x, e_sideways_y, e_sideways_z,
								vec_x, vec_y,vec_z, gm, Sat.ProbTSec[Sat.Elem]);
			    // this will reset mean anomaly to time == dStartJD*24.0*60.0*60.0
			    posn_and_vel( e_epoch, Sat.ProbEcc[nSat], Sat.ProbIncl[nSat], Sat.ProbAscNode[nSat], Sat.ProbArgPer[nSat],Sat.ProbMeanAnom[nSat],
						e_q, e_major_axis,e_t0,e_w0,e_angular_momentum,e_perih_time,
						e_minor_to_major, e_lon_per,
						e_sideways_x, e_sideways_y, e_sideways_z,
						vec_x,vec_y,vec_z,
						loc_x, loc_y, loc_z, loc_r,vel_x,vel_y,vel_z,dStartJD*24.0*60.0*60.0, gm);
						*/
			    long double tVX = tempProbVX;
			    long double tVY = tempProbVY;
			    long double tVZ = tempProbVZ;
			    {
// that is for debugging GPS sattelites only
#if 0

                    //01234567890123456789012345678901234567890123456789012345678901234567890
                    //1 25544U 98067A   04236.56031392  .00020137  00000-0  16538-3 0  5135\
                    //2 25544  51.6335 341.7760 0007976 126.2523 325.9359 15.70406856328903"
                    //    };
                    // "ISS (ZARYA)"    The common name for the object based on information from the SatCat.
                    // "1"              Line Number
                    // "25544"          Object Identification Number
                    // "U"              Elset Classification
                    // "98067A"         International Designator
                    //  98                   - designate the launch year of the object
                    //    067                - launch number, starting from the beginning of the year
                    //       A               - indicates the piece of the launch: "A" is a payload 
                    // "04236.56031392" Element Set Epoch (UTC)
                    //  04                   - year
                    //    236.56031392       - day
                    // "_.00020137"      1st Derivative of the Mean Motion with respect to Time
                    // "_00000-0"        2nd Derivative of the Mean Motion with respect to Time (decimal point assumed)
                    // "_16538-3"        B* Drag Term
                    // "0"              Element Set Type
                    // "_513"            Element Number
                    // "5"              Checksum
                    //                        The checksum is the sum of all of the character in the data line, modulo 10. 
                    //                        In this formula, the following non-numeric characters are assigned the indicated values: 
                    //                        Blanks, periods, letters, '+' signs -> 0
                    //                        '-' signs -> 1
                    // "2"             Line Number
                    // "25544"         Object Identification Number
                    // "_51.6335"       Orbit Inclination (degrees)
                    // "341.7760"      Right Ascension of Ascending Node (degrees)
                    // "0007976"       Eccentricity (decimal point assumed)
                    // "126.2523"      Argument of Perigee (degrees)
                    // "325.9359"      Mean Anomaly (degrees)
                    // "15.70406856"    Mean Motion (revolutions/day)
                    // "328903"        Revolution Number at Epoch

                    // NAVSTAR 54 (USA 177)
                    // 1 28190U 04009A   12103.48551084  .00000032  00000-0  10000-3 0  3101
                    // 2 28190 055.0111 121.3785 0080561 006.7787 353.3670 02.00552748 59113
                    //
                    //         orbit inclanation
                    //                  Right Ascension of Ascending Node



                    // testing GPS sattelite #0x13=19 on UTC= 492765.99999995105 week: 1683 time 16:52:46
                    long double m1ProbTSec, m1ProbEcc, m1ProbIncl, m1ProbAscNode, m1ProbArgPer, m1ProbMeanAnom;
                    long double m2ProbTSec, m2ProbEcc, m2ProbIncl, m2ProbAscNode, m2ProbArgPer, m2ProbMeanAnom;
                    tProbX = -24891582.164414;
                    tProbY = -6179787.716544;
                    tProbZ = 7602393.375224;
                    tProbVX = -760.24581128329362;
                    tProbVY = -633.026865;
                    tProbVZ = -2950.637702;
                    DumpKeplers(m1ProbTSec, // - orbit period in sec
		    		    m1ProbEcc,             // - Eccentricity
                        m1ProbIncl,            // - Inclination
                        m1ProbAscNode,         // - Longitude of ascending node
                        m1ProbArgPer,          // - Argument of perihelion
                        m1ProbMeanAnom,        // - Mean Anomaly (degrees)
                        SolarSystem.M[EARTH],0.0,
                        tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
				    m1ProbIncl = m1ProbIncl / M_PI * 180;
				    m1ProbAscNode = m1ProbAscNode / M_PI * 180;
				    m1ProbArgPer = m1ProbArgPer / M_PI * 180;
				    m1ProbMeanAnom = m1ProbMeanAnom / M_PI * 180;
                    // tProbArgPer	6.0013685480087604	double 343.85308910345697
		            // tProbAscNode	0.28308821186056238	double 16.219759769515520
		            // tProbEcc	0.35006709180942119	double
		            // tProbIncl	1.4375592933975376	double 82.366080311487735
		            // tProbMeanAnom	3.1240698993743474	double 178.99602013800990
		            // tProbTSec	27757.763563012952	double
		        
                    // testing GPS sattelite #13 on UTC= 492841.99999994022 (75.99999998917 sec later)
                    tProbX = -24948509.734849;
                    tProbY = -6227303.238357;
                    tProbZ = 7377692.112397;
                    tProbVX = -737.80818258029888;
                    tProbVY = -617.413576;
                    tProbVZ = -2962.493300;
                    DumpKeplers(m2ProbTSec, // - orbit period in sec
				        m2ProbEcc,             // - Eccentricity
                        m2ProbIncl,            // - Inclination
                        m2ProbAscNode,         // - Longitude of ascending node
                        m2ProbArgPer,          // - Argument of perihelion
                        m2ProbMeanAnom,        // - Mean Anomaly (degrees)
                        SolarSystem.M[EARTH],0.0,
                        tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
				    m2ProbIncl = m2ProbIncl / M_PI * 180;
				    m2ProbAscNode = m2ProbAscNode / M_PI * 180;
				    m2ProbArgPer = m2ProbArgPer / M_PI * 180;
				    m2ProbMeanAnom = m2ProbMeanAnom / M_PI * 180;
                    // tProbArgPer	6.0101929578191724	double 344.35869054228738
		            // tProbAscNode	0.28215312764826711	double 16.166183390661683
		            // tProbEcc	0.34889965261994360	double
		            // tProbIncl	1.4407074570581755	double 82.546456802458735
		            // tProbMeanAnom	3.1243871650157362	double 179.01419812024599
		            // tProbTSec	27795.543500545853	double
                    // position and Keplers does not match == something wrong - also wrong inclanation that sattelite must be 55 degree
#endif
    				// mean amomaly on curent time
	    			DumpKeplers(tProbTSec, // - orbit period in sec
		    		    tProbEcc,             // - Eccentricity
                        tProbIncl,            // - Inclination
                        tProbAscNode,         // - Longitude of ascending node
                        tProbArgPer,          // - Argument of perihelion
                        tProbMeanAnom,        // - Mean Anomaly (degrees)
                        SolarSystem.M[EARTH],0.0,
                        tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
			    	tProbIncl = tProbIncl / M_PI * 180;
				    tProbAscNode = tProbAscNode / M_PI * 180;
    				tProbArgPer = tProbArgPer / M_PI * 180;
	    			tProbMeanAnom = tProbMeanAnom / M_PI * 180;
                    // mean amomaly on curent time
			    	DumpKeplers(tProbTSec, // - orbit period in sec
				        tProbEcc,             // - Eccentricity
                        tProbIncl,            // - Inclination
                        tProbAscNode,         // - Longitude of ascending node
                        tProbArgPer,          // - Argument of perihelion
                        tProbMeanAnom,        // - Mean Anomaly (degrees)
                        SolarSystem.M[EARTH],0.0,
                        tempProbX,tempProbY,tempProbZ,tempProbVX,tempProbVY,tempProbVZ);
	    			tProbIncl = tProbIncl / M_PI * 180;
		    		tProbAscNode = tProbAscNode / M_PI * 180;
			    	tProbArgPer = tProbArgPer / M_PI * 180;
				    tProbMeanAnom = tProbMeanAnom / M_PI * 180;
                    // mean amomaly on curent time
	    			/*DumpKeplers(tProbTSec, // - orbit period in sec
		    		    tProbEcc,             // - Eccentricity
                        tProbIncl,            // - Inclination
                        tProbAscNode,         // - Longitude of ascending node
                        tProbArgPer,          // - Argument of perihelion
                        tProbMeanAnom,        // - Mean Anomaly (degrees)
                        SolarSystem.M[EARTH],0.0,
                        tempProbX,tempProbY,tempProbZ,vel_x,vel_y,vel_z);
				    tProbIncl = tProbIncl / M_PI * 180;
				    tProbAscNode = tProbAscNode / M_PI * 180;
				    tProbArgPer = tProbArgPer / M_PI * 180;
				    tProbMeanAnom = tProbMeanAnom / M_PI * 180;

    				long double VecAngle = AngleBtw(tempProbVX,tempProbVY,tempProbVZ,vel_x,vel_y,vel_z);
	    			printf("\n angle =%f",VecAngle);
		    		*/
			    	//printf(" 
			    }

                //printf("\n Was ProbPer = %f", ProbPer);
                //printf("\n Was ProbAph = %f", ProbAph);

                // adjust Prob + earth position and speed based on a point of a centre of gravity
                // for today only one sattelite
                //Sat.Elem = 1;
                Sat.flInUse[nSat] = 1;
                Sat.Lambda = GreenwichA;
                //Sat.Lambda  = M_PI * 15.0 / 8.0; //3.0 / 8.0; 11.0 / 8.0;
                //Sat.M[0] = ProbM;
                Sat.X[nSat] = tempProbX + SolarSystem.X[EARTH];
                Sat.Y[nSat] = tempProbY + SolarSystem.Y[EARTH];
                Sat.Z[nSat] = tempProbZ + SolarSystem.Z[EARTH];
                Sat.VX[nSat] = tempProbVX + SolarSystem.VX[EARTH];
                Sat.VY[nSat] = tempProbVY + SolarSystem.VY[EARTH];
                Sat.VZ[nSat] = tempProbVZ + SolarSystem.VZ[EARTH];
                // SGP4
                Sat.X[nSat] = tProbX + SolarSystem.X[EARTH];
                Sat.Y[nSat] = tProbY + SolarSystem.Y[EARTH];
                Sat.Z[nSat] = tProbZ + SolarSystem.Z[EARTH];
                Sat.VX[nSat] = tProbVX + SolarSystem.VX[EARTH];
                Sat.VY[nSat] = tProbVY + SolarSystem.VY[EARTH];
                Sat.VZ[nSat] = tProbVZ + SolarSystem.VZ[EARTH];
                // this has to be set to properly calc helper variables
                TimeSl = 1.0 / IterPerSec;
            }    
            // this is done to reduce errors and avoid unnessary 5 mul/div operations
            // temporary X_, VX_ will just added (in paralel calculations can be done actualy faster) 
            for (int i = 0; i <Sat.Elem; i++)
            {
                Sat.VX_[i] = Sat.VX[i]/* Sat.M[i]*/ /TimeSl ;
                Sat.VY_[i] = Sat.VY[i]/* Sat.M[i]*/ /TimeSl;
                Sat.VZ_[i] = Sat.VZ[i]/* Sat.M[i]*/ /TimeSl;

                Sat.X_[i] = Sat.X[i]/* Sat.M[i]*/ /TimeSl/TimeSl ;
                Sat.Y_[i] = Sat.Y[i]/* Sat.M[i]*/ /TimeSl/TimeSl ;
                Sat.Z_[i] = Sat.Z[i]/* Sat.M[i]*/ /TimeSl/TimeSl ;
            }
            for (int n = 0 ; n < MAX_COEF_J; n++)
            {
                
                //l	m	            C                                S
	            //2	0	  -0.10826360229840D-02	                       0.0
                Sat.J[n] = //sqrt(2*(long double)n+1)* // coeff already normalized
                    (Clm[0][n]);
                for (int k = 0; k < MAX_COEF_J; k++)
                {
                    Sat.CNK[n][k] = Clm[k][n];
                    Sat.SNK[n][k] = Slm[k][n];
                    // CNK = sqrt(2*(2*n+1)) * sqrt(((n-k)!/(n+k)!) * Clm
                    // SNK = sqrt(2*(2*n+1) * sqrt((n-k)!/(n+k)!)) * Slm
                }
            }
            // amount of J coeff used in calcualtion
            Sat.iLeg = 8;
            Sat.iLeg_longit = 0; // no longitude in calculation
            Sat.Lambda = -2;
            Sat.LegBody = EARTH;
            
        }

        IF_XML_READ(Targetlongitude) // dolgota
        {
            Targetlongitude = atof(pszQuo);
        }

        IF_XML_READ(Targetlatitude) // shirota
        {
            Targetlatitude = atof(pszQuo);
        }

    XML_SECTION_END;

    XML_SECTION(Engine);
        // this will switch on engine 0,1,2,3 and etc.
        IF_XML_READ(EngineNumber)
        {
            if (++EnginesCount <= MAX_ENGINES)
            {
            }
            else
            {
                EnginesCount--;
                printf("\n Max engines limit reached - .XML file is incorrect");
            }
            Engine[EnginesCount-1].iLine = 0;
            Engine[EnginesCount-1].TotalImpulse = 0.0;
            Engine[EnginesCount-1].EngineOn = 0;
            Engine[EnginesCount-1].EngineDone = 0;
            Engine[EnginesCount-1].iCalculate = 1;
            Engine[EnginesCount-1].XVec = 0;
            Engine[EnginesCount-1].YVec = 0;
            Engine[EnginesCount-1].ZVec = 0;
            Engine[EnginesCount-1].OptimizationFirstDirectionSwitch = 0;
            //Engine[EnginesCount-1].iNextStepNow = 0;
            //Engine[EnginesCount-1].iNextSteps = 0;
            //for (int ii=0; ii< sizeof(Engine[EnginesCount-1].NextEngineToOptimize)/sizeof(int); ii++)
            //{
            //    Engine[EnginesCount-1].NextEngineToOptimize[ii] = -1;
            //    Engine[EnginesCount-1].NextEngineToOptimizationType[ii] = -1;
            //    Engine[EnginesCount-1].NextEngineToOptimizationNearBody[ii] = -1;
            //    Engine[EnginesCount-1].WhatWillBeLastEngine[ii] = -1;
            //    Engine[EnginesCount-1].NextEngineToOptimizeCalculate[ii] = -1;
            //}
            //Engine[EnginesCount-1].iNumberOfTryValues = 0;
            Engine[EnginesCount-1].SeartchForPeriod = 0.0;
            //Engine[EnginesCount-1].Period = 0.0;
            Engine[EnginesCount-1].iCountApogPerig = 0;
            printf("\n getting Engine %d parameters", EnginesCount-1);
        }
        // just for convinience for orbit trajectory calculation
        IF_XML_READ(PropCoeff)
        {
            EngCoeff = atof(pszQuo);
        }
        IF_XML_READ(Calculate)
        {
            Engine[EnginesCount-1].iCalculate = atoi(pszQuo);
        }
        // total weight of propellant
        IF_XML_READ(Weight)
        {
            Engine[EnginesCount-1].Weight = EngCoeff*atof(pszQuo);
            printf("\n Engine %d parameters propelant mass= %f", EnginesCount-1,Engine[EnginesCount-1].Weight);
        }
        IF_XML_READ(TotalWeight)
        {
            Engine[EnginesCount-1].TotalWeight = EngCoeff*atof(pszQuo);
            printf("\n Engine %d parameters total mass= %f", EnginesCount-1,Engine[EnginesCount-1].TotalWeight);
        }
        
        // parameters of an engine telta time from a plot
        IF_XML_READ(DeltaT)
        {
            Engine[EnginesCount-1].IteraPerSec = atoi(pszQuo);
            Engine[EnginesCount-1].DeltaTime = 1.0/((double)(Engine[EnginesCount-1].IteraPerSec));
        }
        // time when engine will fire from begining of a mission
        IF_XML_READ(FireTime)
        {
            Engine[EnginesCount-1].FireTime = atof(pszQuo);
        }
        // starting angle for a firing engine - in plane: centre of a Earth/Moon/ Planet (Y) and vector of velocity (X)
        // angle btw direction to center clockwise
        IF_XML_READ(FireAng1)
        {
            Engine[EnginesCount-1].Ang1 = atof(pszQuo);
        }
        // second starting angle for a firing: in plane prependicular to center of earth/moon/planet
        // angle btw velocity vector (Y) clockwise
        IF_XML_READ(FireAng2)
        {
            Engine[EnginesCount-1].Ang2 = atof(pszQuo);
        }
        // number of satellite in satellite array (based 0)
        IF_XML_READ(EngineOnSatellite)
        {
            Engine[EnginesCount-1].iEngineOnSatellite = atoi(pszQuo);
        }
        // all values from a plot
        IF_XML_READ(ImplVal)
        {
            if (Engine[EnginesCount-1].iLine >= MAX_IMPULSE_LINES)
            {
                printf("\n Error = impulse for the engine %d is out of range.\n Results of trajectory calcualtions will be wrong", EnginesCount-1);
            }
            else
            {
                Engine[EnginesCount-1].TotalImpulse += EngCoeff*atof(pszQuo)*Engine[EnginesCount-1].DeltaTime;
                Engine[EnginesCount-1].ValImpl[Engine[EnginesCount-1].iLine++] = EngCoeff*atof(pszQuo);
            }
        }
		IF_XML_READ(XVector)
        {
			Engine[EnginesCount-1].XVec = atof(pszQuo);
        }
		IF_XML_READ(YVector)
        {
			Engine[EnginesCount-1].YVec = atof(pszQuo);
        }
		IF_XML_READ(ZVector)
        {
			Engine[EnginesCount-1].ZVec = atof(pszQuo);
        }
        IF_XML_READ(NearBody)
        {
            Engine[EnginesCount-1].NearBody = atoi(pszQuo);
        }
        IF_XML_READ(OptimizationInitialStep)
        {
            Engine[EnginesCount-1].OptimizationInitialStep = atof(pszQuo);
            Engine[EnginesCount-1].OptimizationInitialStepCopy = atof(pszQuo);
        }
        IF_XML_READ(OptimizationDecCoef)
        {
            Engine[EnginesCount-1].OptimizationDecCoef = atof(pszQuo);
            Engine[EnginesCount-1].OptimizationDecCoefCopy = atof(pszQuo);
        }
        IF_XML_READ(OptimizationStop)
        {
            Engine[EnginesCount-1].OptimizationStop = atof(pszQuo);
        }
        IF_XML_READ(AngleOnBody)
        {
            Engine[EnginesCount-1].AngleOnBody = atoi(pszQuo);
        }
		IF_XML_READ(AngleType)
        {
            Engine[EnginesCount-1].AngleType = atoi(pszQuo);
        }
        //IF_XML_READ(NextEngineToOptimize)
        //{
        //    Engine[EnginesCount-1].NextEngineToOptimize[Engine[EnginesCount-1].iNextSteps] = atoi(pszQuo);
        //}
        //IF_XML_READ(NextEngineToOptimizationType)
        //{
        //    Engine[EnginesCount-1].NextEngineToOptimizationType[Engine[EnginesCount-1].iNextSteps] = atoi(pszQuo);
        //}
        //IF_XML_READ(NextEngineToOptimizationNearBody)
        //{
        //    Engine[EnginesCount-1].NextEngineToOptimizationNearBody[Engine[EnginesCount-1].iNextSteps] = atoi(pszQuo);
        //}
        //IF_XML_READ(NextEngineToOptimizeCalculate)
        //{
        //    Engine[EnginesCount-1].NextEngineToOptimizeCalculate[Engine[EnginesCount-1].iNextSteps] = atoi(pszQuo);
        //}
        //IF_XML_READ(WhatWillBeLastEngine)
        //{
        //    if (Engine[EnginesCount-1].NextEngineToOptimizationNearBody[Engine[EnginesCount-1].iNextSteps] == -1)
        //        Engine[EnginesCount-1].NextEngineToOptimizationNearBody[Engine[EnginesCount-1].iNextSteps] = Engine[EnginesCount-1].NearBody;
        //
        //    Engine[EnginesCount-1].WhatWillBeLastEngine[Engine[EnginesCount-1].iNextSteps++] = atoi(pszQuo);
        //}
        
        
        //IF_XML_READ(ValTry)
        //{
        //    Engine[EnginesCount-1].dValTryMaxMin[Engine[EnginesCount-1].iNumberOfTryValues] =0.0;
        //    Engine[EnginesCount-1].dValTry[Engine[EnginesCount-1].iNumberOfTryValues++] = atof(pszQuo);
        //}
    
    XML_SECTION_END;

            
            
    XML_SECTION(Optim);
        IF_XML_READ(EngineToOptimize)
        {
            if (++iOptPtr < MAX_OPTIM)
            {
                //iOptPtr++;
                Opt[iOptPtr-1].EngineToOptimize = atoi(pszQuo);
                Opt[iOptPtr-1].OptimizationInitialStep = 0.0;
                Opt[iOptPtr-1].OptimizationDecCoef = 0.0;
                Opt[iOptPtr-1].OptimizationStop = 0.0;
            }


        }
        IF_XML_READ(TrajectoryOptimizationType)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].TrajectoryOptimizationType = atoi(pszQuo);
        }
        IF_XML_READ(NearBody)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].NearBody = atoi(pszQuo);
        }
        IF_XML_READ(Calculate)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].Calculate = atoi(pszQuo);
        }
        IF_XML_READ(OptimizationInitialStep)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].OptimizationInitialStep = atof(pszQuo);
        }
        IF_XML_READ(OptimizationDecCoef)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].OptimizationDecCoef = atof(pszQuo);
        }
        IF_XML_READ(OptimizationStop)
        {
            if (iOptPtr < MAX_OPTIM)
                Opt[iOptPtr-1].OptimizationStop = atof(pszQuo);
        }
        IF_XML_READ(LastEngine)
        {
            if (iOptPtr < MAX_OPTIM)
            {
                Opt[iOptPtr-1].LastEngine = atoi(pszQuo);
                if ((iOptPtr-1) == StartOptim)
                {
                    LastEngine = atoi(pszQuo);
                    TrajectoryOptimizationType = Opt[iOptPtr-1].TrajectoryOptimizationType;
                    EngineToOptimize = Opt[iOptPtr-1].EngineToOptimize;
                    Engine[EngineToOptimize].iCalculate = Opt[iOptPtr-1].Calculate;
                    Engine[EngineToOptimize].NearBody = Opt[iOptPtr-1].NearBody;
                }
            }
        }
        IF_XML_READ(Period)
        {
            if (iOptPtr < MAX_OPTIM)
            {
                Opt[iOptPtr-1].Period = atof(pszQuo);
            }
        }
    XML_SECTION_END;
    XML_END;
}
// TRA.XML processing
void ParamDoAll(FILE *fInput)
{
    char szString[1024];
    char *pszQuo;
    char *pszQuo2;
    while(fgets(szString, sizeof(szString) -1, fInput))
    {
        if (strstr(szString, "TRA:section") != NULL)
        {
            if ((pszQuo = strstr(szString, "name=\"")) != NULL)
            {
                strcpy(szSection, pszQuo + sizeof("name=\"") - 1);
                if ((pszQuo2 = strchr(szSection, '\"')) != 0)
                    *pszQuo2 = 0;
            }
        }
        ParamCommon(szString);
        ParamSun(szString);
        ParamEarth(szString);
        ParamMoon(szString);
        ParamProb(szString);
    }
}
void makeExplanationText(char*szText, int iCalc, int iTraj, int iBody)
{
    switch(iCalc)
    {
    case CALC_PERIGEE: strcpy(szText,"\n<!-- calculates perigee to the ");break;
    case CALC_APOGEE: strcpy(szText,"\n <!-- calculates apogee to the ");break;
    case CALC_TARGET_PRACTICE: strcpy(szText,"\n<!-- calculates error in target practice to the moon");iBody=-1;break;
    case CALC_TARGET_POINT: strcpy(szText,"\n<!-- calculates error in traget prective to the point on the moon");iBody=-1;break;
    case CALC_PERIOD: strcpy(szText,"\n<!-- calculates orbit`s period after engine firing");iTraj=-1;iBody=-1;break;
    case CALC_INIT_PERIOD: strcpy(szText,"\n<!-- calculates orbit`s period ");iTraj=-1;iBody=-1;break;
    case CALC_FIRE_FIRST_ENGINE_TIME: strcpy(szText,"\n<!-- calculates first engine firing time");iTraj=-1;iBody=-1;break;
    case CALC_FIRE_SECOND_ENGINE_TIME:strcpy(szText,"\n<!-- calculates second engine firing time");iTraj=-1;iBody=-1;break;
    case CALC_FIRE_THIRD_ENGINE_TIME_TRY_ONE:strcpy(szText,"\n<!-- calculates third engine firing time(first try)");iTraj=-1;iBody=-1;break;
    case CALC_AT_APOGEE_DIFF_TO_3_4_DIST:strcpy(szText,"\n<!-- calculates difference from 3/4 earth-moon distance and apogee");iBody=-1;break;
    }
    switch(iBody)
    {
    case EARTH:strcat(szText,"Earth");break;
    case MOON:strcat(szText,"Moon");break;
    default:break;
    }
    switch(iTraj)
    {
    case MINIMUM_BY_TIME:strcat(szText,", search for minimum by adjusting time -->");break;
    case MAXIMUM_BY_TIME:strcat(szText,", search for maximum by adjusting time -->");break;
    case MINIMUM_BY_ANGLE:strcat(szText,", search for minimum by adjusting angle -->");break;
    case MINIMUM_BY_WEIGHT:strcat(szText,", search for minimum by adjusting engines weight -->");break;
    default:strcat(szText,".-->");break;

    }

}

void PostXMLToServer(char* URLServer, int urlport, char* URLFileName, char* FileToTransfer)
{
    char szWebServerResp[8096];
    FILE *FileTransfer = fopen(FileToTransfer,"rb");
    if (FileTransfer)
    {
        fseek(FileTransfer, 0L, SEEK_END);
        long iSize = ftell( FileTransfer);
        fseek(FileTransfer, 0L, SEEK_SET);
        char * szFileContent = (char*)malloc((size_t)iSize);
        if (szFileContent)
        {
            fread(szFileContent,iSize,1,FileTransfer);
            fclose(FileTransfer);
            CHttpConnection* m_MainHttpServer = NULL;
            CInternetSession  *m_MainInternetConnection = NULL;

            if (m_MainHttpServer == NULL)
            {
                m_MainInternetConnection = new CInternetSession("SessionToControlServer",12,INTERNET_OPEN_TYPE_DIRECT,NULL, // proxi name
                            NULL, // proxi bypass
				            INTERNET_FLAG_DONT_CACHE|INTERNET_FLAG_TRANSFER_BINARY);
		        try
		        {
		            m_MainHttpServer = 	m_MainInternetConnection->GetHttpConnection( URLServer, 0, urlport, NULL, NULL );
                }
	            catch(CInternetException *e)
		        {
		            m_MainHttpServer = NULL;
		        }
            }
            if (m_MainHttpServer)
            {
                CHttpFile* myCHttpFile = NULL;
                try
                {
                    myCHttpFile = m_MainHttpServer->OpenRequest( CHttpConnection::HTTP_VERB_POST,URLFileName, NULL,NULL, NULL, NULL, INTERNET_FLAG_EXISTING_CONNECT|	INTERNET_FLAG_DONT_CACHE|INTERNET_FLAG_RELOAD );
			    }
			    catch(CInternetException *e)
			    {
				    myCHttpFile = NULL;
			    }

			    if (myCHttpFile !=NULL)
			    {
				    try
				    {
                        CString strHeader = "Accept: text/*\r\n";
                        strHeader += "User-Agent: HttpCall\r\n";
                        strHeader += "Accept-Language: en-us\r\n";

                        //    strHeader += "Content-type: application/x-www-form-urlencoded\r\n";
                        //    strHeader += "REMOTE_USER: "+strUser+"\r\n";
                        //    strHeader += "Accept-Language: en-us\r\n";

                        myCHttpFile->AddRequestHeaders((LPCSTR)strHeader);
                        myCHttpFile->SendRequestEx(iSize,HSR_INITIATE,1);
                        myCHttpFile->WriteString((LPCTSTR)szFileContent);
                        myCHttpFile->EndRequest();

					    memset(szWebServerResp, 0, sizeof(szWebServerResp));
						DWORD dwSize;
						CString strSize;
						myCHttpFile->QueryInfo(HTTP_QUERY_CONTENT_LENGTH,strSize);
						dwSize = atoi(strSize.GetString());
    				    if (dwSize > (sizeof(szWebServerResp)-1))
					    {
						    for (DWORD dwread=0; dwread < dwSize; dwread+= (sizeof(szWebServerResp)-1))
						    {
							    if ((dwSize - dwread) > (sizeof(szWebServerResp)-1))
                                {
	                                if (myCHttpFile->Read(&szWebServerResp,(sizeof(szWebServerResp)-1)))
                                    {
                                    }
                                }
								else
                                {
								    if (myCHttpFile->Read(&szWebServerResp,(dwSize - dwread)))
                                    {
                                    }
                                }
                            }
                        }
						else
                        {
                            if (myCHttpFile->Read(&szWebServerResp,dwSize))
                            {
                            }
                        }
                    }
    				catch(CInternetException *e)
	    			{
		    			//ptrApp->m_MainHttpServer = NULL;
			    	}
    				myCHttpFile->Close();
	    			delete myCHttpFile;
		    	}
                m_MainHttpServer->Close();
                m_MainInternetConnection->Close();
		    }
            free(szFileContent);
        }
    }
}

FILE *VisualFile = NULL;

void dumpTRAvisual(long i)
{
    if (VisualFileSet == FALSE)
        return;
    if (VisualFile == NULL)
    {
		VisualFile = fopen(szTraVisualFileName, "w");
        if (VisualFile)
        {
		    fprintf(VisualFile,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
		    fprintf(VisualFile,"<Universe>\n");
        }
    }
	if (VisualFile)
    {

        for (int iSat = 0; iSat < Sat.Elem; iSat++)
        {
            fprintf(VisualFile,"	<Sat%d>\n",iSat);
            //fprintf(VisualFile,"		<ID>SatTra</ID>\n");
		    //fprintf(VisualFile,"		<time>%.18g</time>\n", dStartJD + ((double)i)/(24.0*60.0*60.0));
			fprintf(VisualFile,"		<X>%.18g</X>\n",Sat.X[iSat]-SolarSystem.X[EARTH]);
			fprintf(VisualFile,"		<Y>%.18g</Y>\n",Sat.Y[iSat]-SolarSystem.Y[EARTH]);
			fprintf(VisualFile,"		<Z>%.18g</Z>\n",Sat.Z[iSat]-SolarSystem.Z[EARTH]);
            if (iSat == 0) // only for a first satellite - check visibility from ground stations
            {
                // check is it from now time till end of the simulation
                if ((iTotalSec - i) <= dMinFromNow* 60)
                {
                    SYSTEMTIME ThatTime; 
                    TIME_ZONE_INFORMATION tmzone;
                    ConvertJulianDayToDateAndTime(dStartJD + ((double)(i))/(24.0*60.0*60.0), &ThatTime);
                    int Iret = GetTimeZoneInformation(&tmzone); 
                    double dEpoch = ConvertDateTimeToTLEEpoch(ThatTime.wDay, ThatTime.wMonth, ThatTime.wYear, ThatTime.wHour, ThatTime.wMinute, ThatTime.wSecond, ThatTime.wMilliseconds);
                    GreenwichA = GreenwichAscension(dEpoch);
                    SUN_08 (ThatTime.wYear,
                                iDayOfTheYearZeroBase(ThatTime.wDay, ThatTime.wMonth, ThatTime.wYear)+1 ,// in that function specified that 1 day is january 1
                                ThatTime.wHour,ThatTime.wMinute,ThatTime.wSecond,
                                GST,SLONG,SRASN,SDEC);
                    // only for red line do the check
                    for (int jGr = 0; jGr < iGr; jGr++)
                    {
                        double PosX;double PosY;double PosZ;
                        GetXYZfromLatLong(GrLong[jGr]+GST*180.0/M_PI-90.0, GrLat[jGr],PosX,PosY,PosZ, EarthR);
                        // check: does Cubesat in proximity of a ground station : angle btw vactor-radius and vector cubesat-ground station > 90 degree
                        if (AngleBtw(-PosX,-PosY,-PosZ,
                            (Sat.X[iSat]-SolarSystem.X[EARTH])-PosX,(Sat.Y[iSat]-SolarSystem.Y[EARTH])-PosY,(Sat.Z[iSat]-SolarSystem.Z[EARTH])-PosZ) > 90.0)
                        {
                            fprintf(VisualFile,"		        <X%d>%.18g</X%d>\n",jGr,PosX,jGr);
                            fprintf(VisualFile,"		        <Y%d>%.18g</Y%d>\n",jGr,PosY,jGr);
                            fprintf(VisualFile,"		        <Z%d>%.18g</Z%d>\n",jGr,PosZ,jGr);
                        }
                    }
                }
            }
            fprintf(VisualFile,"	</Sat%d>\n",iSat);
        }
	    if (OutLast == TRUE)
	    {
            SYSTEMTIME ThatTime; 
            ConvertJulianDayToDateAndTime(dStartJD + ((double)i)/(24.0*60.0*60.0), &ThatTime);
            double dEpoch = ConvertDateTimeToTLEEpoch(ThatTime.wDay, ThatTime.wMonth, ThatTime.wYear, ThatTime.wHour, ThatTime.wMinute, ThatTime.wSecond, ThatTime.wMilliseconds);
            GreenwichA = GreenwichAscension(dEpoch);
            SUN_08 (ThatTime.wYear,
                iDayOfTheYearZeroBase(ThatTime.wDay, ThatTime.wMonth, ThatTime.wYear) + 1 , // in fucntion spec that 1 Jan == 1 
                ThatTime.wHour,ThatTime.wMinute,ThatTime.wSecond,
                GST,SLONG,SRASN,SDEC);
            // GST is a position Greenwich in rad 

            //for moon rotation - moon looks at earth each time
            // that is a angle btw vector (0,-1,0) and moon position vector in geocentrical coordinates
            double MoonRot = AngleBtw(0,-1,0,SolarSystem.X[MOON] - SolarSystem.X[EARTH],SolarSystem.Y[MOON] - SolarSystem.Y[EARTH],SolarSystem.Z[MOON] - SolarSystem.Z[EARTH]);

            fprintf(VisualFile,"	<Object>\n");
            fprintf(VisualFile,"		<type>Earth</type>\n");
			fprintf(VisualFile,"		<X>%.18g</X>\n",SolarSystem.X[EARTH]);
			fprintf(VisualFile,"		<Y>%.18g</Y>\n",SolarSystem.Y[EARTH]);
			fprintf(VisualFile,"		<Z>%.18g</Z>\n",SolarSystem.Z[EARTH]);
			fprintf(VisualFile,"		<R>%.18g</R>\n",EarthR);
            fprintf(VisualFile,"		<Rot>%.18g</Rot>\n",GST*180.0/M_PI - 90.0); // convert to degree
			fprintf(VisualFile,"	</Object>\n");
			fprintf(VisualFile,"	<Object>\n");
			fprintf(VisualFile,"		<type>Sun</type>\n");
			fprintf(VisualFile,"		<X>%.18g</X>\n",SolarSystem.X[SUN]);
			fprintf(VisualFile,"		<Y>%.18g</Y>\n",SolarSystem.Y[SUN]);
			fprintf(VisualFile,"		<Z>%.18g</Z>\n",SolarSystem.Z[SUN]);
			fprintf(VisualFile,"		<R>%.18g</R>\n",SunR);
            fprintf(VisualFile,"		<Rot>%.18g</Rot>\n",0);
			fprintf(VisualFile,"	</Object>\n");
			fprintf(VisualFile,"	<Object>\n");
			fprintf(VisualFile,"		<type>Moon</type>\n");
			fprintf(VisualFile,"		<X>%.18g</X>\n",SolarSystem.X[MOON]);
			fprintf(VisualFile,"		<Y>%.18g</Y>\n",SolarSystem.Y[MOON]);
			fprintf(VisualFile,"		<Z>%.18g</Z>\n",SolarSystem.Z[MOON]);
			fprintf(VisualFile,"		<R>%.18g</R>\n",MoonR);
            fprintf(VisualFile,"		<Rot>%.18g</Rot>\n",MoonRot);
			fprintf(VisualFile,"	</Object>\n");
            fprintf(VisualFile,"	<ObjectTime>\n");
		    fprintf(VisualFile,"		<timeJD>%.18g</timeJD>\n", dStartJD + ((double)i)/(24.0*60.0*60.0));
            
            fprintf(VisualFile,"		<timeYYDDMMHHMMSS>%02d/%02d/%02d %02d:%02d:%02d</timeYYDDMMHHMMSS>\n", ThatTime.wYear-2000,ThatTime.wMonth,ThatTime.wDay,ThatTime.wHour,ThatTime.wMinute,ThatTime.wSecond);
            fprintf(VisualFile,"		<ReloadInSec>00001</ReloadInSec>\n"); // for the best case it is 1 sec refresh == that value has to be  
            fprintf(VisualFile,"		<dMinFromNow>%d</dMinFromNow>\n",(int)dMinFromNow); // how many minutes from now to aproximate
            fprintf(VisualFile,"	</ObjectTime>\n");
            for (int i= 0; i <iGr; i++)
            {
                fprintf(VisualFile,"	<GrSt>\n");
                // now need to calculate coordinates based on latitude and longitude
                double PosX;double PosY;double PosZ;
                GetXYZfromLatLong(GrLong[i], GrLat[i],PosX,PosY,PosZ, EarthR);
                //fprintf(VisualFile,"		<type>GrStn</type>\n");
                fprintf(VisualFile,"		<X>%.18g</X>\n",PosX);
                fprintf(VisualFile,"		<Y>%.18g</Y>\n",PosY);
                fprintf(VisualFile,"		<Z>%.18g</Z>\n",PosZ);
                fprintf(VisualFile,"	</GrSt>\n");
            }
            fprintf(VisualFile,"</Universe>\n");
			fclose(VisualFile);
            VisualFile = NULL;
            if (szTraVisualFileName[0] == '@') // yes! it is agly - that is a case when visualization output must to be submit to some server
            {
                PostXMLToServer(szURLTraVisualServer, UrlTraVisualPort, szURLTraVisualFileName, szTraVisualFileName);
            }
		}
	}
}
// output the same copy of the TRA.XML file
void dumpXMLParam(TRAOBJ *Sat, TRAIMPLOBJ *MyEngine, int iNumbOfEng)
{
    char szText[1024];
    FILE *EnginesFile = fopen("tra_out.xml", "w");
    int i,j;
    if (EnginesFile)
    {
#define XML_DUMPF(XML_PARAM) fprintf(EnginesFile,"\n    <TRA:setting name=\"%s\" value=\"%.18g\" />",#XML_PARAM,XML_PARAM);
#define XML_DUMPI(XML_PARAM) fprintf(EnginesFile,"\n    <TRA:setting name=\"%s\" value=\"%d\" />",#XML_PARAM,XML_PARAM);
        fprintf(EnginesFile,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>");
        fprintf(EnginesFile,"\n<TRA:data version=\"1.00\" xmlns:CT=\"http://www.adobri.com/tra\">");
        fprintf(EnginesFile,"\n<TRA:section name=\"TraInfo\">\n");
        fprintf(EnginesFile,"\n<!-- starting date (1 jan 2000) = 2451544.5JD ");
        fprintf(EnginesFile,"\n     if this is not set then use keplers elements from a satelite 0 -->");
        XML_DUMPF(dStartJD);
        fprintf(EnginesFile,"\n\n    <TRA:setting name=\"RGBImageW\" value=\"%d\" />",bRGBImageW);
        fprintf(EnginesFile,"\n    <TRA:setting name=\"RGBImageH\" value=\"%d\" />",bRGBImageH);

        fprintf(EnginesFile,"\n<!-- int iProfile = 0; ");
        fprintf(EnginesFile,"\n              // 0 == XY , 1 == YZ, 2 == XZ 3 == -YZ 4 == -XZ 5==-XY");
        fprintf(EnginesFile,"\n              // 0 or XY is a view from North to south, 5 (- XY) is a view from south to north");
        fprintf(EnginesFile,"\n              // 1 or YZ is a view to easter -->");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"RGBView\" value=\"%d\" />", iProfile);

        fprintf(EnginesFile,"\n<!--  EARTH 2 MOON  9 -->\n");
        XML_DUMPI(RGBReferenceBody);
        fprintf(EnginesFile,"\n<!-- max amount of pictures -->");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"RGBMaxPictures\" value=\"%d\" />",iMaxSeq);
        fprintf(EnginesFile,"\n<!-- one picture per sec -->");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"RGBSecPerPictures\" value=\"%d\" />", iMaxCounter);
        fprintf(EnginesFile,"\n<!-- scale in m -->");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"RGBScale\" value=\"%.18g\" />\n",dRGBScale);
        XML_DUMPF(IterPerSec);
        XML_DUMPF(StartLandingIteraPerSec);

        fprintf(EnginesFile,"\n<!-- next vaue actualy just a reference it will be used in-->");
        fprintf(EnginesFile,"\n<!-- GMMoon = MoonM * GBig and GMEarth = EarthM*GBig      -->");
        fprintf(EnginesFile,"\n<!-- know GMMoon and GMEarth will give proper value for M -->");
        XML_DUMPF(Gbig);
        XML_DUMPF(TotalDays);
        XML_DUMPF(EarthCurTime);
        XML_DUMPF(AU); 
        for (i = 0; i < Sat->Elem; i++)
        {
            //XML_DUMPF(ProbM);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"ProbM\" value=\"%.18g\" />",Sat->M[i]);

            fprintf(EnginesFile,"\n    <TRA:setting name=\"ProbKeplerLine1\" value=\"%s",Sat->Kepler1);
            fprintf(EnginesFile,"    <TRA:setting name=\"ProbKeplerLine2\" value=\"%s",Sat->Kepler2);
            fprintf(EnginesFile,"    <TRA:setting name=\"ProbKeplerLine3\" value=\"%s",Sat->Kepler3);
        }
        fprintf(EnginesFile,"\n<!-- target point on the Moon -->");
        XML_DUMPF(Targetlongitude);
        XML_DUMPF(Targetlatitude);
        XML_DUMPF(SunM);
        XML_DUMPF(EarthR);
        XML_DUMPF(EarthRP);
        XML_DUMPF(EarthRE);
        XML_DUMPF(EarthM);
        XML_DUMPF(MassRatioSunToEarthPlusMoon);
        XML_DUMPF(EarthTSolSec);
        XML_DUMPF(EarthSmAx);
        XML_DUMPF(EarthSmAxAU);
        XML_DUMPF(EarthTDays);

        XML_DUMPF(MoonR);
        XML_DUMPF(MoonRP);
        XML_DUMPF(MoonRE);
        XML_DUMPF(MoonM);
        fprintf(EnginesFile,"\n\n<!--next line will force recalculation of SunM based on G constant -->");
        XML_DUMPF(GMSun);
        XML_DUMPF(GMEarthMoon);
        XML_DUMPF(GMEarth);
        XML_DUMPF(GMMoon);
        fprintf(EnginesFile,"\n\n<!-- next line will force calculation of a Earth and Moon mass -->");
        fprintf(EnginesFile,"\n<!-- it is not in use - value just for reference and fo trigger -->");
        XML_DUMPF(MassRatioEarthToMoon);
        fprintf(EnginesFile,"\n    <TRA:setting name=\"MoonKeplerLine1\" value=\"%s",szMoonKeplerLine1);
        fprintf(EnginesFile,"    <TRA:setting name=\"MoonKeplerLine2\" value=\"%s",szMoonKeplerLine2);
        fprintf(EnginesFile,"    <TRA:setting name=\"MoonKeplerLine3\" value=\"%s",szMoonKeplerLine3);
        fprintf(EnginesFile,"\n\n<!-- next line (value == \"1.0\") will force calculation"); 
        fprintf(EnginesFile,"\n                         of a Earth kepler position -->");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"EarthCalcKepler\" value=\"0.0\" />");
        fprintf(EnginesFile,"\n\n<!-- next line will force assigning velocities and positions ");
        fprintf(EnginesFile,"\n     based on JPL data -->\n");
        fprintf(EnginesFile,"\n    <TRA:setting name=\"UseJPLxyz\" value=\"1.0\" />");
        XML_DUMPI(StartOptim);
        XML_DUMPI(MaxOptim);

        //fprintf(EnginesFile,"\n    <TRA:setting name=\"UseJPLxyz\" value=\"1.0\" />");
        fprintf(EnginesFile,"\n\n<!-- last used engine (0,1,2,3..) in trajectory optimization or calculations ");
        fprintf(EnginesFile,"\n     setting lastengine == -1 disable optimization -->\n");
        XML_DUMPI(LastEngine);
        fprintf(EnginesFile,"\n\n<!-- do optimization of an engine N... (0,1,2,3...)");
        fprintf(EnginesFile,"\n     it is posible to use 4 engines but optimize engine 3-->\n");
        XML_DUMPI(EngineToOptimize);


        fprintf(EnginesFile,"\n\n<!-- Type of optimization");
        fprintf(EnginesFile,"\n      1 - search for a minimum by adjusting time of firing");
        fprintf(EnginesFile,"\n      2 - search for a maximum by adjusting time of firing");       
        fprintf(EnginesFile,"\n      4 - search fo maximum by adjusting angle of firing");
        fprintf(EnginesFile,"\n      6 - search fo minimum by adjusting weight of fuel");

        fprintf(EnginesFile,"\n  initial step and decrease steps for each search specified");
        fprintf(EnginesFile,"\n  individualy for each engine -->\n");

        XML_DUMPI(TrajectoryOptimizationType);

        fprintf(EnginesFile,"\n</TRA:section>\n\n            <!--     now all engines    -->");

        for (i = 0; i < iNumbOfEng; i++)
        {
            fprintf(EnginesFile,"\n<TRA:section name=\"Engine\" value=\"%d\" />", i);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"EngineNumber\" value=\"%d\" />", i);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"EngineOnSatellite\" value=\"%d\" />",MyEngine[i].iEngineOnSatellite);
            fprintf(EnginesFile,"\n    <!-- convinent coeff - instead of entry real values just assume scaled version-->");
            fprintf(EnginesFile,"\n    <TRA:setting name=\"PropCoeff\" value=\"1.0\" />");
            fprintf(EnginesFile,"\n\n    <TRA:setting name=\"Weight\" value=\"%.18g\" />",MyEngine[i].Weight);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"TotalWeight\" value=\"%.18g\" />",MyEngine[i].TotalWeight);
            fprintf(EnginesFile,"\n\n    <!-- iteration per sec from engine's plot -->");
            fprintf(EnginesFile,"\n    <TRA:setting name=\"DeltaT\" value=\"%.18g\" />",1.0/MyEngine[i].DeltaTime);

            fprintf(EnginesFile,"\n\n     <!-- 2- EARTH 9-MOON for calculation distanses-->");
            fprintf(EnginesFile,"\n           <TRA:setting name=\"NearBody\" value=\"%d\" />",MyEngine[i].NearBody);

            fprintf(EnginesFile,"\n\n     <!-- AngleType 0 - tangent line to orbit (elipse) oposit velocity");
            fprintf(EnginesFile,"\n                    1 - two angles set with reference to NearBody centre direction");
            fprintf(EnginesFile,"\n                    2 - 3 angles set vector fire (constant all fire) ");
            fprintf(EnginesFile,"\n                    3 - oposit vector of velocity");
            fprintf(EnginesFile,"\n                    4 - same direction as vector of velocity -->");
            fprintf(EnginesFile,"\n     <TRA:setting name=\"AngleType\" value=\"%d\" />",MyEngine[i].AngleType);
            fprintf(EnginesFile,"\n\n     <!-- 2- EARTH 9-MOON for firing angle -->");
            fprintf(EnginesFile,"\n     <TRA:setting name=\"AngleOnBody\" value=\"%d\" />",MyEngine[i].AngleOnBody);
            fprintf(EnginesFile,"\n\n    <!-- first angle: in a plane over vector from the center of NearBody and Sat"); 
            fprintf(EnginesFile,"\n         and vector of velocity. Angle: Centre,Sat,Direction ");
            fprintf(EnginesFile,"\n         (aggle == 90 degr is a Tangent line to elipse) -->");
            fprintf(EnginesFile,"\n         <TRA:setting name=\"FireAng1\" value=\"%.18g\" />",MyEngine[i].Ang1);
            fprintf(EnginesFile,"\n\n     <!-- second angle from projection of a velocity vector to a  ");
            fprintf(EnginesFile,"\n         plane perpendicular to direction to centre of nearbody -->");
            fprintf(EnginesFile,"\n    <TRA:setting name=\"FireAng2\" value=\"%.18g\" />",MyEngine[i].Ang2);
            fprintf(EnginesFile,"\n\n    <TRA:setting name=\"XVector\" value=\"%.18g\" />",MyEngine[i].XVec);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"YVector\" value=\"%.18g\" />",MyEngine[i].YVec);
            fprintf(EnginesFile,"\n    <TRA:setting name=\"ZVector\" value=\"%.18g\" />\n",MyEngine[i].ZVec);
            fprintf(EnginesFile,"\n\n     <TRA:setting name=\"OptimizationInitialStep\" value=\"%.18g\" />",MyEngine[i].OptimizationInitialStepCopy);
            fprintf(EnginesFile,"\n     <TRA:setting name=\"OptimizationDecCoef\" value=\"%.18g\" />" ,MyEngine[i].OptimizationDecCoefCopy);
            fprintf(EnginesFile,"\n     <TRA:setting name=\"OptimizationStop\" value=\"%.18g\" />\n",MyEngine[i].OptimizationStop);
            //fprintf(EnginesFile,"\n        <TRA:setting name=\"Period\" value=\"%f\" />",MyEngine[i].Period);
            
            //if (MyEngine[i].iNumberOfTryValues)
            //{
            //    for (j = 0; j < MyEngine[i].iNumberOfTryValues; j++)
            //    {
            //        fprintf(EnginesFile,"\n        <TRA:setting name=\"TryVal\" value=\"%.24g\" /> <!-- %f -->",MyEngine[i].dValTry[j],MyEngine[i].dValTryMaxMin);
            //    }
            //}

            fprintf(EnginesFile,"\n\n    <!-- set impulses in a time -->");
            fprintf(EnginesFile,"\n        <TRA:setting name=\"FireTime\" value=\"%.18g\" />",MyEngine[i].FireTime);

            for (j = 0; j < MyEngine[i].iLine; j++)
            {
                fprintf(EnginesFile,"\n        <TRA:setting name=\"ImplVal\" value=\"%.24g\" /> <!-- %f -->",MyEngine[i].ValImpl[j], j*MyEngine[i].DeltaTime );
            }
            fprintf(EnginesFile,"\n</TRA:section>");
            fprintf(EnginesFile,"\n\n");
        }
        for (j=0; j < iOptPtr; j++)
        {
            fprintf(EnginesFile,"\n<TRA:section name=\"Optim\" value=\"%d\" />", j);
            makeExplanationText(szText, Opt[j].Calculate,
                    Opt[j].TrajectoryOptimizationType, 
                    Opt[j].NearBody);
            fprintf(EnginesFile,szText);
            fprintf(EnginesFile,"\n         <TRA:setting name=\"EngineToOptimize\" value=\"%d\" />",Opt[j].EngineToOptimize);
            if (j == 0)
            {
                fprintf(EnginesFile,"\n\n<!-- Type of optimization");
                fprintf(EnginesFile,"\n      1 - search for a minimum by adjusting time of firing");
                fprintf(EnginesFile,"\n      2 - search for a maximum by adjusting time of firing");       
                fprintf(EnginesFile,"\n      4 - search fo maximum by adjusting angle of firing");
                fprintf(EnginesFile,"\n      6 - search fo minimum by adjusting weight of fuel");

                fprintf(EnginesFile,"\n  initial step and decrease steps for each search specified");
                fprintf(EnginesFile,"\n  individualy for each engine -->\n");
            }

            fprintf(EnginesFile,"\n           <TRA:setting name=\"TrajectoryOptimizationType\" value=\"%d\" />",Opt[j].TrajectoryOptimizationType);
            if (j == 0)
                fprintf(EnginesFile,"\n\n     <!-- 2- EARTH 9-MOON for calculation distanses-->");
            fprintf(EnginesFile,"\n           <TRA:setting name=\"NearBody\" value=\"%d\" />",Opt[j].NearBody);
            if (j == 0)
            {
                fprintf(EnginesFile,"\n\n     <!-- calculates ");
                fprintf(EnginesFile,"\n          (0) Perigee to a center of NearBody");
                fprintf(EnginesFile,"\n          (1) Apogee to a center of NearBody");
                fprintf(EnginesFile,"\n          (3) taget practice not far than distance earth-moon");
                fprintf(EnginesFile,"\n          (5) target practice to apoint on a moon's surface");
                fprintf(EnginesFile,"\n          (6) orbit period around near body");
                fprintf(EnginesFile,"\n          (8) orbit period before engine firing");
                fprintf(EnginesFile,"\n          (9) first engine firing time");
                fprintf(EnginesFile,"\n          (10) second engine firing time");
                fprintf(EnginesFile,"\n          (11) third engine firing time (first try)");
                fprintf(EnginesFile,"\n          (12) at apogee difference from 3/4 of a earth-moon distance");
                fprintf(EnginesFile,"\n               --> ");
            }
            fprintf(EnginesFile,"\n           <TRA:setting name=\"Calculate\" value=\"%d\" />",Opt[j].Calculate);
            fprintf(EnginesFile,"\n           <TRA:setting name=\"OptimizationInitialStep\" value=\"%.18g\" />",Opt[j].OptimizationInitialStep);
            fprintf(EnginesFile,"\n           <TRA:setting name=\"OptimizationDecCoef\" value=\"%.18g\" />",Opt[j].OptimizationDecCoef);
            fprintf(EnginesFile,"\n           <TRA:setting name=\"OptimizationStop\" value=\"%.18g\" />",Opt[j].OptimizationStop);
            fprintf(EnginesFile,"\n           <TRA:setting name=\"LastEngine\" value=\"%d\" />\n",Opt[j].LastEngine);
            fprintf(EnginesFile,"\n        <TRA:setting name=\"Period\" value=\"%f\" />",Opt[j].Period);
            fprintf(EnginesFile,"\n</TRA:section>");
            fprintf(EnginesFile,"\n\n");
        }

        fprintf(EnginesFile,"\n</TRA:section>\n\n");

        fprintf(EnginesFile,"\n</TRA:data>");
        fclose(EnginesFile);
    }
}

// this will be a check doe it need to optimize next engine base on a results of a prev one
int CheckWhatnext(TRAOPTIMOBJ* Opt, TRAOBJ *Sat, TRAIMPLOBJ *Eng, TRAIMPLOBJ *EngStore, 
                  int &CurentEngine, int &CurentEngineOptimizationType,int &WhatWillBeLastEngine, 
                  int iNumbOfEngines, double dCurTime)
{
    int iret = 0;
    int i;
    int iOld;
    int iNewCalculate;
    int iNewNearBody;
    char szText[1024];
    
    //if (Eng[CurentEngine].NextEngineToOptimize[Eng[CurentEngine].iNextStepNow] != -1)
    if (++iOptimizationStep <= MaxOptim)
    {
        //iOptimizationStep++;
        //CurentEngineOptimizationType = Eng[CurentEngine].NextEngineToOptimizationType[Eng[CurentEngine].iNextStepNow];
        CurentEngineOptimizationType = Opt[iOptimizationStep].TrajectoryOptimizationType;
        //WhatWillBeLastEngine = Eng[CurentEngine].WhatWillBeLastEngine[Eng[CurentEngine].iNextStepNow];
        WhatWillBeLastEngine = Opt[iOptimizationStep].LastEngine;
        //iNewCalculate = Eng[CurentEngine].NextEngineToOptimizeCalculate[Eng[CurentEngine].iNextStepNow];
        iNewCalculate = Opt[iOptimizationStep].Calculate;
        //iNewNearBody = Eng[CurentEngine].NextEngineToOptimizationNearBody[Eng[CurentEngine].iNextStepNow];
        iNewNearBody = Opt[iOptimizationStep].NearBody;
        iOld = CurentEngine;
        //CurentEngine = Eng[CurentEngine].NextEngineToOptimize[Eng[CurentEngine].iNextStepNow];
        CurentEngine = Opt[iOptimizationStep].EngineToOptimize;
        Eng[CurentEngine].iCalculate = iNewCalculate;
        Eng[CurentEngine].NearBody = iNewNearBody;
        if (Opt[iOptimizationStep].OptimizationInitialStep != 0.0)
        {
            Eng[CurentEngine].OptimizationInitialStep = Opt[iOptimizationStep].OptimizationInitialStep;
            Eng[CurentEngine].OptimizationStop = Opt[iOptimizationStep].OptimizationStop;
            Eng[CurentEngine].OptimizationDecCoef = Opt[iOptimizationStep].OptimizationDecCoef;
        }
        // for the same engine adjust next step counter
        //if (iOld == CurentEngine)
        //{
        //    Eng[CurentEngine].iNextStepNow++;
        //}
        makeExplanationText(szText, iNewCalculate, CurentEngineOptimizationType,iNewNearBody);
        printf("\n next: %s", szText);

        // chek is it posible to do next step now?
        if (iNewCalculate == CALC_FIRE_FIRST_ENGINE_TIME || iNewCalculate == CALC_FIRE_SECOND_ENGINE_TIME ||
            iNewCalculate == CALC_FIRE_THIRD_ENGINE_TIME_TRY_ONE)
        {
            //Eng[CurentEngine].FireTime = dCurTime + Eng[CurentEngine].Period;
            Eng[CurentEngine].FireTime = dCurTime + Opt[iOptimizationStep-1].Period;
            printf("\n firing time set=%f ",Eng[CurentEngine].FireTime);
            iret = CheckWhatnext(Opt, Sat, Eng, EngStore, CurentEngine, CurentEngineOptimizationType, WhatWillBeLastEngine, iNumbOfEngines, dCurTime);

        }
        else
            iret = 1;
        // may be need to update value in storage for optimization?
        if (EngStore)
        {
            
            for (i = 0; i < MAX_ENGINES; i++)
            {
                EngStore[i].iCalculate = Eng[i].iCalculate;
                EngStore[i].FireTime = Eng[i].FireTime;
                //EngStore[i].iNextStepNow = Eng[i].iNextStepNow;
                EngStore[i].NearBody = Eng[i].NearBody;
                //EngStore[i].Period = Eng[i].Period;
                EngStore[i].OptimizationInitialStep = Eng[i].OptimizationInitialStep;
                EngStore[i].OptimizationStop = Eng[i].OptimizationStop;
                EngStore[i].OptimizationDecCoef = Eng[i].OptimizationDecCoef;

            }
        }
    }
    else
    {
        // that mean it will be exit - just dump curent satelite and engines status
        dumpXMLParam(Sat, Eng, iNumbOfEngines);
    }
    return iret;
}

#define TEST_RUN_CALC_YEAR 1
#define TEST_RUN_ERROR 1
#define TEST_RUN_EARTH_ERROR 1
int main(int argc, char * argv[])
{
    //char szSatelite[] = {"ISS (ZARYA)\
//0         1         2         3         4         5         6         7
//01234567890123456789012345678901234567890123456789012345678901234567890
//1 25544U 98067A   04236.56031392  .00020137  00000-0  16538-3 0  5135\
//2 25544  51.6335 341.7760 0007976 126.2523 325.9359 15.70406856328903"
//    };
    // "ISS (ZARYA)"    The common name for the object based on information from the SatCat.
    // "1"              Line Number
    // "25544"          Object Identification Number
    // "U"              Elset Classification
    // "98067A"         International Designator
    //  98                   - designate the launch year of the object
    //    067                - launch number, starting from the beginning of the year
    //       A               - indicates the piece of the launch: "A" is a payload 
    // "04236.56031392" Element Set Epoch (UTC)
    //  04                   - year
    //    236.56031392       - day
    // "_.00020137"      1st Derivative of the Mean Motion with respect to Time
    // "_00000-0"        2nd Derivative of the Mean Motion with respect to Time (decimal point assumed)
    // "_16538-3"        B* Drag Term
    // "0"              Element Set Type
    // "_513"            Element Number
    // "5"              Checksum
    //                        The checksum is the sum of all of the character in the data line, modulo 10. 
    //                        In this formula, the following non-numeric characters are assigned the indicated values: 
    //                        Blanks, periods, letters, '+' signs -> 0
    //                        '-' signs -> 1
    // "2"             Line Number
    // "25544"         Object Identification Number
    // "_51.6335"       Orbit Inclination (degrees)
    // "341.7760"      Right Ascension of Ascending Node (degrees)
    // "0007976"       Eccentricity (decimal point assumed)
    // "126.2523"      Argument of Perigee (degrees)
    // "325.9359"      Mean Anomaly (degrees)
    // "15.70406856"    Mean Motion (revolutions/day)
    // "328903"        Revolution Number at Epoch
    int iDay;
    int iHour;
    int iSecond;
    int iFlag = 0;
    long i;
    int j;
    int flFindMax;
    int flFindMin;
    int OptimMin = 1;

    stateType  StateEarth;
    stateType  StateMoon;
    
    double minDeltaMinMaxD;
    double maxDeltaMinMaxD;
    long iCountMin;
    int flFindFirst1KmError = 1;
    double dErrorValue = 1000000.0;
    int iApog = 0;
    double Apog = 0.0;
    int iPerig = 0;
    double Perig = 0.0;
    double PerigMoon = 100000000000.0;
    double ApogPergTime = 0;
    double dRE;
    double dRM;
    double dRM0;
    int idRM = 0;
    double dRM1 =0;
    double dRM2 = 0;
    double dRMDelta = 1.0;
    int idRMDelta = 0;
    int StartSequence = 0;
    double dREMV;
    double SCH_Per = 10000000000000.0;
    double SCH_Apg = 0.0;
    double SCH_Dist = 0.0;
    int iSCH_Apg = 0;
    int iSCH_Per = 0;
    int SCH_ApgPerTime = 0;
    

#define MAXTRYANGLESDIR 6
    double FindMin;
    int iFindMin;
    double newTryAnglesDirDelta;

    int iFirstAngleDone = 0;
    int iMaxTryAnglesDir = MAXTRYANGLESDIR;
    int iTryAnglesDir = 0;
    double TryAnglesDirDelta = 0.25;
    double TryAnglesDir[MAXTRYANGLESDIR+1][3] = {0,0,0, 0,0,1, 0,0,-1, 0,1,0, 0,-1,0, 1,0,0, -1,0,0}; 
    double TryAnglesDirValues[MAXTRYANGLESDIR+1];
    double LastStepTryAnglesDirValuesX;
    double LastStepTryAnglesDirValuesY;
    double LastStepTryAnglesDirValuesZ;
#define LASTSTEPHIST 10
    int iLastStepHist;
    double LastStepHistX[LASTSTEPHIST];
    double LastStepHistY[LASTSTEPHIST];
    double LastStepHistZ[LASTSTEPHIST];
    TRAOBJ MyTry = SolarSystem;
    TRAOBJ MyTrySat = Sat;
    TRAIMPLOBJ MyEngine[MAX_ENGINES];

	long double tProbTSec,tProbEcc,tProbIncl,tProbAscNode,tProbArgPer,tProbMeanAnom;
	long double tX,tY,tZ,tVX,tVY,tVZ;
    long double TimeFEpoch=0;
    double ttProbX,ttProbY,ttProbZ,ttProbVX,ttProbVY,ttProbVZ;
    double tttX,tttVX;
    char szXMLFileName[3*_MAX_PATH] = {"tra.xml"};


#ifdef TEST_RUN_ERROR
    flFindFirst1KmError = 1;
#else
    flFindFirst1KmError = 0;
#endif
    //InitTraMap(&CurTraMap);
    Sat.Elem = 0;  // zero asatelites at the begining
    Initialize_Ephemeris("bin410.bin");
    if (argc == 2)
    {
        if (argv[1][1] == '?')
        {
            printf("\n TRA.EXE trajectory calculation software. Adobri Solutions LTD. Team Plan B");
            printf("\n licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.");
            printf("\n");
            printf("\n USAGE:");
            printf("\n  TRA.EXE  -?          - help");
            printf("\n  TRA.EXE              - process TRA.XML file");
            printf("\n  TRA.EXE <FILE NAME>  - process <FILE NAME>");
            printf("\n  TRA.EXE <URL> - process <URL> file");
            exit (0);
        }
        else
        {
            strcpy(szXMLFileName,argv[1]);
        }
    }
    if (strstr(szXMLFileName,"http://"))
    {
        // needs to copy original http file from server to a local with temporary name, than process temporary
        char szURLFileName[3*_MAX_PATH];
        char szURLServer[3*_MAX_PATH];
        char sztempFileName[3*_MAX_PATH];
        char szWebServerResp[8096];
        int UrlPort=80;
       	CHttpConnection* m_MainHttpServer = NULL;
    	CInternetSession  *m_MainInternetConnection = NULL;

        if (ParsURL(szURLServer, &UrlPort, szURLFileName, szXMLFileName))
        {
            strcpy(szXMLFileName, "@tra.xml");
        }
        AfxSocketInit();
        if (m_MainHttpServer == NULL)
	    {
		    m_MainInternetConnection = new CInternetSession("SessionToControlServer",12,INTERNET_OPEN_TYPE_DIRECT,NULL, // proxi name
				            NULL, // proxi bypass
				            INTERNET_FLAG_DONT_CACHE|INTERNET_FLAG_TRANSFER_BINARY);
		    try
		    {
			    m_MainHttpServer = 	m_MainInternetConnection->GetHttpConnection( szURLServer, 0, UrlPort, NULL, NULL );
    		}
	    	catch(CInternetException *e)
		    {
			    m_MainHttpServer = NULL;
		    }
	    }
	    if (m_MainHttpServer)
	    {
			CHttpFile* myCHttpFile = NULL;
			try
			{
				myCHttpFile = m_MainHttpServer->OpenRequest( CHttpConnection::HTTP_VERB_GET,
					szURLFileName,
					NULL,//((CGrStnApp*)AfxGetApp())->szLoginRQ,
					NULL,//12345678,
					NULL, 
					NULL, 
					INTERNET_FLAG_EXISTING_CONNECT|
					INTERNET_FLAG_DONT_CACHE|
					INTERNET_FLAG_RELOAD );
			}
			catch(CInternetException *e)
			{
				myCHttpFile = NULL;
			}

			if (myCHttpFile !=NULL)
			{
				try
				{
					myCHttpFile->SendRequest();
					memset(szWebServerResp, 0, sizeof(szWebServerResp));
					{
						DWORD dwSize;
						CString strSize;
						myCHttpFile->QueryInfo(HTTP_QUERY_CONTENT_LENGTH,strSize);
						dwSize = atoi(strSize.GetString());
                        FILE *TempFile = fopen(szXMLFileName, "wb");
                        if (TempFile)
                        {
						    if (dwSize > (sizeof(szWebServerResp)-1))
						    {
							    for (DWORD dwread=0; dwread < dwSize; dwread+= (sizeof(szWebServerResp)-1))
							    {
								    if ((dwSize - dwread) > (sizeof(szWebServerResp)-1))
                                    {
									    if (myCHttpFile->Read(&szWebServerResp,(sizeof(szWebServerResp)-1)))
                                        {
                                            fwrite(&szWebServerResp,(sizeof(szWebServerResp)-1),1,TempFile);
                                        }
                                    }
								    else
                                    {
									    if (myCHttpFile->Read(&szWebServerResp,(dwSize - dwread)))
                                        {
                                            fwrite(&szWebServerResp,(dwSize - dwread),1,TempFile);
                                        }
                                    }
							    }
						    }
						    else
                            {
							    if (myCHttpFile->Read(&szWebServerResp,dwSize))
                                {
                                    fwrite(&szWebServerResp,dwSize,1,TempFile);
                                }
                            }
                            fclose(TempFile);
                        }
					}
				}
				catch(CInternetException *e)
				{
					//ptrApp->m_MainHttpServer = NULL;
				}
				myCHttpFile->Close();
				delete myCHttpFile;
			}
            m_MainHttpServer->Close();
            m_MainInternetConnection->Close();
		}
	}
    FILE *fInput = fopen(szXMLFileName, "r");
    if (fInput != NULL)
    {
        ParamDoAll(fInput);
        fclose(fInput);
    }
    else
    {
        printf("\n file %s missing", szXMLFileName);
        exit(1);
    }
	{
        Interpolate_State( dStartJD , 2 , &StateEarth );
#ifdef FIND_IMPULSE_TIME
        for (int itry = 0 ; itry < 10000; itry++)
        {
        MyTry = SolarSystem;
        MyTrySat = Sat;
		MyEngine[0]= Engine[0];MyEngine[1]= Engine[1];MyEngine[2]= Engine[2];MyEngine[3]= Engine[3];MyEngine[4]= Engine[4];MyEngine[5]= Engine[5];
        SCH_Per = 10000000000000.0;
        SCH_Apg = 0.0;
        iSCH_Apg = 0;
        iSCH_Per = 0;
        SCH_ApgPerTime = 0;


        iApog = 0;
        Apog = 0.0;
        iPerig = 0;
        Perig = 0.0;
        ApogPergTime = 0;
        PerigMoon = 100000000000.0;
        iStartLandingIteraPerSec = 0;
        idRMDelta = 0;
#endif
        StartSequence = 0;
        iSecond = 0;
        iDay = 0;
        iHour = 0;
        int iPerSec = IterPerSec;//(int)(1.0 /TimeSl);
        printf("\n iterations per sec = %d", iPerSec);
        TimeSl = 1.0 / iPerSec;
        TimeSl_2 = TimeSl*TimeSl;

#ifdef CALC_SOLAR_SYSTEM
//#define FIND_SPEED_BASED_ON_BC 1
        EarthX = SolarSystem.X[EARTH];
        EarthY = SolarSystem.Y[EARTH];
        EarthZ = SolarSystem.Z[EARTH];

        MoonX = SolarSystem.X[MOON];
        MoonY = SolarSystem.Y[MOON];
        MoonZ = SolarSystem.Z[MOON];
#endif
#ifdef TEST_RUN_CALC_YEAR
        {
            MinMaxX = (MoonX * MoonM + EarthX * EarthM) / (MoonM + EarthM);
            MinMaxY = (MoonY * MoonM + EarthY * EarthM) / (MoonM + EarthM);
            MinMaxZ = (MoonZ * MoonM + EarthZ * EarthM) / (MoonM + EarthM);
            double dDi = sqrt((MinMaxX - SolarSystem.X[SUN])* (MinMaxX - SolarSystem.X[SUN]) +
                (MinMaxY - SolarSystem.Y[SUN])* (MinMaxY - SolarSystem.Y[SUN]) +
                (MinMaxZ - SolarSystem.Z[SUN])* (MinMaxZ - SolarSystem.Z[SUN]) 
                );
            MinMaxX = ((MinMaxX - SolarSystem.X[SUN])/ dDi);
            MinMaxY = ((MinMaxY - SolarSystem.Y[SUN])/ dDi);
            MinMaxZ = ((MinMaxZ - SolarSystem.Z[SUN])/ dDi);
        }
#endif
        printf("\nStart \nx=%f y=%f z=%f  \tMx=%f y=%f z= %f\tPx=%f y=%f z= %f\n(run for %ld sec)", 
                 EarthX, 
                 EarthY, 
                 EarthZ,                  
            MoonX - EarthX, MoonY - EarthY, MoonZ - EarthZ, 
            Sat.X[0] - EarthX, Sat.Y[0] - EarthY, Sat.Z[0] - EarthZ, iTotalSec);
        
        double tDeltaEarthN = (( EarthX - StateEarth.Position[0]*1000.0)*( EarthX - StateEarth.Position[0]*1000.0) +
                                             ( EarthY - StateEarth.Position[1]*1000.0)*( EarthY - StateEarth.Position[1]*1000.0) +
                                             ( EarthZ - StateEarth.Position[2]*1000.0)*( EarthZ - StateEarth.Position[2]*1000.0)
                                            );
        printf("\nJPL EPHEMERIDES  \nx=%f y=%f z=%f ", 
            StateEarth.Position[0]*1000.0, 
            StateEarth.Position[1]*1000.0, 
            StateEarth.Position[2]*1000.0);

        for (i = 0; i < iTotalSec; i++, iSecond++)
		{
            for (j = 0; j < iPerSec; j++)
			{
#ifndef CALC_SOLAR_SYSTEM
#endif
#ifdef TEST_RUN_CALC_YEAR
                if (iFirstMinMax == 0)
                {
                    iFirstMinMax = 1;
                    FirstMinMaxX2 = MinMaxX;
                    FirstMinMaxY2 = MinMaxY;
                    FirstMinMaxZ2 = MinMaxZ;
                    maxDeltaMinMaxD = 0.0;
                    flFindMax = 1;
                    flFindMin = 0;
                    iCountMin = 1000;
                }
#endif

                if (RunOrVoidEngine(1, &Engine[0], &SolarSystem, &Sat, i, j, iPerSec, dStartJD))
                {
                    // engine is running
                }
                else
                {
                    if (StartLandingIteraPerSec != 0.0)
                    {
                        if (StartLandingIteraPerSec <= (i + j*TimeSl))
                        {
                            if (iStartLandingIteraPerSec == 0)
                            {
                                TimeSl = 0.1;
                                iPerSec = 10;
                                iStartLandingIteraPerSec = 1;
                            }
                        }
                    }

                }

                IteraSat(1, &SolarSystem, &Sat,dStartJD + ((TimeFEpoch+((long double)j)/((long double)iPerSec))/24.0/60.0/60.0));
                IteraSolarSystem(1, &SolarSystem);
                EarthX = SolarSystem.X[EARTH];
                EarthY = SolarSystem.Y[EARTH];
                EarthZ = SolarSystem.Z[EARTH];
                MoonX = SolarSystem.X[MOON];
                MoonY = SolarSystem.Y[MOON];
                MoonZ = SolarSystem.Z[MOON];
                double ProbX = Sat.X[0];
                double ProbY = Sat.Y[0];
                double ProbZ = Sat.Z[0];
                {
                    dRE = sqrt(
                        (ProbX - EarthX)*(ProbX - EarthX)+
                        (ProbY - EarthY)*(ProbY - EarthY)+
                        (ProbZ - EarthZ)*(ProbZ - EarthZ)
                        );
                    dRM = sqrt(
                        (ProbX - MoonX)*(ProbX - MoonX)+
                        (ProbY - MoonY)*(ProbY - MoonY)+
                        (ProbZ - MoonZ)*(ProbZ - MoonZ)
                        );
                    dRM0 = dRM;
#ifdef FIND_IMPULSE_TIME
                    // apogee and prerigee search each time
                    if (Apog < (dRE - EarthR))
                    {
                        Apog = (dRE - EarthR);
                        iApog = 1;
                    }
                    else
                    {
                        if (iApog)
                        {
                            if (i + j*TimeSl - ApogPergTime >100)
                            {
                                printf("\n Apog = %f km DT = %f at=%d sec", Apog/1000.0, i + j*TimeSl - ApogPergTime, i);
                                if ((EngineToOptimize == 0) || 
                                    (EngineToOptimize > 0 && EngineToOptimize < MAX_ENGINES && (Engine[EngineToOptimize].EngineDone==0) && (Engine[EngineToOptimize-1].EngineDone!=0))
                                   )
                                //if (EngineToOptimize >=0 && EngineToOptimize < MAX_ENGINES)
                                {
                                    if (Engine[EngineToOptimize].iCountApogPerig ==0) // frist apogee skip
                                        Engine[EngineToOptimize].iCountApogPerig = 1;
                                    else
                                    {
                                        // on second apogee store value
                                        Engine[EngineToOptimize].SeartchForPeriod = i + j*TimeSl - ApogPergTime;
                                        Engine[EngineToOptimize].iCountApogPerig = 2;
                                    }
                                }
                                iApog = 0;
                                Perig = (dRE - EarthR);
                                ApogPergTime = i + j*TimeSl;
                                printf(", Disatnce from a Moon %f km",dRM/1000.0);
                            }
                        }
                    }
                    if (Perig > (dRE - EarthR))
                    {
                        Perig = (dRE - EarthR);
                        iPerig = 1;

                    }
                    else
                    {
                        if (iPerig)
                        {
                            if (i + j*TimeSl - ApogPergTime >100)
                            {
                                printf("\n Perig = %f km DT = %f", Perig/1000.0, i + j*TimeSl - ApogPergTime);
                                if ((EngineToOptimize == 0) || 
                                    (EngineToOptimize > 0 && EngineToOptimize < MAX_ENGINES && (Engine[EngineToOptimize].EngineDone==0) && (Engine[EngineToOptimize-1].EngineDone!=0))
                                   )
                                //if (EngineToOptimize >=0 && EngineToOptimize < MAX_ENGINES)
                                {

                                    if (Engine[EngineToOptimize].iCountApogPerig == 2) // this will be second perigee
                                    {
                                        Engine[EngineToOptimize].SeartchForPeriod += i + j*TimeSl - ApogPergTime;
                                        //Engine[EngineToOptimize].Period = Engine[EngineToOptimize].SeartchForPeriod;
                                        Opt[iOptimizationStep].Period = Engine[EngineToOptimize].SeartchForPeriod;
                                        Engine[EngineToOptimize].iCountApogPerig = 3;
                                        if (Engine[EngineToOptimize].iCalculate == CALC_INIT_PERIOD)
                                        {
                                            if (CheckWhatnext(&Opt[0],&Sat, &Engine[0], &MyEngine[0], EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                            {
                                                exit(0);
                                            }
                                        }
                                    }
                                }
                                iPerig = 0;
                                Apog = (dRE - EarthR);
                                ApogPergTime = i + j*TimeSl;
                            }
                        }
                    }
                    if ((LastEngine >= 0) && (LastEngine <= MAX_ENGINES) && (Engine[LastEngine].EngineDone))
                    {
                        SCH_Dist = sqrt(
                            (ProbX - SolarSystem.X[Engine[LastEngine].NearBody])*(ProbX - SolarSystem.X[Engine[LastEngine].NearBody])+
                            (ProbY - SolarSystem.Y[Engine[LastEngine].NearBody])*(ProbY - SolarSystem.Y[Engine[LastEngine].NearBody])+
                            (ProbZ - SolarSystem.Z[Engine[LastEngine].NearBody])*(ProbZ - SolarSystem.Z[Engine[LastEngine].NearBody])
                            );
                        dREMV = sqrt((SolarSystem.VX[Engine[LastEngine].NearBody]-Sat.VX[0])*(SolarSystem.VX[Engine[LastEngine].NearBody]-Sat.VX[0])+
                                     (SolarSystem.VY[Engine[LastEngine].NearBody]-Sat.VY[0])*(SolarSystem.VY[Engine[LastEngine].NearBody]-Sat.VY[0])+
                                     (SolarSystem.VZ[Engine[LastEngine].NearBody]-Sat.VZ[0])*(SolarSystem.VZ[Engine[LastEngine].NearBody]-Sat.VZ[0]));
                        if (Engine[LastEngine].iCalculate == CALC_APOGEE) // value 1
                        {
                            if (SCH_Apg < SCH_Dist)
                            {
                                SCH_Apg = SCH_Dist;
                                iSCH_Apg = 1;
                            }
                            else
                            {
                                if (iSCH_Apg)
                                {
                                    printf("\n SCHApogee = %f km (%f) DT = %f at=%d sec", SCH_Apg/1000.0,
                                        (SCH_Apg - GetRadius(&SolarSystem, Engine[LastEngine].NearBody, &Sat, 0))/1000.0,
                                        i + j*TimeSl - SCH_ApgPerTime, i);
                                    iSCH_Apg = 0;
                                    SCH_Per = SCH_Apg;
                                    SCH_ApgPerTime = (int)(i + j*TimeSl);
                                    printf("====engine==i=%d====<fire at=%f=>===",i,Engine[EngineToOptimize].FireTime);
                                    // SCH_Per is a parameter for that "call"
                                    goto NextTry;
                                }
                            }
                        }
                        else if (Engine[LastEngine].iCalculate == CALC_PERIGEE) // value 0
                        {
                            if (SCH_Per > SCH_Dist)
                            {
                                SCH_Per = SCH_Dist;
                                iSCH_Per = 1;
                            }
                            else
                            {
                                if (iSCH_Per)
                                {
                                    if (SCH_Per < 300000000.0)
                                    {
                                        printf("\n SCHPerigee = %f km (%f)DT = %f at=%d sec", SCH_Per/1000.0, 
                                        (SCH_Per-GetRadius(&SolarSystem, Engine[LastEngine].NearBody, &Sat, 0))/1000.0, 
                                        i + j*TimeSl - SCH_ApgPerTime, i);
                                        iSCH_Apg = 0;
                                        SCH_Apg = SCH_Per;
                                        SCH_ApgPerTime = (int)(i + j*TimeSl);
                                        printf("====engine==i=%d====<fire at=%f=>===",i,Engine[EngineToOptimize].FireTime);
                                        // SCH_Per is a parameter for that "call"
                                        goto NextTry;
                                    }
                                }
                            }
                        }
                        else if (Engine[LastEngine].iCalculate == CALC_TARGET_PRACTICE) // value 3
                        {
                            double dREM = sqrt((MoonX - EarthX)*(MoonX - EarthX)+
                                                (MoonY - EarthY)*(MoonY - EarthY)+
                                                (MoonZ - EarthZ)*(MoonZ - EarthZ));
                            // for target practice next check is fly out of moon-earth distance anyway this is a missing target
                            if (dREM < dRE)
						    {
							    printf("\n Target practice get to Moon = %f km ", PerigMoon/1000.0);
                                printf("===i=%d===========<==%f=>==v=%f=",i,Engine[EngineToOptimize].FireTime-dRMDelta,dREMV);
                                SCH_Per = SCH_Dist;
                                // SCH_Per is a parameter for that "call"
                                goto NextTry;
						    }
                            if (SCH_Per > SCH_Dist )
                            {
                                SCH_Per = SCH_Dist;
                            }
                            else
                            {
                                if (SCH_Per < 300000000.0)
                                {
                                    printf("\n Traget practice Perigee Moon = %f km ", SCH_Per/1000.0);
                                    printf("===i==%d===========<==%f=>==v=%f=",i,Engine[EngineToOptimize].FireTime-dRMDelta,dREMV);
                                    // SCH_Per is a parameter for that "call"
                                    goto NextTry;
                                }
                            }
                        }
                        else if (Engine[LastEngine].iCalculate == CALC_TARGET_POINT)  // value 5
                        {
							
							double dREM = sqrt((MoonX - EarthX)*(MoonX - EarthX)+
                                                (MoonY - EarthY)*(MoonY - EarthY)+
                                                (MoonZ - EarthZ)*(MoonZ - EarthZ))+20.0*MoonR;
                            // for target practice next check is fly out of moon-earth distance anyway this is a missing target
                            if (dREM < dRE)
						    {
							    printf("\n Target practice get to Moon = %f km ", PerigMoon/1000.0);
                                printf("===i=%d===========<==%f=>==v=%f=",i,Engine[EngineToOptimize].FireTime-dRMDelta,dREMV);
                                SCH_Per = SCH_Dist;
                                // SCH_Per is a parameter for that "call"
                                goto NextTry;
						    }
							if (SCH_Dist < MoonR*2)
							{
								double PosXMoon = 0;
								double PosYMoon = 0;
								double PosZMoon = 0;

                                // first Longitude second latitute
                                 // dolgota,  shirota

		                        getXYZMoon(Targetlongitude,Targetlatitude,PosXMoon,PosYMoon,PosZMoon,&SolarSystem,MOON,EARTH,MoonR);
								if (SCH_Dist < MoonR)
								{
									SCH_Dist = sqrt(
				                    (ProbX - PosXMoon)*(ProbX - PosXMoon)+
						            (ProbY - PosYMoon)*(ProbY - PosYMoon)+
								    (ProbZ - PosZMoon)*(ProbZ - PosZMoon)
									);
									SCH_Per = SCH_Dist;
									printf("\n Traget practice Point on a Moon = %f km ", SCH_Per/1000.0);
									printf("===i==%d===========<==%f=>==v=%f=",i,Engine[EngineToOptimize].FireTime-dRMDelta,dREMV);
										// SCH_Per is a parameter for that "call"
										goto NextTry;
								}
		                        SCH_Dist = sqrt(
				                    (ProbX - PosXMoon)*(ProbX - PosXMoon)+
						            (ProbY - PosYMoon)*(ProbY - PosYMoon)+
								    (ProbZ - PosZMoon)*(ProbZ - PosZMoon)
									);
							}
                            if (SCH_Per > SCH_Dist )
                            {
                                SCH_Per = SCH_Dist;
								iSCH_Per = 1;
                            }
                            else
                            {
								if (iSCH_Per)
								{
									if (SCH_Per < 300000000.0)
									{
										printf("\n Traget practice Point on a Moon = %f km ", SCH_Per/1000.0);
										printf("===i==%d===========<==%f=>==v=%f=",i,Engine[EngineToOptimize].FireTime-dRMDelta,dREMV);
										// SCH_Per is a parameter for that "call"
										goto NextTry;
									}
                                }
                            }
                        }
                        else if (Engine[LastEngine].iCalculate == CALC_PERIOD) // value 6
                        {
                            if (Engine[LastEngine].iCountApogPerig == 3)
                            {
                                if (CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                {
                                    exit(0);
                                }
                            }
                        }
                        else if (Engine[LastEngine].iCalculate == CALC_AT_APOGEE_DIFF_TO_3_4_DIST) // 12
                        {
                            if (SCH_Dist > 400000000*1.25)
                                SCH_Dist = 400000000*1.25;
                            if (SCH_Apg < SCH_Dist)
                            {
                                SCH_Apg = SCH_Dist;
                                iSCH_Apg = 1;
                            }
                            else
                            {
                                if (iSCH_Apg)
                                {
                                    printf("\n SCHApogee = %f km (%f) DT = %f at=%d sec", SCH_Apg/1000.0,
                                        (SCH_Apg - GetRadius(&SolarSystem, Engine[LastEngine].NearBody, &Sat, 0))/1000.0,
                                        i + j*TimeSl - SCH_ApgPerTime, i);
                                    iSCH_Apg = 0;
                                    SCH_Per = SCH_Apg;
                                    SCH_ApgPerTime = (int)(i + j*TimeSl);
                                    SCH_Per = abs(SCH_Per - 400000000*0.75);
                                    printf("====engine==i=%d====<weight at=%f=>===",i,Engine[EngineToOptimize].Weight);
                                    // SCH_Per is a parameter for that "call"
                                    goto NextTry;
                                }
                            }
                        }
                    }

#endif
                    if (dRE < EarthR) // TBD Earth is not round!!!
                    {

                        printf("\n Landed on Earth at sec = %d", i);
                        Sat.flInUse[0] = 0;
                    }

                    if (dRM < MoonR) // TBD Moon is not round!!!
                    {
                        double LongOnMoon = 0.0; // dolgota
                        double LatiOnMoon = 0.0; // shirota
                        double PosXMoon = 0;
                        double PosYMoon = 0;
                        double PosZMoon = 0;
                        getLongLatiMoon(LongOnMoon,LatiOnMoon,&SolarSystem,MOON,EARTH,&Sat,0);
                        getXYZMoon(LongOnMoon,LatiOnMoon,PosXMoon,PosYMoon,PosZMoon,&SolarSystem,MOON,EARTH,dRM);
                        dREMV = sqrt((SolarSystem.VX[MOON]-Sat.VX[0])*(SolarSystem.VX[MOON]-Sat.VX[0])+(SolarSystem.VY[MOON]-Sat.VY[0])*(SolarSystem.VY[MOON]-Sat.VY[0])+(SolarSystem.VZ[MOON]-Sat.VZ[0])*(SolarSystem.VZ[MOON]-Sat.VZ[0]));
                        printf("\n Landed on Moon at sec = %f x=%f Y=%f z=%f V=%f", ((double)i) + TimeSl*((double)j), Sat.X[0] -SolarSystem.X[MOON], Sat.Y[0] -SolarSystem.Y[MOON], Sat.Z[0] -SolarSystem.Z[MOON],dREMV);
                        printf("\n Longitute = %f Latitute %f", LongOnMoon, LatiOnMoon);
                        printf("\n Landed weight = %f from initial = %f (%f percent)", Sat.M[0], MyTrySat.M[0],Sat.M[0]/MyTrySat.M[0]);
                        
                        
#ifdef _DO_VISUALIZATION
                        // store last image 
                        //DrawAnimationSequence(&SolarSystem,&Sat, i,"TRA",&SolarSystem, RGBReferenceBody, dRGBScale, StartSequence, 1);
                        DrawFinalBody(&SolarSystem, MOON, &Sat, i,"TRA", &SolarSystem, RGBReferenceBody, dRGBScale, StartSequence);
#endif
                        dumpXMLParam(&MyTrySat, &Engine[0],EnginesCount);
                        Sat.flInUse[0] = 0;
                        exit(0);
                    }
                }



#ifdef TEST_RUN_CALC_YEAR
                // this is for debug only - checking chnages of the orbit points
                MinMaxX = (MoonX * MoonM + EarthX * EarthM) / (MoonM + EarthM);
                MinMaxY = (MoonY * MoonM + EarthY * EarthM) / (MoonM + EarthM);
                MinMaxZ = (MoonZ * MoonM + EarthZ * EarthM) / (MoonM + EarthM);

                double dDi = sqrt(
                    (MinMaxX - SolarSystem.X[SUN])* (MinMaxX - SolarSystem.X[SUN]) +
                    (MinMaxY - SolarSystem.Y[SUN])* (MinMaxY - SolarSystem.Y[SUN]) +
                    (MinMaxZ - SolarSystem.Z[SUN])* (MinMaxZ - SolarSystem.Z[SUN]) );
                MinMaxX = ((MinMaxX - SolarSystem.X[SUN])/ dDi);
                MinMaxY = ((MinMaxY - SolarSystem.Y[SUN])/ dDi);
                MinMaxZ = ((MinMaxZ - SolarSystem.Z[SUN])/ dDi);
                if (iFirstMinMax == 1)
                {
                    double tDeltaMinMaxD = (( MinMaxX - FirstMinMaxX2)*( MinMaxX - FirstMinMaxX2) +
                                            ( MinMaxY - FirstMinMaxY2)*( MinMaxY - FirstMinMaxY2) +
                                            ( MinMaxZ - FirstMinMaxZ2)*( MinMaxZ - FirstMinMaxZ2)
                                                );
                    if (flFindMax)
                    {
                        if (maxDeltaMinMaxD  < tDeltaMinMaxD)
                        {
                            maxDeltaMinMaxD  = tDeltaMinMaxD;
                            iCountMin = 1000;
                        }
                        else
                        {
                            if (iCountMin < 1000)
                            {
                                iCountMin++;
                                maxDeltaMinMaxD  = tDeltaMinMaxD;
                            }
                            else
                            {
                                printf("\nchange(to min)=%f x=%f y=%f z=%f (d=%f) at changed sign at %ld sec + %f - %f sec ", 
                                    sqrt(maxDeltaMinMaxD),
                                    MinMaxX, MinMaxY, MinMaxZ, sqrt(MinMaxX*MinMaxX+MinMaxY*MinMaxY+MinMaxZ*MinMaxZ),
                                    i, 
                                    TimeSl*j, TimeSl*(j+1));

                                flFindMax = 0;
                                flFindMin = 1;
                                minDeltaMinMaxD = tDeltaMinMaxD;
                                iCountMin = 0;
                            }
                        }
                    }
                    else if (flFindMin)
                    {                  
                        if (minDeltaMinMaxD  > tDeltaMinMaxD)
                        {
                            minDeltaMinMaxD  = tDeltaMinMaxD;
                            iCountMin = 1000;
                        }
                        else
                        {
                            if (iCountMin < 1000)
                            {
                                iCountMin++;
                                minDeltaMinMaxD  = tDeltaMinMaxD;
                            }
                            else
                            {
                                printf("\nchange(to max)=%f x=%f y=%f z=%f (d=%f) at changed sign at %ld sec + %f - %f sec ", 
                                    sqrt(minDeltaMinMaxD),
                                    MinMaxX, MinMaxY, MinMaxZ, sqrt(MinMaxX*MinMaxX+MinMaxY*MinMaxY+MinMaxZ*MinMaxZ),
                                    i, 
                                    TimeSl*j, TimeSl*(j+1));
                                Interpolate_State( dStartJD+((double)(i+1))/(24.0*60.0*60.0)+TimeSl*((double)j) , EARTH , &StateEarth );
                                printf("\n JPL=%f %f %f ",StateEarth.Position[0],StateEarth.Position[1],StateEarth.Position[0]);
                                flFindMax = 1;
                                flFindMin = 0;
                                maxDeltaMinMaxD = tDeltaMinMaxD;
                                iCountMin = 0;
                            }
                        }
                    }
                    
                }
#endif
            }

            // this is 1 day position 
            if (iSecond >= 60*60)
            {
                iSecond = 0;
                iHour++;
                //printf("\nh=%d x=%f y=%f z=%f", iHour-1, EarthX, EarthY, EarthZ);
                if (iHour >= 24)
                {
                    iHour = 0;
                    iDay++;
                    printf("\nd=%d x=%f y=%f z=%f \tMx=%f y=%f z= %f", iDay, MinMaxX, MinMaxY, MinMaxZ, 
                        SolarSystem.X[MOON] - SolarSystem.X[EARTH], 
                        SolarSystem.Y[MOON] - SolarSystem.Y[EARTH], 
                        SolarSystem.Z[MOON] - SolarSystem.X[EARTH]);
                }
            }
            // this flag switch on/off comparation of calculated data against JPL 410
           if (flFindFirst1KmError)
           {
                Interpolate_State( dStartJD+((double)(i+1))/(24.0*60.0*60.0) , EARTH , &StateEarth );
                Interpolate_State( dStartJD+((double)(i+1))/(24.0*60.0*60.0) , MOON , &StateMoon );
                MoonX = SolarSystem.X[MOON];
                MoonY = SolarSystem.Y[MOON];
                MoonZ = SolarSystem.Z[MOON];

                EarthX = SolarSystem.X[EARTH];
                EarthY = SolarSystem.Y[EARTH];
                EarthZ = SolarSystem.Z[EARTH];
#ifdef  TEST_RUN_EARTH_ERROR
                // this error checks position ob earth-moon
                // barycentre against JPL
                double EarthBSX = (EarthX*SolarSystem.M[EARTH] + MoonX*SolarSystem.M[MOON])/(SolarSystem.M[EARTH]+SolarSystem.M[MOON]);
                double EarthBSY = (EarthY*SolarSystem.M[EARTH] + MoonY*SolarSystem.M[MOON])/(SolarSystem.M[EARTH]+SolarSystem.M[MOON]);
                double EarthBSZ = (EarthZ*SolarSystem.M[EARTH] + MoonZ*SolarSystem.M[MOON])/(SolarSystem.M[EARTH]+SolarSystem.M[MOON]);

                double EarthBSNX =  StateEarth.Position[0]*1000.0 ;
                double EarthBSNY =  StateEarth.Position[1]*1000.0 ;
                double EarthBSNZ =  StateEarth.Position[2]*1000.0 ;

#else
                // otherwise it will be Moon position
                double EarthBSX = ( MoonX - EarthX);
                double EarthBSY = ( MoonY - EarthY);
                double EarthBSZ = ( MoonZ - EarthZ);

                double EarthBSNX =  StateMoon.Position[0]*1000.0 ;
                double EarthBSNY =  StateMoon.Position[1]*1000.0 ;
                double EarthBSNZ =  StateMoon.Position[2]*1000.0 ;

#endif
                double tDeltaEarthNASA = (( EarthBSX - EarthBSNX)*( EarthBSX - EarthBSNX) +
                                             ( EarthBSY - EarthBSNY)*( EarthBSY - EarthBSNY) +
                                             ( EarthBSZ - EarthBSNZ)*( EarthBSZ - EarthBSNZ)
                                            );
                if (tDeltaEarthNASA> dErrorValue)
                {
                    // error bigger then 1000 M
                    double flX;
                    double flY;
                    double flZ;
                    MoonVX = SolarSystem.VX_[MOON]*TimeSl / SolarSystem.M[MOON];
                    MoonVY = SolarSystem.VY_[MOON]*TimeSl / SolarSystem.M[MOON];
                    MoonVZ = SolarSystem.VZ_[MOON]*TimeSl / SolarSystem.M[MOON];

                    EarthVX = SolarSystem.VX_[EARTH]*TimeSl / SolarSystem.M[EARTH];
                    EarthVX = SolarSystem.VX_[EARTH]*TimeSl / SolarSystem.M[EARTH];
                    EarthVY = SolarSystem.VY_[EARTH]*TimeSl / SolarSystem.M[EARTH];
                    EarthVZ = SolarSystem.VZ_[EARTH]*TimeSl / SolarSystem.M[EARTH];
#ifdef TEST_RUN_EARTH_ERROR
                    double EarthBSVX = EarthVX;
                    double EarthBSVY = EarthVY;
                    double EarthBSVZ = EarthVZ;

                    double EarthBSNVX =  StateEarth.Velocity[0]*1000.0 ;
                    double EarthBSNVY =  StateEarth.Velocity[1]*1000.0 ;
                    double EarthBSNVZ =  StateEarth.Velocity[2]*1000.0 ;
                    printf("\n Error in Earth position bigger then %f M", sqrt(dErrorValue));
#else
                    double EarthBSVX = ( MoonVX - EarthVX);
                    double EarthBSVY = ( MoonVY - EarthVY);
                    double EarthBSVZ = ( MoonVZ - EarthVZ);

                    double EarthBSNVX =  StateMoon.Velocity[0]*1000.0 ;
                    double EarthBSNVY =  StateMoon.Velocity[1]*1000.0 ;
                    double EarthBSNVZ =  StateMoon.Velocity[2]*1000.0 ;
                    printf("\n Error in Moon position bigger then %f M", sqrt(dErrorValue));
#endif

                    
                    dErrorValue*=2.0;
                    printf("\n=%f \nx=%f y=%f z=%f ; JPL EPHEMERIDES:\nx=%f y=%f z=%f %ld sec + %f - %f sec ", 
                                    sqrt(tDeltaEarthNASA),
                                    EarthBSX, EarthBSY, EarthBSZ, EarthBSNX, EarthBSNY, EarthBSNZ,
                                    i, 
                                    TimeSl*j, TimeSl*(j+1));
                    MoonXYZCalc(flX, flY, flZ, (dStartJD+((double)(i+1))/(24.0*60.0*60.0) - 2451544.0)/36525.0);
                    printf("\nMoon position by sin/cos approximation\n x=%f  y=%f  z=%f\nvx=%f vy=%f vz=%f ; JPL EPHEMERIDES:\nvx=%f vy=%f vz=%f ", 
                                    flX,flY,flZ,
                                    EarthBSVX, EarthBSVY, EarthBSVZ, EarthBSNVX, EarthBSNVY, EarthBSNVZ
                                    );

                    //Target0 = tDeltaEarthNASA;
                    //break;
                }
                else
                {
                    //printf("\n tttt");
                }

           }
           if (i == 31558149)
           {
               double tDeltaMinMaxD = (( MinMaxX - FirstMinMaxX2)*( MinMaxX - FirstMinMaxX2) +
                                             ( MinMaxY - FirstMinMaxY2)*( MinMaxY - FirstMinMaxY2) +
                                             ( MinMaxZ - FirstMinMaxZ2)*( MinMaxZ - FirstMinMaxZ2)
                                                );
               printf("\n=======%f x=%f y=%f z=%f (d=%f) at changed sign at %ld sec + %f - %f sec ", 
                                    sqrt(maxDeltaMinMaxD),
                                    MinMaxX, MinMaxY, MinMaxZ, sqrt(MinMaxX*MinMaxX+MinMaxY*MinMaxY+MinMaxZ*MinMaxZ),
                                    i, 
                                    TimeSl*j, TimeSl*(j+1));
           }
#ifdef _DO_VISUALIZATION
           DrawAnimationSequence(&SolarSystem,&Sat, i,"TRA",&SolarSystem, RGBReferenceBody, dRGBScale, StartSequence, 0); 
#endif
           // on first sattelite do compare of the calculated position and SGP4 
           int iCheck = 0;
           TimeFEpoch = (long double) (i+1);
            long double AE = 1.0;
            long double XKMPER = 6378.1350; //XKMPER kilometers/Earth radii 6378.135
			long double XKE = BIG_XKE;//.743669161E-1;

			long double XJ2 = 1.082616E-3;
			long double CK2=.5*XJ2*AE*AE;

			long double XMNPDA = 1440.0; // XMNPDA time units(minutes) /day 1440.0
			long double TEMP=2*M_PI/XMNPDA/XMNPDA; // 2*pi / (1440 **2)
            long double ProbMeanMotion = Sat.ProbMeanMotion[iCheck];
			long double XNO=ProbMeanMotion*TEMP*XMNPDA; // rotation per day * 2*pi /1440 == rotation per day on 1 unit (1 min)
			long double XNDT2O=Sat.ProbFirstDervMeanMotion[iCheck]*TEMP;
			long double XNDD6O=Sat.ProbSecondDervmeanMotion[iCheck]*TEMP/XMNPDA;
			long double BSTAR=Sat.ProbDragterm[iCheck]/AE;
            long double tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ;
            // next lines has to be removed ==> today they are included only to avoid drag effect
#if 1
            BSTAR = 0.0;
            XNDD6O = 0.0;
            XNDT2O =0.0;
#endif

            //(dStartJD - Sat.ProbEpoch[nSat])*XMNPDA
            // TimeFEpoch/60.0
			SGP4((dStartJD + (TimeFEpoch/24.0/60.0/60.0)- Sat.ProbEpoch[iCheck])*XMNPDA, 
                XNDT2O,XNDD6O,BSTAR,Sat.ProbIncl[iCheck], Sat.ProbAscNode[iCheck],Sat.ProbEcc[iCheck], Sat.ProbArgPer[iCheck], Sat.ProbMeanAnom[iCheck],XNO, 
				tProbX,tProbY,tProbZ,tProbVX,tProbVY,tProbVZ);
			tProbX=tProbX*XKMPER/AE*1000.0;                 tProbY=tProbY*XKMPER/AE*1000.0;                 tProbZ=tProbZ*XKMPER/AE*1000.0;
			tProbVX=tProbVX*XKMPER/AE*XMNPDA/86400.*1000.0;	tProbVY=tProbVY*XKMPER/AE*XMNPDA/86400.*1000.0;	tProbVZ=tProbVZ*XKMPER/AE*XMNPDA/86400.*1000.0;

            tX = Sat.X[iCheck] - SolarSystem.X[EARTH];      tY = Sat.Y[iCheck] - SolarSystem.Y[EARTH];      tZ = Sat.Z[iCheck] - SolarSystem.Z[EARTH];
		    tVX = Sat.VX[iCheck] - SolarSystem.VX[EARTH];   tVY = Sat.VY[iCheck] - SolarSystem.VY[EARTH];   tVZ = Sat.VZ[iCheck] - SolarSystem.VZ[EARTH];

            ttProbX	= tX - tProbX;ttProbY= tY - tProbY; ttProbZ	= tZ - tProbZ;
            ttProbVX	= tVX - tProbVX; ttProbVY	= tVY - tProbVY; ttProbVZ	= tVZ - tProbVZ;

            tttX = sqrt(ttProbX*ttProbX + ttProbY*ttProbY + ttProbZ*ttProbZ);
            tttVX = sqrt(ttProbVX*ttProbVX + ttProbVY*ttProbVY + ttProbVZ*ttProbVZ);
            double errorCos = (tVX*tProbVX + tVY*tProbVY + tVZ*tProbVZ)/
                (sqrt(tVX*tVX +tVY*tVY + tVZ*tVZ)*
                sqrt(tProbVX*tProbVX + tProbVY*tProbVY + tProbVZ*tProbVZ));
            double errAngle =  acos(errorCos);
            double errorD = sqrt(tVX*tVX + tVY*tVY + tVZ*tVZ)/sqrt(tProbVX*tProbVX + tProbVY*tProbVY + tProbVZ*tProbVZ);
            double cosCoLAtitude = tZ / sqrt(tX*tX + tY*tY + tZ*tZ);
            if (i%60 == 0)
            {
                //if (abs(cosCoLAtitude) >0.35)
                //    printf("\n%3d=%f err(%f)V=%f angle=%f d=%f ",(int)(acos(cosCoLAtitude)*180/M_PI),cosCoLAtitude,tttX,tttVX,errAngle*1000.0,(errorD-1.0)*1000.0);
                //else
                    printf("\n%3d=%f err(X=%f V=%f) min=%d ",(int)(acos(cosCoLAtitude)*180/M_PI),cosCoLAtitude,tttX,tttVX, i/60);
            }
            //Sat.X[0] = tProbX + SolarSystem.X[EARTH]; Sat.Y[0] = tProbY + SolarSystem.Y[EARTH]; Sat.Z[0] = tProbZ + SolarSystem.Z[EARTH];
            //Sat.VX[0] = tProbVX + SolarSystem.VX[EARTH]; Sat.VY[0] = tProbVY + SolarSystem.VY[EARTH]; Sat.VZ[0] = tProbVZ + SolarSystem.VZ[EARTH];

			// needs to dump data for visualization each XX sec
            if (i%60 == 0) // each min output data to XML file
			    dumpTRAvisual(i);
		}
		OutLast = TRUE;
		dumpTRAvisual(i);
		printf("\n iteration done");
		tProbTSec = 0;
		tProbEcc = 0;
        tProbIncl = 0;
        tProbAscNode = 0;
        tProbArgPer = 0;
        tProbMeanAnom = 0;
            
		tX = Sat.X[0] - SolarSystem.X[EARTH];
		tY = Sat.Y[0] - SolarSystem.Y[EARTH];
		tZ = Sat.Z[0] - SolarSystem.Z[EARTH];
		tVX = Sat.VX[0] - SolarSystem.VX[EARTH];
		tVY = Sat.VY[0] - SolarSystem.VY[EARTH];
		tVZ = Sat.VZ[0] - SolarSystem.VZ[EARTH];
			// mean amomaly on curent time
		DumpKeplers(tProbTSec, // - orbit period in sec
				    tProbEcc,             // - Eccentricity
                    tProbIncl,            // - Inclination
                    tProbAscNode,         // - Longitude of ascending node
                    tProbArgPer,          // - Argument of perihelion
                    tProbMeanAnom,        // - Mean Anomaly (degrees)
                    SolarSystem.M[EARTH],0.0,
                    tX,tY,tZ,tVX,tVY,tVZ);
		tProbIncl = tProbIncl / M_PI * 180;
		tProbAscNode = tProbAscNode / M_PI * 180;
		tProbArgPer = tProbArgPer / M_PI * 180;
		tProbMeanAnom = tProbMeanAnom / M_PI * 180;
        // NO DRAG COMPARATION:
        ttProbX	= tX - (2321904.7964419187);
        ttProbY	= tY - (-6005011.6680821748);
        ttProbZ	= tZ - (1696686.3456158256);
        ttProbVX	= tVX - (2898.9465363029231);
        ttProbVY	= tVY - (-964.7411305511611);
        ttProbVZ	= tVZ - (-7098.2955155539603);

        // no drug 1 day atget:
        //ttProbX	= tX - (2734383.4756962131);
        //ttProbY	= tY - (-6085363.2719152933);
        //ttProbZ	= tZ - (-294685.77944232832);
        //ttProbVX	= tVX - (1964.3006721124766);
        //ttProbVY	= tVY - (1175.6494914836899);
        //ttProbVZ	= tVZ - (-7357.5774951035146);

        tttX = sqrt(ttProbX*ttProbX + ttProbY*ttProbY + ttProbZ*ttProbZ);
        tttVX = sqrt(ttProbVX*ttProbVX + ttProbVY*ttProbVY + ttProbVZ*ttProbVZ);

		printf("\n tProbTSec=%f tProbEcc=%f tProbIncl=%f tProbAscNode%f tProbArgPer=%f tProbMeanAnom=%f",tProbTSec,tProbEcc,tProbIncl,tProbAscNode,tProbArgPer,tProbMeanAnom);
#ifdef FIND_IMPULSE_TIME
NextTry:
#ifdef _DO_VISUALIZATION
            // store last image 
            //DrawAnimationSequence(&SolarSystem,&Sat, i,"TRA",&SolarSystem, RGBReferenceBody, dRGBScale, StartSequence, 1);
            DrawFinalBody(&SolarSystem, MOON, &Sat, i,"TRA", &SolarSystem, RGBReferenceBody, dRGBScale, StartSequence);
            //if (++iProfile > 1)
            //    iProfile = 0;
#endif
            SolarSystem = MyTry;
            Sat = MyTrySat;
            Engine[0] = MyEngine[0];Engine[1] = MyEngine[1];Engine[2] = MyEngine[2];Engine[3] = MyEngine[3];Engine[4] = MyEngine[4];Engine[5] = MyEngine[5]; 
            if ((LastEngine >=0) && (LastEngine < MAX_ENGINES))
            {

                switch(TrajectoryOptimizationType)
                {
                case MAXIMUM_BY_TIME: // value 2
                    // first assign intial values
                    if (idRM == 0)
                    {
                        dRM2 = dRM1;
                        dRM1 = SCH_Per; // this is a parameter for optimization (has to be maximized)
                        Engine[EngineToOptimize].FireTime += Engine[EngineToOptimize].OptimizationInitialStep;
                        idRM = 1;
                        idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                    }
                    else if (idRM == 1)
                    {
                        dRM2 = dRM1;
                        dRM1 = SCH_Per;
                        if (dRM2 < dRM1)
                        {
                            printf("\n Maximum at %d %f ", itry,dRM2);
                            Engine[EngineToOptimize].FireTime += Engine[EngineToOptimize].OptimizationInitialStep;
                            idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                        }
                        else
                        {
                            if (Engine[EngineToOptimize].OptimizationFirstDirectionSwitch++)
                            {
                                Engine[EngineToOptimize].OptimizationInitialStep = Engine[EngineToOptimize].OptimizationInitialStep/Engine[EngineToOptimize].OptimizationDecCoef;
                                if (abs(Engine[EngineToOptimize].OptimizationInitialStep) < Engine[EngineToOptimize].OptimizationStop)
                                {
                                    printf("\n optimization reached max delta value");
                                    if (CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                    {
                                        exit(0);
                                    }
                                    idRM = 0;
                                    break;
                                }
                            }
                            Engine[EngineToOptimize].OptimizationInitialStep = -Engine[EngineToOptimize].OptimizationInitialStep;
                            idRM = 0;
                        }

                    }
                    break;
                case MINIMUM_BY_TIME: // value 1
                    if (idRM == 0)
                    {
                        dRM2 = dRM1;
                        dRM1 = SCH_Per; // this is a parameter for optimization (has to be minimized)
                        Engine[EngineToOptimize].FireTime += Engine[EngineToOptimize].OptimizationInitialStep;
                        idRM = 1;
                        idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                    }
                    
                    else
                    {
                        dRM2 = dRM1;
                        dRM1 = SCH_Per;
                        if (dRM2 > dRM1)
                        {
                            Engine[EngineToOptimize].FireTime += Engine[EngineToOptimize].OptimizationInitialStep;
                            idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                        }
                        else
                        {
                            if (Engine[EngineToOptimize].OptimizationFirstDirectionSwitch++)
                            {
                                Engine[EngineToOptimize].OptimizationInitialStep = Engine[EngineToOptimize].OptimizationInitialStep/Engine[EngineToOptimize].OptimizationDecCoef;
                                if (abs(Engine[EngineToOptimize].OptimizationInitialStep) < Engine[EngineToOptimize].OptimizationStop)
                                {
                                    printf("\n optimization reached min delta value");
                                    if (CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                    {
                                        exit(0);
                                    }
                                    idRM = 0;
                                    break;
                                }
                            }
                            Engine[EngineToOptimize].OptimizationInitialStep = -Engine[EngineToOptimize].OptimizationInitialStep;
							idRM = 0;
                            printf("\n Minimum at %d %f ", itry,dRM2);
                        }
                    }
                    break;
                case MINIMUM_BY_WEIGHT: //  value 6
                    // first assign intial values
                    if (idRM == 0)
                    {
                        int iIm = 0;
                        dRM2 = dRM1;
                        dRM1 = SCH_Per; // this is a parameter for optimization (has to be minimized)
                        for (iIm = 0; iIm <Engine[EngineToOptimize].iLine; iIm++)
                        {
                            Engine[EngineToOptimize].ValImpl[iIm] += Engine[EngineToOptimize].OptimizationInitialStep*
                                Engine[EngineToOptimize].ValImpl[iIm];
                        }
                        Engine[EngineToOptimize].Weight += Engine[EngineToOptimize].OptimizationInitialStep*
                            Engine[EngineToOptimize].Weight;
                        Engine[EngineToOptimize].TotalWeight += Engine[EngineToOptimize].OptimizationInitialStep*
                            Engine[EngineToOptimize].TotalWeight;
                        Engine[EngineToOptimize].TotalImpulse += Engine[EngineToOptimize].OptimizationInitialStep*
                            Engine[EngineToOptimize].TotalImpulse;
                        
                        idRM = 1;
                        idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                    }
                    else if (idRM == 1)
                    {
                        dRM2 = dRM1;
                        dRM1 = SCH_Per;
                        if (dRM2 > dRM1)
                        {
                            int iIm = 0;
                            printf("\n Minimum at %d %f ", itry,dRM2);
                            for (iIm = 0; iIm <Engine[EngineToOptimize].iLine; iIm++)
                            {
                                Engine[EngineToOptimize].ValImpl[iIm] += Engine[EngineToOptimize].OptimizationInitialStep*
                                    Engine[EngineToOptimize].ValImpl[iIm];
                            }
                            Engine[EngineToOptimize].Weight += Engine[EngineToOptimize].OptimizationInitialStep*
                                Engine[EngineToOptimize].Weight;
                            Engine[EngineToOptimize].TotalWeight += Engine[EngineToOptimize].OptimizationInitialStep*
                                Engine[EngineToOptimize].TotalWeight;

                            Engine[EngineToOptimize].TotalImpulse += Engine[EngineToOptimize].OptimizationInitialStep*
                            Engine[EngineToOptimize].TotalImpulse;
                        
                            idRMDelta += (int)Engine[EngineToOptimize].OptimizationInitialStep;
                        }
                        else
                        {
                            if (Engine[EngineToOptimize].OptimizationFirstDirectionSwitch++)
                            {
                                Engine[EngineToOptimize].OptimizationInitialStep = Engine[EngineToOptimize].OptimizationInitialStep/Engine[EngineToOptimize].OptimizationDecCoef;
                                if (abs(Engine[EngineToOptimize].OptimizationInitialStep) < Engine[EngineToOptimize].OptimizationStop)
                                {
                                    printf("\n optimization reached max delta value");
                                    if (CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                    {
                                        exit(0);
                                    }
                                    idRM = 0;
                                    break;
                                }
                            }
                            Engine[EngineToOptimize].OptimizationInitialStep = -Engine[EngineToOptimize].OptimizationInitialStep;
                            idRM = 0;
                        }

                    }
                    break;
                case MINIMUM_BY_ANGLE: //  value 4
					// find closest point to but by adjusting angle
                    if (iFirstAngleDone == 0)
                    {
                        LastStepTryAnglesDirValuesX = Engine[EngineToOptimize].XVec;
                        LastStepTryAnglesDirValuesY = Engine[EngineToOptimize].YVec;
                        LastStepTryAnglesDirValuesZ = Engine[EngineToOptimize].ZVec;
						LastStepHistX[0] = LastStepTryAnglesDirValuesX;
                        LastStepHistY[0] = LastStepTryAnglesDirValuesY;
                        LastStepHistZ[0] = LastStepTryAnglesDirValuesZ;
						//TryAnglesDirValues[0] = dRM;
                        iFirstAngleDone = 1;
                        TryAnglesDirDelta = Engine[EngineToOptimize].OptimizationInitialStep;
                    }
                    else
                    {
                        
                        TryAnglesDirValues[iTryAnglesDir] = SCH_Per; // this is a parameter for optimization (has to be minimized)
                        // attempt to check: may be next step is not nesesary- vector is (0,0,0)!? which is bad!
                        while(iTryAnglesDir + 1 <= iMaxTryAnglesDir)
                        {
                            if ((abs(LastStepTryAnglesDirValuesX + TryAnglesDir[iTryAnglesDir+1][0]*TryAnglesDirDelta) <= 0.0000001) &&
                                (abs(LastStepTryAnglesDirValuesY + TryAnglesDir[iTryAnglesDir+1][1]*TryAnglesDirDelta) <= 0.0000001) &&
                                (abs(LastStepTryAnglesDirValuesZ + TryAnglesDir[iTryAnglesDir+1][2]*TryAnglesDirDelta) <= 0.0000001))
                            {
                                // make it garante big
                                TryAnglesDirValues[iTryAnglesDir+1] = SCH_Per*1000;
                                iTryAnglesDir++;
                                continue;
                            }
                            break;
                        }
                        if (++iTryAnglesDir > iMaxTryAnglesDir)
                        {
                            // find minimum
							iFindMin = 0;
                            FindMin = TryAnglesDirValues[iFindMin];
                            
                            newTryAnglesDirDelta=TryAnglesDirDelta;
                            for(int iSe = iFindMin+1; iSe <= iMaxTryAnglesDir; iSe++)
                            {
                                if (FindMin > TryAnglesDirValues[iSe])
                                {
                                    FindMin = TryAnglesDirValues[iSe];
                                    iFindMin = iSe;
                                }
                            }
							// check to adjust TryAnglesDirDelta to be smaller
							if (iFindMin == 0)
							{
                                newTryAnglesDirDelta = TryAnglesDirDelta/Engine[EngineToOptimize].OptimizationDecCoef;
                                printf("\n<TRA:setting name=\"FireAng1\" value=\"0.0\" />");
                                printf("\n<TRA:setting name=\"FireAng2\" value=\"0.0\" />");
                                printf("\n    <TRA:setting name=\"XVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesX+ TryAnglesDir[iFindMin][0]*TryAnglesDirDelta));
                                printf("\n    <TRA:setting name=\"YVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesY+ TryAnglesDir[iFindMin][1]*TryAnglesDirDelta));
                                printf("\n    <TRA:setting name=\"ZVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesZ+ TryAnglesDir[iFindMin][2]*TryAnglesDirDelta));

                                if (newTryAnglesDirDelta < Engine[EngineToOptimize].OptimizationStop)
                                {
                                    printf("\n optimization minimum by angle done");
                                    if (CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                    {
                                        exit(0);
                                    }
                                    idRM = 0;
                                    break;
                                }
							}
							else
							{
                                int iTemp = 1;
								// check to adjust TryAnglesDirDelta to be smaller
								for (int iHistS = 0; iHistS < 10; iHistS++)
								{
									if ((LastStepHistX[iHistS] == (LastStepTryAnglesDirValuesX+ TryAnglesDir[iFindMin][0]*TryAnglesDirDelta)) &&
										(LastStepHistY[iHistS] == (LastStepTryAnglesDirValuesY+ TryAnglesDir[iFindMin][1]*TryAnglesDirDelta)) &&
										(LastStepHistZ[iHistS] == (LastStepTryAnglesDirValuesZ+ TryAnglesDir[iFindMin][2]*TryAnglesDirDelta)))
									{
                                        newTryAnglesDirDelta = TryAnglesDirDelta/Engine[EngineToOptimize].OptimizationDecCoef;
										printf("\n<TRA:setting name=\"FireAng1\" value=\"0.0\" />");
										printf("\n<TRA:setting name=\"FireAng2\" value=\"0.0\" />");
										printf("\n    <TRA:setting name=\"XVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesX+ TryAnglesDir[iFindMin][0]*TryAnglesDirDelta));
										printf("\n    <TRA:setting name=\"YVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesY+ TryAnglesDir[iFindMin][1]*TryAnglesDirDelta));
										printf("\n    <TRA:setting name=\"ZVector\" value=\"%f\" />",(LastStepTryAnglesDirValuesZ+ TryAnglesDir[iFindMin][2]*TryAnglesDirDelta));

										if (newTryAnglesDirDelta < Engine[EngineToOptimize].OptimizationStop)
										{
                                            printf("\n optimization minimum by angle done");
                                            if (iTemp = CheckWhatnext(&Opt[0], &Sat, &Engine[0], NULL, EngineToOptimize, TrajectoryOptimizationType, LastEngine, EnginesCount,i + j*TimeSl) == 0)
                                            {
                                                exit(0);
                                            }
                                            idRM = 0;
	                                    }
		                                break;
			                        }
				                }
                                if (iTemp)
                                    break;
							}
                            // adjust next step try
                            Engine[EngineToOptimize].XVec = LastStepTryAnglesDirValuesX + TryAnglesDir[iFindMin][0]*TryAnglesDirDelta;
                            Engine[EngineToOptimize].YVec = LastStepTryAnglesDirValuesY + TryAnglesDir[iFindMin][1]*TryAnglesDirDelta;
                            Engine[EngineToOptimize].ZVec = LastStepTryAnglesDirValuesZ + TryAnglesDir[iFindMin][2]*TryAnglesDirDelta;
							TryAnglesDirValues[0] = TryAnglesDirValues[iFindMin];

                            for (int iHist = 10 - 1; iHist > 0; iHist --)
                            {
                                LastStepHistX[iHist] = LastStepHistX[iHist-1];
                                LastStepHistY[iHist] = LastStepHistY[iHist-1];
                                LastStepHistZ[iHist] = LastStepHistZ[iHist-1];
                            }
                            LastStepHistX[0] = LastStepTryAnglesDirValuesX;
                            LastStepHistY[0] = LastStepTryAnglesDirValuesY;
                            LastStepHistZ[0] = LastStepTryAnglesDirValuesZ;

                            iTryAnglesDir = 1;
                                
                            LastStepTryAnglesDirValuesX = Engine[EngineToOptimize].XVec;
                            LastStepTryAnglesDirValuesY = Engine[EngineToOptimize].YVec;
                            LastStepTryAnglesDirValuesZ = Engine[EngineToOptimize].ZVec;

                            TryAnglesDirDelta = newTryAnglesDirDelta;

                            if ((abs(LastStepTryAnglesDirValuesX + TryAnglesDir[iTryAnglesDir][0]*TryAnglesDirDelta) <= 0.0000001) &&
                                (abs(LastStepTryAnglesDirValuesY + TryAnglesDir[iTryAnglesDir][1]*TryAnglesDirDelta) <= 0.0000001) &&
                                (abs(LastStepTryAnglesDirValuesZ + TryAnglesDir[iTryAnglesDir][2]*TryAnglesDirDelta) <= 0.0000001))
                            {
                                // make it garante big
                                TryAnglesDirValues[iTryAnglesDir] = SCH_Per*1000;
                                iTryAnglesDir++;
                                //continue;
                            }
                        }
                        Engine[EngineToOptimize].XVec = LastStepTryAnglesDirValuesX + TryAnglesDir[iTryAnglesDir][0]*TryAnglesDirDelta;
                        Engine[EngineToOptimize].YVec = LastStepTryAnglesDirValuesY + TryAnglesDir[iTryAnglesDir][1]*TryAnglesDirDelta;
                        Engine[EngineToOptimize].ZVec = LastStepTryAnglesDirValuesZ + TryAnglesDir[iTryAnglesDir][2]*TryAnglesDirDelta;
                    
                    }
                    break;
                }
            }

        } // end of attempts to calculate optimum time of impulse to achive min distance from Moon
#endif
        printf("\nEnd x=%f y=%f x=%f", EarthX, EarthY, EarthZ);
        printf("\nCurent second %ld", i);
	}
	return 0;
}

