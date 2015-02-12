#ifdef USE_GLOBAL
#define GLOBAL_VARIABLE
#else
#define GLOBAL_VARIABLE extern
#endif

#define _NORMALIZED_COEF 1

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

#define NPULSARS 150

#define MAX_OPTIM 30

#define IMAGE_W 1280
#define IMAGE_H 720

#define MAX_OUTPUT_TIMES 64

#define MAX_MEASURES 128

#define MAX_COEF_J 18

#define TOTAL_COEF 360

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

//long iCurSec; // current second from begining of the simulation
//int iCurPortionOfTheSecond;
//#define SIMPSON_INTEGRAL 1 
#define ADJUST(__INIT,__FORMULA,__ADDON) ldTemp=(__INIT+(__FORMULA))+__ADDON;ldTemp2=ldTemp-(__INIT +(__FORMULA));__INIT=ldTemp;__ADDON-=ldTemp2;

#define COPYKEPLER(a1,b1,c1) memset(szTempo, 0, sizeof(szTempo)); memcpy(szTempo, b1, c1);  if (szTempo[0] == ' ') {szTempo[0] = '0';if (szTempo[1] == ' ') {szTempo[1] ='0';if (szTempo[2] == ' '){szTempo[2] ='0';}}}a1 = atof(szTempo);

long double GreenwichAscensionFromTLEEpoch(long double EP, long double &preEps, long double &preTetta, long double &preZ, long double &nutEpsilon, long double &nutDFeta);
long double ConverTLEEpochDate2JulianDay(long double KeplerDate);
long double ConvertJulianDayToDateAndTime(double JulianDay, SYSTEMTIME *ThatTime);

void ConvertDateFromXML(char *pszQuo, long double &ld_TotalDays, long double &ld_dStartJD, long double &ld_dStartTLEEpoch);
void ParamSun(char *szString);
void init_tra_XML(void);
void ParamCommon(char *szString);
BOOL ParsURL(char * URLServer, int *port, char* URL,  char * szParsingName);
void ParamDoAll(FILE *fInput);
void ParamEarth(char *szString);
void ParamMoon(char *szString);
void ParamProb(char *szString);

void SUN_08 (int IYEAR,int IDAY,int IHOUR,int MIN,int ISEC,	long double &GST,long double &SLONG,long double &SRASN,long double &SDEC);
long double GreenwichAscensionFromTLEEpoch(long double EP, long double &preEps, long double &preTetta, long double &preZ, long double &nutEpsilon, long double &nutDFeta);



GLOBAL_VARIABLE long double C_S_nk[(TOTAL_COEF+3)*(TOTAL_COEF+3)][2];

GLOBAL_VARIABLE long double __AH[2][7];

GLOBAL_VARIABLE long double __LH[2][7];

GLOBAL_VARIABLE long double __BH[2][7];

GLOBAL_VARIABLE long double __AO[2][7][7];

GLOBAL_VARIABLE long double __L0[2][5][7];

GLOBAL_VARIABLE long double __BO[2][5][7];

GLOBAL_VARIABLE    long double __CH[2][7];
GLOBAL_VARIABLE    long double __CO[2][5][7];
GLOBAL_VARIABLE    long double __NO[2][3][7];
GLOBAL_VARIABLE    long double __FEO[2][1][7];
GLOBAL_VARIABLE    long double __DH[2][7];
GLOBAL_VARIABLE    long double __DO[2][5][7];
GLOBAL_VARIABLE   long double __EH[2][7];
GLOBAL_VARIABLE    long double __EO[2][9][7];
GLOBAL_VARIABLE   long double __ET[2][4][7];



GLOBAL_VARIABLE int CpuCore;

GLOBAL_VARIABLE int nk_lm_Numbers[(TOTAL_COEF+3)*(TOTAL_COEF+3)][2];

GLOBAL_VARIABLE int iItearationsPerSec; // that is "int" == IterPerSec

GLOBAL_VARIABLE long double SimLat[32];

GLOBAL_VARIABLE char EarthModelFile[1024];

GLOBAL_VARIABLE int EarthSmoothCoefStart;

#define USE_MODEL_LOAD

GLOBAL_VARIABLE int EarthModelCoefs;

GLOBAL_VARIABLE long double GM_MODEL;
GLOBAL_VARIABLE long double R0_MODEL;

GLOBAL_VARIABLE long double dStartGreenwichA;

typedef struct Measurement
{
    int NearBody;
    int TypeOfmesaure;

    long double T; // time of measure
    long double TError; // time measurement error (s)
    long double X, Y, Z; // position measure (m)
    long double H,LAT, LON; // H, LAT, LONG of the measure
    long double Err;       // error of the X,Y,Z (m) 
    long double D1, Err1;  // distance measure from XYZ to satelite and the error of the measure (m)
    long double T2, ErrT2; // time processing of the distance (PING) message on the satellite and error for this measure
    long double P1,P2,P3; // period measured by pulsar obrevation
} MEASUREMENT, *PMEASUREMENT;

GLOBAL_VARIABLE MEASUREMENT measures[MAX_MEASURES];

typedef struct tagPulsars
{
    int N;
    char Name[20];
    long double ELONG;
    long double ELAT;
    long double P0;
    long double S400mJy;

} PULSARS, *PPULSARS;

typedef struct MassPointElement
{
    long double X;
    long double Y;
	long double Z;
	long double Mp;
} MASS_POINT_ELEMENT, *PMASS_POINT_ELEMENT;

typedef struct XYZ_Split_Ponter_Var
{
    long double valX;
    long double valY;
    long double valZ;
    long double valXM;
    long double valYM;
    long double valZM;
} XYZ_SPLIT_POINTER_VAR, *PXYZ_SPLIT_POINTER_VAR;

// that is not workinh gor now:
//#define SIMPSON_INTEGRAL


typedef struct Long_Double_Intergal_Var
{
    long long nX0;

    long double Xm;
    long double Ym;
    long double Zm;

    long double X0;
    long double Y0;
    long double Z0;
    long double VX0;
    long double VY0;
    long double VZ0;
#ifndef SIMPSON_INTEGRAL
#else
    long double X6[6];
    long double Y6[6];
    long double Z6[6];

    long double x6[6];
    long double y6[6];
    long double z6[6];
    int i6;

    long double X6_h[6];
    long double Y6_h[6];
    long double Z6_h[6];

    long double X6_hh[6];
    long double Y6_hh[6];
    long double Z6_hh[6];

    long double X6_hhh[6];
    long double Y6_hhh[6];
    long double Z6_hhh[6];
#endif
    long double X;
    long double Y;
    long double Z;

    long double X_h;
    long double Y_h;
    long double Z_h;

    long double X_hh;
    long double Y_hh;
    long double Z_hh;

    long double X_hhh;
    long double Y_hhh;
    long double Z_hhh;


    long double X_temp;
    long double Y_temp;
    long double Z_temp;

    long double VX;
    long double VY;
    long double VZ;

    long double VX_h;
    long double VY_h;
    long double VZ_h;

    long double VX_hh;
    long double VY_hh;
    long double VZ_hh;

    long double VX_hhh;
    long double VY_hhh;
    long double VZ_hhh;


    long double VX_temp;
    long double VY_temp;
    long double VZ_temp;
    long CountNx;
    long CountNy;
    long CountNz;

    long CountNx_h;
    long CountNy_h;
    long CountNz_h;

    long CountNx_hh;
    long CountNy_hh;
    long CountNz_hh;


    long double x() 
    { 
           return X+X_temp+Xm; 
    };
    long double y() 
    { 
            return Y+Y_temp+Ym; 
    };
    long double z() 
    { 
            return Z+Z_temp+Zm; 
    };
    long double vx() 
    { 
           return VX+VX_temp; 
    };
    long double vy() 
    { 
            return VY+VY_temp; 
    };
    long double vz() 
    { 
            return VZ+VZ_temp; 
    };
#ifndef SIMPSON_INTEGRAL
    void Addpos(XYZ_SPLIT_POINTER_VAR *XYZdata,XYZ_SPLIT_POINTER_VAR *XYZdataVel)
    {
        X+=XYZdata->valX;  Y+=XYZdata->valY;   Z+=XYZdata->valZ;
        Xm +=XYZdata->valXM; Ym +=XYZdata->valYM; Zm +=XYZdata->valZM;
        //X0S +=XYZdata->valXS; Y0S +=XYZdata->valYS; Z0S +=XYZdata->valZS;
        //VX+=XYZdataVel->valX;  VY+=XYZdataVel->valY;   VZ+=XYZdataVel->valZ;
        
    }
    void Add(XYZ_SPLIT_POINTER_VAR *XYZdata)
    {
        X+=XYZdata->valX;  Y+=XYZdata->valY;  Z+=XYZdata->valZ;
        Xm +=XYZdata->valXM; Ym +=XYZdata->valYM; Zm +=XYZdata->valZM;
    }
#else
    void Add(XYZ_SPLIT_POINTER_VAR *XYZdata)
    {
        //nX0++;
        x6[i6] =XYZdata->valX; y6[i6] =XYZdata->valY; z6[i6] =XYZdata->valZ;
        if (++i6 == 4)
        {
            X6[1] += x6[1]; X6[2] += x6[2];
            Y6[1] += y6[1]; Y6[2] += y6[2];
            Z6[1] += z6[1]; Z6[2] += z6[2];
            X_temp =  (x6[0] + 4.0* (X6[1]+X6_h[1]+X6_hh[1]+X6_hhh[1]) + 2.0* (X6[2]+X6_h[2]+X6_hh[2]+X6_hhh[2]) + XYZdata->valX)/3.0; X = x6[3];
            Y_temp =  (y6[0] + 4.0* (Y6[1]+Y6_h[1]+Y6_hh[1]+Y6_hhh[1]) + 2.0* (Y6[2]+Y6_h[2]+Y6_hh[2]+Y6_hhh[2]) + XYZdata->valY)/3.0; Y = y6[3];
            Z_temp =  (z6[0] + 4.0* (Z6[1]+Z6_h[1]+Z6_hh[1]+Z6_hhh[1]) + 2.0* (Z6[2]+Z6_h[2]+Z6_hh[2]+Z6_hhh[2]) + XYZdata->valZ)/3.0; Z = z6[3];
            i6 = 1;
            x6[i6] =XYZdata->valX; y6[i6] =XYZdata->valY; z6[i6] =XYZdata->valZ;
            i6 = 2;
        }
        else
        {
            X+=XYZdata->valX;   Y+=XYZdata->valY;   Z+=XYZdata->valZ;
        }
    }
#endif
    void getIntegralpos(long double &valVX, long double &valVY, long double &valVZ)
    {
#ifndef SIMPSON_INTEGRAL
            valVX = (x()+X0+vx()+VX0);    
            valVY = (y()+Y0+vy()+VY0);    
            valVZ = (z()+Z0+vz()+VZ0);
#else
            valVX = (X_temp + X + X0 + X0M + X0S);    
            valVY = (Y_temp + Y + Y0 + Y0M + Y0S);    
            valVZ = (Z_temp + Z + Z0 + Z0M + Z0S);
#endif
    }

    void getIntegral(long double &valVX, long double &valVY, long double &valVZ)
    {
#ifndef SIMPSON_INTEGRAL
            valVX = (x()+X0);    
            valVY = (y()+Y0);    
            valVZ = (z()+Z0);
#else
            valVX = (X_temp+X+X0 + X0M + X0S);    
            valVY = (Y_temp+Y+Y0 + Y0M + Y0S);    
            valVZ = (Z_temp+Z+Z0 + Z0M + Z0S);
#endif
    }
    void adjustAllpos(void)
    {
#ifndef SIMPSON_INTEGRAL
        adjustXpos(10103,10163);      adjustYpos(10463,10559);      adjustZpos(10607,10667);
#else
        adjustXX(10103,10163);      adjustYY(10463,10559);      adjustZZ(10607,10667);
#endif
    }

    void adjustAll(void)
    {
#ifndef SIMPSON_INTEGRAL
        adjustX(10103,10163);      adjustY(10463,10559);      adjustZ(10607,10667);
#else
        adjustXX(10103,10163);      adjustYY(10463,10559);      adjustZZ(10607,10667);
#endif
    }
    void adjustall(void)
    {
#ifndef SIMPSON_INTEGRAL
        adjustX(10799,10883);      adjustY(11279,11423);      adjustZ(11483,11699);
#else
        adjustXX(10799,10883);      adjustYY(11279,11423);      adjustZZ(11483,11699);
#endif
    }
    void adjustX(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNx >= MaxVal)
        {
            CountNx = 0;
            ADJUST(X_h, 0, X);
            X_temp = X_h+X_hh+X_hhh;
            if (++CountNx_h >= MaxVal_h)
            {
                CountNx_h = 0;
                ADJUST(X_hh, 0, X_h);
                X_temp = X_h+X_hh+X_hhh;
                if (++CountNx_hh >= 1019)
                {
                    CountNx_hh = 0;
                    ADJUST(X_hhh, 0, X_hh);
                    X_temp = X_h+X_hh+X_hhh;
                    printf("X");
                }
            }
        }
    };
    void adjustY(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNy >= MaxVal)
        {
            CountNy = 0;
            ADJUST(Y_h, 0, Y);
            Y_temp = Y_h + Y_hh + Y_hhh;
            if (++CountNy_h >= MaxVal_h)
            {
                CountNy_h = 0;
                ADJUST(Y_hh, 0, Y_h);
                Y_temp = Y_h + Y_hh + Y_hhh;
                if (++CountNy_hh >= 1307)
                {
                    CountNy_hh = 0;
                    ADJUST(Y_hhh, 0, Y_hh);
                    Y_temp = Y_h + Y_hh + Y_hhh;
                    printf("Y");
                }
            }
        }
    };
    void adjustZ(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNz >= MaxVal)
        {
            CountNz = 0;
            ADJUST(Z_h, 0, Z);
            Z_temp=Z_h+Z_hh +Z_hhh;
            if (++CountNz_h >= MaxVal_h)
            {
                CountNz_h = 0;
                ADJUST(Z_hh, 0, Z_h);
                Z_temp=Z_h+Z_hh +Z_hhh;
                if (++CountNz_hh >= 983)
                {
                    CountNz_hh = 0;
                    ADJUST(Z_hhh, 0, Z_hh);
                    Z_temp=Z_h+Z_hh +Z_hhh;
                    printf("Z");
                }
            }
        }
    };

    void adjustXpos(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNx >= MaxVal)
        {
            CountNx = 0;
            ADJUST(X_h, 0, X);
            ADJUST(VX_h, 0, VX);
            X_temp = X_h+X_hh+X_hhh;
            VX_temp = VX_h+VX_hh+VX_hhh;
            if (++CountNx_h >= MaxVal_h)
            {
                CountNx_h = 0;
                ADJUST(X_hh, 0, X_h);
                ADJUST(VX_hh, 0, VX_h);
                X_temp = X_h+X_hh+X_hhh;
                VX_temp = VX_h+VX_hh+VX_hhh;
                if (++CountNx_hh >= 1019)
                {
                    CountNx_hh = 0;
                    ADJUST(X_hhh, 0, X_hh);
                    ADJUST(VX_hhh, 0, VX_hh);
                    X_temp = X_h+X_hh+X_hhh;
                    VX_temp = VX_h+VX_hh+VX_hhh;
                    printf("X");
                }
            }
        }
    };
    void adjustYpos(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNy >= MaxVal)
        {
            CountNy = 0;
            ADJUST(Y_h, 0, Y);
            ADJUST(VY_h, 0, VY);
            Y_temp = Y_h + Y_hh + Y_hhh;
            VY_temp = VY_h + VY_hh + VY_hhh;
            if (++CountNy_h >= MaxVal_h)
            {
                CountNy_h = 0;
                ADJUST(Y_hh, 0, Y_h);
                ADJUST(VY_hh, 0, VY_h);
                Y_temp = Y_h + Y_hh + Y_hhh;
                VY_temp = VY_h + VY_hh + VY_hhh;
                if (++CountNy_hh >= 1307)
                {
                    CountNy_hh = 0;
                    ADJUST(Y_hhh, 0, Y_hh);
                    ADJUST(VY_hhh, 0, VY_hh);
                    Y_temp = Y_h + Y_hh + Y_hhh;
                    VY_temp = VY_h + VY_hh + VY_hhh;
                    printf("Y");
                }
            }
        }
    };
    void adjustZpos(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNz >= MaxVal)
        {
            CountNz = 0;
            ADJUST(Z_h, 0, Z);
            ADJUST(VZ_h, 0, VZ);
            Z_temp=Z_h+Z_hh +Z_hhh;
            VZ_temp=VZ_h+VZ_hh +VZ_hhh;
            if (++CountNz_h >= MaxVal_h)
            {
                CountNz_h = 0;
                ADJUST(Z_hh, 0, Z_h);
                ADJUST(VZ_hh, 0, VZ_h);
                Z_temp=Z_h+Z_hh +Z_hhh;
                VZ_temp=VZ_h+VZ_hh +VZ_hhh;
                if (++CountNz_hh >= 983)
                {
                    CountNz_hh = 0;
                    ADJUST(Z_hhh, 0, Z_hh);
                    ADJUST(VZ_hhh, 0, VZ_hh);
                    Z_temp=Z_h+Z_hh +Z_hhh;
                    VZ_temp=VZ_h+VZ_hh +VZ_hhh;
                    printf("Z");
                }
            }
        }
    };

#ifdef SIMPSON_INTEGRAL
    void adjustXX(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNx >= MaxVal)
        {
            CountNx = 0;
            ADJUST(X6_h[1], 0, X6[1]);
            ADJUST(X6_h[2], 0, X6[2]);
            X_temp =  (x6[0] + 4.0* (X6[1]+X6_h[1]+X6_hh[1]+X6_hhh[1]) + 2.0* (X6[2]+X6_h[2]+X6_hh[2]+X6_hhh[2]) + x6[3])/3.0;
            if (++CountNx_h >= MaxVal_h)
            {
                CountNx_h = 0;
                ADJUST(X6_hh[1], 0, X6_h[1]);
                ADJUST(X6_hh[2], 0, X6_h[2]);
                X_temp =  (x6[0] + 4.0* (X6[1]+X6_h[1]+X6_hh[1]+X6_hhh[1]) + 2.0* (X6[2]+X6_h[2]+X6_hh[2]+X6_hhh[2]) + x6[3])/3.0;
                if (++CountNx_hh >= 1019)
                {
                    CountNx_hh = 0;
                    ADJUST(X6_hhh[1], 0, X6_hh[1]);
                    ADJUST(X6_hhh[2], 0, X6_hh[2]);
                    X_temp =  (x6[0] + 4.0* (X6[1]+X6_h[1]+X6_hh[1]+X6_hhh[1]) + 2.0* (X6[2]+X6_h[2]+X6_hh[2]+X6_hhh[2]) + x6[3])/3.0;
                    printf("X");
                }
            }
        }
    };
    void adjustYY(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNy >= MaxVal)
        {
            CountNy = 0;
            ADJUST(Y6_h[1], 0, Y6[1]);
            ADJUST(Y6_h[2], 0, Y6[2]);
            Y_temp =  (y6[0] + 4.0* (Y6[1]+Y6_h[1]+Y6_hh[1]+Y6_hhh[1]) + 2.0* (Y6[2]+Y6_h[2]+Y6_hh[2]+Y6_hhh[2]) + y6[3])/3.0;
            if (++CountNy_h >= MaxVal_h)
            {
                CountNy_h = 0;
                ADJUST(Y6_hh[1], 0, Y6_h[1]);
                ADJUST(Y6_hh[2], 0, Y6_h[2]);
                Y_temp =  (y6[0] + 4.0* (Y6[1]+Y6_h[1]+Y6_hh[1]+Y6_hhh[1]) + 2.0* (Y6[2]+Y6_h[2]+Y6_hh[2]+Y6_hhh[2]) + y6[3])/3.0;
                if (++CountNy_hh >= 1307)
                {
                    CountNy_hh = 0;
                    ADJUST(Y6_hhh[1], 0, Y6_hh[1]);
                    ADJUST(Y6_hhh[2], 0, Y6_hh[2]);
                    Y_temp =  (y6[0] + 4.0* (Y6[1]+Y6_h[1]+Y6_hh[1]+Y6_hhh[1]) + 2.0* (Y6[2]+Y6_h[2]+Y6_hh[2]+Y6_hhh[2]) + y6[3])/3.0;
                    printf("Y");
                }
            }
        }
    };
    void adjustZZ(long MaxVal, long MaxVal_h)
    {
        long double ldTemp, ldTemp2;
        if (++CountNz >= MaxVal)
        {
            CountNz = 0;
            ADJUST(Z6_h[1], 0, Z6[1]);
            ADJUST(Z6_h[2], 0, Z6[2]);
            Z_temp =  (z6[0] + 4.0* (Z6[1]+Z6_h[1]+Z6_hh[1]+Z6_hhh[1]) + 2.0* (Z6[2]+Z6_h[2]+Z6_hh[2]+Z6_hhh[2]) + z6[3])/3.0;
            if (++CountNz_h >= MaxVal_h)
            {
                CountNz_h = 0;
                ADJUST(Z6_hh[1], 0, Z6_h[1]);
                ADJUST(Z6_hh[2], 0, Z6_h[2]);
                Z_temp =  (z6[0] + 4.0* (Z6[1]+Z6_h[1]+Z6_hh[1]+Z6_hhh[1]) + 2.0* (Z6[2]+Z6_h[2]+Z6_hh[2]+Z6_hhh[2]) + z6[3])/3.0;
                if (++CountNz_hh >= 983)
                {
                    CountNz_hh = 0;
                    ADJUST(Z6_hhh[1], 0, Z6_hh[1]);
                    ADJUST(Z6_hhh[2], 0, Z6_hh[2]);
                    Z_temp =  (z6[0] + 4.0* (Z6[1]+Z6_h[1]+Z6_hh[1]+Z6_hhh[1]) + 2.0* (Z6[2]+Z6_h[2]+Z6_hh[2]+Z6_hhh[2]) + z6[3])/3.0;
                    printf("Z");
                }
            }
        }
    };
#endif
    void ZeroIntegral (void)
    {
        Xm = 0.0; Ym = 0.0; Zm = 0.0;
        nX0 = 0;
        X = 0.0;      Y = 0.0;      Z = 0.0;
        X_h =0.0;    Y_h =0.0;    Z_h =0.0;
        X_hh =0.0;   Y_hh =0.0;   Z_hh =0.0;
        X_hhh =0.0;   Y_hhh =0.0;   Z_hhh =0.0;
        VX = 0.0;      VY = 0.0;      VZ = 0.0;
        VX_h =0.0;    VY_h =0.0;    VZ_h =0.0;
        VX_hh =0.0;   VY_hh =0.0;   VZ_hh =0.0;
        VX_hhh =0.0;   VY_hhh =0.0;   VZ_hhh =0.0;

        CountNx = 0; CountNy = 0; CountNz = 0;
        CountNx_h = 0; CountNy_h = 0; CountNz_h = 0;
        CountNx_hh = 0; CountNy_hh = 0; CountNz_hh = 0;

        X0 = 0.0;      Y0 = 0.0;      Z0 = 0.0;
        VX0 = 0.0;      VY0 = 0.0;      VZ0 = 0.0;
#ifndef SIMPSON_INTEGRAL
#else
        for (int i = 0; i < 6; i++)
        {
            x6[i] = 0.0;      x6[i] = 0.0;      x6[i] = 0.0;
            X6[i] = 0.0;      Y6[i] = 0.0;      Z6[i] = 0.0;
            X6_h[i] = 0.0;    Y6_h[i] = 0.0;    Z6_h[i] = 0.0;
            X6_hh[i] = 0.0;   Y6_hh[i] = 0.0;   Z6_hh[i] = 0.0;
            X6_hhh[i] = 0.0;  Y6_hhh[i] = 0.0;  Z6_hhh[i] = 0.0;
        }
        i6 = 0;
#endif
        X_temp = 0.0; Y_temp = 0.0;   Z_temp = 0.0;
        VX_temp = 0.0; VY_temp = 0.0;   VZ_temp = 0.0;
    };
    void SetIterationCoefMult(long long itercoef)
    {
        X *=itercoef;  X_temp  *=itercoef;   X0M  *=itercoef;
        Y *=itercoef;  Y_temp  *=itercoef;   Y0M  *=itercoef;
        Z *=itercoef;  Z_temp  *=itercoef;   Z0M  *=itercoef;
        X_h *=itercoef;    Y_h *=itercoef;    Z_h *=itercoef;
        X_hh *=itercoef;   Y_hh *=itercoef;   Z_hh *=itercoef;
        X_hhh *=itercoef;   Y_hhh *=itercoef;   Z_hhh *=itercoef;
        VX *=itercoef;      VY *=itercoef;      VZ *=itercoef;
        VX_h *=itercoef;    VY_h *=itercoef;    VZ_h *=itercoef;
        VX_hh *=itercoef;   VY_hh *=itercoef;   VZ_hh *=itercoef;
        VX_hhh *=itercoef;   VY_hhh *=itercoef;   VZ_hhh *=itercoef;
#ifndef SIMPSON_INTEGRAL
#else
        for (int i = 0; i < 6; i++)
        {
            x6[i] *=itercoef;      x6[i] *=itercoef;      x6[i] *=itercoef;
            X6[i] *=itercoef;      Y6[i] *=itercoef;      Z6[i] *=itercoef;
            X6_h[i] *=itercoef;    Y6_h[i] *=itercoef;    Z6_h[i] *=itercoef;
            X6_hh[i] *=itercoef;   Y6_hh[i] *=itercoef;   Z6_hh[i] *=itercoef;
            X6_hhh[i] *=itercoef;  Y6_hhh[i] *=itercoef;  Z6_hhh[i] *=itercoef;
        }
#endif
        nX0 *=itercoef;
    }
    void SetIterationCoefDiv(long long itercoef)
    {
        X /=itercoef;  X_temp  /=itercoef;   X0M  /=itercoef;
        Y /=itercoef;  Y_temp  /=itercoef;   Y0M  /=itercoef;
        Z /=itercoef;  Z_temp  /=itercoef;   Z0M  /=itercoef;
        X_h /=itercoef;    Y_h /=itercoef;    Z_h /=itercoef;
        X_hh /=itercoef;   Y_hh /=itercoef;   Z_hh /=itercoef;
        X_hhh /=itercoef;   Y_hhh /=itercoef;   Z_hhh /=itercoef;
        VX /=itercoef;      VY /=itercoef;      VZ /=itercoef;
        VX_h /=itercoef;    VY_h /=itercoef;    VZ_h /=itercoef;
        VX_hh /=itercoef;   VY_hh /=itercoef;   VZ_hh /=itercoef;
        VX_hhh /=itercoef;   VY_hhh /=itercoef;   VZ_hhh /=itercoef;
#ifndef SIMPSON_INTEGRAL
#else
        for (int i = 0; i < 6; i++)
        {
            x6[i] /=itercoef;      x6[i] /=itercoef;      x6[i] /=itercoef;
            X6[i] /=itercoef;      Y6[i] /=itercoef;      Z6[i] /=itercoef;
            X6_h[i] /=itercoef;    Y6_h[i] /=itercoef;    Z6_h[i] /=itercoef;
            X6_hh[i] /=itercoef;   Y6_hh[i] /=itercoef;   Z6_hh[i] /=itercoef;
            X6_hhh[i] /=itercoef;  Y6_hhh[i] /=itercoef;  Z6_hhh[i] /=itercoef;
        }
#endif
        nX0 /=itercoef;
    }

} LONG_DOUBLE_INT_VAR, *PLONG_DOUBLE_INT_VAR;

typedef struct TraObj
{
    long double TimeSl;
    long double TimeSl_div_m[PLANET_COUNT];
    long double TimeSl_2;
    long double TimeSl_2_div_m[PLANET_COUNT];
    long double TimeSlOld;
    int Elem;
    int flInUse[PLANET_COUNT];
    
    LONG_DOUBLE_INT_VAR _position_[PLANET_COUNT];

    LONG_DOUBLE_INT_VAR _velosity_[PLANET_COUNT];

    long double X[PLANET_COUNT];
    long double Y[PLANET_COUNT];
    long double Z[PLANET_COUNT];
    long double VX[PLANET_COUNT];
    long double VY[PLANET_COUNT];
    long double VZ[PLANET_COUNT];

    long double DX[PLANET_COUNT];
    long double DY[PLANET_COUNT];
    long double DZ[PLANET_COUNT];
    long double DVX[PLANET_COUNT];
    long double DVY[PLANET_COUNT];
    long double DVZ[PLANET_COUNT];

    long double SX[PLANET_COUNT];
    long double SY[PLANET_COUNT];
    long double SZ[PLANET_COUNT];
    long double SVX[PLANET_COUNT];
    long double SVY[PLANET_COUNT];
    long double SVZ[PLANET_COUNT];
    long double SFX[PLANET_COUNT];
    long double SFY[PLANET_COUNT];
    long double SFZ[PLANET_COUNT];


    long long CountNx;  // attention only 1 month with 1024 iteration per 1 sec == 2 month with iteration 512 sec
    long long CountNy;
    long long CountNz;
    BOOL RunOne;

    long double FXM[PLANET_COUNT][PLANET_COUNT];
    long double FYM[PLANET_COUNT][PLANET_COUNT];
    long double FZM[PLANET_COUNT][PLANET_COUNT];

    long double Fxm[PLANET_COUNT];
    long double Fym[PLANET_COUNT];
    long double Fzm[PLANET_COUNT];

    long double X_[PLANET_COUNT];
    long double VX_[PLANET_COUNT];
    long double FX[PLANET_COUNT];
    long double X0divDt2[PLANET_COUNT];
    long double VX0divDt[PLANET_COUNT];
    long double FX0[PLANET_COUNT];
    long double FX1[PLANET_COUNT];

    long double Y_[PLANET_COUNT];
    long double VY_[PLANET_COUNT];
    long double FY[PLANET_COUNT];
    long double Y0divDt2[PLANET_COUNT];
    long double VY0divDt[PLANET_COUNT];
    long double FY0[PLANET_COUNT];
    long double FY1[PLANET_COUNT];


    long double Z_[PLANET_COUNT];
    long double VZ_[PLANET_COUNT];
    long double FZ[PLANET_COUNT];
    long double Z0divDt2[PLANET_COUNT];
    long double VZ0divDt[PLANET_COUNT];
    long double FZ0[PLANET_COUNT];
    long double FZ1[PLANET_COUNT];

    long double fx[PLANET_COUNT];
    long double fy[PLANET_COUNT];
    long double fz[PLANET_COUNT];

    long double fsinX[PLANET_COUNT];
    long double fsinY[PLANET_COUNT];
    long double fsinZ[PLANET_COUNT];

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
    //long double ProbEpochOnStart[PLANET_COUNT];
    long double ProbTLEEpoch[PLANET_COUNT];
    long double ProbJD[PLANET_COUNT];
    long double ProbJDSec[PLANET_COUNT];
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
    long double ProbSquare[PLANET_COUNT];
    // satelitte close to body
    // it can be only one body
    int iLeg;
    int iLeg_longit;
    int LegBody;
#ifdef ALL_OLD_CODE
    long double J[MAX_COEF_J];
    long double CNK[MAX_COEF_J][MAX_COEF_J];
    long double SNK[MAX_COEF_J][MAX_COEF_J];
    long double SinTetta;
    long double P[MAX_COEF_J];
    //long double Rn1divR[MAX_COEF_J];
    long double Ptilda[MAX_COEF_J];
    long double Pnk_tilda[MAX_COEF_J][MAX_COEF_J];
    //double Pnk[MAX_COEF_J][MAX_COEF_J];
    long double Qnk[MAX_COEF_J][MAX_COEF_J];
#endif

    long double ForceDD_;
    long double Lambda; // that is direction to a Greenwich
    long double precEps; // that is epsilon in precession
    long double precTet; // that is Tetta in precession
    long double precZ;  // that is Z in precession
    long double nutEpsilon;
    long double nutDFeta;
    long double h[PLANET_COUNT];
    long double ro[PLANET_COUNT];
    long double H;
    long double Ro;
    long double d;// days from the begining of the year

    //long double Sc[PLANET_COUNT];
    //long double Fp[PLANET_COUNT];
    int iAtm[PLANET_COUNT];
    long double DeltaVX[PLANET_COUNT][PLANET_COUNT];
    long double DeltaVY[PLANET_COUNT][PLANET_COUNT];
    long double DeltaVZ[PLANET_COUNT][PLANET_COUNT];

#if 1
    long double Xk[MAX_COEF_J];
    long double Yk[MAX_COEF_J];
    long double XkDxr[MAX_COEF_J];
    long double YkDxr[MAX_COEF_J];
    long double XkDyr[MAX_COEF_J];
    long double YkDyr[MAX_COEF_J];
#endif
#if 0
    long double OneMinusSinTettaInSquare;
    long double OneMinusXdivRInSquare;
    long double OneMinusYdivRInSquare;
#endif
    long double OneMinusSinTettaInSquare2;
    long double OneMinusXdivRInSquare2;
    long double OneMinusYdivRInSquare2;
    long double XdivRval;
    long double YdivRval;
    //long double OneMinusXdivRInSquare_XdivRval[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    //long double OneMinusYdivRInSquare_YdivRval[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    //long double OneMinusZdivRInSquare_ZdivRval[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    //long double OldXSign[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    //long double OldYSign[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    //long double OldZSign[PLANET_COUNT][MAX_COEF_J][MAX_COEF_J];
    long double D_Qnk_Dxr[MAX_COEF_J][MAX_COEF_J];
    long double D_Qnk_Dyr[MAX_COEF_J][MAX_COEF_J];

    long double cos_precEps, sin_precEps;
    long double cos_precTet, sin_precTet;
    long double cos_precZ, sin_precZ;
    long double cos_Lambda, sin_Lambda;
    long double cos_nutEpsilon, sin_nutEpsilon;
    long double cos_nutDFeta, sin_nutDFeta;
    long double _SQRT3;
    long double R0divR[TOTAL_COEF+3];
    long double _p_n_m_1[TOTAL_COEF+3];
    long double _p_n_m_2[TOTAL_COEF+3];
    long double _pt_nk[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tp_nm1_k[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tp_nm2_k[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tpk_n_k[TOTAL_COEF+3];
    long double diagonal[TOTAL_COEF+3];
    long double _p_n_k[TOTAL_COEF+3];
    int i_proc;
typedef struct CpuData {
    int cpuid;
    TraObj *my;
    long double *ParamSinTetta;
    long double xx, xadd, yadd,zadd;
    long double *Xk;
    long double *Yk;
    long double *R0divR;
    volatile int WaitVar;
    volatile int WaitDoVar;
} CPUDATA, *PCPUDATA;

typedef struct CpuMemory {
    long double R0divR[TOTAL_COEF+3];
    long double _p_n_m_1[TOTAL_COEF+3];
    long double _p_n_m_2[TOTAL_COEF+3];
    long double _pt_nk[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tp_nm1_k[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tp_nm2_k[TOTAL_COEF+3][TOTAL_COEF+3];
    long double _tpk_n_k[TOTAL_COEF+3];
    long double diagonal[TOTAL_COEF+3];
    long double _p_n_k[TOTAL_COEF+3];
    long double Xk[TOTAL_COEF+3];
    long double Yk[TOTAL_COEF+3];
    long double C_S_nk[(TOTAL_COEF+3)*(TOTAL_COEF+3)][2];
    long double _SQRT3;
} CPUMEMORY, *PCPUMEMORY;

    HANDLE mainThread;
    CPUMEMORY MainCpu;

    CPUDATA CPUID[32];
    int i_split[32][3];
    DWORD dwServiceStateThreadID[32];

    HANDLE		hWaitForExit[32];
    HANDLE		hWaitCmdDoCalc[32];
    HANDLE		hWaitCmdDoneCalc[32];
    HANDLE		hWaitCmdStop[32];
    HANDLE		Callback_Thread[32];

// uncomment - for thread signalling used SetEvent ; commented - for thread signalling used pooling
#define THREAD_SIGNAL SET_EVENT
#ifdef THREAD_SIGNAL
#define WAIT_THREAD_POINT while((ResWait = WaitForMultipleObjects(2,hList,FALSE,INFINITE)) != WAIT_OBJECT_0 )
#define DONE_THREAD_SIGNAL SetEvent(my->hWaitCmdDoneCalc[cpuid]);
#define CONTINUE_TO_WAIT ;
#else
#define WAIT_THREAD_POINT WAIT_AGAIN:\
        my->CPUID[cpuid].WaitDoVar = 0;\
        SetThreadPriority(my->Callback_Thread[cpuid],THREAD_PRIORITY_IDLE);\
        while(my->CPUID[cpuid].WaitDoVar == 0)\
        {\
            SwitchToThread();\
        }\
        SetThreadPriority(my->Callback_Thread[cpuid],THREAD_PRIORITY_NORMAL);\
        if (my->CPUID[cpuid].WaitDoVar == 1)


#define DONE_THREAD_SIGNAL Param->WaitVar = 1;
#define CONTINUE_TO_WAIT goto WAIT_AGAIN;
#endif

    static DWORD WINAPI CallbackThread_Proc(LPVOID lParm)
    {   
        UINT uResult = 0;
	    DWORD ResWait;
	    HANDLE hList[2];
	    BOOL bFound;
	    BOOL bRes;
	    int iRea;
        int Nb, Ne, Kb, Ke, Nall;
        CPUDATA *Param = (CPUDATA*)lParm;
        TraObj *my = Param->my;
        int cpuid = Param->cpuid;
        long double *pSinTetta = Param->ParamSinTetta;
        long double sin_Tetta;
        long double xx, xadd, yadd,zadd;
        CPUMEMORY ThreadCpuMem;
	
        hList[0] = my->hWaitCmdStop[cpuid];
	    hList[1] = my->hWaitCmdDoCalc[cpuid];

        //bRes = SetThreadAffinityMask(my->Callback_Thread[cpuid],4);

        memcpy(ThreadCpuMem.diagonal, my->diagonal, sizeof(ThreadCpuMem.diagonal));
        memcpy(ThreadCpuMem._pt_nk, my->_pt_nk, sizeof(ThreadCpuMem._pt_nk));
        memcpy(ThreadCpuMem._p_n_m_1, my->_p_n_m_1, sizeof(ThreadCpuMem._p_n_m_1));
        memcpy(ThreadCpuMem._p_n_m_2, my->_p_n_m_2, sizeof(ThreadCpuMem._p_n_m_2));
        memcpy(ThreadCpuMem._tpk_n_k, my->_tpk_n_k, sizeof(ThreadCpuMem._tpk_n_k));
        memcpy(ThreadCpuMem._tp_nm1_k, my->_tp_nm1_k, sizeof(ThreadCpuMem._tp_nm1_k));
        memcpy(ThreadCpuMem._tp_nm2_k, my->_tp_nm2_k, sizeof(ThreadCpuMem._tp_nm2_k));
        memcpy(ThreadCpuMem.C_S_nk, C_S_nk, sizeof(ThreadCpuMem.C_S_nk));
        memcpy(ThreadCpuMem._p_n_k, my->_p_n_k, sizeof(ThreadCpuMem._p_n_k));
        ThreadCpuMem._SQRT3 = my->_SQRT3;
        if (cpuid == 0)
            Nb = 2;
        else
            Nb = my->i_split[cpuid][0];

        Ne = my->iLeg;
        Kb = my->i_split[cpuid][0];
        Ke = my->i_split[cpuid][1];
        WAIT_THREAD_POINT
	    {
            sin_Tetta = *Param->ParamSinTetta;
            memcpy(ThreadCpuMem.Xk, Param->Xk,sizeof(ThreadCpuMem.Xk));
            memcpy(ThreadCpuMem.Yk, Param->Yk,sizeof(ThreadCpuMem.Yk));
            memcpy(ThreadCpuMem.R0divR, Param->R0divR,sizeof(ThreadCpuMem.R0divR));
            my->CpuPartSummXYZ ( &ThreadCpuMem, sin_Tetta,  xx, xadd, yadd, zadd, Nb , Ne, Kb, Ke);
            //my->CpuPartSummXYZ ( &ThreadCpuMem, sin_Tetta,  xx, xadd, yadd, zadd, 2 , Ne, 0, Ne);
            Param->xx = xx; Param->xadd = xadd; Param->yadd = yadd; Param->zadd = zadd;
            DONE_THREAD_SIGNAL
            CONTINUE_TO_WAIT
        }
	    SetEvent(my->hWaitForExit[cpuid]);
	    return(uResult);   
    };
    

    void StartThreads(void)
    {
        if (i_proc > 1)
        {
            mainThread = GetCurrentThread();
#ifndef THREAD_SIGNAL
            DWORD dwProcessAffinityMask;
            DWORD dwSystemAffinityMask;
            DWORD dwThreadAffinityMask = 1;
            BOOL Ret;
            GetProcessAffinityMask(GetCurrentProcess(), &dwProcessAffinityMask, &dwSystemAffinityMask);
            Ret = SetThreadAffinityMask(GetCurrentThread(), dwThreadAffinityMask);
            dwThreadAffinityMask <<= 2;
#endif
            for (int i =1; i <i_proc; i++)
            {
                hWaitForExit[i] = CreateEvent(NULL, FALSE, TRUE, NULL);
                ResetEvent(hWaitForExit[i]);
                hWaitCmdDoCalc[i] = CreateEvent(NULL, FALSE, TRUE, NULL);
                ResetEvent(hWaitCmdDoCalc[i]);
                hWaitCmdDoneCalc[i] = CreateEvent(NULL, FALSE, TRUE, NULL);
                ResetEvent(hWaitCmdDoneCalc[i]);
                hWaitCmdStop[i] = CreateEvent(NULL, FALSE, TRUE, NULL);
                ResetEvent(hWaitCmdStop[i]);

                Callback_Thread[i] = CreateThread(NULL,20980000,(LPTHREAD_START_ROUTINE)CallbackThread_Proc,(LPVOID)&CPUID[i], 0/*STACK_SIZE_PARAM_IS_A_RESERVATION*/,&dwServiceStateThreadID[i]);
#ifndef THREAD_SIGNAL
                SetThreadPriority(Callback_Thread[i],THREAD_PRIORITY_TIME_CRITICAL);
                SwitchToThread();
                if (dwThreadAffinityMask & dwProcessAffinityMask)
                {
                    Ret = SetThreadAffinityMask(Callback_Thread[i],dwThreadAffinityMask);
                    dwThreadAffinityMask <<= 2;
                }
#endif
                //SetThreadPriority(Callback_Thread[i],THREAD_PRIORITY_IDLE);
                DWORD err = GetLastError();
            }
        }
    }
    void StopThreads(void)
    {
        int i;
        if (i_proc > 1)
        {
            for (i =1; i <i_proc; i++)
            {
                SetEvent(hWaitCmdStop[i]);
                CPUID[i].WaitDoVar = 2;
                SetThreadPriority(Callback_Thread[i],THREAD_PRIORITY_TIME_CRITICAL);
            }
            CPUID[0].WaitDoVar = 2;
            
            WaitForMultipleObjects(i_proc-1,&hWaitForExit[1],TRUE,2000);
            for (i =1; i <i_proc; i++)
            {
                CloseHandle(hWaitForExit[i]);
                CloseHandle(hWaitCmdDoCalc[i]);
                CloseHandle(hWaitCmdDoneCalc[i]);
                CloseHandle(hWaitCmdStop[i]);
            }
        }
    }
    void calc_tr_matrix(void)
    {
        cos_precEps= cos(precEps); sin_precEps= sin(precEps);
        cos_precTet = cos(precTet); sin_precTet = sin(precTet);
        cos_precZ =  cos(precZ); sin_precZ = sin(precZ);
        cos_nutEpsilon = cos(nutEpsilon); sin_nutEpsilon = sin(nutEpsilon);
        cos_nutDFeta = cos(nutDFeta); sin_nutDFeta=sin(nutDFeta);
        cos_Lambda = cos(Lambda); sin_Lambda = sin(Lambda);
    }
    void gcrs_2_trs(long double &X, long double &Y, long double &Z)
    {
        long double tempX;
        long double tempY;
        long double tempZ;
        // matrix Deps
            cos_precEps= cos(precEps); sin_precEps= sin(precEps);
            //| cos(eps)   sin(eps) 0 |
            //| -sin(eps)  cos(eps) 0 |
            //|   0           0     1 |
            tempX =  cos_precEps * X + sin_precEps * Y;
            tempY = -sin_precEps * X + cos_precEps * Y;
            X = tempX; Y = tempY;
       // matrix Dtet
            cos_precTet = cos(precTet); sin_precTet = sin(precTet);
            //| cos(tet)   0 sin(tet) | 
            //|   0        1        0 |
            //|- sin(tet)  0 cos(tet) |
            tempX =  cos_precTet * X + sin_precTet * Z;
            tempZ = -sin_precTet * X + cos_precTet * Z;
            X = tempX; Z = tempZ;
       // matrix Dz
            cos_precZ =  cos(precZ); sin_precZ = sin(precZ);
            //| cos(z)  -sin(z) 0 | 
            //| sin(z)   cos(z) 0 |
            //|   0      0     1 |
            tempX =  cos_precZ * X - sin_precZ * Y;
            tempY =  sin_precZ * X + cos_precZ * Y;
            X = tempX; Y = tempY;

       // matrix C nut Epsilon
            cos_nutEpsilon = cos(nutEpsilon); sin_nutEpsilon = sin(nutEpsilon);
            tempY =    cos_nutEpsilon * Y + sin_nutEpsilon * Z;
            tempZ =  - sin_nutEpsilon * Y + cos_nutEpsilon * Z;
            Y = tempY; Z = tempZ;
            
       // matrix C nut dFeta
            cos_nutDFeta = cos(nutDFeta); sin_nutDFeta=sin(nutDFeta);
            tempX =  cos_nutDFeta * X - sin_nutDFeta * Y;
            tempY =  sin_nutDFeta * X + cos_nutDFeta * Y;
            X = tempX; Y = tempY;

       // matrix C nut Epsilon Trans
            tempY =    cos_nutEpsilon * Y - sin_nutEpsilon * Z;
            tempZ =    sin_nutEpsilon * Y + cos_nutEpsilon * Z;
            Y = tempY; Z = tempZ;

            // matrix B:
            cos_Lambda = cos(Lambda); sin_Lambda = sin(Lambda);
            //| cos(h)   sin(h) 0 |  |
            //| -sin(h)  cos(h) 0 |
            //|   0       0     1 |
            tempX = cos_Lambda * X + sin_Lambda * Y;
            tempY = -sin_Lambda * X + cos_Lambda * Y;
            X = tempX; Y = tempY;
    };
    void trs_2_gcrs(long double &X, long double &Y, long double &Z)
    {
        long double tempX;
        long double tempY;
        long double tempZ;
        // matrix (T) B:
        //| cos(h)   -sin(h) 0 |  |
        //| sin(h)   cos(h) 0 |
        //|   0        0     1 |
        tempX = cos_Lambda * X - sin_Lambda * Y;
        tempY = sin_Lambda * X + cos_Lambda * Y;
        X = tempX; Y = tempY;

       // matrix (T) C nut Epsilon Trans
            tempY =    cos_nutEpsilon * Y + sin_nutEpsilon * Z;
            tempZ =  - sin_nutEpsilon * Y + cos_nutEpsilon * Z;
            Y = tempY; Z = tempZ;

       // matrix (T)C nut dFeta
            tempX =  cos_nutDFeta * X + sin_nutDFeta * Y;
            tempY = -sin_nutDFeta * X + cos_nutDFeta * Y;
            X = tempX; Y = tempY;

       // matrix C (T) nut Epsilon
            tempY =    cos_nutEpsilon * Y - sin_nutEpsilon * Z;
            tempZ =    sin_nutEpsilon * Y + cos_nutEpsilon * Z;
            Y = tempY; Z = tempZ;

        // matrix (T) Dz
            //| cos(z)   sin(z) 0 | 
            //|-sin(z)   cos(z) 0 |
            //|   0      0     1 |
            tempX =  cos_precZ * X + sin_precZ * Y;
            tempY = -sin_precZ * X + cos_precZ * Y;
            X = tempX; Y = tempY;
        // matrix(T) Dtet
            //| cos(tet)   0 -sin(tet) | 
            //|   0        1        0  |
            //| sin(tet)   0  cos(tet) |
            tempX =  cos_precTet * X - sin_precTet * Z;
            tempZ =  sin_precTet * X + cos_precTet * Z;
            X = tempX; Z = tempZ;

        // matrix (T) Deps
            //| cos(eps)   -sin(eps) 0 |
            //| sin(eps)    cos(eps) 0 |
            //|   0           0     1 |
            tempX =  cos_precEps * X - sin_precEps * Y;
            tempY =  sin_precEps * X + cos_precEps * Y;
            X = tempX; Y = tempY;
    }
    long double GetDens(void)
    {
#define _E_CONST 2.71828182845904523536028
#if 0
#define _Ro0 1.2250
#define _M 0.0289644
#define _R_STAR 8.31432
#define __G0 9.80665
        long double Tb= 288.15;
        if (H< 11000.0)
        {
#define _LB 288.15
#define _LR -0.0065
            Tb = _LB  + _LR*H;
            Ro = _Ro0 *pow((_LB + _LR* H)/_LB, (long double)-__G0 *_M/_R_STAR/_LR-1.0);
        }
        else if (H < 20000.0)
        {
#define _LB 216.65
#undef _Ro0
#define _Ro0 0.36391802667700957
            Tb = _LB;
            Ro = _Ro0 *pow((long double)_E_CONST,(-__G0 * _M * (H-11000.0))/_R_STAR/_LB);
        }
        else if (H< 32000.0)
        {
#define _LB 216.65
#define _LR 0.001
#undef _Ro0
#define _Ro0 0.088034864330498286
            Tb = _LB  + _LR*H;
            Ro = _Ro0 *pow((_LB + _LR* (H-20000.0))/_LB, (long double)-__G0 *_M/_R_STAR/_LR-1.0);
        }
        else if (H<47000.0)
        {
#define _LB 228.65
#define _LR 0.0028
#undef _Ro0
#define _Ro0 0.013225008760250812
            Tb = _LB  + _LR*H;
            Ro = _Ro0 *pow((_LB + _LR* (H-32000.0))/_LB, (long double)-__G0 *_M/_R_STAR/_LR-1.0);
        }
        else if (H < 51000.0)
        {
#define _LB 270.65
#undef _Ro0
#define _Ro0 0.0014275334960788702
            Tb = _LB;            
            Ro = _Ro0 *pow((long double)_E_CONST,(-__G0 * _M * (H-47000.0))/_R_STAR/_LB);
        }
        else if (H<71000.0)
        {
#define _LB 270.65
#define _LR -0.0028
#undef _Ro0
#define _Ro0 0.00086160550652810688
            Tb = _LB  + _LR*H;
            Ro = _Ro0 *pow((_LB + _LR* (H-51000.00))/_LB, (long double)-__G0 *_M/_R_STAR/_LR-1.0);
        }
        else if (H < 120000.0)
        {
#define _LB 214.65
#define _LR -0.002
#undef _Ro0
#define _Ro0 6.4211030986884037e-005

            Tb = _LB  + _LR*H;
            Ro = _Ro0 *pow((_LB + _LR* (H-71000.0))/_LB, (long double)-__G0 *_M/_R_STAR/_LR-1.0);
        }
#undef _Ro0
#undef _M
#undef _R_STAR
#undef __G0
#undef _LB
#undef _LR
#else
#define _A0 (long double)1.228
#define _K1 -.090764e-3
#define _H0 0
#define _K2 -2.045e-9
#define _H0_TOP 20000.0
        if (H < _H0_TOP)
        {
            Ro = _A0 * pow((long double)_E_CONST, _K1* (H-_H0) + _K2*(H-_H0)*(H-_H0));
        }
        #undef _A0
#undef _K1
#undef _H0
#undef _K2
#undef _H0_TOP
#define _A0 (long double)9.013e-2
#define _K1 -0.16739e-3
#define _H0 20000
#define _K2 6.2669e-10
#define _H0_TOP 60000.0
        else if (H < _H0_TOP)
        {
            Ro = _A0 * pow((long double)_E_CONST, _K1* (H-_H0) + _K2*(H-_H0)*(H-_H0));
        }
#undef _A0
#undef _K1
#undef _H0
#undef _K2
#undef _H0_TOP
#define _A0 (long double)3.104e-4
#define _K1 -0.137e-3
#define _H0 60000
#define _K2 -7.8653e-10
#define _H0_TOP 100000.0
        else if (H < _H0_TOP)
        {
            Ro = _A0 * pow((long double)_E_CONST, _K1* (H-_H0) + _K2*(H-_H0)*(H-_H0));
        }
#undef _A0
#undef _K1
#undef _H0
#undef _K2
#undef _H0_TOP
#define _A0 (long double)3.66e-7
#define _K1 -0.18553e-3
#define _H0 100000
#define _K2 1.5397e-9
#define _H0_TOP 120000.0
        else if (H < _H0_TOP)
        {
            Ro = _A0 * pow((long double)_E_CONST, _K1* (H-_H0) + _K2*(H-_H0)*(H-_H0));
        }
        else if (H < 1500000.0) // 120 km
        {
            long double Hkm = H/1000.0;
#undef _A0
#undef _K1
#undef _H0
#undef _K2
#undef _H0_TOP

            // from GOST P 25645.166-2004
#define _RO_0 1.58868e-8
            long double F10_7 = 115;
            long double F81 = 115;
            long double F0 = 75.0; 
            int SA = 0;
            int Hr_for_AO = 0;
            int Hr_for_L0 = 0;
            int Hr_for_B0 = 0;
            for (SA = 0; SA < 6; SA++)
            {
                if ((F0+25.0) > F81)
                    break;
                F0 += 25.0;
            }
            if ((__AH[Hr_for_AO][SA] < Hkm) && (__AH[Hr_for_AO+1][SA] > Hkm))
                ;
            else
                Hr_for_AO = 1;
            if ((__LH[Hr_for_L0][SA] < Hkm) && (__LH[Hr_for_L0+1][SA] > Hkm))
                ;
            else
                Hr_for_L0 = 1;

            if ((__BH[Hr_for_B0][SA] < Hkm) && (__BH[Hr_for_B0+1][SA] > Hkm))
                ;
            else
                Hr_for_B0 = 1;
            

            long double RoH = _RO_0 * pow((long double)_E_CONST, __AO[Hr_for_AO][0][SA] + (__AO[Hr_for_AO][1][SA] + (__AO[Hr_for_AO][2][SA] + (__AO[Hr_for_AO][3][SA] + (__AO[Hr_for_AO][4][SA]+ (__AO[Hr_for_AO][5][SA]+ __AO[Hr_for_AO][6][SA]*Hkm)*Hkm)*Hkm)*Hkm) *Hkm)*Hkm);
            long double K0 = 1 + (__L0[Hr_for_L0][0][SA] + (__L0[Hr_for_L0][1][SA]+(__L0[Hr_for_L0][2][SA]+(__L0[Hr_for_L0][3][SA]+(__L0[Hr_for_L0][4][SA])*Hkm)*Hkm)*Hkm)*Hkm)*(F81-F0)/F0;
            long double K1 = 0;
            long double K2 = 0;
            long double K3 = 0;
            if (F10_7 != F81)
                K3 = (__BO[Hr_for_B0][0][SA] + (__BO[Hr_for_B0][1][SA]+(__BO[Hr_for_B0][2][SA]+(__BO[Hr_for_B0][3][SA]+(__BO[Hr_for_B0][4][SA])*Hkm)*Hkm)*Hkm)*Hkm)*(F10_7 - F81)/(F81+fabs(F10_7-F81));
            long double K4 = 0;
            Ro = RoH * K0 * (1 + K1 + K2 + K3 + K4);
        }
        else
            Ro = 0;

#endif
        return Ro;
    };
        // LAT  == shirota
        // LON == dolgota
        // (1) sin(LAT) = z/R  => LAT = arcsin(z/R)
        //
        // (2) sin(LON) = - X / (R * cos(LAT))
        //     LON = arcsin(- X / (R * cos(LAT)))
        //     cos(LON) = Y / (R * cos(LAT))
        //                                             AY
        //      0-> PI/2     sin(LON)<0 cos(LON)>0     | sin(LON) >0 cos(LON) >0      0 -> -PI/2
        //      -------------------------------------------------------------------------------> X
        //      PI/2->PI     sin(LON)<0 cos(LON)<0     | sin(LON)>0 cos(LON) <0    -PI/2 -> -PI

    // do not remember == from where code (below) was taken?  Some source code probably
    long double GetH(long double X, long double Y, long double Z, long double a, long double b, long double &ldLAT, long double &ldLON)
    {
        long double r = sqrt(X*X + Y*Y);
        if (r ==0.0)
        {
            ldLAT = 180.0/M_PI* asin(Z/sqrt(X*X + Y*Y + Z*Z));
            ldLON = 0;
            return fabs(Z)-b;
        }
    
        long double f = (a-b)/a;
        long double c= a/(1-f);
        long double e2 = f*(2-f);
        long double e_tilda_2 = e2/(1-e2);
        long double tg_B0 = Z/r*(1+e_tilda_2*b/sqrt(X*X+Y*Y*Z*Z));
        long double sin_;
        long double cos_;
        for (int i = 0; i<2;i++)
        {
            long double tg_tetta = (1-f)*tg_B0;
            long double tetta = atan(tg_tetta);
            sin_ = sin(tetta);
            cos_ = cos(tetta);
            long double sin_3=sin_*sin_*sin_;
            long double cos_3=cos_*cos_*cos_;
            tg_B0 = (Z+e_tilda_2*b*sin_3)/(r-e2*a*cos_3);
        }
        long double lat = atan(tg_B0);
        ldLAT = lat*180/M_PI;
        ldLON = - X / (sqrt(X*X + Y*Y + Z*Z) * cos(lat))*180/M_PI;
        long double n = c /sqrt(1.0 +e_tilda_2* cos(lat)*cos(lat));
        if (fabs(tg_B0) <=1.0)
            return r/cos(lat) -n;
        else
            return Z/sin(lat) - n*(1.0-e_tilda_2);
	};
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // geodesial to trs
    // source: http://www.astronet.ru/db/msg/1190817/node25.html
    // and http://gis-lab.info/qa/geodesic-coords.html (some error has to be fixed)
    void LatLongToTRS(long double Long,long double Lat,long double h, long double a, long double b, long double &PosX,long double &PosY,long double &PosZ)
    {
        //initial calculations
        long double cos_lat = cos(Lat* M_PI/180.0);
        long double sin_lat = sin(Lat* M_PI/180.0);
        long double f = (a-b)/a;
        long double c= a/(1-f);
        long double e2 = f*(2-f);
        long double e_tilda_2 = e2/(1-e2);
        long double n = c / sqrt(1.0 + e_tilda_2 * cos_lat *cos_lat);

        //   LAT = latitude * pi/180    // shirota
        //   LON = longitude * pi/180   // dolgota
        //   Y =  R * cos(LAT) * cos(LON)
        //   Z =  R * sin(LAT) 
        //   X = -R * cos(LAT) * sin(LON)

        //PosZ =   dRadius * sin(Lat* M_PI/180.0);
        //PosX = - dRadius * cos(Lat* M_PI/180.0) * sin(Long* M_PI/180.0); ///??????
        //PosY =   dRadius * cos(Lat* M_PI/180.0) * cos(Long* M_PI/180.0); ///??????

        PosX = (n+h) * cos_lat * cos(Long* M_PI/180.0);
        PosY = (n+h) * cos_lat * sin(Long* M_PI/180.0);
        PosZ = (n + h - e2 * n) * sin_lat;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    // geodesial to crs
    void LatLongToCRS(long double Long,long double Lat,long double h, long double a, long double b, long double curTLETime, long double &PosX,long double &PosY,long double &PosZ)
    {
        Lambda = GreenwichAscensionFromTLEEpoch(curTLETime,precEps,precTet,precZ,nutEpsilon,nutDFeta);
        LatLongToTRS(Long, Lat, h, a, b, PosX, PosY, PosZ);
        calc_tr_matrix();
        trs_2_gcrs(PosX, PosY, PosZ);

    }
#if 1

#if 1
    void CalcSplit(long iCoefs, long iCpu)
    {
        if (iCpu < 0)
        {

            DWORD_PTR dwProcessAffinityMask;
            DWORD_PTR dwSystemAffinityMask;
            DWORD_PTR dwThreadAffinityMask = {1};
            BOOL Ret;
            GetProcessAffinityMask(GetCurrentProcess(), &dwProcessAffinityMask, &dwSystemAffinityMask);
            iCpu = 0;
            while(dwThreadAffinityMask & dwProcessAffinityMask)
            {
                dwThreadAffinityMask <<=2;
                iCpu++;
            }
            CpuCore = iCpu;
            printf("\n detected %d cpu",CpuCore);
            
        }
        i_proc = 0;
        long iAll = 0;
        int i;
        long irow = iCoefs -2;
        for (int i = 0; i< iCoefs; i++)
        {
            iAll += irow;
            if (i>=2)
                irow--;
        }
        long iCount = iAll/iCpu;
        long iEachCount = 0;
        irow = iCoefs -2;
        i_split[i_proc][0] = 0;
        for (int i = 0; i< iCoefs; i++)
        {
            iEachCount += irow;
            if (i>=2)
                irow--;
            if (iEachCount >= (iAll*(i_proc+1))/iCpu)
            {
                CPUID[i_proc].cpuid = i_proc;
                CPUID[i_proc].my = this;
                CPUID[i_proc].WaitVar = 0;
                CPUID[i_proc].WaitDoVar = 0;
                i_split[i_proc][2] = iEachCount;
                i_split[i_proc++][1] = i;
                i_split[i_proc][0] = i+1;
                
            }
        }
        i_split[i_proc-1][1] = i_split[i_proc][0];
    }

    void PowerR(long double *__R0divR)
    {
        long double R0divR_ = __R0divR[1]*__R0divR[1];
        for (int n = 2; n <= iLeg; n++)
        {
            __R0divR[n] = R0divR_;
            R0divR_*= __R0divR[1];
        }
    }
    void CpuPartSummXYZ (CPUMEMORY * CpuMemory, long double sinTetta, long double &MainVal, 
        long double &Xadd, long double &Yadd, long double &Zadd, int Nstart, int Nstop, int Kstart, int Kstop)
    {
        int k;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int n = 0;  //  initial

        long double Ptilda_m_2[TOTAL_COEF+3];
        long double Ptilda_m_1[TOTAL_COEF+3];
        long double Ptilda_[TOTAL_COEF+3];
        long double *ptilda_m_2 = &Ptilda_m_2[Kstart];
        long double *ptilda_m_1 = &Ptilda_m_1[Kstart];
        long double *ptilda_=&Ptilda_[Kstart];
        int cpSize = sizeof(long double) *3;
        long double P_20_x_Q20_ = 0;
        long double Ptilda_20_x_Qnk_ = 0;
        long double C_nk_ip;
        long double S_nk_ip;

        long double P_nk_x_Qnk_;
        long double Ptilda_nk_x_Qnk_;
        long double P_nk_x_K_x_XSumD;
        long double P_nk_x_K_x_YSumD;

        long double p_nk_x_Qnk_ = 0.0;
        long double ptilda_nk_x_Qnk_ = 0.0;
        long double p_nk_x_K_x_XSumD = 0.0;
        long double p_nk_x_K_x_YSumD = 0.0;

        for (k = Kstart; k <= Kstop+1; k++) 
        {
            Ptilda_[k] = 0;  Ptilda_m_1[k] =0;  Ptilda_m_2[k]=0;
        }
        long double P_m_2 = 0;
        long double P_m_1 = 0;
        long double P_ = 1;
        Ptilda_[0]= P_;



        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // next iteration by n
        n = 1;
        P_m_2 = P_m_1; P_m_1 = P_;
        memcpy(ptilda_m_2,ptilda_m_1, cpSize); memcpy(ptilda_m_1,ptilda_, cpSize);
        //P_ = sinTetta;

#ifdef _NORMALIZED_COEF
        P_ = sinTetta*CpuMemory->_SQRT3;
#else
        P_ = sinTetta;
#endif
        Ptilda_[0]= P_;

        //Ptilda_[1] = n * P_m_1 + sinTetta * Ptilda_m_1[1]; // P'[1]  k == '

        // P = sin => d(P)/d(sin) = 1
#ifdef _NORMALIZED_COEF
        Ptilda_[1] =  CpuMemory->_SQRT3;
#else
        Ptilda_[1] =  1;
#endif

        long double R0divR_;// = R0divR[1]*R0divR[1];
        int ip = 0;

        //int iXkYk = 1;
        int Klast;
        if (Nstart != 2)
        {
            for (n= 2; n< Nstart; n++)
            {
                ip +=n+1;
                //R0divR[n] = R0divR_;
                //R0divR_*= R0divR[1];
            }
        }
        
        for (n = Nstart; n <=Nstop; n++)
        {
            R0divR_ = CpuMemory->R0divR[n];
            if (Kstop >= n)
            {
                Klast = n;
                cpSize += sizeof(long double);
            }
            else
                Klast = Kstop;
            ip+=Kstart;
            P_nk_x_Qnk_ = 0;
            Ptilda_nk_x_Qnk_ = 0;
            P_nk_x_K_x_XSumD = 0;
            P_nk_x_K_x_YSumD = 0;

            for (k = Kstart; k <=Klast; k++)
            {
                long double P_nk;
                long double Ptilda_nk;
                long double Qnk_;
                long double XSumD, YSumD;

#if _DEBUG
                // sanity check n:
                if (n != nk_lm_Numbers[ip][0])
                    exit (1);
#endif
                if (k == Kstart)
                {
                    P_m_2 = P_m_1; P_m_1 = P_;
                    memcpy(ptilda_m_2,ptilda_m_1, cpSize); memcpy(ptilda_m_1,ptilda_, cpSize);
                    if (k)
                    {
                        //k=k-1;  // one Pnk left
                        if ((k-1)==n)
                            P_nk =0;
                        else
                        {
#ifdef _NORMALIZED_COEF
                            if ((k-1) == (n-1))
                                P_nk = CpuMemory->_p_n_k[n];
                            else if ((k-1) == (n-2))
                                P_nk = CpuMemory->_tpk_n_k[n]*sinTetta;
                            else
                                P_nk = CpuMemory->_tp_nm1_k [n][(k-1)+1] * Ptilda_m_1[(k-1)+1]*sinTetta - CpuMemory->_tp_nm2_k[n][(k-1)+1] * Ptilda_m_2[(k-1)+1];
#else
                            if ((k-1) == (n-1))
                                P_nk = CpuMemory->diagonal[n];
#ifdef _DO_NOT_SKIP_OBVIOUS
                            else if ((k-1) == (n-2))
                                P_nk = CpuMemory->diagonal[n]*sinTetta;
#endif
                            else
                                P_nk = ((2*n-1) * Ptilda_m_1[(k-1)+1]*sinTetta - (n + ((k-1)+1) -1)*Ptilda_m_2[(k-1)+1])/(n-((k-1)+1));
#endif
                        }
                        Ptilda_[(k-1)+1] = P_nk; // store Pnk and rerstore original K;
                        //k++;
                    }
                    else
                    {
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // next iteration by n
#ifdef _NORMALIZED_COEF
                        P_ = CpuMemory->_p_n_m_1[n] *sinTetta * P_m_1 - CpuMemory->_p_n_m_2[n]*P_m_2;  // P[2]
#else
                        P_ = ((2.0* n-1.0) *sinTetta * P_m_1 - (n-1)*P_m_2)/n;  // P[2]
#endif
                        P_nk = P_;
                        /////////////////////////////////////////////////////////////////////////////  k =================0
                        //k = 0;
                        Ptilda_[k]= P_;
#if _DEBUG
                        // sanity check k:
                        if (k != nk_lm_Numbers[ip][1])
                            exit (1);
#endif
#ifdef _NORMALIZED_COEF
                        if (2 == n) // k==0 & n == 2
                            Ptilda_nk  = CpuMemory->_tpk_n_k[n]*sinTetta;
                        else
                            Ptilda_nk  = CpuMemory->_tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - CpuMemory->_tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
            
#else
#ifdef _DO_NOT_SKIP_OBVIOUS
                        if (2 == n) // k==0 & n == 2
                            Ptilda_nk  = CpuMemory->diagonal[n]*sinTetta;
                        else
#endif
                            Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));   // P'[2]
#endif
                        Ptilda_[k+1] = Ptilda_nk; // store P'[2] for use 
                        Qnk_ = CpuMemory->C_S_nk[ip][0] * CpuMemory->Xk[k] + CpuMemory->C_S_nk[ip][1] * CpuMemory->Yk[k];
                        // J case
                        if (ip == 0)
                        {
                            P_20_x_Q20_ = P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                            Ptilda_20_x_Qnk_ = CpuMemory->_pt_nk[n][k] *Ptilda_nk *  Qnk_;
#else
                            Ptilda_20_x_Qnk_ = Ptilda_nk *  Qnk_;
#endif
                        }
                        else
                        {
                            P_nk_x_Qnk_ += -(n+1) * P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                            Ptilda_nk_x_Qnk_ += - CpuMemory->_pt_nk[n][k] *Ptilda_nk *  Qnk_;
#else
                            Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_ ;
#endif
                        }
                        ++ip;
                        continue;
                    }
                }
                ////////////////////////////////////////////////////////////////////////////////////////
                //for (k = 1; k <=n; k++)
                {
                    ////////////////////////////////////////////////////////////////////////////////////////
                    // next iteration == k ==2
#if _DEBUG
                    // sanity check k:
                    if (k != nk_lm_Numbers[ip][1])
                        exit (1);
#endif
                    C_nk_ip = CpuMemory->C_S_nk[ip][0];
                    S_nk_ip = CpuMemory->C_S_nk[ip][1];
                    Qnk_ = C_nk_ip * CpuMemory->Xk[k] + S_nk_ip * CpuMemory->Yk[k];
                    XSumD = C_nk_ip * CpuMemory->Xk[k-1] + S_nk_ip * CpuMemory->Yk[k-1];
                    YSumD = C_nk_ip * CpuMemory->Yk[k-1] - S_nk_ip * CpuMemory->Xk[k-1];
                    P_nk = Ptilda_[k];
                    if (k==n)
                        Ptilda_nk = 0;
                    else
                    {
#ifdef _NORMALIZED_COEF
                        if (k == (n-1))
                            Ptilda_nk  = CpuMemory->_p_n_k[n];
                        else if (k == (n-2))
                            Ptilda_nk  = CpuMemory->_tpk_n_k[n]*sinTetta;
                        else
                            Ptilda_nk  = CpuMemory->_tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - CpuMemory->_tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
#else
                        if (k == (n-1))
                            Ptilda_nk = CpuMemory->diagonal[n];
#ifdef _DO_NOT_SKIP_OBVIOUS
                        else if (k == (n-2))
                            Ptilda_nk = CpuMemory->diagonal[n]*sinTetta;
#endif
                        else
                            Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));
#endif
                    }
                    Ptilda_[k+1] = Ptilda_nk;
                    P_nk_x_Qnk_ += -(n+k+1) * P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                    Ptilda_nk_x_Qnk_ += - CpuMemory->_pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
                    Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_;
#endif
                    P_nk_x_K_x_XSumD += P_nk * ( k *  XSumD   );
                    P_nk_x_K_x_YSumD += P_nk * ( k * -YSumD   );
                }
                ++ip;
            }
            ////////////////////////////////////////////////////////////////////////////////////////
            // next iteration == k ==2
            ip += n - Klast;
            p_nk_x_Qnk_      += P_nk_x_Qnk_*R0divR_;
            ptilda_nk_x_Qnk_ +=Ptilda_nk_x_Qnk_*R0divR_;
            p_nk_x_K_x_XSumD +=P_nk_x_K_x_XSumD*R0divR_;
            p_nk_x_K_x_YSumD +=P_nk_x_K_x_YSumD*R0divR_;
            //R0divR[n] = R0divR_;
            //R0divR_*= R0divR[1];
        }
        Xadd = (    + p_nk_x_K_x_XSumD ); 
        Yadd = (    + p_nk_x_K_x_YSumD );
        Zadd = (    - ptilda_nk_x_Qnk_ );
        Zadd+= (  + Ptilda_20_x_Qnk_ * (1.0))* CpuMemory->R0divR[2];
        MainVal = (p_nk_x_Qnk_ + ptilda_nk_x_Qnk_ * sinTetta ) + (-(2+1) * P_20_x_Q20_ - Ptilda_20_x_Qnk_ * sinTetta)* CpuMemory->R0divR[2] ; 
    }
    void FillXkYk(long double XdivR, long double YdivR, long double *Xk, long double *Yk)
    {
        Xk[0] = 1.0;
        Yk[0] = 0.0;
        Xk[1] = Xk[0]*XdivR - Yk[0]*YdivR;
        Yk[1] = Yk[0]*XdivR + Xk[0]*YdivR;
        for (int k = 2; k <= iLeg; k++)
        {
            Xk[k] = Xk[k-1]*XdivR - Yk[k-1]*YdivR;
            Yk[k] = Yk[k-1]*XdivR + Xk[k-1]*YdivR;
        }
    }

    void PartSummXYZ ( long double *Xk, long double *Yk, long double sinTetta, long double &MainVal, 
        long double &Xadd, long double &Yadd, long double &Zadd, int Nstart, int Nstop, int Kstart, int Kstop)
    {
        int k;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int n = 0;  //  initial

        long double Ptilda_m_2[TOTAL_COEF+3];
        long double Ptilda_m_1[TOTAL_COEF+3];
        long double Ptilda_[TOTAL_COEF+3];
        long double *ptilda_m_2 = &Ptilda_m_2[Kstart];
        long double *ptilda_m_1 = &Ptilda_m_1[Kstart];
        long double *ptilda_=&Ptilda_[Kstart];
        int cpSize = sizeof(long double) *3;
        long double P_20_x_Q20_ = 0;
        long double Ptilda_20_x_Qnk_ = 0;
        long double C_nk_ip;
        long double S_nk_ip;

        long double P_nk_x_Qnk_;
        long double Ptilda_nk_x_Qnk_;
        long double P_nk_x_K_x_XSumD;
        long double P_nk_x_K_x_YSumD;

        long double p_nk_x_Qnk_ = 0.0;
        long double ptilda_nk_x_Qnk_ = 0.0;
        long double p_nk_x_K_x_XSumD = 0.0;
        long double p_nk_x_K_x_YSumD = 0.0;

        //long double Xk[TOTAL_COEF+3];
        //long double Yk[TOTAL_COEF+3];
        //memcpy(Xk, xk, sizeof(Xk));
        //memcpy(Yk, yk, sizeof(Yk));


        for (k = Kstart; k <= Kstop+1; k++) 
        {
            Ptilda_[k] = 0;  Ptilda_m_1[k] =0;  Ptilda_m_2[k]=0;
        }
        long double P_m_2 = 0;
        long double P_m_1 = 0;
        long double P_ = 1;
        Ptilda_[0]= P_;



        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // next iteration by n
        n = 1;
        P_m_2 = P_m_1; P_m_1 = P_;
        memcpy(ptilda_m_2,ptilda_m_1, cpSize); memcpy(ptilda_m_1,ptilda_, cpSize);
        //P_ = sinTetta;

#ifdef _NORMALIZED_COEF
        P_ = sinTetta*_SQRT3;
#else
        P_ = sinTetta;
#endif
        Ptilda_[0]= P_;

        //Ptilda_[1] = n * P_m_1 + sinTetta * Ptilda_m_1[1]; // P'[1]  k == '

        // P = sin => d(P)/d(sin) = 1
#ifdef _NORMALIZED_COEF
        Ptilda_[1] =  _SQRT3;
#else
        Ptilda_[1] =  1;
#endif

        long double R0divR_;// = R0divR[1]*R0divR[1];
        int ip = 0;

        //int iXkYk = 1;
        int Klast;
        if (Nstart != 2)
        {
            for (n= 2; n< Nstart; n++)
            {
                ip +=n+1;
                //R0divR[n] = R0divR_;
                //R0divR_*= R0divR[1];
            }
        }
        
        for (n = Nstart; n <=Nstop; n++)
        {
            R0divR_ = R0divR[n];
            if (Kstop >= n)
            {
                Klast = n;
                cpSize += sizeof(long double);
            }
            else
                Klast = Kstop;
            ip+=Kstart;
            P_nk_x_Qnk_ = 0;
            Ptilda_nk_x_Qnk_ = 0;
            P_nk_x_K_x_XSumD = 0;
            P_nk_x_K_x_YSumD = 0;

            for (k = Kstart; k <=Klast; k++)
            {
                long double P_nk;
                long double Ptilda_nk;
                long double Qnk_;
                long double XSumD, YSumD;

#if _DEBUG
                // sanity check n:
                if (n != nk_lm_Numbers[ip][0])
                    exit (1);
#endif
                if (k == Kstart)
                {
                    P_m_2 = P_m_1; P_m_1 = P_;
                    memcpy(ptilda_m_2,ptilda_m_1, cpSize); memcpy(ptilda_m_1,ptilda_, cpSize);
                    if (k)
                    {
                        //k=k-1;  // one Pnk left
                        if ((k-1)==n)
                            P_nk =0;
                        else
                        {
#ifdef _NORMALIZED_COEF
                            if ((k-1) == (n-1))
                                P_nk = _p_n_k[n];
                            else if ((k-1) == (n-2))
                                P_nk = _tpk_n_k[n]*sinTetta;
                            else
                                P_nk = _tp_nm1_k [n][(k-1)+1] * Ptilda_m_1[(k-1)+1]*sinTetta - _tp_nm2_k[n][(k-1)+1] * Ptilda_m_2[(k-1)+1];
#else
                            if ((k-1) == (n-1))
                                P_nk = diagonal[n];
#ifdef _DO_NOT_SKIP_OBVIOUS
                            else if ((k-1) == (n-2))
                                P_nk = diagonal[n]*sinTetta;
#endif
                            else
                                P_nk = ((2*n-1) * Ptilda_m_1[(k-1)+1]*sinTetta - (n + ((k-1)+1) -1)*Ptilda_m_2[(k-1)+1])/(n-((k-1)+1));
#endif
                        }
                        Ptilda_[(k-1)+1] = P_nk; // store Pnk and rerstore original K;
                        //k++;
                    }
                    else
                    {
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // next iteration by n
#ifdef _NORMALIZED_COEF
                        P_ = _p_n_m_1[n] *sinTetta * P_m_1 - _p_n_m_2[n]*P_m_2;  // P[2]
#else
                        P_ = ((2.0* n-1.0) *sinTetta * P_m_1 - (n-1)*P_m_2)/n;  // P[2]
#endif
                        P_nk = P_;
                        /////////////////////////////////////////////////////////////////////////////  k =================0
                        //k = 0;
                        Ptilda_[k]= P_;
#if _DEBUG
                        // sanity check k:
                        if (k != nk_lm_Numbers[ip][1])
                            exit (1);
#endif
#ifdef _NORMALIZED_COEF
                        if (2 == n) // k==0 & n == 2
                            Ptilda_nk  = _tpk_n_k[n]*sinTetta;
                        else
                            Ptilda_nk  = _tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - _tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
            
#else
#ifdef _DO_NOT_SKIP_OBVIOUS
                        if (2 == n) // k==0 & n == 2
                            Ptilda_nk  = diagonal[n]*sinTetta;
                        else
#endif
                            Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));   // P'[2]
#endif
                        Ptilda_[k+1] = Ptilda_nk; // store P'[2] for use 
                        Qnk_ = C_S_nk[ip][0] * Xk[k] + C_S_nk[ip][1] * Yk[k];
                        // J case
                        if (ip == 0)
                        {
                            P_20_x_Q20_ = P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                            Ptilda_20_x_Qnk_ = _pt_nk[n][k] *Ptilda_nk *  Qnk_;
#else
                            Ptilda_20_x_Qnk_ = Ptilda_nk *  Qnk_;
#endif
                        }
                        else
                        {
                            P_nk_x_Qnk_ += -(n+1) * P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                            Ptilda_nk_x_Qnk_ += - _pt_nk[n][k] *Ptilda_nk *  Qnk_;
#else
                            Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_ ;
#endif
                        }
                        ++ip;
                        continue;
                    }
                }
                ////////////////////////////////////////////////////////////////////////////////////////
                //for (k = 1; k <=n; k++)
                {
                    ////////////////////////////////////////////////////////////////////////////////////////
                    // next iteration == k ==2
#if _DEBUG
                    // sanity check k:
                    if (k != nk_lm_Numbers[ip][1])
                        exit (1);
#endif
                    C_nk_ip = C_S_nk[ip][0];
                    S_nk_ip = C_S_nk[ip][1];
                    Qnk_ = C_nk_ip * Xk[k] + S_nk_ip * Yk[k];
                    XSumD = C_nk_ip * Xk[k-1] + S_nk_ip * Yk[k-1];
                    YSumD = C_nk_ip * Yk[k-1] - S_nk_ip * Xk[k-1];
                    P_nk = Ptilda_[k];
                    if (k==n)
                        Ptilda_nk = 0;
                    else
                    {
#ifdef _NORMALIZED_COEF
                        if (k == (n-1))
                            Ptilda_nk  = _p_n_k[n];
                        else if (k == (n-2))
                            Ptilda_nk  = _tpk_n_k[n]*sinTetta;
                        else
                            Ptilda_nk  = _tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - _tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
#else
                        if (k == (n-1))
                            Ptilda_nk = diagonal[n];
#ifdef _DO_NOT_SKIP_OBVIOUS
                        else if (k == (n-2))
                            Ptilda_nk = diagonal[n]*sinTetta;
#endif
                        else
                            Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));
#endif
                    }
                    Ptilda_[k+1] = Ptilda_nk;
                    P_nk_x_Qnk_ += -(n+k+1) * P_nk * Qnk_;
#ifdef _NORMALIZED_COEF
                    Ptilda_nk_x_Qnk_ += - _pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
                    Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_;
#endif
                    P_nk_x_K_x_XSumD += P_nk * ( k *  XSumD   );
                    P_nk_x_K_x_YSumD += P_nk * ( k * -YSumD   );
                }
                ++ip;
            }
            ////////////////////////////////////////////////////////////////////////////////////////
            // next iteration == k ==2
            ip += n - Klast;
            p_nk_x_Qnk_      += P_nk_x_Qnk_*R0divR_;
            ptilda_nk_x_Qnk_ +=Ptilda_nk_x_Qnk_*R0divR_;
            p_nk_x_K_x_XSumD +=P_nk_x_K_x_XSumD*R0divR_;
            p_nk_x_K_x_YSumD +=P_nk_x_K_x_YSumD*R0divR_;
            //R0divR[n] = R0divR_;
            //R0divR_*= R0divR[1];
        }
        Xadd = (    + p_nk_x_K_x_XSumD ); 
        Yadd = (    + p_nk_x_K_x_YSumD );
        Zadd = (    - ptilda_nk_x_Qnk_ );
        Zadd+= (  + Ptilda_20_x_Qnk_ * (1.0))*R0divR[2];
        MainVal = (p_nk_x_Qnk_ + ptilda_nk_x_Qnk_ * sinTetta ) + (-(2+1) * P_20_x_Q20_ - Ptilda_20_x_Qnk_ * sinTetta)*R0divR[2] ; 
    }
#if 1
    void FastSummXYZ( long double ValX, long double ValY, long double ValZ, long double ValR, long double &X, long double &Y, long double &Z, 
            long double &Xadd, long double &Yadd, long double &Zadd,
            int iCurSat)
    {
        int n,k;
        long double tempX;
        long double tempY;
        long double tempZ;
        long double sinTetta, XdivR, YdivR;
        X = 0; Y = 0; Z = 0;

        //long double _x[TOTAL_COEF][3];
        //long double _y[TOTAL_COEF][3];
        //long double _z[TOTAL_COEF][3];
        //long double _x20,_y20,_z20;
        long double P_20_x_Q20_;
        long double Ptilda_20_x_Qnk_;

        long double P_nk_x_Qnk_;//[TOTAL_COEF];
        long double Ptilda_nk_x_Qnk_;//[TOTAL_COEF];
        long double P_nk_x_K_x_XSumD;//[TOTAL_COEF];
        long double P_nk_x_K_x_YSumD;//[TOTAL_COEF];

        long double p_nk_x_Qnk_ = 0.0;
        long double ptilda_nk_x_Qnk_ = 0.0;
        long double p_nk_x_K_x_XSumD = 0.0;
        long double p_nk_x_K_x_YSumD = 0.0;

        tempX = ValX; tempY = ValY; tempZ = ValZ;
        gcrs_2_trs(tempX, tempY, tempZ);
        // now earth in Terra Ref System
        // it is possible to calculate H to get air drag
        if (--iAtm[iCurSat] == 0)
        {
            long double dlLAT, dlLON;
            iAtm[iCurSat] = iItearationsPerSec;//479; // onc per 1000 iteration == 1 per sec
            h[iCurSat] = H =GetH(tempX, tempY, tempZ, 6378245.000, 6356863.019,dlLAT,dlLON);
            ro[iCurSat] = Ro=GetDens(); // 15C
        }

        sinTetta =tempZ/ValR;

        XdivR =   tempX/ValR;
        YdivR =   tempY/ValR;
        XdivRval = XdivR;
        YdivRval = YdivR;
#ifdef ALL_OLD_CODE
        SinTetta = sinTetta;
#endif

        // loop iteration starts from n=2 k = 0
        // formula 8 on page 92
        long double Xk[TOTAL_COEF+3];
        long double Yk[TOTAL_COEF+3];
#if 0
        FillXkYk(XdivR, YdivR, Xk, Yk);
        PowerR(R0divR);
        if (i_proc == 0)
        {
            PartSummXYZ (Xk,Yk, sinTetta,  X, Xadd, Yadd, Zadd, 2, iLeg, 0, iLeg);

                           Y=X;            Z=X;
            X=1-X;         Y=1-Y;          Z=1-Z;
            Xadd = -Xadd;  Yadd = -Yadd;   Zadd = -Zadd;
        }

#else
        FillXkYk(XdivR, YdivR, MainCpu.Xk, MainCpu.Yk);
        MainCpu.R0divR[0] = R0divR[0];
        MainCpu.R0divR[1] = R0divR[1];
        PowerR(MainCpu.R0divR);
        if (i_proc <= 1)
        {
            CpuPartSummXYZ (&MainCpu, sinTetta,  X, Xadd, Yadd, Zadd, 2, iLeg, 0, iLeg);
                           Y=X;            Z=X;
            X=-X;         Y=-Y;          Z=-Z;
            Xadd = -Xadd;  Yadd = -Yadd;   Zadd = -Zadd;
        }
#endif 
        else
        {
#define CALC_VIA_THREADS
#ifdef CALC_VIA_THREADS

#ifdef THREAD_SIGNAL
#define KICK_THREAD_SIGNAL SetEvent(hWaitCmdDoCalc[ipr]);
#define KICK_ALL_THREAD_SIGNAL ;
#define WAIT_ALL_THREAD_DONE WaitForMultipleObjects(i_proc-1,&hWaitCmdDoneCalc[1],TRUE,INFINITE);
#else
#define KICK_THREAD_SIGNAL CPUID[ipr].WaitDoVar = 1;\
                SetThreadPriority(Callback_Thread[ipr],THREAD_PRIORITY_TIME_CRITICAL);
#define KICK_ALL_THREAD_SIGNAL CPUID[0].WaitDoVar = 1;
#define WAIT_ALL_THREAD_DONE    int iDoneCound = 0;\
            for (ipr = 1; ipr< i_proc; ipr++)\
            {\
                iDoneCound += CPUID[ipr].WaitVar;\
            }\
            SetThreadPriority(mainThread,THREAD_PRIORITY_IDLE);\
            while(iDoneCound != (i_proc-1))\
            {\
                SwitchToThread(); \
                iDoneCound = 0;\
                for (ipr = 1; ipr< i_proc; ipr++)\
                {\
                    iDoneCound += CPUID[ipr].WaitVar;\
                }\
            }\
            SetThreadPriority(mainThread,THREAD_PRIORITY_NORMAL);\
            for (ipr = 0; ipr< i_proc; ipr++)\
                CPUID[ipr].WaitVar = 0; 
#endif
            int ipr;
            
            
            for (ipr = 1; ipr< i_proc; ipr++)
            {
                CPUID[ipr].ParamSinTetta = &sinTetta; CPUID[ipr].Xk =  MainCpu.Xk; CPUID[ipr].Yk = MainCpu.Yk;
                CPUID[ipr].R0divR = MainCpu.R0divR;
                KICK_THREAD_SIGNAL
            }
            //SetThreadPriority(mainThread,THREAD_PRIORITY_TIME_CRITICAL);
            CpuPartSummXYZ (&MainCpu, sinTetta,  X, Xadd, Yadd, Zadd, 2, iLeg, 0, i_split[0][1]);
            
            WAIT_ALL_THREAD_DONE
            
            for (ipr = 1; ipr< i_proc; ipr++)
            {
                X+= CPUID[ipr].xx;
                Xadd += CPUID[ipr].xadd;  Yadd += CPUID[ipr].yadd;  Zadd += CPUID[ipr].zadd;
            }
                           Y=X;            Z=X;
            X=-X;         Y=-Y;          Z=-Z;
#else
            long double xx[4], xadd[4], yadd[4],zadd[4];
            PowerR(R0divR);
            // i_split[0][0] - i_split[0][1]
            PartSummXYZ ( MainCpu.Xk, MainCpu.Yk, sinTetta,  xx[0], xadd[0], yadd[0], zadd[0], 2            , iLeg, i_split[0][0], i_split[0][1]);
            PartSummXYZ ( MainCpu.Xk, MainCpu.Yk, sinTetta,  xx[1], xadd[1], yadd[1], zadd[1], i_split[1][0], iLeg, i_split[1][0], i_split[1][1]);
            PartSummXYZ ( MainCpu.Xk, MainCpu.Yk, sinTetta,  xx[2], xadd[2], yadd[2], zadd[2], i_split[2][0], iLeg, i_split[2][0], i_split[2][1]);
            PartSummXYZ ( MainCpu.Xk, MainCpu.Yk, sinTetta,  xx[3], xadd[3], yadd[3], zadd[3], i_split[3][0], iLeg, i_split[3][0], iLeg         );
            X = xx[0]+xx[1]+xx[2]+xx[3];
                           Y=X;            Z=X;
            X=-X;         Y=-Y;          Z=-Z;
            Xadd = xadd[0]+xadd[1]+xadd[2]+xadd[3];
            Yadd = yadd[0]+yadd[1]+yadd[2]+yadd[3];
            Zadd = zadd[0]+zadd[1]+zadd[2]+zadd[3];
#endif
            Xadd = -Xadd;  Yadd = -Yadd;   Zadd = -Zadd;
        }
        trs_2_gcrs(Xadd, Yadd, Zadd);
    };
#else
    // outdated : just for historic reference (to keep in source code instead in VSS)
    void FastSummXYZ( long double ValX, long double ValY, long double ValZ, long double ValR, long double &X, long double &Y, long double &Z, 
            long double &Xadd, long double &Yadd, long double &Zadd,
            int iCurSat)
    {
        int n,k;
        long double tempX;
        long double tempY;
        long double tempZ;
        long double sinTetta, XdivR, YdivR;
        X = 0; Y = 0; Z = 0;
        //long double _x[TOTAL_COEF][3];
        //long double _y[TOTAL_COEF][3];
        //long double _z[TOTAL_COEF][3];
        //long double _x20,_y20,_z20;
        long double P_20_x_Q20_;
        long double Ptilda_20_x_Qnk_;

        long double P_nk_x_Qnk_;//[TOTAL_COEF];
        long double Ptilda_nk_x_Qnk_;//[TOTAL_COEF];
        long double P_nk_x_K_x_XSumD;//[TOTAL_COEF];
        long double P_nk_x_K_x_YSumD;//[TOTAL_COEF];

        long double p_nk_x_Qnk_ = 0.0;
        long double ptilda_nk_x_Qnk_ = 0.0;
        long double p_nk_x_K_x_XSumD = 0.0;
        long double p_nk_x_K_x_YSumD = 0.0;

        tempX = ValX; tempY = ValY; tempZ = ValZ;
        gcrs_2_trs(tempX, tempY, tempZ);
        // now earth in Terra Ref System
        // it is possible to calculate H to get air drag
        if (--iAtm[iCurSat] == 0)
        {
            long double dlLAT, dlLON;
            iAtm[iCurSat] = iItearationsPerSec;//479; // onc per 1000 iteration == 1 per sec
            h[iCurSat] = H =GetH(tempX, tempY, tempZ, 6378245.000, 6356863.019,dlLAT,dlLON);
            ro[iCurSat] = Ro=GetDens(); // 15C
        }

        sinTetta =tempZ/ValR;

        XdivR =   tempX/ValR;
        YdivR =   tempY/ValR;
        XdivRval = XdivR;
        YdivRval = YdivR;

        SinTetta = sinTetta;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        n = 0;  //  initial

        long double Ptilda_m_2[TOTAL_COEF+3];
        long double Ptilda_m_1[TOTAL_COEF+3];
        long double Ptilda_[TOTAL_COEF+3];
        for (k = 0; k < TOTAL_COEF; k++) 
        {
            Ptilda_[k] = 0;  Ptilda_m_1[k] =0;  Ptilda_m_2[k]=0;
        }
        long double P_m_2 = 0;
        long double P_m_1 = 0;
        long double P_ = 1;
        Ptilda_[0]= P_;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // next iteration by n
        n = 1;
        P_m_2 = P_m_1; P_m_1 = P_;
        memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
        //P_ = sinTetta;

#ifdef _NORMALIZED_COEF
        P_ = sinTetta*_SQRT3;
#else
        P_ = sinTetta;
#endif
        Ptilda_[0]= P_;

        //Ptilda_[1] = n * P_m_1 + sinTetta * Ptilda_m_1[1]; // P'[1]  k == '

        // P = sin => d(P)/d(sin) = 1
#ifdef _NORMALIZED_COEF
        Ptilda_[1] =  _SQRT3;
#else
        Ptilda_[1] =  1;
#endif

        long double R0divR_ = R0divR[1]*R0divR[1];
        int ip = 0;
        // loop iteration starts from n=2 k = 0
        // formula 8 on page 92
        long double Xk[TOTAL_COEF+3];
        long double Yk[TOTAL_COEF+3];
        
        Xk[0] = 1.0;
        Yk[0] = 0.0;
        Xk[1] = Xk[0]*XdivR - Yk[0]*YdivR;
        Yk[1] = Yk[0]*XdivR + Xk[0]*YdivR;
        for (k = 2; k <= iLeg; k++)
        {
            Xk[k] = Xk[k-1]*XdivR - Yk[k-1]*YdivR;
            Yk[k] = Yk[k-1]*XdivR + Xk[k-1]*YdivR;
        }

        //int iXkYk = 1;
        int nb = 2;
        int ne = iLeg;
        int kb = 0;
        for (n = nb; n <=ne; n++)
        {
            long double x[3],y[3],z[3];
#if _DEBUG
            // sanity check n:
            if (n != nk_lm_Numbers[ip][0])
                exit (1);
#endif
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // next iteration by n
            P_m_2 = P_m_1; P_m_1 = P_;
            memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
#ifdef _NORMALIZED_COEF
            P_ = _p_n_m_1[n] *sinTetta * P_m_1 - _p_n_m_2[n]*P_m_2;  // P[2]
#else
            P_ = ((2.0* n-1.0) *sinTetta * P_m_1 - (n-1)*P_m_2)/n;  // P[2]
#endif
            long double P_nk = P_;
            long double XSumD, YSumD;
            /////////////////////////////////////////////////////////////////////////////  k =================0
            k = 0;
            Ptilda_[k]= P_;

#if _DEBUG
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);
#endif

            long double Ptilda_nk;
#ifdef _NORMALIZED_COEF
            if (2 == n) // k==0 & n == 2
                Ptilda_nk  = _tpk_n_k[n]*sinTetta;
            else
                Ptilda_nk  = _tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - _tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
            
#else
            //long double Ptilda_nk = n * P_m_1 + sinTetta * Ptilda_m_1[1];                        // P'[2]
            // this is equivalent, but better for iterations
            //long double Ptilda_nk  = (2*n-1) * Ptilda_m_1[k] + Ptilda_m_2[k+1];                      // P'[2]
            // this is equivalent
#ifdef _DO_NOT_SKIP_OBVIOUS
            if (2 == n) // k==0 & n == 2
                Ptilda_nk  = diagonal[n]*sinTetta;
            else
#endif
                Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));   // P'[2]
#endif
            Ptilda_[k+1] = Ptilda_nk; // store P'[2] for use 

            long double Qnk_ = C_S_nk[ip][0] * Xk[k] + C_S_nk[ip][1] * Yk[k];
            // on k=0 iteration!! i.e. n=2, k=0
            // Qnk = Cnk=0*Xk=0 +Snk=0*Yk=0
            // Xk=0 = 1; and Yk=0 = 0;
            // Qn0 = Cn0  => D_Qnk_Dxr =0; D_Qnk_Dyr=0
            //long double D_Qnk_Dxr_ = 0;
            //long double D_Qnk_Dyr_ = 0;
            // k is derivative

            // J case
            P_nk_x_Qnk_ = 0;       // sumv_n => sumgam_n
            Ptilda_nk_x_Qnk_ = 0;  // sumh_n
            P_nk_x_K_x_XSumD = 0;  // sumj_n
            P_nk_x_K_x_YSumD = 0;  // sumk_n

            if (ip == 0)
            {           //Sumgam_N := Pn[0]*Cn[O]*(n + 1)
                                                       // Sumh_N := Pn[1]* Cn[0];
                //_x20 = -(n+1) *XdivR    * P_nk * Qnk_  - Ptilda_nk *  Qnk_ * XdivR * SinTetta ;
                //_y20 = -(n+1) *YdivR    * P_nk * Qnk_  - Ptilda_nk *  Qnk_ * YdivR * SinTetta     ;
                //_z20 = -(n+1) *SinTetta * P_nk * Qnk_  +  Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta);
                P_20_x_Q20_ = P_nk * Qnk_;              // sumv_n => sumgam_n
#ifdef _NORMALIZED_COEF
                Ptilda_20_x_Qnk_ = _pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
                Ptilda_20_x_Qnk_ = Ptilda_nk *  Qnk_;   // sumh_n 
#endif

            }
            else
            {
                //x = (-(n+1) *XdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * XdivR * SinTetta   );
                //y = (-(n+1) *YdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * YdivR * SinTetta   );
                //z = (-(n+1) *SinTetta * P_ * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) );

                //x[0] = -(n+1) *XdivR    * P_nk * Qnk_;
                //y[0] = -(n+1) *YdivR    * P_nk * Qnk_;
                //z[0] = -(n+1) *SinTetta * P_nk * Qnk_;
                P_nk_x_Qnk_ += -(n+1) * P_nk * Qnk_;
                //x[1] = - Ptilda_nk *  Qnk_ * XdivR * SinTetta     ;
                //y[1] = - Ptilda_nk *  Qnk_ * YdivR * SinTetta     ;
                //z[1] =   Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta);
#ifdef _NORMALIZED_COEF
                Ptilda_nk_x_Qnk_ += - _pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
                Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_ ;
#endif
                //x[2] =0; y[2]=0;z[2]=0;
            }
            
            //////////////////////////////////////////////////////////////////////////   k ==================1
            // next iteration by k
            ip++;
            k = 1;
#if _DEBUG
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);
#endif
            //if (iXkYk < k)
            //{
            //    Xk[k] = Xk[k-1]*XdivR - Yk[k-1]*YdivR;
            //    Yk[k] = Yk[k-1]*XdivR + Xk[k-1]*YdivR;
            //    iXkYk = k;
            //}
            Qnk_ = C_S_nk[ip][0] * Xk[k] + C_S_nk[ip][1] * Yk[k];  //Bnmtil := Cnm*ctll[M] + Snm*stil[M];
            XSumD = C_S_nk[ip][0] * Xk[k-1] + C_S_nk[ip][1] * Yk[k-1];
            YSumD = C_S_nk[ip][0] * Yk[k-1] - C_S_nk[ip][1] * Xk[k-1];
            //XkPrev = Xk; YkPrev = Yk;

            P_nk = Ptilda_[k]; // P'[2] == (k= 1)
#ifdef _NORMALIZED_COEF
            if (2 == n) // k==1 && n==2
                Ptilda_nk  = _p_n_k[n];// + Ptilda_m_2[k+1];
            else
                Ptilda_nk  = _tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - _tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
#else
            if (2 == n) // k==1 && n==2
            {
                //Ptilda_nk  = (2*n-1) * Ptilda_m_1[k];// + Ptilda_m_2[k+1];
                Ptilda_nk  = diagonal[n];
            }
#ifdef _DO_NOT_SKIP_OBVIOUS
            else if (3 == n) // k==1 && n == 3
            {
                Ptilda_nk  = diagonal[n]*sinTetta;
            }
#endif
            else
            {
                Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));
            }
#endif
            Ptilda_[k+1] = Ptilda_nk; // store P"[2] for next use
            //x += (-(n+1+1) *XdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivR    * SinTetta     + P_nk * ( 1 *  XSumD   ));
            //y += (-(n+1+1) *YdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivR    * SinTetta     + P_nk * ( 1 * -YSumD   ));
            //z += (-(n+1+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta));


            //x[0] += -(n+1+1) *XdivR    * P_nk * Qnk_;
            //y[0] += -(n+1+1) *YdivR    * P_nk * Qnk_;
            //z[0] += -(n+1+1) *SinTetta * P_nk * Qnk_;
            P_nk_x_Qnk_ += -(n+1+1) * P_nk * Qnk_;      // sumgam_n (normalized == n+k+10 P_nk * Qnk_)

            //x[1] +=  - Ptilda_nk *  Qnk_ * XdivR    * SinTetta ;
            //y[1] +=  - Ptilda_nk *  Qnk_ * YdivR    * SinTetta;
            //z[1] +=    Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta);

#ifdef _NORMALIZED_COEF
            Ptilda_nk_x_Qnk_ += - _pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
            Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#endif

            //x[2] += P_nk * ( 1 *  XSumD   );
            //y[2] += P_nk * ( 1 * -YSumD   );
            //z[2] += 0;
            P_nk_x_K_x_XSumD += P_nk * ( 1 *  XSumD   );   // sumj_n (normalized == k * P_nk * XSumD)
            P_nk_x_K_x_YSumD += P_nk * ( 1 * -YSumD   );   // sumk_n (normakized ==-k * P_nk * YSumD

            ////////////////////////////////////////////////////////////////////////////////////////
            for (k = 2; k <=n; k++)
            {
                ////////////////////////////////////////////////////////////////////////////////////////
                // next iteration == k ==2
                ip++;
#if _DEBUG
                // sanity check k:
                if (k != nk_lm_Numbers[ip][1])
                    exit (1);
#endif
                P_nk = Ptilda_[k];
                if (k==n)
                {
                    //Xk[k] = Xk[k-1]*XdivR - Yk[k-1]*YdivR;
                    //Yk[k] = Yk[k-1]*XdivR + Xk[k-1]*YdivR;
                    //iXkYk = k;
                    Ptilda_nk = 0;//(2*n-1) * Ptilda_m_1[k] + Ptilda_m_2[k+1];
                }
                else
                {
#ifdef _NORMALIZED_COEF
                    if (k == (n-1))
                        Ptilda_nk  = _p_n_k[n];
                    else if (k == (n-2))
                        Ptilda_nk  = _tpk_n_k[n]*sinTetta;// + Ptilda_m_2[k+1];
                    else
                        Ptilda_nk  = _tp_nm1_k [n][k+1] * Ptilda_m_1[k+1]*sinTetta - _tp_nm2_k[n][k+1] * Ptilda_m_2[k+1];
#else
                    if (k == (n-1))
                        Ptilda_nk = diagonal[n];// + Ptilda_m_2[k+1];
#ifdef _DO_NOT_SKIP_OBVIOUS
                    else if (k == (n-2))
                        Ptilda_nk = diagonal[n]*sinTetta;
#endif
                    else
                        Ptilda_nk  = ((2*n-1) * Ptilda_m_1[k+1]*sinTetta - (n + (k+1) -1)*Ptilda_m_2[k+1])/(n-(k+1));
#endif
                }
                Ptilda_[k+1] = Ptilda_nk; // store P'"[2] (third derivative) for next use

                Qnk_ = C_S_nk[ip][0] * Xk[k] + C_S_nk[ip][1] * Yk[k];
                XSumD = C_S_nk[ip][0] * Xk[k-1] + C_S_nk[ip][1] * Yk[k-1];
                YSumD = C_S_nk[ip][0] * Yk[k-1] - C_S_nk[ip][1] * Xk[k-1];

                //XkPrev = Xk; YkPrev = Yk;

                

                //x += (-(n+k+1) *XdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivR    * SinTetta   + P_nk * ( k *  XSumD   ));
                //y += (-(n+k+1) *YdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivR    * SinTetta   + P_nk * ( k * -YSumD   ));
                //z += (-(n+k+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta));


                //x[0] += -(n+k+1) *XdivR    * P_nk * Qnk_;
                //y[0] += -(n+k+1) *YdivR    * P_nk * Qnk_;
                //z[0] += -(n+k+1) *SinTetta * P_nk * Qnk_;
                P_nk_x_Qnk_ += -(n+k+1) * P_nk * Qnk_;

                //x[1] += - Ptilda_nk *  Qnk_ * XdivR    * SinTetta;
                //y[1] += - Ptilda_nk *  Qnk_ * YdivR    * SinTetta;
                //z[1] +=   Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta);
#ifdef _NORMALIZED_COEF
                Ptilda_nk_x_Qnk_ += - _pt_nk[n][k] *Ptilda_nk *  Qnk_;   // sumh_n    (normalized == z[n][k] * Ptilda_nk *  Qnk_
#else
                Ptilda_nk_x_Qnk_ += - Ptilda_nk *  Qnk_;
#endif

                //x[2] +=  P_nk * ( k *  XSumD   );
                //y[2] +=  P_nk * ( k * -YSumD   );
                //z[2] += 0;
                P_nk_x_K_x_XSumD += P_nk * ( k *  XSumD   );
                P_nk_x_K_x_YSumD += P_nk * ( k * -YSumD   );
            }
            p_nk_x_Qnk_      += P_nk_x_Qnk_*R0divR_;      // sumgam
            ptilda_nk_x_Qnk_ +=Ptilda_nk_x_Qnk_*R0divR_;  // sumh
            p_nk_x_K_x_XSumD +=P_nk_x_K_x_XSumD*R0divR_;  // sumj
            p_nk_x_K_x_YSumD +=P_nk_x_K_x_YSumD*R0divR_;  // sumk
            //_x[n][0] = x[0];_x[n][1] = x[1];_x[n][2] = x[2];
            //_y[n][0] = y[0];_y[n][1] = y[1];_y[n][2] = y[2];
            //_z[n][0] = z[0];_z[n][1] = z[1];_z[n][2] = z[2];
            R0divR[n] = R0divR_;


            //X += x* R0divR_; Y += y *R0divR_; Z += z * R0divR_;
            R0divR_*= R0divR[1];
            ////////////////////////////////////////////////////////////////////////////////////////
            // next iteration == k ==2
            ip++;
        }
        //for (n=2; n <= iLeg; n++)
        //{
        //    _x[n] *= R0divR[n]; _y[n] *= R0divR[n]; _z[n] *= R0divR[n];
        //}
        //for (n=iLeg; n >=2; n--)
        //{
        //    //X += (_x[n][0]+_x[n][1]+_x[n][2])*R0divR[n]; 
        //    //Y += (_y[n][0]+_y[n][1]+_y[n][2])*R0divR[n];
        //    //Z += (_z[n][0]+_z[n][1]+_z[n][2])*R0divR[n];
        //    X += (P_nk_x_Qnk_[n] * XdivR    + Ptilda_nk_x_Qnk_[n] * XdivR * SinTetta         + P_nk_x_K_x_XSumD[n] )*R0divR[n]; 
        //    Y += (P_nk_x_Qnk_[n] * YdivR    + Ptilda_nk_x_Qnk_[n] * YdivR * SinTetta         + P_nk_x_K_x_YSumD[n] )*R0divR[n];
        //    Z += (P_nk_x_Qnk_[n] * SinTetta - Ptilda_nk_x_Qnk_[n] * (1- SinTetta * SinTetta)                       )*R0divR[n];
        //}
        //for (n=iLeg; n >=2; n--)
        //for (n=2; n <=iLeg; n++)
        //{
        //    p_nk_x_Qnk_      += P_nk_x_Qnk_[n]*R0divR[n];
        //    ptilda_nk_x_Qnk_ +=Ptilda_nk_x_Qnk_[n]*R0divR[n];
        //    p_nk_x_K_x_XSumD +=P_nk_x_K_x_XSumD[n]*R0divR[n];
        //    p_nk_x_K_x_YSumD +=P_nk_x_K_x_YSumD[n]*R0divR[n];
        //}

        //    lambda = sumgam + ep*sumh    
        //        sumgam==-p_nk_x_Qnk_   sumh==-ptilda_nk_x_Qnk_
        //        sumj == p_nk_x_K_x_XSumD   sumk == p_nk_x_K_x_YSumD
        //    -(lambda * XdivR - sumj)  
        // =>   -sumgam   * XdivR          -sumh        * XdivR *  ep              + sumj 
#if 0
        Xadd = (p_nk_x_Qnk_ * XdivR    + ptilda_nk_x_Qnk_ * XdivR * SinTetta         + p_nk_x_K_x_XSumD ); 

        //    -(lamda * YdivR - sumk)
        //      -sumgam   * YdivR          -sumh        * YdivR *  ep              + sunk
        Yadd = (p_nk_x_Qnk_ * YdivR    + ptilda_nk_x_Qnk_ * YdivR * SinTetta         + p_nk_x_K_x_YSumD );
        //    lamda * SinTetta - sumh
        //     - simgam   * sintetta        -sumh       * SinTetta * SinTetta      + sumh
        Zadd = (p_nk_x_Qnk_ * SinTetta - ptilda_nk_x_Qnk_ * (1- SinTetta * SinTetta)                    );

        //X += _x20*R0divR[2];  Y += _y20*R0divR[2];  Z += _z20*R0divR[2];

        Xadd+= (-(2+1) *XdivR    * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * XdivR * SinTetta)*R0divR[2] ;
        Yadd+= (-(2+1) *YdivR    * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * YdivR * SinTetta)*R0divR[2]     ;
        Zadd+= (-(2+1) *SinTetta * P_20_x_Q20_  + Ptilda_20_x_Qnk_ * (1-SinTetta*SinTetta))*R0divR[2];
        X =1; Y= 1; Z =1;
        trs_2_gcrs(Xadd, Yadd, Zadd);
        //trs_2_gcrs(X, Y, Z);
#else
        // Xadd = (p_nk_x_Qnk_ * XdivR    + ptilda_nk_x_Qnk_ * XdivR * SinTetta         + p_nk_x_K_x_XSumD );
        // Yadd = (p_nk_x_Qnk_ * YdivR    + ptilda_nk_x_Qnk_ * YdivR * SinTetta         + p_nk_x_K_x_YSumD );
        // Zadd = (p_nk_x_Qnk_ * SinTetta - ptilda_nk_x_Qnk_ * (1- SinTetta * SinTetta)                    );
        // Zadd = (p_nk_x_Qnk_ * SinTetta + ptilda_nk_x_Qnk_ * SinTetta * SinTetta                         - ptilda_nk_x_Qnk_ );
        // xadd = XdivR *(p_nk_x_Qnk_ + ptilda_nk_x_Qnk_ * SinTetta      + p_nk_x_K_x_XSumD/XdivR); 
        // Yadd = YdivR *(p_nk_x_Qnk_ + ptilda_nk_x_Qnk_ * SinTetta      + p_nk_x_K_x_YSumD/YdivR );
        // Zadd = SinTetta *(p_nk_x_Qnk_  + ptilda_nk_x_Qnk_ * SinTetta  - ptilda_nk_x_Qnk_/SinTetta );

        // Xadd+= (-(2+1) *XdivR    * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * XdivR * SinTetta)*R0divR[2] ;
        // Yadd+= (-(2+1) *YdivR    * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * YdivR * SinTetta)*R0divR[2]     ;
        // Zadd+= (-(2+1) *SinTetta * P_20_x_Q20_  + Ptilda_20_x_Qnk_ * (1-SinTetta*SinTetta))*R0divR[2];
        // Xadd+= XdivR  * (-(2+1) * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * SinTetta)*R0divR[2] ;
        // Yadd+= YdivR  * (-(2+1) * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * SinTetta)*R0divR[2] ;
        // Zadd+= SinTetta*(-(2+1) * P_20_x_Q20_  - Ptilda_20_x_Qnk_ * SinTetta)  + Ptilda_20_x_Qnk_/SinTetta)*R0divR[2];
        {
            Xadd = (    + p_nk_x_K_x_XSumD ); 
            Yadd = (    + p_nk_x_K_x_YSumD );
            Zadd = (    - ptilda_nk_x_Qnk_ );
            //Xadd+= (  - Ptilda_20_x_Qnk_ * XdivR * SinTetta)*R0divR[2] ;
            //Yadd+= (  - Ptilda_20_x_Qnk_ * YdivR * SinTetta)*R0divR[2]     ;
            Zadd+= (  + Ptilda_20_x_Qnk_ * (1.0))*R0divR[2];

            //X =XdivR; Y= YdivR; Z =SinTetta;

            //trs_2_gcrs(X, Y, Z); // that will be original X0divR Y0divR Z0divR
             // prove :
            //X =1.1*XdivR;       Y= 1.1*YdivR; Z =1.1*SinTetta;

            //trs_2_gcrs(X, Y, Z);
            // such way reduce error

            X = (p_nk_x_Qnk_ + ptilda_nk_x_Qnk_ * SinTetta ) + (-(2+1) * P_20_x_Q20_ - Ptilda_20_x_Qnk_ * SinTetta)*R0divR[2] ; 
                           Y=X;            Z=X;
            X=1-X;         Y=1-Y;          Z=1-Z;
            Xadd = -Xadd;  Yadd = -Yadd;   Zadd = -Zadd;
            trs_2_gcrs(Xadd, Yadd, Zadd);
        }
#endif
    };
#endif
#else
    // outdated : just for historic reference (to keep in source code instead in VSS)
    void FastSummXYZ( long double ValX, long double ValY, long double ValZ, long double ValR, long double &X, long double &Y, long double &Z, int iCurSat)
    {
        int n,k;
        long double tempX;
        long double tempY;
        long double tempZ;
        long double sinTetta, XdivR, YdivR;
        X = 0; Y = 0; Z = 0;
        long double _x[TOTAL_COEF];
        long double _y[TOTAL_COEF];
        long double _z[TOTAL_COEF];
        long double _x20,_y20,_z20;

        gcrs_2_trs(ValX, ValY, ValZ);
        tempX = ValX; tempY = ValY; tempZ = ValZ;
       // now earth in Terra Ref System
       // it is possible to calculate H to get air drag
       if (--iAtm[iCurSat] == 0)
       {
           long double dlLAT, dlLON;
           iAtm[iCurSat] = iItearationsPerSec;//479; // onc per 1000 iteration == 1 per sec
           h[iCurSat] = H =GetH(tempX, tempY, tempZ, 6378245.000, 6356863.019,dlLAT,dlLON);
           ro[iCurSat] = Ro=GetDens(); // 15C
       }

       sinTetta =ValZ/ValR;

        XdivR =   tempX/ValR;
        YdivR =   tempY/ValR;
        XdivRval = XdivR;
        YdivRval = YdivR;

        SinTetta = sinTetta;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        n = 0;  //  initial

        long double Ptilda_m_2[TOTAL_COEF];
        long double Ptilda_m_1[TOTAL_COEF];
        long double Ptilda_[TOTAL_COEF];
        for (k = 0; k < TOTAL_COEF; k++) 
        {
            Ptilda_[k] = 0;  Ptilda_m_1[k] =0;  Ptilda_m_2[k]=0;
        }
        long double P_m_2 = 0;
        long double P_m_1 = 0;
        long double P_ = 1;
        Ptilda_[0]= P_;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // next iteration by n
        n = 1;
        P_m_2 = P_m_1; P_m_1 = P_;
        memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
        //P_ = sinTetta;
        P_ = sinTetta;
        Ptilda_[0]= P_;

        //Ptilda_[1] = n * P_m_1 + sinTetta * Ptilda_m_1[1]; // P'[1]  k == '

        // P = sin => d(P)/d(sin) = 1
        Ptilda_[1] =  1;

        long double R0divR_ = R0divR[1]*R0divR[1];
        int ip = 0;
        // loop iteration starts from n=2 k = 0
        // formula 8 on page 92
        for (n = 2; n <=iLeg; n++)
        {
            long double x,y,z;
            x = 0;  y = 0;  z = 0;
            // sanity check n:
            if (n != nk_lm_Numbers[ip][0])
                exit (1);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // next iteration by n
            P_m_2 = P_m_1; P_m_1 = P_;
            memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
            P_ = ((2.0* n-1.0) *sinTetta * P_m_1 - (n-1)*P_m_2)/n;  // P[2]
            Ptilda_[0]= P_;
            long double XkDxrPrev =0;
            long double XkDyrPrev =0;
            long double YkDxrPrev =0;
            long double YkDyrPrev =0;
            long double XkPrev =1;
            long double YkPrev =0;
            //long double XkDxr, XkDyr, YkDxr, YkDyr;
            long double XSumD, YSumD;
            long double Xk =1;
            long double Yk =0;
            /////////////////////////////////////////////////////////////////////////////  k =================0
            k = 0;
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);

            long double Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;
            // on k=0 iteration!! i.e. n=2, k=0
            // Qnk = Cnk=0*Xk=0 +Snk=0*Yk=0
            // Xk=0 = 1; and Yk=0 = 0;
            // Qn0 = Cn0  => D_Qnk_Dxr =0; D_Qnk_Dyr=0
            //long double D_Qnk_Dxr_ = 0;
            //long double D_Qnk_Dyr_ = 0;
            // k is derivative
            long double Ptilda_nk = n * P_m_1 + sinTetta * Ptilda_m_1[1];                        // P'[2]
            Ptilda_[1] = Ptilda_nk; // store P'[2] for use 
            // J case
            if (ip == 0)
            {           //Sumgam_N := Pn[0]*Cn[O]*(n + 1)
                                                       // Sumh_N := Pn[1]* Cn[0];
                _x20 = (-(n+1) *XdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * XdivR * SinTetta     );
                _y20 = (-(n+1) *YdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * YdivR * SinTetta     );
                _z20 = (-(n+1) *SinTetta * P_ * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta));
            }
            else
            {
                x = (-(n+1) *XdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * XdivR * SinTetta   );
                y = (-(n+1) *YdivR    * P_ * Qnk_ - Ptilda_nk *  Qnk_ * YdivR * SinTetta   );
                z = (-(n+1) *SinTetta * P_ * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) );

            }
            
            //////////////////////////////////////////////////////////////////////////   k ==================1
            // next iteration by k
            ip++;
            k = 1;
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);
            Xk = XkPrev*XdivR - YkPrev*YdivR;
            Yk = YkPrev*XdivR + XkPrev*YdivR;
            Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;  //Bnmtil := Cnm*ctll[M] + Snm*stil[M];
            XSumD = C_S_nk[ip][0] * XkPrev + C_S_nk[ip][1] * YkPrev;
            YSumD = C_S_nk[ip][0] * YkPrev - C_S_nk[ip][1] * XkPrev;
            //XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR; // only XkPrev and YkPrev in use => Sumj_N := Sumj_N + M*Pnm (Cnm*ctil[M-1] + Snm*stil[M-1]);   
            //XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;          // only YlPrev and XkPrev in use => Sumk_N := Sumk_N - M*Pnm (Cnm*stil[M-1l - Snm*ctil[M-1]):
            //YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR; 
            //YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;
                                                                         // at the end multiplied Sumj := Sumj + Reorn * Sumj_N; and Sumk := Sumk + Reorn * Sumk_N
                                                                         // Cnk*(XkDxrPrev*XdivR + YkDxrPrev*XdivR - YkDxrPrev*YdivR + XkDxrPrev*YdivR)
                                                                         // Snk*(XkDyrPrev*XdivR + YkDyrPrev*XdivR - YkDyrPrev*YdivR + XkDyrPrev*YdivR)
                                                                         // Cnk*((XkDxrPrev+YkDxrPrev)*XdivR + (- YkDxrPrev + XkDxrPrev)*YdivR)
                                                                         // Snk*((XkDyrPrev+YkDyrPrev)*XdivR + (- YkDyrPrev + XkDyrPrev)*YdivR)
                                                                         // somehow comes to Qnk_ * k = ((XkDxrPrev*XdivR*XdivR  - YkDxrPrev*YdivR*XdivR + XkDyrPrev*XdivR* YdivR - YkDyrPrev*YdivR* YdivR ))*C_S_nk[ip][0]+((YkDxrPrev*XdivR*XdivR  + XkDxrPrev*YdivR*XdivR + YkDyrPrev*XdivR* YdivR + XkDyrPrev*YdivR* YdivR ))*C_S_nk[ip][1]
            //D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)
            //D_Qnk_Dxr_ = C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr;
            //D_Qnk_Dyr_ = C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr;
            //XkDxrPrev =XkDxr; XkDyrPrev =XkDyr; YkDxrPrev =YkDxr; YkDyrPrev =YkDyr;
            XkPrev = Xk; YkPrev = Yk;

            long double P_nk = Ptilda_[1]; // P'[2] == (k= 1)
            //Ptilda_nk  = n * Ptilda_[1] + sinTetta * Ptilda_m_1[2]; // P"[2] 
            Ptilda_nk  = (2*n-1) * Ptilda_m_1[1] + Ptilda_m_2[2];
            Ptilda_[2] = Ptilda_nk; // store P"[2] for next use
                  
                  // Sumgam_N := Sumgam_N + (N + m + 1) * Pnm * Bnmtil;
                                                   // Sumh_N += Pn(m+l)*Bnmtil;
            x += (-(n+1+1) *XdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivR    * SinTetta     + P_nk * ( 1 *  XSumD   ));
            y += (-(n+1+1) *YdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivR    * SinTetta     + P_nk * ( 1 * -YSumD   ));
            z += (-(n+1+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta));
                                                                                                  //  (-(C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr)*XdivR*SinTetta +
                                                                                                  //                                    - (C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr) * YdivR*SinTetta)
                                                                                                  //  (-(C_S_nk[ip][0]*(+ XkPrev*XdivR*SinTetta - YkPrev* YdivR*SinTetta ) + C_S_nk[ip][0]*((XkDxrPrev*XdivR*XdivR*SinTetta  - YkDxrPrev*YdivR*XdivR*SinTetta + XkDyrPrev*XdivR* YdivR*SinTetta - YkDyrPrev*YdivR* YdivR*SinTetta )) + 
                                                                                                  //     C_S_nk[ip][1]*(+ YkPrev*XdivR*SinTetta + XkPrev* YdivR*SinTetta ) + C_S_nk[ip][1]*((YkDxrPrev*XdivR*XdivR*SinTetta  + XkDxrPrev*YdivR*XdivR*SinTetta + YkDyrPrev*XdivR* YdivR*SinTetta + XkDyrPrev*YdivR* YdivR*SinTetta )))
                                                                                                  //  (-((C_S_nk[ip][0]*XkPrev*XdivR*SinTetta      + C_S_nk[ip][1] *YkPrev*XdivR*SinTetta)  + C_S_nk[ip][0]*((XkDxrPrev*XdivR*XdivR*SinTetta  - YkDxrPrev*YdivR*XdivR*SinTetta + XkDyrPrev*XdivR* YdivR*SinTetta - YkDyrPrev*YdivR* YdivR*SinTetta )) + 
                                                                                                  //     (C_S_nk[ip][0]*(- YkPrev* YdivR*SinTetta) + C_S_nk[ip][1] *XkPrev* YdivR*SinTetta) + C_S_nk[ip][1]*((YkDxrPrev*XdivR*XdivR*SinTetta  + XkDxrPrev*YdivR*XdivR*SinTetta + YkDyrPrev*XdivR* YdivR*SinTetta + XkDyrPrev*YdivR* YdivR*SinTetta )))

            ////////////////////////////////////////////////////////////////////////////////////////
            for (k = 2; k <=n; k++)
            {
                ////////////////////////////////////////////////////////////////////////////////////////
                // next iteration == k ==2
                ip++;

                // sanity check k:
                if (k != nk_lm_Numbers[ip][1])
                    exit (1);
                Xk = XkPrev*XdivR - YkPrev*YdivR;
                Yk = YkPrev*XdivR + XkPrev*YdivR;
                Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;
                XSumD = C_S_nk[ip][0] * XkPrev + C_S_nk[ip][1] * YkPrev;
                YSumD = C_S_nk[ip][0] * YkPrev - C_S_nk[ip][1] * XkPrev;

                //XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR;
                //XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;
                //YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR;
                //YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;

                //D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)
                //D_Qnk_Dxr_ = C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr;
                //D_Qnk_Dyr_ = C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr;
                //XkDxrPrev =XkDxr; XkDyrPrev =XkDyr; YkDxrPrev =YkDxr; YkDyrPrev =YkDyr;
                XkPrev = Xk; YkPrev = Yk;

                P_nk = Ptilda_[k];
                Ptilda_nk = (2*n-1) * Ptilda_m_1[k] + Ptilda_m_2[k+1];
                Ptilda_[k+1] = Ptilda_nk; // store P'"[2] (third derivative) for next use
                x += (-(n+k+1) *XdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivR    * SinTetta   + P_nk * ( k *  XSumD   ));
                y += (-(n+k+1) *YdivR    * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivR    * SinTetta   + P_nk * ( k * -YSumD   ));
                z += (-(n+k+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1- SinTetta * SinTetta));
                
            }
            _x[n] = x;_y[n] = y;_z[n] = z;
            R0divR[n] = R0divR_;

            //X += x* R0divR_; Y += y *R0divR_; Z += z * R0divR_;
            R0divR_*= R0divR[1];
            ////////////////////////////////////////////////////////////////////////////////////////
            // next iteration == k ==2
            ip++;
        }
        for (n=2; n <= iLeg; n++)
        {
            _x[n] *= R0divR[n]; _y[n] *= R0divR[n]; _z[n] *= R0divR[n];
        }
        for (n=iLeg; n >=2; n--)
        {
            X += _x[n]; Y += _y[n]; Z += _z[n];
        }
        X += _x20*R0divR[2];  Y += _y20*R0divR[2];  Z += _z20*R0divR[2];

        trs_2_gcrs(X, Y, Z);
    };
#endif
#else
    // outdated : just for historic reference (to keep in source code instead in VSS)
    void FastSummXYZ( long double ValX, long double ValY, long double ValZ, long double ValR, long double &X, long double &Y, long double &Z, int iCurSat)
    {
        int n,k;
        long double tempX;
        long double tempY;
        long double sinTetta, XdivR, YdivR, tempValX, tempValY;
        X = 0; Y = 0; Z = 0;
        long double _x[TOTAL_COEF];
        long double _y[TOTAL_COEF];
        long double _z[TOTAL_COEF];
        long double _x20,_y20,_z20;

#if 0
            tempX = cos(Lambda) * ValX - sin(Lambda) * ValY;
            tempY = sin(Lambda) * ValX + cos(Lambda) * ValY;
#else
            tempX = cos(Lambda) * ValX + sin(Lambda) * ValY;
            tempY = -sin(Lambda) * ValX + cos(Lambda) * ValY;

#endif
        sinTetta =ValZ/ValR;

        XdivR =   tempX/ValR;
        YdivR =   tempY/ValR;
        XdivRval = XdivR;
        YdivRval = YdivR;

        SinTetta = sinTetta;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        n = 0;  //  initial

        long double Ptilda_m_2[TOTAL_COEF];
        long double Ptilda_m_1[TOTAL_COEF];
        long double Ptilda_[TOTAL_COEF];
        for (k = 0; k < TOTAL_COEF; k++) 
        {
            Ptilda_[k] = 0;  Ptilda_m_1[k] =0;  Ptilda_m_2[k]=0;
        }
        long double P_m_2 = 0;
        long double P_m_1 = 0;
        long double P_ = 1;
        Ptilda_[0]= P_;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // next iteration by n
        n = 1;
        P_m_2 = P_m_1; P_m_1 = P_;
        memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
        //P_ = sinTetta;
        P_ = sinTetta;
        Ptilda_[0]= P_;

        //Ptilda_[1] = n * P_m_1 + sinTetta * Ptilda_m_1[1]; // P'[1]  k == '

        // P = sin => d(P)/d(sin) = 1
        Ptilda_[1] =  1;

        long double R0divR_ = R0divR[1]*R0divR[1];
        int ip = 0;
        // loop iteration starts from n=2 k = 0
        // formula 8 on page 92
        for (n = 2; n <=iLeg; n++)
        {
            long double x,y,z;
            x = 0;  y = 0;  z = 0;
            // sanity check n:
            if (n != nk_lm_Numbers[ip][0])
                exit (1);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // next iteration by n
            P_m_2 = P_m_1; P_m_1 = P_;
            memcpy(Ptilda_m_2,Ptilda_m_1, sizeof(Ptilda_m_2)); memcpy(Ptilda_m_1,Ptilda_, sizeof(Ptilda_m_1));
            P_ = ((2.0* n-1.0) *sinTetta * P_m_1 - (n-1)*P_m_2)/n;  // P[2]
            Ptilda_[0]= P_;
            long double XkDxrPrev =0;
            long double XkDyrPrev =0;
            long double YkDxrPrev =0;
            long double YkDyrPrev =0;
            long double XkPrev =1;
            long double YkPrev =0;
            long double XkDxr, XkDyr, YkDxr, YkDyr;
            long double Xk =1;
            long double Yk =0;
            /////////////////////////////////////////////////////////////////////////////  k =================0
            k = 0;
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);

            long double Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;
            // on k=0 iteration!! i.e. n=2, k=0
            // Qnk = Cnk=0*Xk=0 +Snk=0*Yk=0
            // Xk=0 = 1; and Yk=0 = 0;
            // Qn0 = Cn0  => D_Qnk_Dxr =0; D_Qnk_Dyr=0
            long double D_Qnk_Dxr_ = 0;
            long double D_Qnk_Dyr_ = 0;
            // k is derivative
            long double Ptilda_nk = n * P_m_1 + sinTetta * Ptilda_m_1[1];                        // P'[2]
            Ptilda_[1] = Ptilda_nk; // store P'[2] for use 
            // J case
            if (ip == 0)
            {           //Sumgam_N := Pn[0]*Cn[O]*(n + 1)
                                                       // Sumh_N := Pn[1]* Cn[0];
                _x20 = (-(n+1) *XdivRval * P_ * Qnk_ - Ptilda_nk *  Qnk_ * XdivRval * SinTetta   + P_ * ( D_Qnk_Dxr_*(1-XdivRval*XdivRval)- D_Qnk_Dyr_ * YdivRval*XdivRval     ));
                _y20 = (-(n+1) *YdivRval * P_ * Qnk_ - Ptilda_nk *  Qnk_ * YdivRval * SinTetta   + P_ * (-D_Qnk_Dxr_*XdivRval*YdivRval    + D_Qnk_Dyr_ * (1-YdivRval*YdivRval) ));
                _z20 = (-(n+1) *SinTetta * P_ * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) + P_ * (-D_Qnk_Dxr_*XdivRval*SinTetta    - D_Qnk_Dyr_ * YdivRval*SinTetta     ));
            }
            else
            {
                x = (-(n+1) *XdivRval * P_ * Qnk_ - Ptilda_nk *  Qnk_ * XdivRval * SinTetta   + P_ * ( D_Qnk_Dxr_*(1-XdivRval*XdivRval)- D_Qnk_Dyr_ * YdivRval*XdivRval     ));
                y = (-(n+1) *YdivRval * P_ * Qnk_ - Ptilda_nk *  Qnk_ * YdivRval * SinTetta   + P_ * (-D_Qnk_Dxr_*XdivRval*YdivRval    + D_Qnk_Dyr_ * (1-YdivRval*YdivRval) ));
                z = (-(n+1) *SinTetta * P_ * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) + P_ * (-D_Qnk_Dxr_*XdivRval*SinTetta    - D_Qnk_Dyr_ * YdivRval*SinTetta     ));

            }
            
            //////////////////////////////////////////////////////////////////////////   k ==================1
            // next iteration by k
            ip++;
            k = 1;
            // sanity check k:
            if (k != nk_lm_Numbers[ip][1])
                exit (1);
            Xk = XkPrev*XdivR - YkPrev*YdivR;
            Yk = YkPrev*XdivR + XkPrev*YdivR;
            Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;  //Bnmtil := Cnm*ctll[M] + Snm*stil[M];
            XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR; // only XkPrev and YkPrev in use => Sumj_N := Sumj_N + M*Pnm (Cnm*ctil[M-1] + Snm*stil[M-1]);   
            XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;          // only YlPrev and XkPrev in use => Sumk_N := Sumk_N - M*Pnm (Cnm*stil[M-1l - Snm*ctil[M-1]):
            YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR; 
            YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;
                                                                         // at the end multiplied Sumj := Sumj + Reorn * Sumj_N; and Sumk := Sumk + Reorn * Sumk_N
                                                                         // Cnk*(XkDxrPrev*XdivR + YkDxrPrev*XdivR - YkDxrPrev*YdivR + XkDxrPrev*YdivR)
                                                                         // Snk*(XkDyrPrev*XdivR + YkDyrPrev*XdivR - YkDyrPrev*YdivR + XkDyrPrev*YdivR)
                                                                         // Cnk*((XkDxrPrev+YkDxrPrev)*XdivR + (- YkDxrPrev + XkDxrPrev)*YdivR)
                                                                         // Snk*((XkDyrPrev+YkDyrPrev)*XdivR + (- YkDyrPrev + XkDyrPrev)*YdivR)
                                                                         // somehow comes to Qnk_ * k = ((XkDxrPrev*XdivR*XdivR  - YkDxrPrev*YdivR*XdivR + XkDyrPrev*XdivR* YdivR - YkDyrPrev*YdivR* YdivR ))*C_S_nk[ip][0]+((YkDxrPrev*XdivR*XdivR  + XkDxrPrev*YdivR*XdivR + YkDyrPrev*XdivR* YdivR + XkDyrPrev*YdivR* YdivR ))*C_S_nk[ip][1]
            //D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)
            D_Qnk_Dxr_ = C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr;
            D_Qnk_Dyr_ = C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr;
            XkDxrPrev =XkDxr; XkDyrPrev =XkDyr; YkDxrPrev =YkDxr; YkDyrPrev =YkDyr;
            XkPrev = Xk; YkPrev = Yk;

            long double P_nk = Ptilda_[1]; // P'[2] == (k= 1)
            //Ptilda_nk  = n * Ptilda_[1] + sinTetta * Ptilda_m_1[2]; // P"[2] 
            Ptilda_nk  = (2*n-1) * Ptilda_m_1[1] + Ptilda_m_2[2];
            Ptilda_[2] = Ptilda_nk; // store P"[2] for next use
                  
                  // Sumgam_N := Sumgam_N + (N + m + 1) * Pnm * Bnmtil;
                                                   // Sumh_N += Pn(m+l)*Bnmtil;
            x += (-(n+1) *XdivRval * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivRval * SinTetta   + P_nk * ( D_Qnk_Dxr_*(1-XdivRval*XdivRval)- D_Qnk_Dyr_ * YdivRval*XdivRval     ));
            y += (-(n+1) *YdivRval * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivRval * SinTetta   + P_nk * (-D_Qnk_Dxr_*XdivRval*YdivRval    + D_Qnk_Dyr_ * (1-YdivRval*YdivRval) ));
            z += (-(n+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) + P_nk * (-D_Qnk_Dxr_*XdivRval*SinTetta    - D_Qnk_Dyr_ * YdivRval*SinTetta     ));
                                                                                                  //  (-(C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr)*XdivR*SinTetta +
                                                                                                  //                                    - (C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr) * YdivR*SinTetta)
                                                                                                  //  (-(C_S_nk[ip][0]*(+ XkPrev*XdivR*SinTetta - YkPrev* YdivR*SinTetta ) + C_S_nk[ip][0]*((XkDxrPrev*XdivR*XdivR*SinTetta  - YkDxrPrev*YdivR*XdivR*SinTetta + XkDyrPrev*XdivR* YdivR*SinTetta - YkDyrPrev*YdivR* YdivR*SinTetta )) + 
                                                                                                  //     C_S_nk[ip][1]*(+ YkPrev*XdivR*SinTetta + XkPrev* YdivR*SinTetta ) + C_S_nk[ip][1]*((YkDxrPrev*XdivR*XdivR*SinTetta  + XkDxrPrev*YdivR*XdivR*SinTetta + YkDyrPrev*XdivR* YdivR*SinTetta + XkDyrPrev*YdivR* YdivR*SinTetta )))
                                                                                                  //  (-((C_S_nk[ip][0]*XkPrev*XdivR*SinTetta      + C_S_nk[ip][1] *YkPrev*XdivR*SinTetta)  + C_S_nk[ip][0]*((XkDxrPrev*XdivR*XdivR*SinTetta  - YkDxrPrev*YdivR*XdivR*SinTetta + XkDyrPrev*XdivR* YdivR*SinTetta - YkDyrPrev*YdivR* YdivR*SinTetta )) + 
                                                                                                  //     (C_S_nk[ip][0]*(- YkPrev* YdivR*SinTetta) + C_S_nk[ip][1] *XkPrev* YdivR*SinTetta) + C_S_nk[ip][1]*((YkDxrPrev*XdivR*XdivR*SinTetta  + XkDxrPrev*YdivR*XdivR*SinTetta + YkDyrPrev*XdivR* YdivR*SinTetta + XkDyrPrev*YdivR* YdivR*SinTetta )))

            ////////////////////////////////////////////////////////////////////////////////////////
            for (k = 2; k <=n; k++)
            {
                ////////////////////////////////////////////////////////////////////////////////////////
                // next iteration == k ==2
                ip++;

                // sanity check k:
                if (k != nk_lm_Numbers[ip][1])
                    exit (1);
                Xk = XkPrev*XdivR - YkPrev*YdivR;
                Yk = YkPrev*XdivR + XkPrev*YdivR;
                Qnk_ = C_S_nk[ip][0] * Xk + C_S_nk[ip][1] * Yk;
                XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR;
                XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;
                YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR;
                YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;

                //D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)
                D_Qnk_Dxr_ = C_S_nk[ip][0]*XkDxr + C_S_nk[ip][1]*YkDxr;
                D_Qnk_Dyr_ = C_S_nk[ip][0]*XkDyr + C_S_nk[ip][1]*YkDyr;
                XkDxrPrev =XkDxr; XkDyrPrev =XkDyr; YkDxrPrev =YkDxr; YkDyrPrev =YkDyr;
                XkPrev = Xk; YkPrev = Yk;

                P_nk = Ptilda_[k];
                Ptilda_nk = (2*n-1) * Ptilda_m_1[k] + Ptilda_m_2[k+1];
                Ptilda_[k+1] = Ptilda_nk; // store P'"[2] (third derivative) for next use
                x += (-(n+1) *XdivRval * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * XdivRval * SinTetta   + P_nk * ( D_Qnk_Dxr_*(1-XdivRval*XdivRval)- D_Qnk_Dyr_ * YdivRval*XdivRval     ));
                y += (-(n+1) *YdivRval * P_nk * Qnk_ - Ptilda_nk *  Qnk_ * YdivRval * SinTetta   + P_nk * (-D_Qnk_Dxr_*XdivRval*YdivRval    + D_Qnk_Dyr_ * (1-YdivRval*YdivRval) ));
                z += (-(n+1) *SinTetta * P_nk * Qnk_ + Ptilda_nk *  Qnk_ * (1-SinTetta*SinTetta) + P_nk * (-D_Qnk_Dxr_*XdivRval*SinTetta    - D_Qnk_Dyr_ * YdivRval*SinTetta     ));
                // on last (z) P_nk * Qnk_ * k == - P_nk * (-D_Qnk_Dxr_*XdivRval*SinTetta    - D_Qnk_Dyr_ * YdivRval*SinTetta     )
                
            }
            _x[n] = x;_y[n] = y;_z[n] = z;
            R0divR[n] = R0divR_;

            //X += x* R0divR_; Y += y *R0divR_; Z += z * R0divR_;
            R0divR_*= R0divR[1];
            ////////////////////////////////////////////////////////////////////////////////////////
            // next iteration == k ==2
            ip++;
        }
        for (n=2; n <= iLeg; n++)
        {
            _x[n] *= R0divR[n]; _y[n] *= R0divR[n]; _z[n] *= R0divR[n];
        }
        for (n=iLeg; n >=2; n--)
        {
            X += _x[n]; Y += _y[n]; Z += _z[n];
        }
        X += _x20*R0divR[2];  Y += _y20*R0divR[2];  Z += _z20*R0divR[2];
#if 0
        tempX = cos(-Lambda) * X - sin(-Lambda) * Y;
        tempY = sin(-Lambda) * X + cos(-Lambda) * Y;
#else
        tempX = cos(-Lambda) * X + sin(-Lambda) * Y;
        tempY = -sin(-Lambda) * X + cos(-Lambda) * Y;
#endif
        X = tempX;
        Y = tempY;

    };
#endif    
#ifdef ALL_OLD_CODE
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int CalcP( long double ValX, long double ValY, long double ValZ, long double ValR)
    {
        int n,k;
        long double tempX;
        long double tempY;
        long double sinTetta, XdivR, YdivR;
        // Lambda                               // 0 - 143 2 - 165 4 - 049
        //Lambda =  -Lambda -M_PI/2;            // 0 - 162 2 - 152 4 - 101
        //Lambda =  Lambda -M_PI/2;//           // 0 - 164 2 - 143 4 - 147
        //Lambda = -Lambda +M_PI;               // 0 - 082 2 - 106 4 - 048
        //Lambda =  Lambda -3*M_PI/2;           // 0 - 052 2 - 059 4 - 149 5 - 038 6 - 103
        //Lambda = -Lambda +M_PI/2;             // 0 - 047 2 - 149 4 - 159
        //Lambda =  Lambda +M_PI/2;             // 0 - 052 2 - 059 4 - 149
        //Lambda = - Lambda;                    // 0 - 161 2 - 134 4 - 161
        //Lambda = - Lambda + M_PI; // 0 min - 082 2 - 059 4 - 149
        //Lambda =0 ;
        //Lambda = 0.1;
        //Lambda = 0.2;
        //Lambda = 0.3;
        //Lambda = 0.4;
        //Lambda = 0.5;
        //Lambda = 0.6;
        //Lambda = 0.7;
        //Lambda = 0.8;
        //Lambda = 0.9;
        //Lambda = 1.0;
        //Lambda = 1.1;
        //Lambda = 1.2;
        //Lambda = 1.3;
        //Lambda = 1.4;
        //Lambda = 1.5;
        //Lambda = 1.6;
        //Lambda = 1.7;
        //Lambda = 1.8;
        //Lambda = 1.9;
        //Lambda = 2.0;
        //Lambda = 2.08;
        //Lambda = 2.1;
        //Lambda = 2.2;
        //Lambda = 2.3;
        //Lambda = 2.4;
        //Lambda = 2.45;
        //Lambda = 2.5;
        //Lambda = 2.6;
        //Lambda = 2.7;
        //Lambda = 2.8;
        //Lambda = 2.9;
        //Lambda = 3.0;
        //Lambda = 3.1;
        //Lambda = 3.2;
        //Lambda = 3.3;
        //Lambda = 3.4;
        //Lambda = 3.5;
        // Lambda = 3.6;
        //Lambda = 3.7;
        //Lambda = 3.8;
        //Lambda = 3.9;
        //Lambda = 4.0;
        //Lambda = 4.1;
        //Lambda = 4.2;
        //Lambda = 4.3;
        //Lambda = 4.4;
        //Lambda = 4.5;
        //Lambda = 4.6;
        //Lambda = 4.7;
        //Lambda = 4.8;
        //Lambda = 4.9;
        //Lambda = 5.0;
        //Lambda = 5.1;
        //Lambda = 5.2;
        //Lambda = 5.3;
        //Lambda = 5.4;
        //Lambda = 5.5; 
        //Lambda = 5.6;  
        //Lambda = 5.7;     
        //Lambda = 5.8;
        //Lambda = 5.9; 
        //Lambda = 6.0;
        //Lambda = 6.1;
        //Lambda = 6.2;

         {
            //tempX = cos(Lambda) * XdivR - sin(Lambda) * YdivR;
            //tempY = sin(Lambda) * XdivR + cos(Lambda) * YdivR;
            //XdivR = tempX;YdivR = tempY;
#if 0
            tempX = cos(Lambda) * ValX - sin(Lambda) * ValY;
            tempY = sin(Lambda) * ValX + cos(Lambda) * ValY;
#else
            tempX = cos(Lambda) * ValX + sin(Lambda) * ValY;
            tempY = -sin(Lambda) * ValX + cos(Lambda) * ValY;

#endif
            //ValX = tempX;
            //ValY = tempY;

            sinTetta =ValZ/ValR;
            //sinTetta =sqrt(tempX*tempX+ tempY*tempY)/ValR;
            XdivR =   tempX/ValR;
            YdivR =   tempY/ValR;

        }

        XdivRval = XdivR;
        YdivRval = YdivR;
#ifdef ALL_OLD_CODE
        SinTetta = sinTetta;
#endif
#define CPV 0.0000005
#if 0
        if ((XdivRval > -CPV) && (XdivRval < CPV))
            return 1;
        if ((YdivRval > -CPV) && (YdivRval < CPV))
            return 2;
#endif
#define _ACCOUNT_SIN

#ifndef _ACCOUNT_SIN
        if ((sinTetta > -CPV) && (sinTetta < CPV))
            return 3;
#endif


#if 0
        // power of a cos
        OneMinusSinTettaInSquare = 1.0 - sinTetta*sinTetta;
        OneMinusXdivRInSquare =    1.0 - XdivR * XdivR;
        OneMinusYdivRInSquare =    1.0 - YdivR * YdivR;
#endif

        //OneMinusSinTettaInSquare2 = (ValR*ValR - ValZ*ValZ)/(ValR*ValZ);
        //OneMinusXdivRInSquare2 =    (ValR*ValR - tempX*tempX)/(ValR*tempX);
        //OneMinusYdivRInSquare2 =    (ValR*ValR - tempY*tempY)/(ValR*tempY);

        // clean all array
        //for (n = 0; n<=iLeg+1;n++)
        //{
        //    for (k=0; k <=iLeg+1;k++)
        //    {
        //        Pnk_tilda[n][k] = 0;
        //    }
        //}
        // legandr functions from sinTetta
        P[0] = 1.0;
//#define CPV 0.0000005

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
#if 1
        Ptilda[1] = 1;//sqrt(ValX*ValX + ValY*ValY)/ValR;; // P1 was a X == P1' == 1

        Pnk_tilda[0][1] = Ptilda[0];
        Pnk_tilda[1][1] = Ptilda[1];
        for (n=2;n<=iLeg;n++)
        {   // all derivatives
            // for P2' == 2* P1 + sin(tetta) * 1
            Ptilda[n] = n * P[n-1] + sinTetta * Ptilda[n-1];
            Pnk_tilda[n][1] = Ptilda[n];
        }
#else
        Pnk_tilda[0][1] = Ptilda[0];
        for (n=1;n<=iLeg;n++)
        {   // all derivatives
            // for P2' == 2* P1 + sin(tetta) * 1
            Ptilda[n] = n * P[n-1] + sinTetta * Ptilda[n-1];
            Pnk_tilda[n][1] = Ptilda[n];
        }

#endif
        Pnk_tilda[0][2] = 0;
        Pnk_tilda[1][2] = 0;
        for (k= 3;k <=iLeg+1;k++)
        {
            Pnk_tilda[0][k] = 0;
            Pnk_tilda[1][k] =0;
        }
        // derivatives for dK(Pn(sinTetta))/d(sinTetta)**K
        for (n= 2;n <=iLeg;n++)
        {
            for (k = 2; k<=n+1;k++)
            {
                Pnk_tilda[n][k] = (2*n-1) * Pnk_tilda[n-1][k-1] + Pnk_tilda[n-2][k];
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
#if 0
        Xk[0] = 1; Yk[0] =0; 
        Xk[1] = XdivR; Yk[1] = YdivR;
        for (k = 2; k<=(iLeg+1); k++)
        {
            Xk[k] = Xk[k-1]*XdivR - Yk[k-1] * YdivR;
            Yk[k] = Yk[k-1]*XdivR + Xk[k-1] * YdivR;
        }
        XkDxr[0] = 0; XkDyr[0] = 0; YkDxr[0] = 0; YkDyr[0] = 0;
        for (k = 1; k <=iLeg+1; k++)
        {
            XkDxr[k] = XkDxr[k-1]*XdivR + Xk[k-1]          - YkDxr[k-1]*YdivR;
            XkDyr[k] = XkDyr[k-1]*XdivR - YkDyr[k-1]*YdivR - Yk[k-1];
            YkDxr[k] = YkDxr[k-1]*XdivR + Yk[k-1]          + XkDxr[k-1]*YdivR;
            YkDyr[k] = YkDyr[k-1]*XdivR + XkDyr[k-1]*YdivR + Xk[k-1];
        }
        {
            // formula 8 on page 92
             for (n = 2; n <=iLeg+1; n++)
             {
                 for (k = 0; k <=iLeg+1; k++)
                 {
                     Qnk[n][k] = CNK[n][k] * Xk[k] + SNK[n][k] * Yk[k];
                     D_Qnk_Dxr[n][k] = CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k];
                     D_Qnk_Dyr[n][k] = CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k];
                 }
             }
        }
#else
        // formula 8 on page 92
        for (n = 2; n <=iLeg; n++)
        {
            long double XkDxrPrev =0;
            long double XkDyrPrev =0;
            long double YkDxrPrev =0;
            long double YkDyrPrev =0;
            long double XkPrev =1;
            long double YkPrev =0;
            long double XkDxr, XkDyr, YkDxr, YkDyr;
            long double Xk =1;
            long double Yk =0;

            Qnk[n][0] = CNK[n][0] * Xk + SNK[n][0] * Yk;
            // on k=0 iteration!! i.e. n=2, k=0
            // Qnk = Cnk=0*Xk=0 +Snk=0*Yk=0
            // Xk=0 = 1; and Yk=0 = 0;
            // Qn0 = Cn0  => D_Qnk_Dxr =0; D_Qnk_Dyr=0


            D_Qnk_Dxr[n][0] = 0;
            D_Qnk_Dyr[n][0] = 0;
            //XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR;
            //XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;
            //YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR;
            //YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;

            //D_Qnk_Dxr[n][0] = CNK[n][k]*XkDxr + SNK[n][k]*YkDxr;
            //D_Qnk_Dyr[n][0] = CNK[n][k]*XkDyr + SNK[n][k]*YkDyr;
            for (k = 1; k <=n; k++)
            {
                Xk = XkPrev*XdivR - YkPrev*YdivR;
                Yk = YkPrev*XdivR + XkPrev*YdivR;
                Qnk[n][k] = CNK[n][k] * Xk + SNK[n][k] * Yk;
                XkDxr = XkDxrPrev*XdivR + XkPrev          - YkDxrPrev*YdivR;
                XkDyr = XkDyrPrev*XdivR - YkDyrPrev*YdivR - YkPrev;
                YkDxr = YkDxrPrev*XdivR + YkPrev          + XkDxrPrev*YdivR;
                YkDyr = YkDyrPrev*XdivR + XkDyrPrev*YdivR + XkPrev;

                //D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)
                D_Qnk_Dxr[n][k] = CNK[n][k]*XkDxr + SNK[n][k]*YkDxr;
                D_Qnk_Dyr[n][k] = CNK[n][k]*XkDyr + SNK[n][k]*YkDyr;
                XkDxrPrev =XkDxr; XkDyrPrev =XkDyr; YkDxrPrev =YkDxr; YkDyrPrev =YkDyr;
                XkPrev = Xk; YkPrev = Yk;
            }
        }
#endif
        return 0;
    };

    // outdated : just for historic reference (to keep in source code instead in VSS)
    void SummXYZ(int SatCalc, long double &X, long double &Y, long double &Z)
    {
        int n,k;
        X=0; Y=0; Z=0;
        /*
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
            //
            //            = (n+1)* fm * (r0/r)**n * (-x/r**3) * Znk * Qnk +
            //              Rn * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (-x*z/r**3) * Qnk +
            //              Rn * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*D(Yk)/D(x/r)) * (1/r - x**2/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*D(Yk)/D(y/r))*(-x*y/r**3)]
            // as a result for example n=2 and k = 0
            // D(U20)/D(X) = (2+1)* fm * (r0/r)**2 * (-x/r**3) * d(P2(sinTetta))/d(sinTetta) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(0+1)(P2(sinTetta)/d(sintetta)**(0+1) *(-x*z/r**3) * J2 +
            //              fm * r0**2 * (1/r)**3 * d(Pn(sinTetta))/d(sinTetta) * [D(Q20)/D(x/r) * (1/r - x**2/r**3) + D(Q20)/D(y/r)*(-x*y/r**3)]
            //            = - 3 * fm * (ro/r)**2 * (x/r**3) * Pnk_tilda[2][0] * J2 
            //              -     fm * (r0/r)**2 * (1/r) * Pnk_tilda[2][1] * x*z/r**3 *J2 
            //                    fm * (r0/r)**2 * (1/r) Pnk_tilda[2][0] * [0*(1/r-x**2/r**3) + 0*(-x*y/r**3)]

            //X += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]//CNK[n][0] 
            //     - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            X += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]
                 - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]
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
            //Y += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]
            //     - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            Y += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]
                 - R0divR[n] * SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]
            ;

            // Rn = fm /r * (r0/r)**n
            // Znk =   d(K)  (Pn(sintetta))/d(sintetta)**(k)
            // Znk+1 = d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1)

            // D(Unk)/D(z) = d(Rn) / d(1/r)    * D(1/r)/D(z)    * Znk    * Qnk +
            //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(z) * Qnk +
            //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(z) + D(Qnk)/D(y/r)*D(y/r)/D(z)]
            //            = (n+1)* fm * (r0/r)**n * (-z/r**3) * Znk * Qnk +
            //                     fm  * (r0/r)**n *1/r * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * (1/r-z**2/r**3) * Qnk +
            //                     fm  * (r0/r)**n *1//r * Znk * [D(Qnk)/D(x/r) * (- x*z/r**3) + D(Qnk)/D(y/r)*(-y*z/r**3)]
            //            = (n+1)* fm/r**2  * (r0/r)**n * (-z/r) * Znk * Qnk +
            //                     fm  * (r0/r)**n *1/r * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * 1/r* (1-z**2/r**2) * Qnk +
            //                     fm  * (r0/r)**n *1/r * Znk * 1/r * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x*z/r**2) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-y*z/r**2)]
            //            = (n+1)* fm/r**2  * (r0/r)**n * (-z/r) * Znk * Qnk +
            //                     fm/r**2  * (r0/r)**n * Znk+1 * (1-(z/r)**2) * Qnk +
            //                     fm/r**2  * (r0/r)**n * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x/r * z/r) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-y/r *z/r)]
            //            = fm/r**2  * (r0/r)**n ( (n+1)  * (-z/r) * Znk * Qnk +
            //                                     Znk+1 * (1-(z/r)**2) * Qnk +
            //                                     Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x/r * z/r) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-y/r *z/r)]

            //Z += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0] 
            //+ R0divR[n] * CosTetta[2]/ SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0];//CNK[n][0]; 
            Z += -(n+1) * R0divR[n] * Pnk_tilda[n][0] * Qnk[n][0]
            + R0divR[n] * CosTetta[2]/ SinTetta[1] * Pnk_tilda[n][1] * Qnk[n][0]
            ;
        }
        // second round
        */
        for (n= 2;n <=iLeg;n++)
        {
            // k=0 already done
            for (k = 0; k<=n; k++)
            {
                // D(Unk)/D(x) = d(Rn) / d(1/r)    * D(1/r)/D(x)    * Znk    * Qnk +
                //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(x) * Qnk +
                //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(x) + D(Qnk)/D(y/r)*D(y/r)/D(x)]

                //            = fm */r D((r0/r)**n)/D(1/r) * D(1/r)/D(x) * Znk * Qnk +
                //              fm/r (r0/r)**n * D(K+1)(Pn(sinTetta)/D(sintetta) * D(Z/r)/D(x) * Qnk +
                //              fm/r (r0/r)**n * Znk * [ (Cnk*D(Xk)/D(x/r)+Snk*D(Yk)/D(x/r)) * D(x/r)/D(x) + (Cnk*D(Xk)/D(y/r)+Snk*D(Yk)/D(y/r)) *D(y/r)/D(x) ]

                //            = fm * r0**n * (n+1)* (1/r)**n * (-x/r**3) * Znk * Qnk +
                //              fm/r (r0/r)**n * d(K+1)(Pn(sinTetta)/d(sintetta)**(k+1) * 1/r *(-x*z/r**2) * Qnk +
                //              fm/r (r0/r)**n * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*D(Yk)/D(x/r)) * (1/r - x**2/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-x*y/r**3)]

                //            = fm /r**2   *   -(n+1) * (r0/r)**n * x/r    Znk * Qnk +
                //              fm /r**2  *(r0/r)**n *   d(K+1)(Pn(sinTetta)/d(sintetta)        * -x/r * z/r *    Qnk +
                //              fm/r (r0/r)**n * Znk * 1/r * [(Cnk*XkDxr+Snk*YkDxr) * (1 - x**2/r**2) + (Cnk*XkDyr+Snk*YkDyr)*(-x*y/r**2)]

                //            = fm /r**2   * (r0/r)**n [  -(n+1) * Znk * Qnk x/r +
                //                                           Znk+1        * -x/r * SinTetta *    Qnk +
                //                                             Znk * ((Cnk*XkDxr+Snk*YkDxr) * (1 - (x/r)**2) + (Cnk*XkDyr+Snk*YkDyr)*(-x/r * y/r))
#if 0
                long double tempVal;
#if 1
                if (((XdivRval > CPV ) || (XdivRval < -CPV )))// && (((k+n)&1) == 1))
                {
#endif
                    //OneMinusXdivRInSquare_XdivRval[SatCalc][n][k] = tempVal = R0divR[n] * Pnk_tilda[n][k]*(D_Qnk_Dxr[n][k]*OneMinusXdivRInSquare2- D_Qnk_Dyr[n][k]*YdivRval);// (CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * OneMinusXdivRInSquare2;///XdivRval;
                    OneMinusXdivRInSquare_XdivRval[SatCalc][n][k] = tempVal = (R0divR[n] * Pnk_tilda[n][k]*(D_Qnk_Dxr[n][k]*(1-XdivRval*XdivRval)- D_Qnk_Dyr[n][k]*YdivRval*XdivRval))/XdivRval;
                    OldXSign[SatCalc][n][k] = XdivRval;
#if 1
                }
                else
                {
                    
                    if (((XdivRval >0) && (OldXSign[SatCalc][n][k] < 0)) || ((XdivRval <0) && (OldXSign[SatCalc][n][k] > 0)))
                    {
                        OneMinusXdivRInSquare_XdivRval[SatCalc][n][k] = tempVal = -OneMinusXdivRInSquare_XdivRval[SatCalc][n][k];
                        OldXSign[SatCalc][n][k] = - OldXSign[SatCalc][n][k];
                    }
                    else
                        tempVal = OneMinusXdivRInSquare_XdivRval[SatCalc][n][k];
                    //tempVal = (R0divR[n] * Pnk_tilda[n][k]*(D_Qnk_Dxr[n][k]*(1-XdivRval*XdivRval)- D_Qnk_Dyr[n][k]*YdivRval*XdivRval))/XdivRval;
                    //if ((OldXSign[SatCalc][n][k]/XdivRval) > 4)
                    {
                        printf("zero x ");
                    }
                        
                    //tempVal = OneMinusXdivRInSquare_XdivRval[SatCalc][n][k] * (OldXSign[SatCalc][n][k]/XdivRval);
                }
#endif
#endif
#if 0
                    X += R0divR[n] *
                                  (-(n+1)  /**XdivRval*/ *  Pnk_tilda[n][k] * Qnk[n][k]
                                   - Pnk_tilda[n][k+1] *  Qnk[n][k] /** XdivRval*/ * SinTetta
                                   //+ Pnk_tilda[n][k] * (/*(CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * OneMinusXdivRInSquare/XdivRval*/ //OneMinusXdivRInSquare_XdivRval[SatCalc][n][k]
                                   //                     - D_Qnk_Dyr[n][k]*YdivRval)//(CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])/**XdivRval*/*YdivRval)
                                                        
                                  )/*/XdivRval*/
                                  +tempVal
                    ;
#else
                    X += R0divR[n] *
                                  (-(n+1)  *XdivRval *  Pnk_tilda[n][k] * Qnk[n][k]
                                   - Pnk_tilda[n][k+1] *  Qnk[n][k] * XdivRval * SinTetta
                                   + Pnk_tilda[n][k]*(D_Qnk_Dxr[n][k]*(1-XdivRval*XdivRval)- D_Qnk_Dyr[n][k]*YdivRval*XdivRval)
                                   //+ Pnk_tilda[n][k] * (/*(CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * OneMinusXdivRInSquare/XdivRval*/ //OneMinusXdivRInSquare_XdivRval[SatCalc][n][k]
                                   //                     - D_Qnk_Dyr[n][k]*YdivRval)//(CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])/**XdivRval*/*YdivRval)
                                                        
                                  )/*/XdivRval*/

                    ;
#endif

                // D(Unk)/D(y) = d(Rn) / d(1/r)    * D(1/r)/D(y)    * Znk    * Qnk +
                //                 Rn * d(Znk)/d(z/r) * D(z/r)/D(y) * Qnk +
                //                 Rn * Znk * [D(Qnk)/D(x/r) * D(x/r)/D(y) + D(Qnk)/D(y/r)*D(y/r)/D(y)]

                //             = fm/r D((r0/r)**n)/D(1/r) * D(1/r)/D(y) * Znk * Qnk +
                //               fm/r (r0/r)**n * D(K+1)(Pn(sinTetta)/D(sintetta) * D(Z/r)/D(y) * Qnk +
                //               fm/r (r0/r)**n * Znk * [ (Cnk*D(Xk)/D(x/r)+Snk*D(Yk)/D(x/r)) * D(x/r)/D(y) + (Cnk*D(Xk)/D(y/r)+Snk*D(Yk)/D(y/r)) *D(y/r)/D(y) ]

                //            =  fm * r0**n * (n+1)* (1/r)**n * (-y/r**3) * Znk * Qnk +
                //              fm/r (r0/r)**n * Znk+1 * (-y*z/r**3) * Qnk +
                //              fm/r (r0/r)**n * Znk * [(Cnk*D(Xk)/D(x/r)+Snk*D(Yk)/D(x/r)) * (-xy/r**3) + (Cnk*D(Xk)/D(y/r)+Snk*D(Yk)/D(y/r))*(1/r-y**2/r**3)]

                //            = fm /r**2   * (r0/r)**n [  -(n+1) * Znk * *y/r * Qnk +
                //                                           Znk+1       * y/r * -SinTetta *    Qnk +
                //                                             Znk * ((Cnk*XkDxr+Snk*YkDxr) * (-x/r * y/r) + (Cnk*XkDyr+Snk*YkDyr)*(1-(y/r)**2))
#if 0
#if 1
                    if (((YdivRval > CPV ) || (YdivRval < -CPV )))// && (((k+n)&1) == 1))
                    {
#endif
                        //OneMinusYdivRInSquare_YdivRval[SatCalc][n][k] = tempVal= R0divR[n] *Pnk_tilda[n][k]*(-D_Qnk_Dxr[n][k]*XdivRval+D_Qnk_Dyr[n][k] * OneMinusYdivRInSquare2);//(CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])*OneMinusYdivRInSquare2;///YdivRval;
                        OneMinusYdivRInSquare_YdivRval[SatCalc][n][k] = tempVal= (R0divR[n] *Pnk_tilda[n][k]*(-D_Qnk_Dxr[n][k]*XdivRval*YdivRval +D_Qnk_Dyr[n][k] * (1-YdivRval*YdivRval)))/YdivRval;
                        OldYSign[SatCalc][n][k] = YdivRval;
#if 1
                    }
                    else
                    {
                        if (((YdivRval >0) && (OldYSign[SatCalc][n][k] < 0)) || ((YdivRval <0) && (OldYSign[SatCalc][n][k] > 0)))
                        {
                            OneMinusYdivRInSquare_YdivRval[SatCalc][n][k] = tempVal= - OneMinusYdivRInSquare_YdivRval[SatCalc][n][k];
                            OldYSign[SatCalc][n][k] = -OldYSign[SatCalc][n][k];
                        }
                        else
                            tempVal= OneMinusYdivRInSquare_YdivRval[SatCalc][n][k];
                        //tempVal= (R0divR[n] *Pnk_tilda[n][k]*(-D_Qnk_Dxr[n][k]*XdivRval*YdivRval +D_Qnk_Dyr[n][k] * (1-YdivRval*YdivRval)))/YdivRval;
                        //if ((OldYSign[SatCalc][n][k]/YdivRval)>4)
                            printf("zero y ");
                        //tempVal = OneMinusYdivRInSquare_YdivRval[SatCalc][n][k] * (OldYSign[SatCalc][n][k]/YdivRval);
                    }
#endif
#endif
#if 0
                    Y += R0divR[n] *
                                  (-(n+1) /**YdivRval*/ * Pnk_tilda[n][k] * Qnk[n][k]
                                     - Pnk_tilda[n][k+1] *  Qnk[n][k] /** YdivRval*/ * SinTetta
                                     //+ Pnk_tilda[n][k] * (-D_Qnk_Dxr[n][k]*XdivRval// + //(CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * XdivRval /**YdivRval*/ + 
                                         /*(CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])*OneMinusYdivRInSquare/YdivRval*/ //OneMinusYdivRInSquare_YdivRval[SatCalc][n][k]
                                     
                                     )/*/YdivRval*/
                            + tempVal
                    ;
#else
                    Y += R0divR[n] *
                                  (-(n+1) *YdivRval * Pnk_tilda[n][k] * Qnk[n][k]
                                     - Pnk_tilda[n][k+1] *  Qnk[n][k] * YdivRval * SinTetta
                                     + Pnk_tilda[n][k]*(-D_Qnk_Dxr[n][k]*XdivRval*YdivRval +D_Qnk_Dyr[n][k] * (1-YdivRval*YdivRval))
                                     
                                     )
                    ;

#endif
               
            //            = fm/r**2  * (r0/r)**n ( (n+1)  * (-z/r) * Znk * Qnk +
            //                                     Znk+1 * (1-(z/r)**2) * Qnk +
            //                                     Znk * [(Cnk*D(Xk)/D(x/r)+Snk*DYk)/D(x/r)) * (- x/r * z/r) + (Cnk*D(Xk)/D(y/r)+Snk*DYk)/D(y/r))*(-y/r *z/r)]
#if 0
#if 1
                    if (((SinTetta > CPV ) || (SinTetta < -CPV )))// && (((k+n)&1) == 1))
                    {
#endif
                        //OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k] = tempVal=R0divR[n] *OneMinusSinTettaInSquare2 * Pnk_tilda[n][k+1] * Qnk[n][k];// /SinTetta[1];
                        OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k] = tempVal =(Pnk_tilda[n][k+1] * Qnk[n][k]*R0divR[n]*(1-SinTetta*SinTetta))/SinTetta;
                        OldZSign[SatCalc][n][k] = SinTetta;
#if 1
                    }
                    else
                    {
                        if (((SinTetta >0) && (OldZSign[SatCalc][n][k] < 0)) || ((SinTetta <0) && (OldZSign[SatCalc][n][k] > 0)))
                        {
                            OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k] = tempVal = - OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k];
                            OldZSign[SatCalc][n][k] = -OldZSign[SatCalc][n][k];
                        }
                        else
                            tempVal = - OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k];
                        //tempVal = OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k] *(OldZSign[SatCalc][n][k]/SinTetta);
                        //tempVal =(Pnk_tilda[n][k+1] * Qnk[n][k]*R0divR[n]*(1-SinTetta*SinTetta))/SinTetta;
                        //if ((OldZSign[SatCalc][n][k]/SinTetta) > 4)
                        {
                            printf("zeroZ");
                        }
                    }
#endif
#endif
#if 0
                    Z += R0divR[n] *
                                (-(n+1)  /** SinTetta[1]*/ * Pnk_tilda[n][k] * Qnk[n][k] 
                                 //+ /*OneMinusSinTettaInSquare * Pnk_tilda[n][k+1] * Qnk[n][k] /SinTetta[1]*/ tempVal//OneMinusZdivRInSquare_ZdivRval[SatCalc][n][k]
                                 + Pnk_tilda[n][k] * (-D_Qnk_Dxr[n][k]*XdivRval - D_Qnk_Dyr[n][k]*YdivRval)//(CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * XdivRval/**SinTetta[1]*/ - (CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])*YdivRval/**SinTetta[1]*/)
                                 )/*/SinTetta[1]*/
                                 +tempVal
                    ;
#else

                    Z += R0divR[n] *
                                (-(n+1)  * SinTetta * Pnk_tilda[n][k] * Qnk[n][k] 
                                 + Pnk_tilda[n][k+1] * Qnk[n][k]*(1-SinTetta*SinTetta)
                                 + Pnk_tilda[n][k] * (-D_Qnk_Dxr[n][k]*XdivRval*SinTetta - D_Qnk_Dyr[n][k]*YdivRval*SinTetta)//(CNK[n][k]*XkDxr[k] + SNK[n][k]*YkDxr[k]) * XdivRval/**SinTetta[1]*/ - (CNK[n][k]*XkDyr[k] + SNK[n][k]*YkDyr[k])*YdivRval/**SinTetta[1]*/)
                                 )/*/SinTetta[1]*/
                    ;

#endif
            }
        }
        long double tempX = cos(-Lambda) * X - sin(-Lambda) * Y;
        long double tempY = sin(-Lambda) * X + cos(-Lambda) * Y;
        X = tempX;
        Y = tempY;
#if 0
        if (((XdivRval > CPV ) || (XdivRval < -CPV )))
        {
            fx[SatCalc] = X;
            X /= XdivRval;
            fsinX[SatCalc] =XdivRval;
        }
        else
        {
            printf("x");
            if (((fsinX[SatCalc]>0) && (XdivRval > 0)) || ((fsinX[SatCalc]<=0) && (XdivRval <= 0)))
            {
                X = fx[SatCalc]/fsinX[SatCalc];// Y = fy[SatCalc]; Z = fz[SatCalc];
            }
            else
            {
                fx[SatCalc]=-fx[SatCalc]; fsinX[SatCalc] = -fsinX[SatCalc];
                X = fx[SatCalc]/fsinX[SatCalc];//  Y = fy[SatCalc]; Z = fz[SatCalc];
            }
        }
        if (((YdivRval > CPV ) || (YdivRval < -CPV )))
        {
            fy[SatCalc] = Y;
            Y /= YdivRval;
            fsinY[SatCalc] =YdivRval;
        }
        else
        {
            printf("y");
            if (((fsinY[SatCalc]>0) && (YdivRval > 0)) || ((fsinY[SatCalc]<=0) && (YdivRval <= 0)))
            {
                //X = fx[SatCalc]; 
                Y = fy[SatCalc]/fsinY[SatCalc]; //Z = fz[SatCalc];
            }
            else
            {
                fy[SatCalc]=-fy[SatCalc]; fsinY[SatCalc] = -fsinY[SatCalc];
                //X = fx[SatCalc];  
                Y = fy[SatCalc]/fsinY[SatCalc];// Z = fz[SatCalc];
            }
        }
#endif
#ifndef _ACCOUNT_SIN
        if (((SinTetta > CPV ) || (SinTetta < -CPV )))
        {
            fz[SatCalc] = Z;
            Z /= SinTetta;
            fsinZ[SatCalc] =SinTetta;
        }
        else
        {
            printf("z");
            if (((fsinZ[SatCalc]>0) && (SinTetta > 0)) || ((fsinZ[SatCalc]<=0) && (SinTetta <= 0)))
            {
                //X = fx[SatCalc]; Y = fy[SatCalc]; 
                Z = fz[SatCalc]/fsinZ[SatCalc];
            }
            else
            {
                fz[SatCalc]=-fz[SatCalc]; fsinZ[SatCalc] = -fsinZ[SatCalc];
                //X = fx[SatCalc];  Y = fy[SatCalc]; 
                Z = fz[SatCalc]/fsinZ[SatCalc];
            }
        }
#endif
    };
    long double SummJ(void)
    {
        // this is a formula for gravitation potential (NOT a force)
        // to get a force needs to get analiticaly derivitive from Gravitation potencial
        // page 90 on Aksenov lectures : http://vadimchazov.narod.ru/lepa_zov/lesat.pdf
        //	- J2 * (r0/r)**2 * P2(sinPHI)  
        //	- J3 * (r0/r)**3 * P3(sinPHI) 
        //  - J4 * (r0/r)**4 * P4(sinPHI)
    
        long double Summ = 1.0;
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
                Summ += R0divR[n] * Ptilda[n] *(CNK[n][k]*cos((long double)k*Lambda) +SNK[n][k]*sin((long double)k*Lambda));
            }
        }

        return Summ;
    };
#endif
} TRAOBJ, *PTRAOBJ;

void AssignFromNASAData(TRAOBJ * SlS, double JDSec);

GLOBAL_VARIABLE TRAOBJ SolarSystem;
GLOBAL_VARIABLE TRAOBJ Sat;

GLOBAL_VARIABLE long double SunX;
GLOBAL_VARIABLE long double SunY;
GLOBAL_VARIABLE long double SunZ;
GLOBAL_VARIABLE long double SunR;
GLOBAL_VARIABLE long double dStartJD;//2451544.5; // if value dStartJD not set (==0.0) then use value from keplers elements of a satelite 0

GLOBAL_VARIABLE long double dMinFromNow;

GLOBAL_VARIABLE int JustFlySimulation;

GLOBAL_VARIABLE long double TotalDays;

GLOBAL_VARIABLE int iGr;

// ground stations
GLOBAL_VARIABLE long double GrLat[10];
GLOBAL_VARIABLE long double GrLong[10];


// if it will be more then one satellite needs to set this value to a last epoch of all satellites
GLOBAL_VARIABLE long double dStartTLEEpoch;

GLOBAL_VARIABLE long double TimeSl;
GLOBAL_VARIABLE long double TimeSl_2;

GLOBAL_VARIABLE long double Gbig;
GLOBAL_VARIABLE long double IterPerSec;

GLOBAL_VARIABLE long iTotalSec;
GLOBAL_VARIABLE long double StepsValInDay; // step's value in day measurement

GLOBAL_VARIABLE long double EarthCurTime;
GLOBAL_VARIABLE long double EarthCurTimeS;

GLOBAL_VARIABLE long double StartLandingIteraPerSec;

GLOBAL_VARIABLE char szTraVisualFileName[_MAX_PATH*3];

GLOBAL_VARIABLE BOOL VisualFileSet;

GLOBAL_VARIABLE char szURLTraVisualFileName[3*_MAX_PATH];
GLOBAL_VARIABLE char szURLTraVisualServer[3*_MAX_PATH];

GLOBAL_VARIABLE int UrlTraVisualPort;

GLOBAL_VARIABLE int bRGBImageW;
GLOBAL_VARIABLE int bRGBImageH;

GLOBAL_VARIABLE int iProfile;

GLOBAL_VARIABLE int iMaxSeq;
GLOBAL_VARIABLE int iMaxCounter;

GLOBAL_VARIABLE double dRGBScale;
GLOBAL_VARIABLE int RGBReferenceBody;

GLOBAL_VARIABLE int EngineToOptimize;
GLOBAL_VARIABLE int TrajectoryOptimizationType;

GLOBAL_VARIABLE int LastEngine;

GLOBAL_VARIABLE int MaxOptim;

GLOBAL_VARIABLE int StartOptim;

GLOBAL_VARIABLE int iOptimizationStep;

GLOBAL_VARIABLE PULSARS Pulsars[150];

GLOBAL_VARIABLE int nPulsars;

GLOBAL_VARIABLE long double EarthX;
GLOBAL_VARIABLE long double EarthY;
GLOBAL_VARIABLE long double EarthZ;

GLOBAL_VARIABLE long double EarthVX;
GLOBAL_VARIABLE long double EarthVY;
GLOBAL_VARIABLE long double EarthVZ;

GLOBAL_VARIABLE long double EarthR;
GLOBAL_VARIABLE long double EarthM;
GLOBAL_VARIABLE long double EarthSmAx;

GLOBAL_VARIABLE long double AU;

GLOBAL_VARIABLE long double EarthCalcKepler;

GLOBAL_VARIABLE long double GMSun;

GLOBAL_VARIABLE long double GMEarth;

GLOBAL_VARIABLE long double GMMoon;

GLOBAL_VARIABLE long double MoonM;

GLOBAL_VARIABLE long double MoonX;
GLOBAL_VARIABLE long double MoonY;
GLOBAL_VARIABLE long double MoonZ;
GLOBAL_VARIABLE long double MoonR;

GLOBAL_VARIABLE long double MoonVX;
GLOBAL_VARIABLE long double MoonVY;
GLOBAL_VARIABLE long double MoonVZ;

GLOBAL_VARIABLE char szTraCalcFileName[1024];

GLOBAL_VARIABLE char szURLServerCalulationOutputFile[1024]; // in CALC mode the output file with 

GLOBAL_VARIABLE int UrlTraCalcPort;

GLOBAL_VARIABLE char szURLTraCalcFileName[1024];
GLOBAL_VARIABLE char CalcinfoFile[_MAX_PATH];
GLOBAL_VARIABLE char szURLServerCalcinfoFilename[256];
GLOBAL_VARIABLE int UrlCalcinfoFilePort;
GLOBAL_VARIABLE char szURLCalcinfoFilename[256];

GLOBAL_VARIABLE int iMaxMeasures;

GLOBAL_VARIABLE char SimulationType[1024];
                          //    aloowed simulation modes
                          //    TLE - data == read TLE and calculate position on the specific time (TLE_EARTH_CRS,TLE_EARTH_TRS,TLE_SOLAR_CRS == 3 type of data
                          //    PING - data == simulate ping messages from ground station to satellite (at spesific time from all ground stations)
                          //    GPS  - data == simulate GPS's raw data from GPS satellites at specific time from GPS satellites
                          //    PULSAR - data == simulate PULSAR receving signal from all pulsars at specific time
                          // TLE and PING output data is the same: for PING:
                          //     Time=<time>, TimeErr=<time(d)>, PosX=<X(m)>, PosY=<Y(m)>, PosZ=<Z(m)>, PosErr=<error(m)> pingD1=<time(d) from GS->sat>, PingD11Err=<time(d)>
                          //         PingDel=<time(d) of processing data on sat>, PingDelErr=<time(d)>, PingD2=<time(d) from Sat to GS>, PingD2Err=<time(d)>
                          // for TLE:
                          //   Time=<time-from-simulation-output-time>, TimeErr=0, 
                          //           PosX=<X(m) position of the sat related to the earth>, PosY=<Y(m)>, PosZ=<Z(m)>, PosErr=<error(m) of the position>
                          //           pingD1=0, PingD11Err=0, PingDel=0, PingDelErr=0, PingD2=0, PingD2Err=0
                          // for GPS:
                          //  Time=<time-from-simulation-output-time>, GPS=<ID>, UTC=<time UTC from SAT>,
                          //       RawX=<X(m)> RawY=<Y(m)> RawZ=<Z(m)> RawVX=<VX(m)> RawVY=<VY(m)> RawVZ=<VZ(m)>
                          // for PULSAR:
                          //  Time=<time-from-simulation-output-time>, PULSAR=<ID>, DaltaT=<time btw two pulsars signals>

GLOBAL_VARIABLE char Mode[1024]; // allowed modes:
                          //    PROP == (default) proparation of the orbit
                          //    SIM ==  simulation of the data from satellite
                          //    CALC == calculation of the orbit based on any of the data : PING, GPS, PULSAR 
                          //    OPTIM = calculation of the engines firing to reach the moon

GLOBAL_VARIABLE char szTraSimFileName[1024];

GLOBAL_VARIABLE char szURLServerSimulationOutputFile[1024];
GLOBAL_VARIABLE int UrlTraSimPort; 
GLOBAL_VARIABLE char szURLTraSimFileName[1024];

GLOBAL_VARIABLE int SimulationOutputCount;
GLOBAL_VARIABLE long double SimulationOutputTime[MAX_OUTPUT_TIMES];
GLOBAL_VARIABLE char SimNAME[32][80];
GLOBAL_VARIABLE int iSimTotalLocations;
GLOBAL_VARIABLE long double SimLong[32];
GLOBAL_VARIABLE long double SimH[32];
GLOBAL_VARIABLE char szMassPointsModelFile[1024];

GLOBAL_VARIABLE char MidRandPointsFile[1024];
GLOBAL_VARIABLE long double MinH;
GLOBAL_VARIABLE long double MaxH;

GLOBAL_VARIABLE long double GST,SLONG,SRASN, SDEC;

void AssignAllSatelites(TRAOBJ * SlS, int iBody, TRAOBJ * sat, double JDSec);
GLOBAL_VARIABLE  long double ModelCoef;

GLOBAL_VARIABLE int iCounter_nk_lm_Numbers;
GLOBAL_VARIABLE char UseSatData[1024];         // allowed initial data:
                                   // SGP - use SGP from Space Track Report 3 to calculate position and velocity based on TLE
                                   // SGP4 - use SPG4 to calulate initial position and velocity from TLE
                                   // SGP8 - use SGP8
                                   // KEPLER - use internal "kepler" function to calculate position
                                   // INTERNAL use of internal data
                                   // ProbTime=<time>, ProbC=<count of integral iterations>, ProbT=<delta time of iterations>
                                   // ProbX0=<X0>, ProbY0=<Y0>, ProbZ0=<Z0>, ProbVX0=<VX0>, ProbVY0=<VY0>, ProbVZ0=<VZ0>, 
                                   //   ProbIFX1=<integral Fx>, ProbIFY1=<integral Fy>, ProbIFZ1=<integral Fz>
                                   //   ProbIFX2=<integral Fx>, ProbIFY2=<integral Fy>, ProbIFZ2=<integral Fz>
                                   //   ProbIFX3=<integral Fx>, ProbIFY3=<integral Fy>, ProbIFZ3=<integral Fz>
                                   //   ProbIVX1=<integral VX>, ProbIVY1=<integral VY>, ProbIVZ1=<integral VZ>
                                   //   ProbIVX2=<integral VX>, ProbIVY2=<integral VY>, ProbIVZ2=<integral VZ>
                                   //   ProbIVX3=<integral VX>, ProbIVY3=<integral VY>, ProbIVZ3=<integral VZ>
GLOBAL_VARIABLE long double Targetlongitude; // dolgota
GLOBAL_VARIABLE long double Targetlatitude; // shirota
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

GLOBAL_VARIABLE int EnginesCount;
GLOBAL_VARIABLE double EngCoeff;


#define MAX_ENGINES 6
#define MAX_IMPULSE_LINES 100


typedef struct TraImplObj
{
    int iLine;
    int iEngineOnSatellite;
    int EngineOn;
    int EngineDone;
    int ImplsPointer;
    int iCalculate;
    long double TotalImpulse;
    long double Weight;
    long double DeltaTime;
    int IteraPerSec;
    long double FireTime;
    long double Ang1;
    long double Ang2;
    long double XVec;
    long double YVec;
    long double ZVec;
    long double ValImpl[MAX_IMPULSE_LINES];
    int NearBody;
    long double TotalWeight;
    int AngleType;
    int AngleOnBody;
    long double OptimizationInitialStep;
    long double OptimizationDecCoef;
    long double OptimizationInitialStepCopy;
    long double OptimizationDecCoefCopy;
    long double OptimizationStop;
    int OptimizationFirstDirectionSwitch;
    long double SeartchForPeriod;
    int iCountApogPerig;
} TRAIMPLOBJ, *PTRAIMPLOBJ;

typedef struct TraOptimObj
{
    long double FireTime; // firing time of an engine
    long double Ang1; // first angle
    long double Ang2; // second angle
    long double XVec; // direction firing vector (X component)
    long double YVec; // direction (Y component)
    long double ZVec; // direction (Z component)
    int NearBody; 
    int AngleType;
    int AngleOnBody;
    //int OptimizationFirstDirectionSwitch;
    int EngineToOptimize;
    int TrajectoryOptimizationType;
    int LastEngine;
    int Calculate;
    long double OptimizationInitialStep;
    long double OptimizationDecCoef;
    long double OptimizationStop;
    int iNumberOfTryValues;
#define _VAL_TRY 30*24
    long double dValTry[_VAL_TRY]; 
    long double dValTryMaxMin[_VAL_TRY]; 
    //long double SeartchForPeriod;
    long double Period;
}TRAOPTIMOBJ, *PTRAOPTIMOBJ;

GLOBAL_VARIABLE TRAIMPLOBJ Engine[MAX_ENGINES];
GLOBAL_VARIABLE int iItaration;

GLOBAL_VARIABLE TRAOPTIMOBJ Opt[MAX_OPTIM];
GLOBAL_VARIABLE int iOptPtr;


