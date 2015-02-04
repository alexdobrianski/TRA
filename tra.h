#ifdef USE_GLOBAL
#define GLOBAL_VARIABLE
#else
#define GLOBAL_VARIABLE extern
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

#define NPULSARS 150

#define MAX_OPTIM 30

#define IMAGE_W 1280
#define IMAGE_H 720

typedef struct tagPulsars
{
    int N;
    char Name[20];
    long double ELONG;
    long double ELAT;
    long double P0;
    long double S400mJy;

} PULSARS, *PPULSARS;

void ConvertDateFromXML(char *pszQuo, long double &ld_TotalDays, long double &ld_dStartJD, long double &ld_dStartTLEEpoch);
void ParamSun(char *szString);
void init_tra_XML(void);
void ParamCommon(char *szString);
BOOL ParsURL(char * URLServer, int *port, char* URL,  char * szParsingName);

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
GLOBAL_VARIABLE int iItearationsPerSec; // that is "int" == IterPerSec
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
