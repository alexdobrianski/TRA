#include "stdafx.h"
#include "afxinet.h"
#include "afxsock.h"
#include <string.h>
#define _USE_MATH_DEFINES 1
#include <math.h>
#include <stdio.h>

#include "procXML.h"
#include "tra.h"

long double __AO_[2][7][7] = {     26.8629,     27.4598,     28.6395,      29.6418,     30.1671,     29.7578,     30.7854,
                                -0.451674,   -0.463668,   -0.490987,    -0.514957,   -0.527837,   -0.517915,   -0.545695,
                                0.00290397,    0.002974,  0.00320649,   0.00341926,  0.00353211,  0.00342699,  0.00370328,
                                -1.06953e-5,  -1.0753e-5,  -1.1681e-5,  -1.25785e-5, -1.30227e-5, -1.24137e-5, -1.37072e-5,
                                2.21598e-8,  2.17059e-8,  2.36847e-8,    2.5727e-8,  2.66455e-8,  2.48209e-8,  2.80614e-8,
                                -2.42941e-11,-2.30249e-11,-2.51809e-11, -2.75874e-11,-2.85432e-11,-2.58413e-11,-3.00184e-11,
                                1.09926e-14, 1.00123e-14, 1.09536e-14,  1.21091e-14, 1.25009e-14, 1.09383e-14, 1.31142e-14,

                                    17.8781,    -2.54909,    -13.9599,     -23.3079,    -14.7264,      -4.912,    -5.40952,
                                -0.132025,   0.0140064,   0.0844951,     0.135141,   0.0713256,   0.0108326,  0.00550749,
                                0.000227717, -0.00016946,-0.000328875, -0.000420802,-0.000228015, -8.10546e-5, -3.78851e-5,
                                -2.2543e-7,  3.27196e-7,  5.05918e-7,   5.73717e-7,   2.8487e-7,  1.15712e-7,   2.4808e-8,
                                1.33574e-10, -2.8763e-10,-3.92299e-10, -4.03238e-10,-1.74383e-10, -8.13296e-11,4.92183e-12,
                                -4.50458e-14, 1.22625e-13, 1.52279e-13,  1.42846e-13, 5.08071e-14, 3.04913e-14,-8.65011e-15,
                                6.72086e-18,-2.05736e-17,-2.35576e-17, -2.01726e-17,-5.34955e-18,-4.94989e-18,  1.9849e-18};

long double __L0_[2][5][7] =  { -0.407768,   -0.902739,    -0.73303,     -1.31444,    -1.20026,    -1.52158,     -1.67664,
                                0.00148506,  0.00826803,  0.00523396,    0.0133124,   0.0114087,    0.015704,   1.77194e-2,
                                1.25357e-5, -1.25448e-5,  6.35667e-6,  -2.55585e-5, -1.47324e-5, -3.02859e-5,  -3.69498e-5,
                                3.77311e-8,  6.12853e-8,  1.09065e-8,   5.43981e-8,   2.7804e-8,  4.57668e-8,   5.09134e-8,
                            -7.78953e-11,-7.07966e-11,-2.61427e-11, -4.33784e-11, -2.2632e-11,-2.82926e-11, -2.82878e-11,

                                    48.6536,     54.4867,     60.1267,      47.0996,     50.6174,     8.01942,     -15.5728,
                                -0.170291,   -0.178298,   -0.183144,     -0.12526,   -0.129047,   0.0185302,   9.36704e-2,
                                2.26242e-4,  2.22725e-4,  2.12481e-4,   1.26352e-4,  1.24842e-4, -6.14733e-5,  -1.49036e-4,
                                -1.32032e-7,   -1.227e-7, -1.08497e-7,  -5.51584e-8, -5.24993e-8,  4.97674e-8,   9.42151e-8,
                                2.85193e-11, 2.51316e-11,  2.0571e-11,  8.75272e-12, 8.08272e-12,-1.26162e-11,  -2.0961e-11};

long double __BO_[2][5][7] = {   0.0687894,     0.15073,   0.0479451,    0.0223448, -0.00326391,  -0.0514749,   -0.107255,
                                -0.00284077, -0.00400889, -0.00239453,   -0.0019798, -0.00159869,-0.000921059,-0.000174343,
                                1.83922e-5,  2.43937e-5,  1.70335e-5,   1.54101e-5,  1.40443e-5,  1.15147e-5,  9.02759e-6,
                                9.19605e-9, -9.92772e-9, -1.31626e-9,   -2.3543e-9, -3.02287e-9,- 1.22901e-9,-3.16512e-10,
                                -4.16873e-11,-1.82239e-11,-1.74032e-11, -1.24994e-11, -9.2016e-12,-8.13104e-12,   -6.14e-12,

                                    23.1584,     33.2732,     39.1961,     43.2469,      49.5738,      11.278,    -52.6184,
                                -0.0802147,   -0.111099,    -0.12352,   -0.126973,    -0.138613,  0.00143478,    0.214689,
                                0.000105824, 0.000141421, 0.000149015, 0.000142637,  0.000147851, -3.69846e-5,-0.000294882,
                                -6.15036e-8, -7.94952e-8,  -7.9705e-8, -7.09985e-8,  -6.96361e-8,  3.58318e-8,  1.71171e-7,
                                1.32453e-11, 1.65836e-11, 1.58772e-11, 1.31646e-11,  1.21595e-11,-9.91225e-12,-3.60582e-11};

void init_tra_XML(void)
{
    SunX = .0;
    SunY = .0;
    SunZ = .0;
    dStartJD = 0.0;//2451544.5; // if value dStartJD not set (==0.0) then use value from keplers elements of a satelite 0

    //NOTE, SunR does not get initialized

    dMinFromNow = 3.0;

    JustFlySimulation = 0;
    iGr = 0;

    TimeSl = 0;//0.01;
    TimeSl_2 = 0;

    Gbig = 0;//6.6725E-11;

    StartLandingIteraPerSec = 0.0;

    strcpy_s(szTraVisualFileName, "travisual.xml");

    VisualFileSet = FALSE;

    UrlTraVisualPort = 80;

    bRGBImageW = IMAGE_W;
    bRGBImageH = IMAGE_H;
    
    iProfile = 0; // 0 == XY , 1 == YZ, 2 == XZ 3 == -YZ 4 == -XZ 5==-XY
    // 0 or XY is a view from North to south, 5 (- XY) is a view from south to north
    // 1 or YZ is a view to easter

    iMaxSeq = 128;
    iMaxCounter = 24*60*60;
    dRGBScale = 1000000000.0;

    RGBReferenceBody = EARTH;

    EngineToOptimize = 4; // 0 == first engine firing - search for apogee 
                       // 1 == second engine firing - search for perigee
                       // 2 == third engine firing 
                       // 3 == 4th impulse
                       // 4 == 5th impulse
    TrajectoryOptimizationType = 0; 
          //1 - search for a minimum by adjusting time of firing
          //2 - search for a maximum by adjusting time of firing
          //3 - search for minimum by adjusting time of firing
          //4 - search fo minimum by adjusting angle of firing 

    LastEngine = 0;
    
    MaxOptim = MAX_OPTIM;

    StartOptim = 0;

    iOptimizationStep = 0;

    nPulsars = 0;

    EarthX = .0;
    EarthY = .0;
    EarthZ = .0;

    EarthVX = .0;
    EarthVY = .0;
    EarthVZ = .0;

    EarthR = 6371000.0;
    EarthM = 5.9736E24;

    MoonM = 7.3477E22;
    MoonX = 363104000.0;
    MoonY = .0;
    MoonZ = .0;
    MoonR = 1737100.0;

    MoonVX = .0;
    MoonVY = 1022.0;
    MoonVZ = .0;
    szURLServerCalulationOutputFile[0] = 0;

    UrlTraCalcPort = 0;
    CalcinfoFile[0] = 0;
    szURLServerCalcinfoFilename[0] = 0;
    UrlCalcinfoFilePort= 0;
    szURLCalcinfoFilename[0] = 0;

    szURLServerSimulationOutputFile[0] = 0;
    
    UrlTraSimPort = 0; 
    SimulationOutputCount = 0;
    iSimTotalLocations = 0;

    strcpy(Mode, "PROP");

//    __AH[2][7] = {            120,         120,         120,          120,         120,         120,         120,
//                              500,         500,         500,          500,         500,         500,         500};
    __AH[0][0] = 120;
    __AH[0][1] = 120;
    __AH[0][2] = 120;
    __AH[0][3] = 120;
    __AH[0][4] = 120;
    __AH[0][5] = 120;
    __AH[0][6] = 120;

    __AH[1][0] = 120;
    __AH[1][1] = 500;
    __AH[1][2] = 500;
    __AH[1][3] = 500;
    __AH[1][4] = 500;
    __AH[1][5] = 500;
    __AH[1][6] = 500;

//    __LH[2][7] = {            120,         120,         120,          120,        120,         120,         120,  
//                              640,         660,         740,          800,        860,         900,         900};

    __LH[0][0] = 120;
    __LH[0][1] = 120;
    __LH[0][2] = 120;
    __LH[0][3] = 120;
    __LH[0][4] = 120;
    __LH[0][5] = 120;
    __LH[0][6] = 120;

    __LH[1][0] = 640;
    __LH[1][1] = 660;
    __LH[1][2] = 740;
    __LH[1][3] = 800;
    __LH[1][4] = 860;
    __LH[1][5] = 900;
    __LH[1][6] = 900;

//    __BH[2][7] = {            120,         120,         120,          120,         120,         120,         120,
//                              600,         660,         760,          800,         860,         920,         980};

    __BH[0][0] = 120;
    __BH[0][1] = 120;
    __BH[0][2] = 120;
    __BH[0][3] = 120;
    __BH[0][4] = 120;
    __BH[0][5] = 120;
    __BH[0][6] = 120;

    __BH[1][0] = 600;
    __BH[1][1] = 660;
    __BH[1][2] = 760;
    __BH[1][3] = 800;
    __BH[1][4] = 860;
    __BH[1][5] = 920;
    __BH[1][6] = 980;



    memcpy(__AO, __AO_, sizeof(__AO));
    memcpy(__L0, __L0_, sizeof(__L0));
    memcpy(__BO, __BO_, sizeof(__BO));

    CpuCore = 0;

    //    EarthModelFile[1024]={"egm96"};
    strcpy(EarthModelFile, "egm96");

    EarthSmoothCoefStart = 0;
    
#ifdef USE_MODEL_LOAD

    EarthModelCoefs = 16;

    GM_MODEL = 398600.4418e9;
    R0_MODEL =  6378137.00;
#endif
#ifdef USE_MODEL_0
    GM_MODEL = 398600.4415E9;
    R0_MODEL =  6378137.0;
#endif

#ifdef USE_MODEL_1
    GM_MODEL = 398600.4415E9;
    R0_MODEL= 6378136.30;
#endif
#ifdef USE_MODEL_2
    GM_MODEL = 398600.4415E9;
    R0_MODEL= 6378136.30;
#endif
#ifdef USE_MODEL_3
    GM_MODEL = 398600.4418E95;
    R0_MODEL= 6378137;
#endif 

}

void ParamSun(char *szString)
{
    XML_BEGIN;
    XML_SECTION(TraInfo);
        
        XML_READ(SunX);
        XML_READ(SunY);
        XML_READ(SunZ);
		XML_READ(SunR);
        //XML_READ(SunM);
        //IF_XML_READ(GMSun) 
        //{
        //    GMSun = atof(pszQuo);
        //    printf("\n Was  SunM= %f", SunM);
        //    SunM = GMSun / Gbig;
        //    printf("\n Calc SunM= %f", SunM);
        //}
    XML_SECTION_END;
    XML_END;
}

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
        IF_XML_READ(JustFlySimulation)
        {
            JustFlySimulation = atoi(pszQuo); // if 0 then do enegine firing ; if 1 no engines at all
        }
        IF_XML_READ(dStartJD) 
        {
            ConvertDateFromXML(pszQuo, TotalDays, dStartJD, dStartTLEEpoch);
        }
        XML_READ(TimeSl);
        XML_READ(Gbig);
        //XML_READ(IterPerSec);
        IF_XML_READ(IterPerSec)
        {
             IterPerSec = atol(pszQuo);
             iItearationsPerSec = (int)(IterPerSec);
             TimeSl = 1.0 / IterPerSec;
             TimeSl_2 = 1.0 / ((long double)IterPerSec*(long double)IterPerSec);
             StepsValInDay = (1.0/((long double)IterPerSec))/24.0/60.0/60.0;
             printf("\n IterPerSec =%d ", (int)IterPerSec);
        }
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
        //IF_XML_READ(EarthSmAxAU)
        //{
        //    EarthSmAxAU = atof(pszQuo);
        //    AUcalc = EarthSmAx / EarthSmAxAU;
        //}
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
            Pulsars[nPulsars].ELONG = atof(pszQuo);
        }
        IF_XML_READ(ELAT)
        {
            Pulsars[nPulsars].ELAT = atof(pszQuo);
        }
        IF_XML_READ(P0)
        {
            Pulsars[nPulsars].P0 = atof(pszQuo);
        }
        IF_XML_READ(S400mJy)
        {
            Pulsars[nPulsars].S400mJy = atof(pszQuo);
            if (++nPulsars >= NPULSARS)
                nPulsars = NPULSARS-1;
        }
    XML_SECTION_END;
    XML_END;
}

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
        //XML_READ(EarthRP);
        //XML_READ(EarthRE);
        //XML_READ(EarthM);
        XML_READ(AU);
        //XML_READ(EarthTSolSec);
        //XML_READ(EarthSmAx);
        //double Temp = 0.0;
        //IF_XML_READ(EarthTDays)  
        //{
        //    
        //    EarthTDays  = atof(pszQuo);
        //    printf("\n calc EarthTDays=%f ", EarthTDays);
        //    EarthTSec = EarthTDays * 24.0*60.0*60.0;
        //    printf("\n calc EarthTSec =%f ", EarthTSec);
        //    Temp = EarthTSec / EarthTSolSec;

        //}

        //XML_READ(GMEarth);
        //XML_READ(MassRatioSunToEarthPlusMoon);
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
        //XML_READ(MoonRP);
        //XML_READ(MoonRE);
        //XML_READ(MoonM);
        //XML_READ(GMMoon);
        //XML_READ(GMEarthMoon);


        //IF_XML_READ(MassRatioEarthToMoon)
        //{
        //    MassRatioEarthToMoon = atof(pszQuo);
        //    //if (MassRatioSunToEarthPlusMoon)
        //    {
        //        double Temp = SunM / MassRatioSunToEarthPlusMoon;
        //        if (GMEarth != 0.0 && GMMoon != 0.0 && GMSun != 0.0)
        //        {
        //            printf("\n calculation mass based on GMSun GMMoon and GMEarth");
        //            printf("\n Was  GMSun = %f", GMSun);
        //            printf("\n Was  GMEarth = %f", GMEarth);
        //            printf("\n Was  GMMoon = %f", GMMoon);
        //            printf("\n Was  SunM = %f", SunM);
        //            SunM = GMSun /Gbig;
        //            printf("\n Calc SunM = %f", SunM);                    
        //            printf("\n Was  EarthM = %f", EarthM);
        //            EarthM = GMEarth / Gbig;
        //            printf("\n Calc EarthM = %f", EarthM);
        //            printf("\n Was  MoonM = %f", MoonM);
        //            MoonM = GMMoon / Gbig;
        //            printf("\n Calc MoonM = %f", MoonM);
        //        }
        //        else
        //        {
        //            printf("\n calculation mass based on mass ratio earth to moon");
        //            printf("\n Was  MoonM = %f", MoonM);
        //            MoonM = Temp / (1.0 + MassRatioEarthToMoon);
        //            printf("\n Calc MoonM = %f", MoonM);
        //            printf("\n Was  EarthM = %f", EarthM);
        //            EarthM = MoonM * MassRatioEarthToMoon;
        //            printf("\n Calc EarthM = %f", EarthM);
        //        }
        //    }
        //}
    XML_SECTION_END;
    XML_END;
}

void ParamProb(char *szString)
{
    
    char szTemp[128] = {"0."};;
    char valX, valY, valZ, valT;
    XML_BEGIN;
    ///////////////////////////////////////////////////////////////////////////////////////Calcinfo
    XML_SECTION(CalcInfo)
    IF_XML_READ(CalculationOutputFile)
    {
        strcpy(szTraCalcFileName, pszQuo);
        if (strchr(szTraCalcFileName, '\"'))
            *strchr(szTraCalcFileName, '\"')=0;
        if (ParsURL(szURLServerCalulationOutputFile, &UrlTraCalcPort, szURLTraCalcFileName,  szTraCalcFileName))
        {
        }
        else
            UrlTraCalcPort = 0;
    }
    IF_XML_READ(CalcinfoFile)
    {
        strcpy(CalcinfoFile, pszQuo);
        if (strchr(CalcinfoFile, '\"'))
            *strchr(CalcinfoFile, '\"')=0;
        if (ParsURL(szURLServerCalcinfoFilename, &UrlCalcinfoFilePort, szURLCalcinfoFilename,  CalcinfoFile))
        {
        }
        else
            UrlCalcinfoFilePort = 0;
        
    }
    XML_SECTION_GROUP_SEPARATOR
    //===============================================================gCRSmeasure Geocent Cel ref Sys
    XML_GROUP(gTRSmeasure)
        IF_XML_ELEMENT(T)
        {
            if (++iMaxMeasures>MAX_MEASURES)
                iMaxMeasures--;
            long double ld1=0,ld2=0;
            ConvertDateFromXML(pszQuo, ld1, measures[iMaxMeasures-1].T, ld2);
            measures[iMaxMeasures-1].TypeOfmesaure = 1;
            measures[iMaxMeasures-1].X =0;
            measures[iMaxMeasures-1].Y =0;
            measures[iMaxMeasures-1].Z =0;
        }
        IF_XML_ELEMENT(M)
            measures[iMaxMeasures-1].NearBody =atoi(pszQuo);
        IF_XML_ELEMENT(X)
            measures[iMaxMeasures-1].X =atof(pszQuo);
        IF_XML_ELEMENT(Y)
            measures[iMaxMeasures-1].Y =atof(pszQuo);
        IF_XML_ELEMENT(Z)
            measures[iMaxMeasures-1].Z =atof(pszQuo);
        IF_XML_ELEMENT(H)
            measures[iMaxMeasures-1].H =atof(pszQuo);
        IF_XML_ELEMENT(LAT)
            measures[iMaxMeasures-1].LAT =atof(pszQuo);
        IF_XML_ELEMENT(LON)
            measures[iMaxMeasures-1].LON =atof(pszQuo);
        IF_XML_ELEMENT(E)
            measures[iMaxMeasures-1].Err =atof(pszQuo);
        IF_XML_ELEMENT(D1)
            measures[iMaxMeasures-1].D1 =atof(pszQuo);
        IF_XML_ELEMENT(E1)
            measures[iMaxMeasures-1].Err1 =atof(pszQuo);
        IF_XML_ELEMENT(T2)
            measures[iMaxMeasures-1].T2 =atof(pszQuo);
        IF_XML_ELEMENT(E2)
            measures[iMaxMeasures-1].ErrT2 =atof(pszQuo);
    XML_GROUP_END

    XML_GROUP(gCRSmeasure)

        IF_XML_ELEMENT(T)
        {
            if (++iMaxMeasures>MAX_MEASURES)
                iMaxMeasures--;
            long double ld1=0,ld2=0;
            ConvertDateFromXML(pszQuo, ld1, measures[iMaxMeasures-1].T, ld2);
            measures[iMaxMeasures-1].TypeOfmesaure = 2;
        }
        IF_XML_ELEMENT(M)
            measures[iMaxMeasures-1].NearBody =atoi(pszQuo);
        IF_XML_ELEMENT(X)
            measures[iMaxMeasures-1].X =atof(pszQuo);
        IF_XML_ELEMENT(Y)
            measures[iMaxMeasures-1].Y =atof(pszQuo);
        IF_XML_ELEMENT(Z)
        {
            measures[iMaxMeasures-1].Z =atof(pszQuo);
            measures[iMaxMeasures-1].H =0;
            measures[iMaxMeasures-1].LAT =0;
            measures[iMaxMeasures-1].LON =0;
        }
        IF_XML_ELEMENT(E)
            measures[iMaxMeasures-1].Err =atof(pszQuo);
        IF_XML_ELEMENT(D1)
            measures[iMaxMeasures-1].D1 =atof(pszQuo);
        IF_XML_ELEMENT(E1)
            measures[iMaxMeasures-1].Err1 =atof(pszQuo);
        IF_XML_ELEMENT(T2)
            measures[iMaxMeasures-1].T2 =atof(pszQuo);
        IF_XML_ELEMENT(E2)
            measures[iMaxMeasures-1].ErrT2 =atof(pszQuo);
    XML_GROUP_END

    XML_GROUP(hCRSmeasure)
        IF_XML_ELEMENT(T)
        {
            if (++iMaxMeasures>MAX_MEASURES)
                iMaxMeasures--;
            long double ld1=0,ld2=0;
            ConvertDateFromXML(pszQuo, ld1, measures[iMaxMeasures-1].T, ld2);
            measures[iMaxMeasures-1].TypeOfmesaure = 3;
            measures[iMaxMeasures-1].X =0;
            measures[iMaxMeasures-1].Y =0;
            measures[iMaxMeasures-1].Z =0;
            measures[iMaxMeasures-1].H =0;
            measures[iMaxMeasures-1].LAT =0;
            measures[iMaxMeasures-1].LON =0;
            measures[iMaxMeasures-1].Err =0;
            measures[iMaxMeasures-1].D1 =0;
            measures[iMaxMeasures-1].Err1 =0;
            measures[iMaxMeasures-1].T2 =0;
            measures[iMaxMeasures-1].ErrT2 =0;
        }
        IF_XML_ELEMENT(M)
            measures[iMaxMeasures-1].NearBody =atoi(pszQuo);
        IF_XML_ELEMENT(X)
            measures[iMaxMeasures-1].X =atof(pszQuo);
        IF_XML_ELEMENT(Y)
            measures[iMaxMeasures-1].Y =atof(pszQuo);
        IF_XML_ELEMENT(Z)
            measures[iMaxMeasures-1].Z =atof(pszQuo);
        IF_XML_ELEMENT(H)
            measures[iMaxMeasures-1].H =atof(pszQuo);
        IF_XML_ELEMENT(LAT)
            measures[iMaxMeasures-1].LAT =atof(pszQuo);
        IF_XML_ELEMENT(LON)
            measures[iMaxMeasures-1].LON =atof(pszQuo);
        IF_XML_ELEMENT(E)
            measures[iMaxMeasures-1].Err =atof(pszQuo);
        IF_XML_ELEMENT(D1)
            measures[iMaxMeasures-1].D1 =atof(pszQuo);
        IF_XML_ELEMENT(E1)
            measures[iMaxMeasures-1].Err1 =atof(pszQuo);
        IF_XML_ELEMENT(T2)
            measures[iMaxMeasures-1].T2 =atof(pszQuo);
        IF_XML_ELEMENT(E2)
            measures[iMaxMeasures-1].ErrT2 =atof(pszQuo);
    XML_GROUP_END

    XML_GROUP(hPULSARmeasure)
        IF_XML_ELEMENT(M)
        {
            if (++iMaxMeasures>MAX_MEASURES)
                iMaxMeasures--;
            measures[iMaxMeasures-1].NearBody =atoi(pszQuo);
            measures[iMaxMeasures-1].TypeOfmesaure = 4; // pulsar type
            measures[iMaxMeasures-1].X =0;
            measures[iMaxMeasures-1].Y =0;
            measures[iMaxMeasures-1].Z =0;
            measures[iMaxMeasures-1].H =0;
            measures[iMaxMeasures-1].LAT =0;
            measures[iMaxMeasures-1].LON =0;
            measures[iMaxMeasures-1].Err =0;
            measures[iMaxMeasures-1].D1 =0;
            measures[iMaxMeasures-1].Err1 =0;
            measures[iMaxMeasures-1].T2 =0;
            measures[iMaxMeasures-1].ErrT2 =0;
        }

        IF_XML_ELEMENT(T)
        {
            long double ld1=0,ld2=0;
            ConvertDateFromXML(pszQuo, ld1, measures[iMaxMeasures-1].T, ld2);
        }
        IF_XML_ELEMENT(P1)
            measures[iMaxMeasures-1].P1 =atof(pszQuo); // first pulsar Measurement period
        IF_XML_ELEMENT(P2)
            measures[iMaxMeasures-1].P2 =atof(pszQuo); // second pulsar Measurement period
        IF_XML_ELEMENT(P3)
            measures[iMaxMeasures-1].P3 =atof(pszQuo); // third pulsar Measurement period
        IF_XML_ELEMENT(E)
            measures[iMaxMeasures-1].Err =atof(pszQuo); // error in time of measurements
    XML_GROUP_END

    XML_SECTION_END

    ///////////////////////////////////////////////////////////////////////////////////////SimInfo
    XML_SECTION(SimInfo)
    IF_XML_READ(SimulationType)
    {
        strcpy(SimulationType, pszQuo);
        if (strchr(SimulationType, '\"'))
            *strchr(SimulationType, '\"')=0;
    }
    IF_XML_READ(szTraSimFileName)
    {
        strcpy(szTraSimFileName, pszQuo);
        if (strchr(szTraSimFileName, '\"'))
            *strchr(szTraSimFileName, '\"')=0;

        if (ParsURL(szURLServerSimulationOutputFile, &UrlTraSimPort, szURLTraSimFileName,  szTraSimFileName))
        {
        }
        else
            UrlTraSimPort = 0;

    }
    IF_XML_READ(SimulationOutputTime)
    {
        if (++SimulationOutputCount <= MAX_OUTPUT_TIMES)
        {
            long double ld1=0,ld2=0;
            ConvertDateFromXML(pszQuo, ld1, SimulationOutputTime[SimulationOutputCount-1], ld2);
        }
        else
        {
            SimulationOutputCount--;
            printf("\n Max output times for simulation limit reached - .XML file is incorrect");
        }
    }
    IF_XML_READ(SimNAME)
    {
        strcpy(SimNAME[iSimTotalLocations], pszQuo);
        if (strchr(SimNAME[iSimTotalLocations], '\"'))
            *strchr(SimNAME[iSimTotalLocations], '\"')=0;
    }
    IF_XML_READ(SimLat)
    {
        SimLat[iSimTotalLocations]=  atof(pszQuo);
    }
    IF_XML_READ(SimLong)
    {
        SimLong[iSimTotalLocations]=  atof(pszQuo);
    }
    IF_XML_READ(SimH)
    {
        SimH[iSimTotalLocations++]=  atof(pszQuo);
    }


    XML_SECTION_END

	///////////////////////////////////////////////////////////////////////////////////////MassInfo
    XML_SECTION(MassInfo)
    IF_XML_READ(ModelOutput)
    {
        strcpy(szMassPointsModelFile, pszQuo);
        if (strchr(szMassPointsModelFile, '\"'))
            *strchr(szMassPointsModelFile, '\"')=0;
    }
    IF_XML_READ(MidRandPointsFile)
    {
        strcpy(MidRandPointsFile, pszQuo);
        if (strchr(MidRandPointsFile, '\"'))
            *strchr(MidRandPointsFile, '\"')=0;
    }
    IF_XML_READ(MinH)
    {
        MinH = atof(pszQuo);
    }

    IF_XML_READ(MaxH)
    {
        MaxH = atof(pszQuo);
    }
    XML_SECTION_END

    ///////////////////////////////////////////////////////////////////////////////////////ModeInfo
    XML_SECTION(ModeInfo)
    IF_XML_READ(Mode)
    {
        strcpy(Mode, pszQuo);
        if (strchr(Mode, '\"'))
            *strchr(Mode, '\"')=0;
    }

    XML_SECTION_END
    ////////////////////////////////////////////////////////////////////////////////////////probs
    XML_SECTION(probs)
    // pisition of the prob (active by firing engines) can be set by X,Y,Z   VX,VY,VZ
        IF_XML_READ(ProbDX)
        {
            Sat.DX[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbDY)
        {
            Sat.DY[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbDZ)
        {
            Sat.DZ[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbDVX)
        {
            Sat.DVX[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbDVY)
        {
            Sat.DVY[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbDVZ)
        {
            Sat.DVZ[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbM)
        {
            // has to be first!!!!
            Sat.M[Sat.Elem] = atof(pszQuo);
        }
        // or by 3 punch card
        IF_XML_READ(ProbKeplerLine1)
        {
            strcpy(Sat.Kepler1[Sat.Elem], pszQuo);
            Sat.DX[Sat.Elem] = 0;Sat.DY[Sat.Elem] = 0;Sat.DZ[Sat.Elem] = 0;Sat.DVX[Sat.Elem] = 0;Sat.DVY[Sat.Elem] = 0;Sat.DVZ[Sat.Elem] = 0;
        }
        IF_XML_READ(ProbKeplerLine2)
        {
            strcpy(Sat.Kepler2[Sat.Elem], pszQuo);
        }
        IF_XML_READ(ProbSquare)
        {
            Sat.ProbSquare[Sat.Elem] = atof(pszQuo);
        }
        IF_XML_READ(ProbKeplerLine3)
        {
            char szTempo[1024];
            //int iYear;
            //int iDays;
            //double dflTemp;

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
            Sat.ProbTLEEpoch[Sat.Elem] = atof(&Sat.Kepler2[Sat.Elem][18]);
            Sat.ProbJDSec[Sat.Elem] = Sat.ProbJD[Sat.Elem] = ConverTLEEpochDate2JulianDay(Sat.ProbTLEEpoch[Sat.Elem]);
			Sat.ProbJDSec[Sat.Elem] *=  60*60*24;
            SYSTEMTIME ThatTime;
            // now for initial step asign emulation starting time:
			// if nothing is set then starting point is a last satellite epoch
            if (dStartJD == 0.0)
            {
                dStartJD = ConverTLEEpochDate2JulianDay(Sat.ProbTLEEpoch[Sat.Elem]);
                dStartTLEEpoch = ConvertJulianDayToDateAndTime(dStartJD, &ThatTime);
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
            if (Sat.ProbSquare[Sat.Elem] <0) // if munus value then needs to calculate square area:
            {
                Sat.ProbSquare[Sat.Elem] = Sat.ProbDragterm[Sat.Elem] / 0.15696615 *1000000;
                printf("\n Sat %d square area was calculated as =%f",Sat.Elem,Sat.ProbSquare[Sat.Elem]);
            }
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
    XML_SECTION_END

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////TraInfo
    XML_SECTION(TraInfo)
        IF_XML_READ(EarthModelFile)
        {
            strcpy(EarthModelFile, pszQuo);
            if (strchr(EarthModelFile, '\"'))
            {
                *strchr(EarthModelFile, '\"')=0;
            }
        }
        IF_XML_READ(EarthModelCoefs)
        {
            EarthModelCoefs = atoi(pszQuo);
        }
		IF_XML_READ(EarthSmoothCoefStart)
		{
			EarthSmoothCoefStart = atoi(pszQuo);
			//EarthSmoothCoefStart = 0;
		}
        IF_XML_READ(CpuCore)
        {
            CpuCore = atoi(pszQuo);
        }
        IF_XML_READ(GM_MODEL)
        {
            GM_MODEL = atof(pszQuo);
        }
        IF_XML_READ(R0_MODEL)
        {
            R0_MODEL = atof(pszQuo);
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
			
			// just cheking GST must be eq GreenwichA
			SUN_08 (1950,1,0,0,0,GST,SLONG,SRASN,SDEC);
			dStartGreenwichA = GreenwichAscensionFromTLEEpoch(50000.0,Sat.precEps,Sat.precTet,Sat.precZ,Sat.nutEpsilon,Sat.nutDFeta); // 50 000.0 == 1950, day 1( first of january), 00:00:00 
            // error is:
            //   GST =  1.7466460286076149
            // dStartGreenwichA = -0.0000017029127985
            // error in seconds== -0.01561115459725791718482933112546
            SYSTEMTIME ThatTime;
            int YY_Days = ConvertJulianDayToDateAndTime(dStartJD, &ThatTime);
            int Days_in_Year = YY_Days % 1000;
            // now for <TRA:setting name="dStartJD" value="07/05/14 08:28:00:000" />
            // dStartJD = 2456784.8527777777
            // dStartEpoch = 14127.352777777705
            SUN_08 (ThatTime.wYear,Days_in_Year,ThatTime.wHour,ThatTime.wMinute,ThatTime.wSecond,GST,SLONG,SRASN,SDEC);
            dStartGreenwichA = GreenwichAscensionFromTLEEpoch(dStartTLEEpoch,Sat.precEps,Sat.precTet,Sat.precZ,Sat.nutEpsilon,Sat.nutDFeta);
            // error is:
            // GST = 6.1453506996084499
            // dStartGreenwichA 6.1454312968999147
            // error in sec 0.7388615425790099904549342040562

            // all sat positions will be in absolut coordinates 
            // for first sat TLE (at Epoch time) earth was in one position
            // at dStartEpoch it will be on another
            // to account in graviational potential the shape of the earth 
            // needs to know this position
            //dStartGreenwichA = GreenwichAscensionFromTLEEpoch(Sat.ProbEpoch[0],Sat->precEps,Sat->precTet,Sat->precZ,Sat->nutEpsilon,Sat->nutDFeta);
            AssignAllSatelites(&SolarSystem, EARTH, &Sat, dStartJD);
            // the amoun of yhe core can be grabed from kernel - but it will be nice to have a control


            ModelCoef = 1;
            
#ifdef USE_MODEL_LOAD
            Sat.iLeg = EarthModelCoefs;
            iCounter_nk_lm_Numbers =0;
            FILE *FileC_S = fopen(EarthModelFile,"r");
            //FILE *FileC_S = fopen("JGM3.txt","r");
            if (FileC_S == NULL)
            {
                printf("\n file with C and S dose not exsists");
                exit(555);
            }
            char szTemp[512];
            Sat._SQRT3= sqrt((long double)3.0);
            Sat._p_n_k[1]= Sat._SQRT3;
            long double _P_N_K = 3.0;

            for (int n = 2 ; n <= Sat.iLeg; n++)
            {
                //l	m	            C                                S
	            //2	0	  -0.10826360229840D-02	                       0.0
                long double CNK;
                long double SNK;
                long double Factor1;
                long double Factor2;
                long double Factor3;
                long double Betta;
                for (int k = 0; k <= n; k++)
                {
                    memset(szTemp, 0, sizeof(szTemp));
                    fgets(szTemp, sizeof(szTemp), FileC_S);
                    char *ptrD = szTemp;
                    int iDataCount = 0;
                    for (int ic= 0; ic <sizeof(szTemp); ic++,ptrD++)
                    {
                        if (*ptrD != ' ')
                        {
                            switch(iDataCount)
                            {
                            case 0:iDataCount++; nk_lm_Numbers[iCounter_nk_lm_Numbers][0] = atoi(ptrD);  if (nk_lm_Numbers[iCounter_nk_lm_Numbers][0] != n) { printf("\n worng C S file"); exit(777); }
                                break;
                            case 2:iDataCount++; nk_lm_Numbers[iCounter_nk_lm_Numbers][1] = atoi(ptrD);  if (nk_lm_Numbers[iCounter_nk_lm_Numbers][1] != k) { printf("\n worng C S file"); exit(777); }
                                break;
                            case 4:iDataCount++; C_S_nk[iCounter_nk_lm_Numbers][0] = atof(ptrD);  
                                break;
                            case 6:iDataCount++; C_S_nk[iCounter_nk_lm_Numbers][1] = atof(ptrD);  
                                goto DONE_WITH_LINE;
                                break;
                            }
                        }
                        else
                        {
                            switch(iDataCount)
                            {
                            case 1: iDataCount++;break;
                            case 3: iDataCount++;break;
                            case 5: iDataCount++;break;
                            case 7: iDataCount++;break;
                            }
                        }
                    }
DONE_WITH_LINE:
                    if(EarthSmoothCoefStart > 0 && EarthSmoothCoefStart < EarthModelCoefs) 
                    {
					    if(nk_lm_Numbers[iCounter_nk_lm_Numbers][0] >= EarthSmoothCoefStart && nk_lm_Numbers[iCounter_nk_lm_Numbers][1] >= EarthSmoothCoefStart) 
					    {
						    //method : Low Pass Filtering of Gravity Field Models by Gently Cutting the Spherical Harmonic Coe±cients of Higher Degrees
						    //C_S_nk[iCounter_nk_lm_Numbers][0]
						    long double coefL = pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][0]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)4.0 ) - 
                                                (long double)(2) * (pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][0]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)2.0 )) + 1;
						


					    	long double coefM = pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][1]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)4.0 ) - 
                                                (long double)(2) * (pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][1]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)2.0 )) + 1;
						    C_S_nk[iCounter_nk_lm_Numbers][0] *= coefL * coefM;
                            C_S_nk[iCounter_nk_lm_Numbers][1] *= coefL * coefM;
						

						    //C_S_nk[iCounter_nk_lm_Numbers][1]
					    }
					    else if(nk_lm_Numbers[iCounter_nk_lm_Numbers][0] >= EarthSmoothCoefStart)
					    {
                            long double coefL = pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][0]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)4.0 ) - 
                                                (long double)(2) * (pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][0]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)2.0 )) + 1;
						    C_S_nk[iCounter_nk_lm_Numbers][0] *= coefL;
                            C_S_nk[iCounter_nk_lm_Numbers][1] *= coefL;
					    }
					    else if(nk_lm_Numbers[iCounter_nk_lm_Numbers][1] >= EarthSmoothCoefStart) 
					    {
						    long double coefM = pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][1]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)4.0 ) - 
                                                (long double)(2) * (pow( (long double)(((long double)(nk_lm_Numbers[iCounter_nk_lm_Numbers][1]) - (long double)(EarthSmoothCoefStart)) / ((long double)(EarthModelCoefs) - (long double)(EarthSmoothCoefStart))), (long double)2.0 )) + 1;
						    C_S_nk[iCounter_nk_lm_Numbers][0] *= coefM;
                            C_S_nk[iCounter_nk_lm_Numbers][1] *= coefM;
					    }
                    }
                    Factor1 =1.0;
                    Factor2 =1.0;
                    Factor3 =1.0;
                    int m,m1;
                    for (m= n-k; m >=1; m--)
                    {
                        Factor1 *= (long double)m;
                    }
                    for (m= n+k; m >=1; m--)
                    {
                        Factor2 *= (long double)m;
                    }

                    
                    for (m = n-k, m1= n+k; (m >=1) || (m1 >=1); m--,m1--)
                    {
                        if ((m >= 1) && (m1 >=1))
                        {
                            Factor3 *= (long double)m/(long double)m1;
                        }
                        else if (m >= 1)
                        {
                            Factor3 *= (long double)m;
                        }
                        else // (m1 >=1)
                        {
                            Factor3 /= (long double)m1;
                        }

                    }

                    // CNK = sqrt(2*(2*n+1)) * sqrt(((n-k)!/(n+k)!) * Clm
                    // SNK = sqrt(2*(2*n+1) * sqrt((n-k)!/(n+k)!)) * Slm
                    if (k == 0)
                        Betta = 1.0;
                    else
                        Betta = 2.0;
                    
                     if (n==2 && k ==2)
                     {
                         //C_S_nk[iCounter_nk_lm_Numbers][0]+=1.3*1.39e-8;
                         //Sat.CNK[n][k]+=1.3*1.39e-8;
                     }
#ifndef _NORMALIZED_COEF
                    CNK = sqrt(Betta*(2*(long double)n+1) * Factor1/Factor2) * C_S_nk[iCounter_nk_lm_Numbers][0];

                     CNK = sqrt(Betta*(2*(long double)n+1) * Factor3) * C_S_nk[iCounter_nk_lm_Numbers][0];
                     
                    SNK = sqrt(Betta*(2*(long double)n+1) * Factor1/Factor2) * C_S_nk[iCounter_nk_lm_Numbers][1];

                     SNK = sqrt(Betta*(2*(long double)n+1) * Factor3) * C_S_nk[iCounter_nk_lm_Numbers][1];
#endif
                     if (k==0)
                        Sat._pt_nk[n][k] = sqrt((long double)n*(n + 1)/2.0); // z
                     else
                     {
                         if ((n-k)*(n+k+1) == 0)
                            Sat._pt_nk[n][k] = 0.0;
                         else
                            Sat._pt_nk[n][k] = sqrt((long double)(n-k)*(n+k+1));
                     }
                     if (k!=n)
                     {
                        Sat._tp_nm1_k[n][k] = sqrt((long double)((2*n-1)*(2*n+1))/((long double)((n+k)*(n-k))));  // xin
                        if (k==0)
                             Sat._tp_nm2_k[n][k]= 0;  //eta
                        else
                             Sat._tp_nm2_k[n][k]= sqrt((long double)((2*n+1)*(n+k-1)*(n-k-1))/(long double)((n+k)*(n-k)*(2*n-3)));
                     }
                     
                     if (iCounter_nk_lm_Numbers == 0 && n < 2)
                         continue;
                     nk_lm_Numbers[iCounter_nk_lm_Numbers][0] = n;
                     nk_lm_Numbers[iCounter_nk_lm_Numbers][1] = k;
#ifndef _NORMALIZED_COEF
                     C_S_nk[iCounter_nk_lm_Numbers][0] =CNK;
                     C_S_nk[iCounter_nk_lm_Numbers][1] =SNK;
#endif
                     /* some starange corrections
                     if (n==2 && k == 0)
                        C_S_nk[iCounter_nk_lm_Numbers][0] =0.108262982131e-2;
                     if (n==4 && k == 0)
                        C_S_nk[iCounter_nk_lm_Numbers][0] =-.237091120053e-05;
                     if (n==6 && k == 0)
                        C_S_nk[iCounter_nk_lm_Numbers][0] =0.608346498882e-8;
                     if (n==8 && k == 0)
                        C_S_nk[iCounter_nk_lm_Numbers][0] =-0.142681087920e-10;
                     if (n==10 && k == 0)
                        C_S_nk[iCounter_nk_lm_Numbers][0] =0.121439275882e-13;
                    */
                     if (iCounter_nk_lm_Numbers++ >= TOTAL_COEF*TOTAL_COEF)
                     {
                         printf("\n something wrong with coef");
                         exit(555);
                     }
                }
                Sat._p_n_m_1[n] = sqrt((long double)(2*n +1)*(2*n -1))/long double(n);  // alfa
                Sat._p_n_m_2[n] = sqrt((long double)(2*n +1)/(long double)(2*n -3)) * (long double)(n -1)/long double(n); // betta
                Sat.diagonal[n] = _P_N_K;
                _P_N_K *= (long double)(2*n+1);
                Sat._tpk_n_k[n] = sqrt((long double)(2*n+1))*Sat._p_n_k[n-1];
                Sat._p_n_k[n] = sqrt((long double)(2*n+1)/(long double)(2*n)) * Sat._p_n_k[n-1];
            }

#else
            // amount of J coeff used in calcualtion
#ifdef USE_MODEL_0
            Sat.iLeg = 6;
#endif
#ifdef USE_MODEL_1
            Sat.iLeg = 8;
#endif
#ifdef USE_MODEL_2
            Sat.iLeg = 16;
#endif
#ifdef USE_MODEL_3
            Sat.iLeg = 16;
#endif

            iCounter_nk_lm_Numbers =0;
            for (int n = 0 ; n <= Sat.iLeg; n++)
            {
                
                //l	m	            C                                S
	            //2	0	  -0.10826360229840D-02	                       0.0
                Sat.J[n] = //sqrt(2*(long double)n+1)* // coeff already normalized
                    (Clm[0][n]);
                long double CNK;
                long double SNK;
                long double Factor1;
                long double Factor2;
                long double Factor3;
                long double Betta;
                for (int k = 0; k <= n; k++)
                {
                    Sat.CNK[n][k] = Clm[k][n];
                    Sat.SNK[n][k] = Slm[k][n];
                    int m,m1;
                    Factor3 =1;
                    for (m = n-k, m1= n+k; (m >=1) || (m1 >=1); m--,m1--)
                    {
                        if ((m >= 1) && (m1 >=1))
                        {
                            Factor3 *= (long double)m/(long double)m1;
                        }
                        else if (m >= 1)
                        {
                            Factor3 *= (long double)m;
                        }
                        else // (m1 >=1)
                        {
                            Factor3 /= (long double)m1;
                        }
                    }

                    // CNK = sqrt(2*(2*n+1)) * sqrt(((n-k)!/(n+k)!) * Clm
                    // SNK = sqrt(2*(2*n+1) * sqrt((n-k)!/(n+k)!)) * Slm
                    if (k == 0)
                        Betta = 1.0;
                    else
                        Betta = 2.0;
                     CNK = sqrt(Betta*(2*(long double)n+1) * Factor3) * ClmNN[k][n];
                     
                     SNK = sqrt(2*(Betta*(long double)n+1) * Factor3) * SlmNN[k][n];
#ifdef USE_MODEL_2
                     Sat.CNK[n][k] = CNK;
                     Sat.SNK[n][k] = SNK;
#endif
#ifdef USE_MODEL_3
                     Sat.CNK[n][k] = CNK;
                     Sat.SNK[n][k] = SNK;
#endif
                     if (iCounter_nk_lm_Numbers == 0 && n < 2)
                         continue;
                     nk_lm_Numbers[iCounter_nk_lm_Numbers][0] = n;
                     nk_lm_Numbers[iCounter_nk_lm_Numbers][1] = k;
                     C_S_nk[iCounter_nk_lm_Numbers][0] =Sat.CNK[n][k];
                     C_S_nk[iCounter_nk_lm_Numbers][1] =Sat.SNK[n][k];
                     //if (n==2 && k ==2)
                     //{
                     //    C_S_nk[iCounter_nk_lm_Numbers][0]+=1.3*1.39e-8;
                     //    Sat.CNK[n][k]+=1.3*1.39e-8;
                     //}
                     if (iCounter_nk_lm_Numbers++ >= TOTAL_COEF*TOTAL_COEF/2)
                     {
                         printf("\n something wrong with coef");
                         exit(555);
                     }
                }
                Sat.J[n] = (Clm[0][n]);
            }
#endif

            Sat.iLeg_longit = 0; // no longitude in calculation
            Sat.Lambda = -2;
            Sat.LegBody = EARTH;

            memcpy(Sat.MainCpu.diagonal, Sat.diagonal, sizeof(Sat.MainCpu.diagonal));
            memcpy(Sat.MainCpu._pt_nk, Sat._pt_nk, sizeof(Sat.MainCpu._pt_nk));
            memcpy(Sat.MainCpu._p_n_m_1, Sat._p_n_m_1, sizeof(Sat.MainCpu._p_n_m_1));
            memcpy(Sat.MainCpu._p_n_m_2, Sat._p_n_m_2, sizeof(Sat.MainCpu._p_n_m_2));
            memcpy(Sat.MainCpu._tpk_n_k, Sat._tpk_n_k, sizeof(Sat.MainCpu._tpk_n_k));
            memcpy(Sat.MainCpu._tp_nm1_k, Sat._tp_nm1_k, sizeof(Sat.MainCpu._tp_nm1_k));
            memcpy(Sat.MainCpu._tp_nm2_k, Sat._tp_nm2_k, sizeof(Sat.MainCpu._tp_nm2_k));
            memcpy(Sat.MainCpu.C_S_nk, C_S_nk, sizeof(Sat.MainCpu.C_S_nk));
            memcpy(Sat.MainCpu._p_n_k, Sat._p_n_k, sizeof(Sat.MainCpu._p_n_k));
            Sat.MainCpu._SQRT3 = Sat._SQRT3;
            if (CpuCore)
            {
                Sat.CalcSplit(EarthModelCoefs,CpuCore);
                Sat.StartThreads();
            }
            else
                Sat.i_proc = 0;

        }
        IF_XML_READ(UseSatData)
        {
            strcpy(UseSatData, pszQuo);
            if (strchr(UseSatData, '\"'))
                *strchr(UseSatData, '\"')=0;
        }
        IF_XML_READ(Targetlongitude) // dolgota
        {
            Targetlongitude = atof(pszQuo);
        }

        IF_XML_READ(Targetlatitude) // shirota
        {
            Targetlatitude = atof(pszQuo);
        }

    XML_SECTION_END
    ///////////////////////////////////////////////////////////////////////////////////////////////Engine
    XML_SECTION(Engine)
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
    
    XML_SECTION_END

            
    ///////////////////////////////////////////////////////////////////////////////////////////////////Optim
    XML_SECTION(Optim)
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
    XML_SECTION_END
    XML_END;
}