#include "stdafx.h"
#include "afxinet.h"
#include "afxsock.h"
#include <string.h>

#include <stdio.h>

#include "procXML.h"
#include "tra.h"

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