#include "stdafx.h"
#include "procXML.h"


///////////////////////////////////////////////////////////////////////////
// quick XML parser
///

void init_proc_XML(void) 
{
    szSection[0] =  0;
    szGroup[0] =  0;
}

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

char *CallElementPars(char *szString, char *XML_Params)
{
    char szFullComapre[1024] = {"<"};
    char szFullComapre2[1024] = {"</"};
    strcat(szFullComapre, XML_Params);
    strcat(szFullComapre2, XML_Params);
    strcat(szFullComapre, ">");
    strcat(szFullComapre2, ">");
    if (strstr(szString, szFullComapre) != NULL)   
    {
        if (strstr(szString, szFullComapre2) != NULL)
        {
            char *szPrt = strstr(szString, szFullComapre);
            szPrt += strlen(szFullComapre);
            return szPrt;
        }
        return NULL;
    }
    else
        return NULL;
}

char *CallGroupPars(char *szString, char *XML_Params)
{
    char szFullComapre[1024] = {"<"};
    strcat(szFullComapre, XML_Params);
    strcat(szFullComapre, ">");
    if (strstr(szString, szFullComapre) != NULL)   
    {
        char *szPrt = strstr(szString, szFullComapre);
        return XML_Params;
    }
    else
        return NULL;
}