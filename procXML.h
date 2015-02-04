///////////////////////////////////////////////////////////////////////////
// quick XML parser
///
#define XML_READ(XML_PARAM) if (CallXMLPars(szString, #XML_PARAM)) XML_PARAM = atof(pszQuo);
#define IF_XML_READ(XML_PARAM) if (CallXMLPars(szString, #XML_PARAM))
#define XML_BEGIN char *pszQuo;
#define XML_END ;
#define XML_SECTION(XML_SEC_NAME) if (strcmp(szSection, #XML_SEC_NAME)==0){pszQuo = strstr(szString, "value=\"");if (pszQuo != NULL){pszQuo += sizeof("value=\"") -1;
#define XML_SECTION_END }}
#define XML_SECTION_GROUP_SEPARATOR } else {

#define XML_GROUP(XML_ELEMENT) if (pszQuo=CallGroupPars(szString,#XML_ELEMENT)) strcpy(szGroup, pszQuo); if (strcmp(szGroup, #XML_ELEMENT) ==0) {
#define XML_GROUP_END }

#define IF_XML_ELEMENT(XML_ELEMENT) if (pszQuo=CallElementPars(szString,#XML_ELEMENT))

void init_proc_XML(void);
int CallXMLPars(char *szString, char *XML_Params);
char *CallElementPars(char *szString, char *XML_Params);
char *CallGroupPars(char *szString, char *XML_Params);
void ParamSun(char *szString);

#ifdef USE_GLOBAL
#define GLOBAL_VARIABLE
#else
#define GLOBAL_VARIABLE extern
#endif

GLOBAL_VARIABLE char szSection[1024];
GLOBAL_VARIABLE char szGroup[1024];