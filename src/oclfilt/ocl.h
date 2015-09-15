/* Include file for program oclfilt and function getOCLStationData()         */

/* Exit Statuses */
#define SUCCESSFUL 0
#define UNSPECIFIED_PROBLEM 1
#define ZERO_LENGTH_FIELD 2
#define SKIPPED 3
#define HELP_LISTING 4

/* Constants */
#define MAX_VARS 10
#define MAX_LEVELS 6000  /* 100 okay for stdlvl data, 6000 needed for obslvl */


/* Structures */
typedef struct OCLStation {

      /* actual data-file contents */
      long int stationNumber;
      long int bytesLeftInStation;
      long int bytesInStation;
      long int oclStationNumber;
      long int countryCode;
      long int cruiseNumber;
      long int year;
      long int month;
      long int day;
      long int numberOfLevels;
      long int stationType;
      long int numberOfVarCodes;
      long int varCode[MAX_VARS];
      long int errCodeForVarCode[MAX_VARS];
      long int bytesInCharPI;
      long int bytesInSecHdr;
      long int numberOfSecHdrEntries;
      long int secHdrCode[MAX_VARS];
      long int bytesInBioHdr;
      double time;
      double lat;
      double lon;
      double secHdrValue[MAX_VARS];
      double depthValue[MAX_LEVELS];
      double varValue[MAX_VARS][MAX_LEVELS];
      long int errCodeForDepthValue[MAX_LEVELS];
      long int errCodeForVarValue[MAX_VARS][MAX_LEVELS];

      /* added info for this station */
      double *bottomDepthPtr;  /* points to either secondary hdr, last profile
                                  depth, or bathy db value, depending on
                                  bottomDepthSource.
                                  This value (and bottomDepthSource) will
                                  reflect the best bottomDepth value available
                                  given input flags, as per these priorities:
                                   1.) if wantDatabaseBathy flag true, then
                                       return this value - the most dependable.
                                   2.) otherwise, if wantProfile flag true,
                                       then return the best value between sec
                                       header depth value and last profile
                                       depth value :
                                        a.) if bottomDepth value in secondary
                                            header exists and is not shallower
                                            than deepest legit profile depth
                                            (if there is a profile), then
                                            return sec hdr bottomDepth.
                                        b.) otherwise, return deepest legit
                                            profile depth
                                   3.) otherwise, return sec header depth value
                                       if it exists
                                   4.) otherwise, return NULL pointer.       */
                                            
      char bottomDepthSource;  /* 'h'=secondary hdr, 'p'=last profile depth,
                                  'd'=bathy database                         */
      double dbBathy;          /* bathy value for this lat/lon from database */
      int varListChecksOut;    /* flag, specifies whether stn includes all
                                  variables on varList                       */
      int badLatLon;           /* flag specifying that lat & lon values of zero
                                  were checked against wmo-square location -
                                  this flag is true if wmo-square is not along
                                  equator when lat=0, or if wmo-square is not
                                  along prime meridian when lon=0.           */
      int latlonInRange;       /* flag saying stn's latlon values are in the
                                  range specified by latlonRange[] array */
      int monthInRange;        /* flag saying stn's month value is in the
                                  range specified by monthRange[] array */
      int yearInRange;         /* flag saying stn's year value is in the
                                  range specified by yearRange[] array */
      int enoughProfileLevels; /* flag saying minimum number of profile levels
                                  for reporting this profile was okay */

}  OCLStationType;



/* Function Prototypes -
   (not all these functions are globally used, most only within one other
   function, but declaring them here keeps them out of the way and makes for
   more readable code in the modules, with neglibable performance loss)      */
int outputAllStationData(FILE *fp_out, long int i, OCLStationType *stnData);
int getOCLStationData( FILE *fp_in, long int stn, OCLStationType *stnData,
   int wantProfileFlag, int skipFlag, long int stnToSkipTo,
   int varListFlag, long int *varList, long int numVarsOnVarList,
   int minLevelsFlag, long int minLevels,
   int latlonRegionFlag, double *latlonRegion,
   int yearRangeFlag, long int *yearRange,
   int monthRangeFlag, long int *monthRange,
   int dbBathyFlag, FILE *fp_dbBathy, int zeroLatLonFlag, char *wmoSquare );
int parse_commandline( int argc, char **argv, FILE **fpIn, FILE **fpOut,
   int *botDepthFiltFlag, double *shallowerDLimit, double *deeperDLimit,
   int *varListFlag, long int *numVarsOnVarList, long int *varList,
   int *debugFlag, int *endStatsFlag, int *titlesFlag, int *queryFlag,
   int *databaseBathyFlag, char *databaseBathyFilename,
   int *numStnsFlag, long int *numStnsToOutput,
   int *skipFlag, long int *stnToSkipTo,
   int *zeroLatLonFlag, char *wmoSquare,
   int *minLevelsFlag, long int *minLevels,
   int *latlonRegionFlag, double *latlonRegion,
   int *yearRangeFlag, long int *yearRange,
   int *monthRangeFlag, long int *monthRange,
   int *includeErrorFlaggedData );
int getIntDigits(FILE *fp, int numDigits, long int *value);
int getVarlenIntField(FILE *fp, long int *value, long int *bytesLeftInStation);
int getVarlenFloatField(FILE *fp, double *value, long int *bytesLeftInStation);
int skipToNextStation(FILE *fp, long int bytesLeftInStation);
char *varCodeLabel(long int oneVarCode);
char *varCodeUnits(long int oneVarCode);
int checkVarsInclAndNoErrors(long int *varRequestedList, long int *varCodeList,
   long int *errCodeList, long int numVarsRequested, long int numVarCodes);
int zeroLatLonOkay( char *wmoSquare, char *latlon );
double nan();
