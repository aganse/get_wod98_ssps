/* getOCLStationData() -
 *             Reads one station of data from an OCL-formatted data file into
 *             a data structure for use in the calling routine.  If filters
 *             are specified in args, getOCLStationData will stop reading and
 *             skip to the next station in the input file (and return to the
 *             calling function) as soon as a filter applies, in order to
 *             reduce computation time (mostly I/O time).
 *             Note that getOCLStationData skips right over the station's
 *             character data, PI data, and bio/taxo data (since I'm not
 *             interested in those), but the framework set up in this code
 *             makes it easy enough to add those if desired.  Comments in the
 *             code label the places to change.  (seach for PI, bio, taxo...)
 *
 * other required sources/files: ocl.h
 *
 * language:   ANSI C
 *
 * author:     Andy Ganse, Applied Physics Laboratory, University of Washington
 *             (aganse@apl.washington.edu)
 *
 * last update:2/23/00
 *
 * notes:
 *             OCL is the internal data format of the NODC Ocean Climate Lab.
 *             Assumes input files have \r's stripped (ie., in UNIX not DOS
 *             text format)
 *
 *             The data files consist of "stations" (unique instances of
 *             measurement-location and -time), generally read in a loop.
 *             The files use a proprietary format of byte conditions ("this
 *             byte tells how many bytes are in the following byte-field",
 *             "if this byte is a 3 the following 2 bytes mean such-and-such")
 *             They rely strictly on sequential access - one misread byte
 *             can potentially ruin the reading of an ensuing 20MB of data.
 *
 *             getOCLStationData() takes in the data stream-wise and includes
 *             logic to both skip over whole unwanted stations and skip over
 *             unwanted parts of station data in order to decrease computation
 *             time.  This becomes very significant when getOCLStationData() is
 *             used in a loop to read a 10 or 20 MB data file (typical).
 *             The extra logic and count values in the function args are for
 *             this purpose.
 *
 *             The OCL format used was specified in the format.txt file that
 *             accompanies the WOD98 CD set.  However, note that there is an
 *             error at the end of that spec - in the profile section the spec
 *             says (in so many words):
 *               loop over levels:
 *                  get depth value.
 *                  get variable value.
 *             when in fact it should say (in so many words):
 *               loop over levels:
 *                  get depth value.
 *                  loop over variable columns:
 *                     get variable value.
 *             The latter is of course what this program does (you don't get
 *             much data otherwise).
 *
 *             The function prototypes for getOCLStationData and its
 *             subfunctions are found in ocl.h, as well as the data structure
 *             OCLStationType, which consists of all the data values taken
 *             from one "station" of the data file.  Some symbolic constants
 *             are in there too.
 *
 * Usage:      int getOCLStationData( FILE *fp_in, long int stn, 
 *               OCLStationType *stnData, int wantProfileFlag,
 *               int skipFlag, long int stnToSkipTo,
 *               int varListFlag, long int *varList, long int numVarsOnVarList,
 *               int minLevelsFlag, long int minLevels,
 *               int latlonRegionFlag, double *latlonRegion,
 *               int yearRangeFlag, long int *yearRange,
 *               int monthRangeFlag, long int *monthRange,
 *               int dbBathyFlag, FILE *fp_dbBathy,
 *               int zeroLatLonFlag, char *wmoSquare );
 *
 *             Explanation of function args:
 *
 *              FILE *fp_in - OCL-formatted input file we read data from.
 *                 File opened & closed in the calling function.
 *
 *              long int stn - current station # (from loop in calling funct)
 *
 *              OCLStationType *stnData - pointer to structure which will hold
 *                 all the station data we read; note that since this is a
 *                 pointer to struct we use stnData->element in the code below,
 *                 rathen than stnData.element.  Referencing from a pointer
 *                 allows us to pass all our data back to the calling routine.
 *
 *              int wantProfileFlag - true (1) = yes, we want to read the
 *                                       profile data for this station.
 *                                    false (0) = no, skip the profile data for
 *                                       this station - only want its header
 *                                       header data.  Saves lots of I/O time
 *                                       when filtering or doing stats
 *                                       (like report station bytecounts
 *                                       and which vars they have)
 *
 *              int skipFlag - true (1) = yes, we're skipping the station
 *                                specified by stnToSkipTo, so if we're
 *                                not there yet, only read enough station
 *                                data so we know how many bytes to skip
 *                             false (0) = no, we're not skipping any stations
 *                                so don't worry about stnToSkipTo value
 *
 *              long int stnToSkipTo - the station number we want to skip to,
 *                 if skipFlag true (count begins at station 0)
 *
 *              int varListFlag - true (1) = yes, there's a variable list that
 *                                   is specified in *varList, and we
 *                                   want to skip stations that do not
 *                                   contain those variables.  But we
 *                                   had to read a little bit of the
 *                                   station to check that.
 *                                false (0) = no, we're not skipping any
 *                                   stations based on which variables they
 *                                   have, so don't bother with *varList
 *
 *              long int *varList - a 1-D array with the numerical var codes of
 *                 the variables that we want to make sure are in our data
 *                 stations.  See NODC's format.txt from WOD98 for these var
 *                 codes.  For example if we want to skip stations that don't
 *                 have at least temp and pressure vars, then varListFlag
 *                 would be 1 and varList[] here would be {2,25}.
 *                 (1=temp & 25=pres).
 *
 *              long int numVarsOnVarList - # of variables on varList above.
 *                 so for that little example above, numVarsOnVarList=2.
 *
 *              int minLevelsFlag - true (1) = yes, we want to filter out
 *                                     stations having fewer profile levels
 *                                     than minLevels
 *                                  false (0) = no, don't filter stations by
 *                                     how many profile points (levels) they
 *                                     have, and ignore the minLevels arg.
 *              long int minLevels - minimum number of profile levels station
 *                                   must have to pass the filter
 *              int latlonRegionFlag - true (1) = yes, we want to filter out
 *                                        stations whose lat/lon values are
 *                                        not within the range specified
 *                                        in *latlonRegion.
 *                                     false (0) = no, don't filter stations
 *                                        by their lat/lon region, and ignore
 *                                        the *latlonRegion argument.
 *                          
 *              double *latlonRegion - four-element array holding lat/lon
 *                                  region in the following order:
 *                                  westbound,eastbound,southbound,northbound.
 *                                  lat/lon values are in decimal degrees, with
 *                                  lons as -90 to +90 degrees, and lats may be
 *                                  either -180 to +180 or 0 to 360 degrees.
 *
 *              int yearRangeFlag - true (1) = yes, we want to filter out
 *                                     stations whose recorded year is not
 *                                     within the range specified by yearRange.
 *                                  false (0) = no, don't filter stations by
 *                                     their year value, and ignore the
 *                                     *yearRange argument.
 *
 *              long int *yearRange - two-element array holding minyear,maxyear
 *                                  to filter between, inclusive of both.
 *
 *              int monthRangeFlag - true (1) = yes, we want to filter out
 *                                    stations whose recorded month is not
 *                                    within the range specified by monthRange.
 *                                  false (0) = no, don't filter stations by
 *                                    their month value, and ignore the
 *                                    *monthRange argument.
 *
 *              long int *monthRange - two-element array holding
 *                                  minmonth,maxmonth to filter between,
 *                                  inclusive of both.
 *
 *              int dbBathyFlag - true (1) = yes, we want to look up the
 *                                   bathy value for this station's
 *                                   location from a database (more data
 *                                   integrity but much more comp time)
 *                                false (0) = no, don't look up the bathy
 *                                   value for this station's location
 *                                   in a database.  (Ignore fp_dbBathy.)
 *
 *              FILE *fp_dbBathy - pointer to database-bathy file (if used).
 *                 File is a premade ASCII lat-lon-depth file made specifically
 *                 for use with this OCL input-file, by the shell script 
 *                 bathyForThisOCLfile.  File opened/closed in calling function
 *
 *              int zeroLatLonFlag - true (1) = yes, we want to filter out
 *                                      stations with an invalid value of zero
 *                                      for latitude or longitude, as loosely
 *                                      determined by looking at the wmoSquare
 *                                      value and seeing if this station's
 *                                      square lies along equator or prime
 *                                      meridian.
 *                                  false (0) = no, don't filter out stations
 *                                      by their zeroed lats & lons, and ignore
 *                                      wmoSquare arg.
 *
 *              char *wmoSquare - four char string specifying station's WMO
 *                 Square (10-deg lat-lon square, as defined by UN WMO).
 *                 Generally the calling function would gleen this from the
 *                 name of the OCL file.
 *
 *
 * Example Implementation:
 *   Most basic usage.  No filtering, returns all station data.
 *   Calling function must include "ocl.h" and open the input file
 *   (if not using stdin).
 *
 *
 * #include <stdio.h>
 * #include "ocl.h"
 * int main() {
 *   OCLStationType stnData;
 *   int stnNum;
 *   for( stnNum=0; !feof(stdin) && stnNum<10; stnNum++ ) {
 *     if( getOCLStationData( stdin, stnNum, &stnData, 1, 0, 0, 0, 0, 0, 0, 0,
 *                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ) != SUCCESSFUL ) {
 *       printf("error at station %ld\n", stnNum);
 *       exit(1);
 *     }
 *     else
 *       printf("Stn %ld Lat=%f Lon=%f\n", stnNum, stnData.lat, stnData.lon);
 *   }
 *   return SUCCESSFUL;
 * }
 *
 *
 * Compile as :
 *   gcc -O -o exampleprog examplecode.c getOCLStation.c -lm
 * Then run as (note stripping of DOS linefeeds with tr) :
 *   gunzip -c nbds1106.gz | tr -d '\r' | exampleprog | more
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>  /* <-- remember need to specify "-lm" on compile cmd line */
#include "ocl.h"


int getOCLStationData( FILE *fp_in, long int stn, OCLStationType *stnData,
   int wantProfileFlag, int skipFlag, long int stnToSkipTo,
   int varListFlag, long int *varList, long int numVarsOnVarList,
   int minLevelsFlag, long int minLevels,
   int latlonRegionFlag, double *latlonRegion,
   int yearRangeFlag, long int *yearRange,
   int monthRangeFlag, long int *monthRange,
   int dbBathyFlag, FILE *fp_dbBathy, int zeroLatLonFlag, char *wmoSquare ) {

   long int j, k, ld_dummy;
   double lf_dummy;
   int status=SUCCESSFUL, reallyWantProfile, assignLastProfileDepth=0;
   static double databaseBathyValue;
   /* array of standard-level depths : */
   double stdLevelDepth[] = { 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250,
      300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500,
      1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000,
      7500, 8000, 8500, 9000 };


   /* first two fields of station tell how many bytes in station,
      so for now must flag bytesLeftInStation as unusable */
   stnData->bytesLeftInStation=-1;


   /* The byte-reading sequence below can be understood by comparing it
      line-for-line with the NODC OCL format.txt (included with WOD98),
      which defines the byte sequence in OCL formatted files */


   /* Read station header : */

   if( getVarlenIntField( fp_in, &(stnData->bytesInStation),
      &(stnData->bytesLeftInStation) ) == ZERO_LENGTH_FIELD )
      stnData->bytesInStation=0;

   if( getVarlenIntField( fp_in, &(stnData->oclStationNumber),
      &(stnData->bytesLeftInStation) ) == ZERO_LENGTH_FIELD )
      stnData->oclStationNumber=-1;

   if( skipFlag && stn<stnToSkipTo ) {
     skipToNextStation( fp_in, stnData->bytesLeftInStation);
     /* if using bathy file, skip past line in there, too */
     if( dbBathyFlag ) {
       fscanf( fp_dbBathy, "%lf %lf %ld %lf\n", &lf_dummy, &lf_dummy,
            &ld_dummy, &lf_dummy );
     }
     return SKIPPED;
   }

   if( getIntDigits( fp_in, 2, &(stnData->countryCode) ) == ZERO_LENGTH_FIELD)
     stnData->countryCode=-1;
   stnData->bytesLeftInStation-=2;

   if( getVarlenIntField( fp_in, &(stnData->cruiseNumber),
      &(stnData->bytesLeftInStation) ) == ZERO_LENGTH_FIELD)
     stnData->cruiseNumber=-1;

   getIntDigits( fp_in, 4, &(stnData->year) );
   stnData->bytesLeftInStation-=4;

   getIntDigits( fp_in, 2, &(stnData->month) );
   stnData->bytesLeftInStation-=2;

   getIntDigits( fp_in, 2, &(stnData->day) );
   stnData->bytesLeftInStation-=2;

   getVarlenFloatField( fp_in, &(stnData->time),
      &(stnData->bytesLeftInStation) );

   getVarlenFloatField(fp_in, &(stnData->lat), &(stnData->bytesLeftInStation));

   getVarlenFloatField(fp_in, &(stnData->lon), &(stnData->bytesLeftInStation));

   if( getVarlenIntField( fp_in, &(stnData->numberOfLevels),
      &(stnData->bytesLeftInStation) ) == ZERO_LENGTH_FIELD)
     stnData->numberOfLevels=0;

   getIntDigits( fp_in, 1, &(stnData->stationType) );
   stnData->bytesLeftInStation-=1;

   getIntDigits( fp_in, 2, &(stnData->numberOfVarCodes) );
   stnData->bytesLeftInStation-=2;

   for(j=0; j<stnData->numberOfVarCodes; j++) {
     /* varCodes are the data column descriptors, eg 1=temp, 2=sal, etc.
        as defined in Table 4 in the NODC OCL readmev1 doc */
     getVarlenIntField( fp_in, &(stnData->varCode[j]),
        &(stnData->bytesLeftInStation) );
     getIntDigits( fp_in, 1, &(stnData->errCodeForVarCode[j]) );
     stnData->bytesLeftInStation-=1;
   }



   /* Read character data & principal investigator data : */

   /* (note this section has some char data - would need to write a
      getCharField function if I decided to actually take data from here) */

   status = getVarlenIntField( fp_in, &(stnData->bytesInCharPI),
      &(stnData->bytesLeftInStation) );

   if( status != ZERO_LENGTH_FIELD ) {
      /* skipping rest of this section for now */
      for(j=0; j<stnData->bytesInCharPI; j++) {
         getIntDigits( fp_in, 1, &ld_dummy );
         stnData->bytesLeftInStation-=1;
      }
   }
   else stnData->bytesInCharPI=0;



   /* Read secondary header : */

   status = getVarlenIntField( fp_in, &(stnData->bytesInSecHdr),
      &(stnData->bytesLeftInStation) );

   if( status != ZERO_LENGTH_FIELD ) {

      getVarlenIntField( fp_in, &(stnData->numberOfSecHdrEntries),
         &(stnData->bytesLeftInStation) );

      for(j=0; j<stnData->numberOfSecHdrEntries; j++) {

         /* secHdr Codes are the secondary header data descriptors,
            eg 10=bottom depth, 18=sea state, etc.
            as defined in Table 6 in the NODC OCL readmev1 doc */

         getVarlenIntField( fp_in, &(stnData->secHdrCode[j]),
            &(stnData->bytesLeftInStation) );

         getVarlenFloatField( fp_in, &(stnData->secHdrValue[j]),
            &(stnData->bytesLeftInStation) );
      }
   }
   else { /* ie, if it IS a zero-length field... */
      stnData->bytesInSecHdr=0;
      stnData->numberOfSecHdrEntries=0;
   }



   /* Interlude before reading profile to assign and calculate bottomDepth
    * values, and to decide whether or not to read rest of station for profile
    */


   /* Get bottomDepth values now, in case we don't read rest of station : */

      /* initializing bottomDepth elements */
         stnData->bottomDepthPtr = NULL;
         stnData->bottomDepthSource = '-';

      /* if there's a bottomDepth in the secondary hdr, assign it */
      for(j=0; j<stnData->numberOfSecHdrEntries; j++)
         /* search for header code 10, which is depth */
         if( stnData->secHdrCode[j]==10 ) {
            stnData->bottomDepthPtr = &(stnData->secHdrValue[j]);
            stnData->bottomDepthSource='h';
         }

      /* if we're using the dbBathy file, read/increment the file pointer, and
         reassign bottomDepth from the dbBathy file if we want it.  (otherwise
         the bottomDepth from header will remain) */
      if( dbBathyFlag ) {
         /* Get the bathy value for current lat/lon from bathy file.  Note
            there should one line in bathy file for each station in the
            input file - for my use, the file is created automatically in the
            shell script processOcean, using outputAllLatsLons & grdtrack.  */
         fscanf( fp_dbBathy, "%lf %lf %ld %lf\n", &lf_dummy, &lf_dummy,
            &ld_dummy, &databaseBathyValue );
         databaseBathyValue *= -1;        /* (converting neg depths to pos) */

         /* set bottomDepthPtr as databaseBathyValue if:
            we're in domain of dbBathy, and:
               no hdrDepth available, or
               hdrDepth available, but diff between hdrDepth value and dbBathy
                  value is to big. */
         if( stnData->lat<=72. && stnData->lat>=-72. ) {
            if( stnData->bottomDepthSource!='h' ||
                (-80 > *(stnData->bottomDepthPtr)-databaseBathyValue) ||
                (*(stnData->bottomDepthPtr)-databaseBathyValue > 80 )    ) {
               stnData->bottomDepthPtr = &databaseBathyValue;
               stnData->bottomDepthSource = 'd';
            }
         }
      }



   /* Deciding whether to read rest of station */

      /* Reading the rest of the station (including profile) is conditional
         upon whether we need that data - getOCLStationData arg list might
         specify that we aren't interested in the profile, or perhaps we've
         found that the list of required variables isn't covered, so we don't
         want this station. (Skipping saves a bunch of comp & I/O time...) */
   
      /* Won't need rest of station if varList isn't covered and error-free */
      if( varListFlag ) {
         /* checking if station vars include those on varList, and have no
            errors covering the whole variable column */
         stnData->varListChecksOut = checkVarsInclAndNoErrors( varList, 
            stnData->varCode, stnData->errCodeForVarCode, numVarsOnVarList,
            stnData->numberOfVarCodes);
      }
      else stnData->varListChecksOut = 1;
   
      /* Won't need rest of station if we specified we don't want the profile.
         But even if we don't want profile, if bottomDepth still has not been
         set from either the Secondary Hdr or from the bathy database, we must
         read profile in order to assign last profile depth to bottomDepth */
      if(!wantProfileFlag && stnData->bottomDepthPtr!=NULL)
         reallyWantProfile=0;
      else reallyWantProfile=1;

      /* Won't need rest of station if we're filtering out Lat/Lon values and
         the current ones aren't in latlonRegion */
      if( latlonRegionFlag && 
          (stnData->lon < latlonRegion[0] ||
           stnData->lon > latlonRegion[1] ||
           stnData->lat < latlonRegion[2] ||
           stnData->lat > latlonRegion[3] )  )
         stnData->latlonInRange=0;
      else stnData->latlonInRange=1;

      /* Won't need rest of station if we're filtering out year values and the
         current one isn't in yearRange */
      if( yearRangeFlag &&
          (stnData->year < yearRange[0] || stnData->year > yearRange[1]) )
         stnData->yearInRange=0;
      else stnData->yearInRange=1;

      /* Won't need rest of station if we're filtering out month values and the
         current one isn't in monthRange */
      if( monthRangeFlag &&
          (stnData->month < monthRange[0] || stnData->month > monthRange[1]) )
         stnData->monthInRange=0;
      else stnData->monthInRange=1;

      /* Won't need rest of station if we're filtering out stations based on
         minimum number of profile levels, and current one doesn't have that
         many */
      if( minLevelsFlag && stnData->numberOfLevels<minLevels )
         stnData->enoughProfileLevels=0;
      else stnData->enoughProfileLevels=1;

      /* Won't need rest of station if we're filtering out bad Lat/Lon values
         and we find zero-values for lat or lon when we're not on the equator
         or prime meridian (respectively).  We check this merely by looking at
         the wmo-square value that was specified on the command line. */
      stnData->badLatLon=0;
      if( zeroLatLonFlag ) {
         /* (remember can't rely on doubles being exactly equal...) */
         if( stnData->lat<0.0000001 && stnData->lat>-0.0000001 ) {  /*if zero*/
            stnData->badLatLon += !zeroLatLonOkay( wmoSquare, "lat" );
         }
         if( stnData->lon<0.0000001 && stnData->lon>-0.0000001 ) {  /*if zero*/
            stnData->badLatLon += !zeroLatLonOkay( wmoSquare, "lon" );
         }
      }




   /* Go ahead and read rest of station if our criteria is met */
   if( reallyWantProfile && stnData->varListChecksOut && !stnData->badLatLon &&
      stnData->latlonInRange && stnData->monthInRange && stnData->yearInRange
      && stnData->enoughProfileLevels ) {

      /* Read biological header : */
      status = getVarlenIntField( fp_in, &(stnData->bytesInBioHdr),
         &(stnData->bytesLeftInStation) );

      if ( status != ZERO_LENGTH_FIELD ) {
         /* skipping rest of this section for now */
         for(j=0; j<stnData->bytesInBioHdr; j++) {
            getIntDigits( fp_in, 1, &ld_dummy);
            stnData->bytesLeftInStation-=1;
         }
      }
      else stnData->bytesInBioHdr=0;



      /* Read taxonomic and biomass data : */
      /* (skipped already by virtue of the fact that all the bio data is
          skipped from within the above bio-hdr section - but here's a place-
          holder/reminder for future expansion) */



      /* Read profile data : */

      /* make sure profile will fit in memory set aside for it */
      /* (this should really be variably-allocated... I'll get to it someday)*/
      if(stnData->numberOfLevels>MAX_LEVELS) {
         fprintf(stderr, "%% getOCLStationData: warning: stn %ld numLevels=%ld"
            " > MAX_LEVELS=%d.\n", stn, stnData->numberOfLevels,
            MAX_LEVELS);
         fprintf(stderr, "%%   botDepth might be wrong and not all the profile"
            " levels will be output.\n");
      }

      /* loop over profile levels */
      for(j=0; j<stnData->numberOfLevels && j<MAX_LEVELS; j++) {

         /* depth values */
         if( stnData->stationType==0 /*observed level*/ ) {
            status = getVarlenFloatField( fp_in, &(stnData->depthValue[j]),
               &(stnData->bytesLeftInStation) );
            if( status != ZERO_LENGTH_FIELD ) {
               getIntDigits( fp_in, 1, &(stnData->errCodeForDepthValue[j]) );
               stnData->bytesLeftInStation-=1;
            }
         }
         else stnData->depthValue[j]=stdLevelDepth[j];

         /* the values for each varCode (temp, sal, etc) */
         for(k=0; k<stnData->numberOfVarCodes; k++) {
            status = getVarlenFloatField( fp_in, &(stnData->varValue[k][j]),
               &(stnData->bytesLeftInStation) );
            if( status != ZERO_LENGTH_FIELD ) {
               getIntDigits( fp_in, 1, &(stnData->errCodeForVarValue[k][j]) );
               stnData->bytesLeftInStation-=1;
            }
         }
      }




      /* recheck bottomDepth value against lowest profile depth: if lowest
         profile depth is deeper than bottomDepth (=hdrdepth or =dbdepth)
         we rechoose which value we use for bottomDepth, or use lowest profile
         depth for bottomDepth if they're both bad */
      if( stnData->numberOfLevels>0 ) {
         if( stnData->bottomDepthPtr==NULL) assignLastProfileDepth=1;/*if none*/
         else if( stnData->bottomDepthSource=='h' &&     /* if hdrdepth bad */
                  ( *(stnData->bottomDepthPtr) < 
                  stnData->depthValue[stnData->numberOfLevels-1] ) ) {
            if( dbBathyFlag ) {  /* maybe substitute db value if using db */
               if( databaseBathyValue < 
                  stnData->depthValue[stnData->numberOfLevels-1] ) {
                  assignLastProfileDepth=1;
               }
               else {
                  stnData->bottomDepthPtr = &databaseBathyValue;
                  stnData->bottomDepthSource = 'd';
               }
            }
            else assignLastProfileDepth=1;  /*no db value so just use profile*/
         }
         else if( stnData->bottomDepthSource=='d' &&      /* if dbdepth bad */
                  ( *(stnData->bottomDepthPtr) < 
                  stnData->depthValue[stnData->numberOfLevels-1] ) ) {
            assignLastProfileDepth=1;
         }
         else assignLastProfileDepth=0;  /* we have a good bottomDepth value */

         if( assignLastProfileDepth ) {
            stnData->bottomDepthPtr =
               &(stnData->depthValue[stnData->numberOfLevels-1]);
            stnData->bottomDepthSource='p';
         }
      }


   } /* done reading rest of profile */


   /* even at end of station, still need to do this to get past any remaining
      white space till end of line */
   skipToNextStation( fp_in, stnData->bytesLeftInStation);

   return SUCCESSFUL;

} /* end of getOCLStationData */







/* "Get integer digits" - number of digits is specified, take from fp and place
   in integer */
int getIntDigits(FILE *fp, int numDigits, long int *value) {
   char nextch=0, valueStr[10]="", *p;
   int i, status;
   double badValue=nan();  /* assigning NaN */

   p=valueStr;

   for(i=0; i<numDigits; i++) {
      if( (nextch=fgetc(fp))==EOF ) {  /* read & check for eof at same time */
         fprintf(stderr,"%%oclfilt: unexpected EOF - empty or truncated input file?\n");
         exit(1);
      }
      if(nextch=='\n' || nextch=='\r') i--; /*don't inc j for newlines or CRs*/
      else *(p++)=nextch;  /* tack chars onto valueStr, keeping place with p */
   }

   /*check if at least rightmost char is a digit (the others may be spaces)*/
   if(numDigits>0 && isdigit((int)nextch) ) {
      sscanf(valueStr,"%ld", value);
      status=SUCCESSFUL;
   }
   else if(numDigits==0 || (numDigits==1 && nextch=='-') ) {
      *value=(long int)badValue;
      status=ZERO_LENGTH_FIELD;
   }
   else {
      status=UNSPECIFIED_PROBLEM;
   }

  return status;

}






/* "Get variable-length integer field" */
int getVarlenIntField(FILE *fp, long int *value, long int *bytesLeftInStation) {
   int status;
   long int bytesInNextField;
   double badValue=nan();  /* assigning NaN */

   status=getIntDigits(fp,1, &bytesInNextField);
   *bytesLeftInStation-=1;

   /* note 'value' is set inside the if statement here */
   if( getIntDigits(fp, bytesInNextField, value)==ZERO_LENGTH_FIELD ) {
      *value=(long int)badValue;
      status=ZERO_LENGTH_FIELD;
   }
   else {
     /* if bytesLeftInStation not set yet (ie beginning of station), set it */
     if(*bytesLeftInStation<0)
        *bytesLeftInStation=*value-(long int)bytesInNextField-1;
     /* otherwise decrement the actual number of bytesLeftInStation */
     else *bytesLeftInStation-=(long int)bytesInNextField;
     status=SUCCESSFUL;
   }

   return status;
}






/* "Get variable-length floating-point field" (data type is actually double) */
int getVarlenFloatField(FILE *fp, double *value, long int *bytesLeftInStation) {
   int status;
   long int sigDigits, totalDigits, precision, intvalue;
   double badValue=nan();  /* assigning NaN */


   status = getIntDigits(fp,1, &sigDigits);
   *bytesLeftInStation-=1;

   if (status==SUCCESSFUL) {
      status = getIntDigits(fp,1, &totalDigits);
      *bytesLeftInStation-=1;

      status = getIntDigits(fp,1, &precision);
      *bytesLeftInStation-=1;

      status = getIntDigits(fp, totalDigits, &intvalue);
      *bytesLeftInStation-=totalDigits;

      *value=(double)intvalue/pow(10.,(double)precision);
      status=SUCCESSFUL;
   }

   else if (status==ZERO_LENGTH_FIELD) {
      *value=badValue;
      status=ZERO_LENGTH_FIELD;
   }

   return status;
}






/* "Skip to next station" - skip over the remaining bytes left in the station
   (which are kept track of in the conversion process) to get to next station*/
int skipToNextStation(FILE *fp, long int bytesLeftInStation) {
   long int i;
   int nextch=0, status=SUCCESSFUL;
   char dummy[80];

   /* if we're not done with station bytes yet, skip over those */
   for(i=0; i<(bytesLeftInStation); i++) {
      nextch=fgetc(fp);
      if(nextch=='\n' || nextch=='\r') i--; /*don't inc j for newlines or CRs*/
   }

   /* skip past blank space and newline/CR to new station (or eof) */
   fgets(dummy,80,fp);

   /* sorry, dumb hack to expose eof character so next loop stops on it */
   nextch=fgetc(fp); if (nextch!=EOF) ungetc(nextch, fp);

   return status;
}






/* "Variable code label" - returns variable name given variable code number */
char *varCodeLabel(long int oneVarCode) {
   char *label[26];

   if(oneVarCode!=1 && oneVarCode!=2 && oneVarCode!=3 && oneVarCode!=4 &&
      oneVarCode!=6 && oneVarCode!=7 && oneVarCode!=8 && oneVarCode!=9 &&
      oneVarCode!=11 && oneVarCode!=17 && oneVarCode!=25) {

      fprintf(stderr,"%% oclfilt: error in varCodeLabel: invalid varCode number\n");
      exit(1);
   }

   label[1]=  "Temp";
   label[2]=  "Sal";
   label[3]=  "Oxy";
   label[4]=  "Phos";
   label[6]=  "Silic";
   label[7]=  "Nitri";
   label[8]=  "Nitra";
   label[9]=  "pH";
   label[11]= "Chlor";
   label[17]= "Alka";
   label[25]= "Pres";

   return label[oneVarCode];
}




/* "Variable code units" - returns units name given variable code number */
char *varCodeUnits(long int oneVarCode) {
   char *units[26];

   if(oneVarCode!=1 && oneVarCode!=2 && oneVarCode!=3 && oneVarCode!=4 &&
      oneVarCode!=6 && oneVarCode!=7 && oneVarCode!=8 && oneVarCode!=9 &&
      oneVarCode!=11 && oneVarCode!=17 && oneVarCode!=25) {

      fprintf(stderr,"%% oclfilt: error in varCodeUnits: invalid varCode number\n");
     exit(1);
   }

   units[1]=  "deg C";
   units[2]=  "ppt";
   units[3]=  "ml/l";
   units[4]=  "micromolar";
   units[6]=  "micromolar";
   units[7]=  "micromolar";
   units[8]=  "micromolar";
   units[9]=  "unitless";
   units[11]= "ug/l";
   units[17]= "meq/l";
   units[25]= "dbars";

   return units[oneVarCode];
}




/* Checking if station vars include those on varList, and have no errors
   relating to the whole variable column */
int checkVarsInclAndNoErrors(long int *varRequestedList, long int *varCodeList,
   long int *errCodeList, long int numVarsRequested, long int numVarCodes) {

   long int j, k, numMatch=0, applicableErrorFlags=0;

   /* See if elements of varRequestedList are a subset of elements of
      varCodeList (incrementing numMatch for each element match).
      Also keep track if one of the errCodes for a variable on requested var
      list shows an errorFlag */
   for( j=0; j<numVarsRequested; j++ ) {
      for( k=0; k<numVarCodes; k++ ) {
         if( varCodeList[k]==varRequestedList[j] ) {
            numMatch++;
            if( errCodeList[k]>0 ) applicableErrorFlags++;
         }
      }
   }

   if( numMatch>=numVarsRequested && !applicableErrorFlags ) return 1;
   else return 0;

}




int zeroLatLonOkay( char *wmoSquare, char *latlon ) {
   if( !strcmp(latlon,"lat") && wmoSquare[1]=='0' )
      return 1;
   else if( !strcmp(latlon,"lon") && wmoSquare[2]=='0' && wmoSquare[3]=='0' )
      return 1;
   else
      return 0;

}




/* little function to assign NaN */
double nan() {
   double x=0;
   return sqrt(-1/x);
}

