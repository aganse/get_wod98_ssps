/* oclfilt.c - version 1.85
 *             Converts data from an ocl-formatted data file to ascii columns
 *             but also can considerably filter that data as per cmdline
 *             params, in ways try to that minimize computation/search times.
 *
 *             ("ocl" is a data format created by the NODC Ocean Climate Lab)
 *             Assumes input files have \r's stripped (ie., UNIX not DOS text
 *             format)
 * 
 * required sources/libs: getOCLStationData.c, ocl.h, Makefile;
 *
 * required input files for use: NODC/OCL-formatted data as input files (I'm
 *                              using files from NODC/OCL WOD98).
 *
 *                              If using "bathy database" option, requires a
 *                              3-column ASCII file of lat lon depth (depth
 *                              is negative in this file).  Lat-lons must
 *                              match up with lat-lons in OCL input data file.
 *                              I create mine like:
 *                                 set badPoint = ( 70, 30 )  # land location
 *                                 gunzip -c ncts1206.gz | tr -d '\r' | \
 *                                 outputAllLatsLons 1206 $badPoint | \
 *                                 grdtrack -Gtopo.nind.grd >! my.bathyfile,
 *                              which extracts subsets of the Sandwell
 *                              bathymetry database (see
 *                              http://topex.ucsd.edu/marine_topo/mar_topo.html)
 *                              using the GMT3.3 geographical pacakge
 *                              (see http://imina.soest.hawaii.edu/gmt)
 *                              for each station from an OCL file.
 *
 * language:   ANSI C
 *
 * author:     Andy Ganse, Applied Physics Laboratory, University of Washington
 *             (aganse@apl.washington.edu)
 *             As of this writing, this oclfilt-1.85 package is available at:
 *             http://staff.washington.edu/aganse/src
 *             Please email me if you use this code so I can alert you to 
 *             updates and have an estimate for how many folks use it.
 *
 * notes:      The OCL format used was specified in the format.txt file that
 *             accompanies the WOD98 CD set.  However, note that there is an
 *             error at the end of the spec - in the profile section the spec
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
 * usage:      oclfilt [ optional params -bdefhilmnopqrstvwy]
 *             (so note that its default is to use stdin and stdout)
 *
 * where the optional parameters are:
 *             -b <shallower_dlimit>,<deeper_dlimit>
 *                bottom depth filter : only output data for the stations
 *                with bottom depths between <shallower_dlimit>
 *                and <deeper_dlimit>.  These bottom depths are generally
 *                taken from the secondary header info in the stations
 *                (secondary header code #10).
 *                For those stations with no secondary header bottom depth
 *                info, the station is reported if its deepest profile data
 *                depth is within that bottom depth range, but is flagged
 *                as such.  (default yields all stations)
 *             -d <bathy-database filename>
 *                use the bathy values within the named lat-lon-depth file
 *                as the basis for the bottom depth filtering.  This file
 *                is generally created with shell script bathyForThisOCLfile.
 *                (default uses either the bottom depth from secndry header or
 *                from deepest profile depth if sec hdr value is missing/bad)
 *             -e
 *                output "end" statistics for file, instead of full output:
 *                gives number of stations and bytecount for file, according to
 *                the filtering criteria on cmdline
 *                This option superceeds the formatted profile data output.
 *             -f 
 *                full/debugging output : each and every field of data from the
 *                input is listed.  Very lengthy and ugly; just meant for
 *                debugging or looking at example data during development.
 *                This option superceeds the formatted profile data output.
 *             -h
 *                lists brief help/description screen
 *             -i <infilename>
 *                specifies filename of input (default uses stdin)
 *             -l <westbound>/<eastbound>/<southbound>/<northbound>
 *                specifies a lat-lon subregion of interest within the file to
 *                select from the rest.  bounds are in decimal degrees, using
 *                the "/" char as a separator.  Region boundaries are included
 *                in the selected data.  Currently lon bounds don't wrap around
 *                the 0 and 360 marks, they must match the lon style of the
 *                data; this will be fixed at some point.
 *             -m <minmonth>,<maxmonth>
 *                specifies a month range to select data by; eg. -m 1,3
 *                filter is inclusive of both max and min months.
 *                (default does not filter by months - ie returns all months)
 *             -n <numberofstations>
 *                only output <numberofstations> stations:  if
 *                <numberofstations> is greater than number of stations in
 *                file, no error but of course will stop at end of file.
 *                (default does not limit number of stations in this way)
 *             -o <outfilename>
 *                specifies filename of output (default uses stdout)
 *             -p <profilepts>
 *                specify minimum number of profile points (levels) to filter
 *                by.  (default returns station regardless of how few profile
 *                levels it has)
 *             -q
 *                query for variables and data quantity: for each station,
 *                only ouput a list of which variables have data in the
 *                profiles (eg temp, sal, pressure, etc), and the number
 *                of bytes in each profile.
 *                This option superceeds the usual profile data output.
 *             -r
 *                include error-flagged data in output (if it's not already
 *                included) and append error-flags next to each var in output.
 *                (default: if varlist specified with -v option, each profile
 *                level will be output unless it has a nonzero error code
 *                in one of its vars that are in varlist.
 *                if no varlist specified, default outputs all profile levels
 *                for all vars)
 *             -s <stationnumber>
 *                skip to specified station number and start from there.
 *                (default starts with first station in file)
 *             -t
 *                do *NOT* output title header that labels the station number, 
 *                position, time, and data column descriptors in comment lines
 *                which begin with a % character.
 *                This option only works with the usual profile data output.
 *                (default outputs that title header)
 *             -v <var_list>
 *                variables filter : only output profile data for variables
 *                included in <var_list>.  And if any variables in <var_list>
 *                are not in profile or have error codes for the whole
 *                variable, that station is skipped completely.
 *                <var_list> is specified as comma-separated list of the
 *                variable code numbers.  (1=temp, 2=sal, etc. as listed
 *                in Table 4 in NODC OCL readmev1 doc)  eg: -v 1,2
 *                (default outputs data for all variables for all profiles)
 *             -w <wmo_square>
 *                do not output stations which have zero for lat or lon value
 *                when the station is not in a wmo-square that's along the 
 *                equator or prime meridian (respectively).  <wmo_square> is
 *                generally parsed out of the OCL filenames...
 *                (default outputs station even if there are invalid zero
 *                values for lat & lon)
 *             -y <minyear>,<maxyear>
 *                specifies a year range to select data by; eg. -y 1976,1980
 *                filter is inclusive of both max and min years.
 *                (default does not filter by year - ie returns all years)
 *
 * history:
 *     5/05/99-AG-initial program functioning
 *                program does not parse command line yet for above options.
 *                for now I'm hardwiring the ones I'm using just to quickly
 *                get the project done, command line parsing will come next.
 *     6/14/99-AG-separated out file-reading functionality to function
 *                getOCLStationData(), to try to modularize things more and
 *                give others a function to read OCL files without having to
 *                cut & paste & rewrite from original oclfilt.c code,
 *                added the command line parsing/control.
 *     7/20/99-AG-added capability of reading bottom depth values from
 *                bathy-database if desired.  Involved altering stnData struct
 *                slightly in addition to the code changes.
 *     7/27/99-AG-added -w filter for catching lat/lon bad zero values
 *     2/17/00-AG-increased %.2f to %.4f in lat/lon output formatting
 *     2/22/00-AG-added -r flag for including error codes in output, and added
 *                related default of not outputing profile levels that have
 *                errors in individual data within the columns specified by -v
 *     2/23/00-AG-added -p, -l, -m, & -y flags (see above for description)
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ocl.h"



int main (int argc, char **argv) {

   /* these vars are set according to cmdline params, but defaults are here */
   int botDepthFiltFlag=0, debugFlag=0, endStatsFlag=0, numStnsFlag=0;
   int queryFlag=0, skipFlag=0, titlesFlag=1, varListFlag=0;
   int databaseBathyFlag=0, zeroLatLonFlag=0;
   char dbBathyFilename[256], wmoSquare[5];
   int latlonRegionFlag=0, yearRangeFlag=0, monthRangeFlag=0, minLevelsFlag=0;
   double latlonRegion[4];
   long int yearRange[2], monthRange[2];
   long int minLevels;
   int wantProfileFlag, outputThisStation;
   int includeErrorFlaggedData=0;
   long int numVarsOnVarList=0, varList[MAX_VARS];
   long int stnToSkipTo, numStnsToOutput;
   double shallowerDLimit, deeperDLimit;
   FILE *fp_in, *fp_out, *fp_dbBathy=NULL;
 
   /* other vars for just internal bookeeping in main() */
   long int i, j, k, l, totalStationBytes=0, numLevelsWithErrorFlags=0;
   long int stationOutputCount=0, totalStationOutputBytes=0;
   char vars[150], botDepthStr[10], tmp[10];
   int status, errorFlaggedDataExists=0;
 
   OCLStationType stnData;  /* (one station's worth of data in a big struct) */

 
 
   /* Get values from the command line: */
   if( parse_commandline( argc, argv, &fp_in, &fp_out, 
      &botDepthFiltFlag, &shallowerDLimit, &deeperDLimit,
      &varListFlag, &numVarsOnVarList, varList,
      &debugFlag, &endStatsFlag, &titlesFlag, &queryFlag,
      &databaseBathyFlag, dbBathyFilename,
      &numStnsFlag, &numStnsToOutput,
      &skipFlag, &stnToSkipTo,
      &zeroLatLonFlag, wmoSquare,
      &minLevelsFlag, &minLevels,
      &latlonRegionFlag, latlonRegion,
      &yearRangeFlag, yearRange, &monthRangeFlag, monthRange,
      &includeErrorFlaggedData) != SUCCESSFUL )
      exit(1);

      /* Note above that by sending the _addresses_ of the filepointers I made
      it so I can get the filepointers returned to main after they're set in
      the function - that's the reason for the 'FILE **' declarations within
      the parse_commandline function itself, rather than just 'FILE *'. */


   /* If databaseBathy specified on cmdline, then open the bathy file */
   if( databaseBathyFlag ) {
      if ((fp_dbBathy = fopen(dbBathyFilename,"r")) == NULL) {
         fprintf(stderr, "Unable to open bathy file %s.\n", dbBathyFilename);
         exit(1);
      }
   }

 
   /* Set flag - we'll want the profile data if we specified the query or
      formatted output mode (not endStats), or if we're in "spew-everything"
      full output (ie debug) mode */
   wantProfileFlag = !endStatsFlag || debugFlag;


   /* Need to output file header before loop if using query mode (& want hdr)*/
   if(queryFlag && titlesFlag) {
      fprintf(fp_out, "%%  stn year mo dy  time       lat       lon   bytes "
         "numlvls botdepth  vars\n");
      fprintf(fp_out, "%%----- ---- -- -- ----- --------- --------- ------- "
         "------- --------  ----------\n");
   }


   /* loop over stations in this file */
   for (i=0; !feof(fp_in); i++) {

     numLevelsWithErrorFlags=0; /* resetting for new profile */

      /* read in one station of data */
      status = getOCLStationData( fp_in, i, &stnData, wantProfileFlag,
         skipFlag, stnToSkipTo, varListFlag, varList, numVarsOnVarList, 
         minLevelsFlag, minLevels, latlonRegionFlag, latlonRegion,
         yearRangeFlag, yearRange, monthRangeFlag, monthRange,
         databaseBathyFlag, fp_dbBathy, zeroLatLonFlag, wmoSquare );
      if( status==SKIPPED ) continue;    /* skip to next i loop (station) */
      else if( status!=SUCCESSFUL ) {
         fprintf(stderr, 
            "oclfilt: error: failure in getOCLStationData at stn#%ld.\n",i);
         exit(1);
      }

 
      /* keep track of this to compare with totalStationOutputBytes later */
      totalStationBytes += stnData.bytesInStation;

 
      /* Setting up output conditional :
         We want to output/count the data from this station UNLESS -
         1.) a bottomDepth filter was specified and this station was cut by it
          or
         2.) a varList filter was specified and this station was cut by it
          or
         3.) bad-lat/lon filter was specified and a lat or lon is zero when
             our filedata is not near the equator or prime meridian
         4.) latlonRegion was specified and this station was cut by it
         5.) yearRange was specified and this station was cut by it
         6.) monthRange was specified and this station was cut by it
         7.) minLevels was specified and this station was cut by it
         (sorry about using the convoluted "if" statement below, rather than
         logical ops, but it protects against referencing a null pointer...)
      */
      outputThisStation=0;
      if( botDepthFiltFlag && stnData.bottomDepthPtr!=NULL ) {
         if( *(stnData.bottomDepthPtr)>=shallowerDLimit &&
             *(stnData.bottomDepthPtr)<=deeperDLimit ) outputThisStation=1;
         else outputThisStation=0;
      }       
      else outputThisStation=1;
      /* Note below that the *= acts as an appending "and" operator */
      outputThisStation*=( !varListFlag || stnData.varListChecksOut );
      outputThisStation*=( !zeroLatLonFlag || !(stnData.badLatLon) );
      outputThisStation*=( !latlonRegionFlag || stnData.latlonInRange );
      outputThisStation*=( !yearRangeFlag || stnData.yearInRange );
      outputThisStation*=( !monthRangeFlag || stnData.monthInRange );
      outputThisStation*=( !minLevelsFlag || stnData.enoughProfileLevels );


  
      /* If we're going to output the station... */
      if( outputThisStation ) {

         /* need to keep track of these for stats later */
         stationOutputCount++;
         totalStationOutputBytes += stnData.bytesInStation;


         /* full debugging (lengthy & sloppy) output */
         if( debugFlag )
            outputAllStationData( fp_out, i, &stnData );


         /* Query output - one line summary from station's header */
         else if( queryFlag ) {

            /* set up depth string output */
            if(stnData.bottomDepthPtr!=NULL)
              sprintf(botDepthStr, "%6.1f %c", *(stnData.bottomDepthPtr),
                      stnData.bottomDepthSource);
            else
              sprintf(botDepthStr, "   --  -");
        
            /* set up variables string output */
            strcpy(vars,"");
            for(j=0; j<stnData.numberOfVarCodes; j++) {
              sprintf(tmp,"%ld",stnData.varCode[j]);
              strcat(vars, tmp);
              if(stnData.errCodeForVarCode[j]>0) strcat(vars, "*");
              if(j<stnData.numberOfVarCodes-1) strcat(vars, ",");
            }
            if(!strcmp(vars,"")) strcpy(vars,"  --  ");
        
            fprintf(fp_out,
               "%6ld %4ld %2ld %2ld %5.2f %9.4f %9.4f %7ld %7ld %8s  %-9s\n",
               i, stnData.year, stnData.month, stnData.day, stnData.time,
               stnData.lat, stnData.lon, stnData.bytesInStation,
               stnData.numberOfLevels, botDepthStr, vars );
         }
   

         /* Not doing the endStats (or one of the above possibilities) means
            we want the regular formatted output of the profile data. */
         else if( !endStatsFlag ) {
            /* Output title header first if needed */
            if( titlesFlag ) {
               if(stnData.bottomDepthPtr!=NULL)
                  sprintf(botDepthStr,"%.2f m", *(stnData.bottomDepthPtr));
               else strcpy(botDepthStr,"[no data]");
               fprintf(fp_out, "%%\n%%Station #%ld, bottom depth %9s (from %c),"
                  "  %s level data\n",i, botDepthStr, stnData.bottomDepthSource,
                  (stnData.stationType==0) ? "observed" : "standard" );
               fprintf(fp_out, "%%Columns: Lat, Lon, Year, Month, Day, Time, "
                  "Depth");
               for(j=0; j<stnData.numberOfVarCodes; j++)
                  fprintf(fp_out, ", %s", varCodeLabel(stnData.varCode[j]));
               fprintf(fp_out, "\n");
               fprintf(fp_out, "%%Units:   deg, deg, yyyy, mm, dd, hrs, m");
               for(j=0; j<stnData.numberOfVarCodes; j++)
                  fprintf(fp_out, ", %s", varCodeUnits(stnData.varCode[j]));
               fprintf(fp_out, "\n");
            }
            /* Now output the profile data itself */
            for(j=0; j<stnData.numberOfLevels; j++) { /* loop over prof lvls */

 	       /* Conditionals to find whether 'errorFlaggedDataExists' on
		  this profile level, in any of the variables specified as
		  required for this profile (ie in varList) : */
 	       if( varListFlag ) {
		  errorFlaggedDataExists=0;
		  for(k=0; k<stnData.numberOfVarCodes; k++) {
 		     for(l=0; l<numVarsOnVarList; l++) {
		        if( varList[l]==stnData.varCode[k] )
			  if( stnData.errCodeForVarValue[k][j]!=0 /*bad data*/
			      || 
			      !(stnData.varValue[k][j]>0 ||
			        stnData.varValue[k][j]<=0) /*=NaN,ie missing*/)
			      errorFlaggedDataExists=1;
		     }
		  }
	       }

	       /* print out one line = one profile level of output */
 	       if( !errorFlaggedDataExists || includeErrorFlaggedData ) {
		 fprintf(fp_out, "%.4f  %.4f  %4ld %2ld %2ld %.2f  %.2f",
			 stnData.lat, stnData.lon, stnData.year, stnData.month,
			 stnData.day, stnData.time, stnData.depthValue[j] );
		 /* if we want to include error codes in output, append this */
		 if( includeErrorFlaggedData )
		    fprintf(fp_out, " (%ld)", stnData.errCodeForDepthValue[j]);
		 for(k=0; k<stnData.numberOfVarCodes; k++) {
		    fprintf(fp_out, "  %.3f", stnData.varValue[k][j]);
		    /* if we want to include error codes in output, append: */
		    if( includeErrorFlaggedData ) fprintf(fp_out, " (%ld)",
		       stnData.errCodeForVarValue[k][j]);
		 }
		 fprintf(fp_out, "\n");
	       }
	       if( errorFlaggedDataExists ) numLevelsWithErrorFlags++;
            }
         }
      }
 
 
      /* if there was only a specified number of stations we were to output,
         and we've reached that number, break out of station loop to end of
         the program. */
      if( numStnsFlag && stationOutputCount>=numStnsToOutput ) break;
 
   }  /* end of stations loop (i) */
 
 
   /* Output the final statistics if needed */
   if( endStatsFlag || queryFlag ) {
      fprintf(fp_out,"%% summary value units: #Stns / total#Stns, Bytes / totalBytes\n");
      fprintf(fp_out,"%% summary:  %ld / %ld , %ld / %ld\n", stationOutputCount,
              i, totalStationOutputBytes, totalStationBytes);
   }

     
   return SUCCESSFUL;


} /* end of main() */
 







/* "Output all station data" - full, messy output of everything in station for
   debugging purposes */
int outputAllStationData(FILE *fp_out, long int i, OCLStationType *stnData) {

  int status = SUCCESSFUL;
  long int j, k;

  fprintf(fp_out, "bytesInStation(%ld)=%ld\n", i, stnData->bytesInStation);
  fprintf(fp_out, "oclStationNumber(%ld)=%ld\n", i, stnData->oclStationNumber);
  fprintf(fp_out, "countryCode(%ld)=%ld\n", i, stnData->countryCode);
  fprintf(fp_out, "cruiseNumber(%ld)=%ld\n", i, stnData->cruiseNumber);
  fprintf(fp_out, "date(%ld)=%ld-%ld-%ld\n", i, stnData->year, stnData->month,
     stnData->day);
  fprintf(fp_out, "time(%ld)=%f\n", i, stnData->time);
  fprintf(fp_out, "lat(%ld)=%f\n", i, stnData->lat);
  fprintf(fp_out, "lon(%ld)=%f\n", i, stnData->lon);
  fprintf(fp_out, "numberOfLevels(%ld)=%ld\n", i, stnData->numberOfLevels);
  fprintf(fp_out, "stationType(%ld)=%ld\n", i, stnData->stationType);
  fprintf(fp_out, "numberOfVarCodes(%ld)=%ld\n", i, stnData->numberOfVarCodes);
  for(j=0; j<stnData->numberOfVarCodes; j++) {
    fprintf(fp_out, "  varCode(%2ld)=%3ld     errCodeForVarCode(%2ld)=%ld\n",
              j, stnData->varCode[j], j, stnData->errCodeForVarCode[j]);
  }
  fprintf(fp_out, "bytesInCharPI(%ld)=%ld\n", i, stnData->bytesInCharPI);
  fprintf(fp_out, "bytesInSecHdr(%ld)=%ld\n", i, stnData->bytesInSecHdr);
  fprintf(fp_out, "bytesInBioHdr(%ld)=%ld\n", i, stnData->bytesInBioHdr);
  fprintf(fp_out, "numberOfSecHdrEntries(%ld)=%ld\n", i,
     stnData->numberOfSecHdrEntries);
  for(j=0; j<stnData->numberOfSecHdrEntries; j++) {
    fprintf(fp_out, "  secHdrCode(%2ld)=%3ld     secHdrValue(%2ld)=%f\n",
              j, stnData->secHdrCode[j], j,  stnData->secHdrValue[j]);
  }
  fprintf(fp_out, "depth, var1, var2, etc:\n");
  for(j=0; j<stnData->numberOfLevels; j++) {
    fprintf(fp_out, "%f (%ld)     ", 
      stnData->depthValue[j], stnData->errCodeForDepthValue[j]);
    for(k=0; k<stnData->numberOfVarCodes; k++)
      fprintf(fp_out, "%f (%ld)     ", 
        stnData->varValue[k][j], stnData->errCodeForVarValue[k][j]);
    fprintf(fp_out, "\n");
  }
  fprintf(fp_out, "bytesLeftInStation(%ld)=%ld\n", i,
     stnData->bytesLeftInStation);
  fprintf(fp_out, "bottomDepth(%ld)=%f\n", i, *(stnData->bottomDepthPtr));

  return status;

}






/* "Parse Command Line" - get the appropriate command line info for oclfilt */
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
   int *includeErrorFlaggedData ){

  /* note that by using pointers to the filepointers, I made it so I can
     access the filepointers from main after they're set in the function -
     that's the reason for the FILE ** declarations, and why *fp... is used
     below instead of fp... */

  int c, i_flag=0, o_flag=0, status=SUCCESSFUL;
  long int *vp;  /* tmp pointer for filling in varList */
  char filename[256];
  char *tmp;  /* tmp pointer for searching thru *argv for commas */

  /* Loop thru and parse the command line options */
  while (--argc > 0 && (*++argv)[0] == '-') {
    c = *++argv[0];
    switch (c) {
      case 'b': /* bottom depth filter */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *botDepthFiltFlag=1;
          if( sscanf(*argv, "%lf,%lf", shallowerDLimit, deeperDLimit) !=2 ) {
             fprintf(stderr, "The -b param requires an argument of "
                "<shallowerDLimit>,<deeperDLimit>.\n");
             status=UNSPECIFIED_PROBLEM;
          }
        }
        else {
          fprintf(stderr, "The -b param requires an argument of "
                "<shallowerDLimit>,<deeperDLimit>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'd': /* database-bathy file*/
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(databaseBathyFilename,"%s",*argv);
          *databaseBathyFlag=1;
        }
        else {
          fprintf(stderr, "The -d param requires an argument of <filename>\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'e':  /* end-stats flag */
        *endStatsFlag=1;
        break;
      case 'f':  /* full (debug) output flag */
        *debugFlag=1;
        break;
      case 'i': /* input file*/
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(filename,"%s",*argv);
          ++i_flag;
        }
        else {
          fprintf(stderr, "The -i param requires an argument of <infilename>."
                  "\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'l':  /* lat/lon range */
	++argv;
	--argc;
	/* if(*argv!=NULL && *argv[0] != '-') {   this prevents neg lon! */
	if(*argv!=NULL ) {
	  *latlonRegionFlag=1;
          sscanf(*argv, "%lf/%lf/%lf/%lf", &latlonRegion[0], &latlonRegion[1],
	                                   &latlonRegion[2], &latlonRegion[3]);
        }
        else {
          fprintf(stderr, "The -l param requires an argument of <latlonRange>,"
                  " as w/e/s/n in decimal degrees.\n");
          status=UNSPECIFIED_PROBLEM;
        }
	break;
      case 'm':  /* month range */
	++argv;
	--argc;
	if(*argv!=NULL && *argv[0] != '-') {
	  *monthRangeFlag=1;
          sscanf( *argv, "%ld,%ld", &monthRange[0], &monthRange[1] );
        }
        else {
          fprintf(stderr, "The -m param requires an argument of <monthRange>,"
                  " as month1,month2, eg 1,3.\n");
          status=UNSPECIFIED_PROBLEM;
        }
	break;
      case 'n':  /* limit number of stations */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *numStnsFlag=1;
          *numStnsToOutput=atoi(*argv);
        }
        else {
          fprintf(stderr, "The -n param requires an argument of <numStations>."
                  "\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'o': /* output file*/
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(filename,"%s",*argv);
          ++o_flag;
        }
        else {
          fprintf(stderr, "The -o param requires an argument of <outfilename>."
                  "\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'p':  /* require a minimum number of profile points */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *minLevelsFlag=1;
          *minLevels=atoi(*argv);
        }
        else {
          fprintf(stderr, "The -p param requires an argument of <minPts>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'q':  /* query-output flag */
        *queryFlag=1;
        break;
      case 'r':  /* include error-flagged data */
	*includeErrorFlaggedData=1;
	break;
      case 's':  /* skip to this station number */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *skipFlag=1;
          *stnToSkipTo=atoi(*argv);
        }
        else {
          fprintf(stderr,"The -s param requires an argument of <stationnumber>"
                  "\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 't':  /* do NOT print out data-column title headers */
        *titlesFlag=0;
        break;
      case 'v': /* filter by data-variables - only output if these vars good */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *varListFlag=1;
	  tmp=*argv;
	  vp=varList;
	  do {  /* keep allocating varList values as long as there are
                   commas left to delimit them */
	     sscanf(tmp, "%ld", vp++);
	     (*numVarsOnVarList)++;
	     tmp=strchr(tmp,',');
	     if( tmp!=NULL ) tmp++;
	  } while(tmp!=NULL);
        }
        else {
          fprintf(stderr, "The -v param requires an argument of <varlist>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'w': /* WMO square to filter out bad lat/lon zero values */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(wmoSquare,"%s",*argv);
          *zeroLatLonFlag=1;
        }
        else {
          fprintf(stderr,"The -w param requires an argument of <wmo_square>\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'y':  /* year range */
	++argv;
	--argc;
	if(*argv!=NULL && *argv[0] != '-') {
	  *yearRangeFlag=1;
          sscanf( *argv, "%ld,%ld", &yearRange[0], &yearRange[1] );
        }
        else {
          fprintf(stderr, "The -y param requires an argument of <yearRange>,"
                  " eg 1976,1980 (must be 4 digits).\n");
          status=UNSPECIFIED_PROBLEM;
        }
	break;
      case 'h':  /* "help" - show usage listing */
        fprintf(stderr, "\n");
        fprintf(stderr, "oclfilt: Reads an OCL-formatted datafile and "
           "converts to ascii output,\n");
        fprintf(stderr, "         and can filter out certain data from "
           "that input file.\n");
        fprintf(stderr, "         (last compiled: %s, %s)\n\n", __DATE__,
           __TIME__);
        fprintf(stderr, "usage:   oclfilt [optional params -bdefhilmnopqrstvwy]\n");
	fprintf(stderr, "         See oclfilt.manpage for details.\n");
        fprintf(stderr, "         Note that no args assumes stdin & stdout.\n");
        fprintf(stderr, "\n");
        status=HELP_LISTING;
        break;
      default:
        fprintf(stderr, "Illegal Option:  -%c\n", c);
        status=UNSPECIFIED_PROBLEM;
        break;
    }
  }
  /* catch any remaining parsing errors */
  if (argc!=0) {
    fprintf(stderr, "There was some kind of parsing error, probably a\n");
    fprintf(stderr, "missing dash or missing parameter value...\n");
    status=UNSPECIFIED_PROBLEM;
  }

  /* show error messsage and exit if trouble */
  if( status==UNSPECIFIED_PROBLEM ) {
    fprintf(stderr, "For usage list, type oclfilt -h\n\n");
  }
  if( status!=SUCCESSFUL ) return status;

  /* done parsing command line */



  /* assign stdin or open file depending on args */
  if( i_flag ) {
    if ((*fpIn = fopen(filename,"r")) == NULL) {
      fprintf(stderr, "Unable to open file %s.\n", filename);
      status=UNSPECIFIED_PROBLEM;
    }
  }
  else {
    *fpIn = stdin;
  }

  /* assign stdout or open file depending on args */
  if( o_flag ) {
    if ((*fpOut = fopen(filename,"w")) == NULL) {
      fprintf(stderr, "Unable to open file %s.\n", filename);
      status=UNSPECIFIED_PROBLEM;
    }
  }
  else {
    *fpOut = stdout;
  }

  return status;

} /* end of parsing cmd line */

