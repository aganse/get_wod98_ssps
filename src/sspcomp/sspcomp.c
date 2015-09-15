/* sspcomp.c - version 1.11
 *             "Sound Speed Comparison" - meant to mainly take output from 
 *             oclfilt via stdin, to help show us whether there is a strong
 *             effect from salinity variations in the calculation of ssp in
 *             shallow water depths.  sspcomp tablulates ssp values which are
 *             calculated both from the actual salinity and from a constant
 *             comparison salinity (specified on the command line)
 *             Values are binned according to profile depths for consistancy.
 *             Output columns are:
 *              lat, lon, year, month, day, time, depthbin, temp, sal, sndspd,
 *               [comp_sal, comp_sndspd, diff_sndspd, [stdev_ssp, numbins] ]
 *             sspcomp assumes input data is grouped by station, each within
 *             profile depth order (as oclfilt outputs).
 * 
 * required sources/files: sspcomp.c, sspcm2.c, Makefile
 *
 * language:   ANSI C
 *
 * author:     Andy Ganse, Applied Physics Laboratory, University of Washington
 *             (aganse@apl.washington.edu)
 *             As of this writing, this sspcomp-1.11 package is available at:
 *             http://staff.washington.edu/aganse/src
 *             Please email me if you use this code so I can alert you to 
 *             updates and have an estimate for how many folks use it.
 *
 * notes:      sspcm2 and depth2pres were from Mike Boyd & Kristin Kulman
 *             (I translated the FORTRAN to to C).  See credits in sspcm2.c
 *             for details.
 *             (the little formula in depth2pres was actually just gleaned out
 *             of tsspcm2.f - "test sspcm2")
 * 
 * usage:      sspcomp [optional params -dhilsASt]
 *             (so note that its default is to use stdin and stdout)
 *
 * where the optional parameters are:
 *             -d <depthbinsize>
 *                specify depth bin size in meters.  profile depths and their
 *                data will be binned and statistically reported within bins
 *                of: 0-depthbinsize meters, depthbinsize-2depthbinsize meters, 
 *                2depthbinsize-3depthbinsize meters, etc.  Output depthbin
 *                values are the lower end of the bin (ie start with 0).
 *                So for a depthbinsize of 10.00 meters, the bins would be
 *                0-10.00, 10.00-20.00, 20.00-30.00, etc.  The output
 *                listing then will list 0.00, 10.00, 20.00,...
 *                (default uses no depth binning at all - report on every depth)
 *             -h 
 *                show help/usage listing
 *             -i <infilename>
 *                specify input file (default assumes stdin)
 *             -l <labelstring>
 *                specify an extra header label line to add to top of output.
 *                <labelstring> must be in quotes, and may not be more than
 *                77 chars in length.  (default is no extra label line)
 *             -s <compsal>
 *                specify a constant comparison salinity value
 *                May not be used with -S or -A.
 *                (default without this or -S or -A param outputs no comparison
 *                soundspeeds or salinities at all)
 *             -S [win_sal_file,spr_sal_file,sum_sal_file,fall_sal_file]
 *                specify filenames with a databases of seasonal salinity
 *                values to use for comparison.  Files are based on 5deg
 *                convention from WOA94.  (default with -S but no sal_filenames
 *                specified is to use the filenames:
 *                  sal13m.5d, sal14m.5d, sal15m.5d, sal16m.5d;
 *                May not be used with -s or -A.
 *                (default without this or -s or -A param outputs no comparison
 *                soundspeeds or salinities at all)
 *             -A [compsal_filename]
 *                specify a filename with a database of annual salinity
 *                values to use for comparison.  File is based on 5deg
 *                convention from WOA94.  (default with -A but no
 *                compsal_filename specified is to use filename sal00m.5d;
 *                May not be used with -s or -S.
 *                (default without this or -s or -S param outputs no comparison
 *                soundspeeds or salinities at all)
 *             -t
 *                DON'T show title header (default shows header)
 *
 * history:
 *     5/09/99-AG-initial program functioning
 *     5/14/99-AG-modified to bin data by depth (-d option)
 *     6/16/99-AG-modified to use compsals from salinity db file (-S option)
 *     6/29/99-AG-fixed binning so that last line is output.
 *     9/07/99-AG-added capability to use seasonal salfiles for comparison,
 *                where before could only use annual salfile.
 *     2/28/99-AG-modified so that no -s, -S, or -A option means that just 
 *                measured sndspds are calculated, no comparison ones.
 *                Also, in this mod, the Calc'dSSP and CompSal columns swapped
 *                places to accomodate the "dropping off" of the extra columns
 *     9/07/01-AG-modified so that if input profile from oclfilt has no salin
 *                column, then 35ppt (worldwide salinity average) is assumed
 *                and used in sndspd calculation.  Before, to acheive this
 *                I cheated and used "-s 35" on cmdline and separately parsed
 *                out "comparison sndspeed" from output; now the substitution
 *                is done automatically when input has no salinity column, and
 *                the output lists a comment when this happens.
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* function return statuses */
#define SUCCESSFUL 0
#define FAILED 1
#define HELP_LISTING 2
#define UNSPECIFIED_PROBLEM 3

/* array sizes */
#define MAX_SDEPTHS 33
#define MAX_LAT_INDS 36
#define MAX_LON_INDS 72
#define MAX_BIN_ARRAY 100

/* compSalType's: */
#define NONE 0
#define ANNUAL 1
#define SEASONAL 2
#define CONST 3


/* function prototypes */
int parse_commandline(int argc, char **argv, FILE **fp_In, FILE **fp_Out,
  double *compSal, double *depthBinSize, int *depthBinsUsed,
  int *compSalType, double salArray[][MAX_SDEPTHS][MAX_LAT_INDS][MAX_LON_INDS],
  int *showTitleHeader, char *labelString);
int read5DegData( FILE *fpIn,
   double data5Deg[MAX_SDEPTHS][MAX_LAT_INDS][MAX_LON_INDS] );
int sspcm2(double pres, double temp, double sal, double *sndspd);
double depth2pres(double depth);
double stdev(double *cumDiffSsp, double avg, long int N);
double nan();
int getStdLevelInd(double depth);
int getLatInd(double lat);
int getLonInd(double lon);



main(int argc, char *argv[]) {

  long int i, j, N=0;
  int depthBinsUsed=0, showTitleHeader=1, firstLine=1, newDepthBin, newStation;
  int compSalType=0, statusActual, statusComp, lastLinePassed=0;
  int status=SUCCESSFUL, salPresent=1;
  int year, month, day, oldyear=-1, oldmonth=-1, oldday=-1, season;
  static double salArray[4][MAX_SDEPTHS][MAX_LAT_INDS][MAX_LON_INDS];
  double time, oldtime=-1.;
  double badValue=nan();  /* just assigns NaN */
  double lat, oldLat=361., lon, oldLon=361.;
  double depth, depthBin=0., depthBinSize=10.00, pres, temp,sal,compSal=35.000;
  double sspActual, sspComp, diffSsp;
  double cumTemp[MAX_BIN_ARRAY],cumSal[MAX_BIN_ARRAY],cumSspComp[MAX_BIN_ARRAY];
  double cumSspActual[MAX_BIN_ARRAY],cumDiffSsp[MAX_BIN_ARRAY];
  double cumCompSal[MAX_BIN_ARRAY];
  double avgTemp=0., avgSal=0., avgSspComp=0., avgSspActual=0., avgDiffSsp=0.;
  double avgCompSal=0., stdevDiffSsp;
  char inputLine[256], labelString[78]="";
  FILE *fpIn, *fpOut;



  /* Get params from the command line: */
  status=parse_commandline( argc, argv, &fpIn, &fpOut, &compSal, &depthBinSize,
     &depthBinsUsed, &compSalType, salArray, &showTitleHeader, labelString);
  if( status!=SUCCESSFUL ) {
    if( status!=HELP_LISTING )
      fprintf(stderr, "sspcomp: parse_commandline() failed: \n");
    exit(FAILED);
  }
  /* Note above that by sending the addresses of the filepointers I made it so
     I can get the filepointers returned to main after they're set in the
     function - that's the reason for the FILE ** declarations (rather than
     just FILE * ) within the function itself. */



  /* Output title header if specified in cmdline */
  if(showTitleHeader) {
    /* if want the extra label line in header, add it: */
    if( strcmp(labelString,"") ) printf("%% %s\n", labelString);
    /* data column labels */
    printf("%%%7s %8s %4s %2s %2s %5s %8s %8s %8s %9s", "Lat  ", "Lon  ",
       "Year", "Mo", "Dy", " Time", "Depth ", "Temp  ", "Saln ", "Calcd_SSP" );
    if(compSalType!=NONE) {
      printf(" %8s %9s %7s", "CompSaln", "CompSSP", "DiffSSP");
      if(depthBinsUsed) printf(" %7s %2s", "StdvDif", "N");
    }
    printf("\n");
    /* and units */
    printf("%%%7s %8s %4s %2s %2s %5s %8s %8s %8s %9s", "deg  ", "deg  ",
	   "yyyy", "mm", "dd", "hrs", "meters ", "deg C ", "ppt ", "m/s   ");
    if(compSalType!=0) {
      printf(" %8s %9s %7s", "ppt ", "m/s ", "m/s ");
      if(depthBinsUsed) printf(" %7s %2s", "m/s ", "#");
    }
    printf("\n");
    /* and the header line */
    printf("%%--------------------------------------------------------"
           "----------------");
    if(compSalType!=0) {
      printf("-----------------------------");
      if(depthBinsUsed) printf("----------");
      }
    printf("\n");
  }



  /* Loop over the lines in the input stream */
  for (i=0; !lastLinePassed; i++) {

    /* get a line of data from the input file */
    if(fgets(inputLine, 255, fpIn)==NULL) lastLinePassed=1;

    /* comment lines - we want to keep the station-info line,
       check for existence of salinity column in input, and toss the other
       comment lines; and afterwards skip to next line-reading */
    if( !strncmp(inputLine, "%Station", 8) ) {
      printf("%s", inputLine);
      continue;
    }
    else if( !strncmp(inputLine, "%Columns", 8) ) {
      if( strstr(inputLine,"Sal")==NULL) {
        salPresent=0;
        sal=35.;
        printf("%%(salinity data not present in input profile - assuming 35ppt.)\n");
      } else salPresent=1;
      continue;
    }
    else if(inputLine[0]=='%' && !lastLinePassed) continue;

    /* If last line of input file not passed already, read and compute data
       for current line */
    if(!lastLinePassed) {

      /* Read the line's data into vars */
	  if(salPresent)
      sscanf(inputLine,"%lf %lf %d %d %d %lf %lf %lf %lf", &lat, &lon, &year,
        &month, &day, &time, &depth, &temp, &sal);
	  else
      sscanf(inputLine,"%lf %lf %d %d %d %lf %lf %lf", &lat, &lon, &year,
        &month, &day, &time, &depth, &temp);
  
      /* If using salfile for salinities, look up sal for this region/depth */
      if( compSalType == ANNUAL ) {
        compSal =
          salArray[0][getStdLevelInd(depth)][getLatInd(lat)][getLonInd(lon)];
      }
      else if( compSalType == SEASONAL ) {
        if( month>=1 && month<=12 ) {
          season = (month-1)/3;  /* note data seasons were only defined via
                                    month, not down to day. */
          compSal = salArray[season][getStdLevelInd(depth)][getLatInd(lat)][getLonInd(lon)];
        }
	else compSal=badValue; /* ie if bad month value can't find db value. */

      }
      /* (else the salFile wasn't used - either there's a const salinity of
	 compSal or no comparisons will be outputted) */


      /* catching the bad-value of -99.999999, which is what NODC used,
         and changing it to a more clear one */
      if( compSal<-99. && compSal>-101. ) compSal=badValue;

   
      /* Calculate actual ssp value from input data */
      pres=depth2pres(depth);
      statusActual = sspcm2(pres, temp, sal, &sspActual);
      if(statusActual!=0) sspActual=badValue;  /* set bad flag if error */

      if(compSalType!=0) {  
	/* Calculate comparison (const-sal based) ssp value */
	statusComp = sspcm2(pres, temp, compSal, &sspComp);
	if(statusComp!=0) sspComp=badValue;  /* set bad flag if error */
  
	/* Calculate ssp diff values */
	if(!statusActual && !statusComp /* both good values */)
	  diffSsp = sspActual - sspComp;
	else diffSsp = badValue; /* set bad flag if sspActual or sspComp bad */
      }

    }


      
    /* Outputting the data :
       If binning data, there's some data processing before outputting data;
       if not binning, just output the single resulting line (at the "else") */
    if(depthBinsUsed) {

      /* (setting flags for the conditional that follows) */
      newDepthBin = depth>=(depthBin+depthBinSize);
      newStation  = lat!=oldLat || lon!=oldLon || year!=oldyear ||
                    month!=oldmonth || day!=oldday ||
                    (int)(time*100)!=(int)(oldtime*100); /* <-- (since can't
                                                             reliably compare
                                                            floating points) */

      /* If the depth bin or station changes, or if it's the last line of the
         input file, output data & reinit arrays */
      if( ( !firstLine && (newDepthBin || newStation) ) || lastLinePassed ) {

        /* calculate avgs from arrays */
         for(j=0;j<N;j++) {
           avgTemp+=cumTemp[j];
           avgSal+=cumSal[j];
           avgSspActual+=cumSspActual[j];
	   if(compSalType!=0) {
	     avgCompSal+=cumCompSal[j];
	     avgSspComp+=cumSspComp[j];
	     avgDiffSsp+=cumDiffSsp[j];
	   }
         }
         avgTemp/=(double)N;
         avgSal/=(double)N;
         avgSspActual/=(double)N;
	 if(compSalType!=0) {
	   avgCompSal/=(double)N;
	   avgSspComp/=(double)N;
	   avgDiffSsp/=(double)N;
	   
	   /* calculate stdev from avg and array */
	   stdevDiffSsp=stdev(cumDiffSsp, avgDiffSsp, N);
	 }

         /* output bin data line */
         printf("%7.4lf %7.4lf %4d %2d %2d %5.2lf %8.3lf %8.3lf %8.3lf %9.3lf",
           oldLat, oldLon, oldyear, oldmonth, oldday, oldtime, depthBin,
	   avgTemp, avgSal, avgSspActual);
	 if(compSalType!=0) {
	   printf(" %8.3lf %9.3lf %7.3lf %7.3lf %2ld", avgCompSal, avgSspComp,
		  avgDiffSsp, stdevDiffSsp, N);
	 }
	 printf("\n");

        /* reinitialize arrays and avgs */
        for(j=0; j<MAX_BIN_ARRAY; j++) {
          cumTemp[j]=0.;
          cumSal[j]=0.;
          cumCompSal[j]=0.;
          cumSspComp[j]=0.;
          cumSspActual[j]=0.;
          cumDiffSsp[j]=0.0;
        }
        avgTemp=0.;
        avgSal=0.;
        avgCompSal=0.;
        avgSspComp=0.;
        avgSspActual=0.;
        avgDiffSsp=0.;
        N=0;

        /* if lat-lon-datetime changed, reset bins */
        newStation  = lat!=oldLat || lon!=oldLon || year!=oldyear ||
                      month!=oldmonth || day!=oldday ||
                      (int)(time*100)!=(int)(oldtime*100); 
        if( newStation ) depthBin=0.;   /* reset bins */
        else depthBin+=depthBinSize;    /* increment to next depth bin */

      }

      /* if there's no data within new bin, skip to next appropriate bin */
      newDepthBin = depth>=depthBin+depthBinSize;
      if ( newDepthBin )
        for(; depth>=depthBin+depthBinSize; depthBin+=depthBinSize);

      /* copy lat-lon-datetime to old-lat-lon-datetime vars for next round */
      oldLat=lat;
      oldLon=lon;
      oldyear=year;
      oldmonth=month;
      oldday=day;
      oldtime=time;

      /* accumulate data from current line into array */
      cumTemp[N]=temp;
      cumSal[N]=sal;
      cumSspActual[N]=sspActual;
      if(compSalType!=0) {
	cumCompSal[N]=compSal;
	cumSspComp[N]=sspComp;
	cumDiffSsp[N]=diffSsp;
      }
      N++;

      firstLine=0;
    }

    else if (!depthBinsUsed && !lastLinePassed) {
      /* just output the single resulting line of data */
      printf("%7.4lf %7.4lf %4d %2d %2d %5.2lf %8.3lf %8.3lf %8.3lf %9.3lf",
             lat, lon, year, month, day, time, depth, temp, sal, sspActual);
      if(compSalType!=0)
	printf(" %8.3lf %9.3lf %7.3lf", compSal, sspComp, diffSsp);
      printf("\n");
  }


  }  /* end of loop over lines in input stream */

  return SUCCESSFUL;

}  /* end of main */




/* "Depth to Pressure" - quick conversion function */
double depth2pres(double depth) {
  return .1 * depth / .99;
}





/* "Parse Command Line" - get the appropriate command line info for sspcomp */
int parse_commandline(int argc, char **argv, FILE **fp_In, FILE **fp_Out,
  double *compSal, double *depthBinSize, int *depthBinsUsed,
  int *compSalType, double salArray[][MAX_SDEPTHS][MAX_LAT_INDS][MAX_LON_INDS],
  int *showTitleHeader, char *labelString) {
  /* (note that by using pointers to the filepointers, I can access the
     filepointers from main after they're set in this function - that's of
     course the reason for the FILE ** declarations, and why *fp... is used
     below instead of fp...) */

  int i, i_flag=0, c, status=SUCCESSFUL;
  char inFileName[256], salFileName[4][256];
  FILE *fpSal;

  /* in case none specified from cmdline options below: */
  strcpy(labelString,"");

  /* Loop thru and parse the command line options */
  while (--argc > 0 && (*++argv)[0] == '-') {
    c = *++argv[0];
    switch (c) {
      case 'd': /* depth bin size */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *depthBinSize=atof(*argv);
          *depthBinsUsed=1;
        }
        else {
          printf("The -d param requires an argument of <depthBinSize>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'i': /* input file*/
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(inFileName,"%s",*argv);
          ++i_flag;
        }
        else {
          printf("The -i param requires an argument of <infilename>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'l': /* label string */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          sprintf(labelString,"%s",*argv);
        }
        else {
          printf("The -l param requires an argument of <labelString>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 's':  /* comparison-salinity value */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          *compSal=atof(*argv);
	  *compSalType=CONST;
        }
        else {
          printf("The -s param requires an argument of <compsal>.\n");
          status=UNSPECIFIED_PROBLEM;
        }
        break;
      case 'A': /* use annual salinity db file, default name is "sal00m.5d" */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') sprintf(salFileName[0],"%s",*argv);
        else strcpy(salFileName[0],"sal00m.5d");
        *compSalType=ANNUAL;
        break;
      case 'S': /* use seasonal salinity db files, default names from WOA94 */
        ++argv;
        --argc;
        if(*argv!=NULL && *argv[0] != '-') {
          if( sscanf( *argv,"%s,%s,%s,%s", salFileName[0], salFileName[1],
              salFileName[2], salFileName[3] ) != 4 ) {
            printf("The optional arg after the -S param must be of form:\n");
            printf("[wint_filename,spr_filename,sum_filename,fall_filename]\n");
            status=UNSPECIFIED_PROBLEM;
          }
        }
        else {
          strcpy(salFileName[0],"sal13m.5d");
          strcpy(salFileName[1],"sal14m.5d");
          strcpy(salFileName[2],"sal15m.5d");
          strcpy(salFileName[3],"sal16m.5d");
        }
        *compSalType=SEASONAL;
        break;
      case 'h':  /* "help" - show usage listing */
        printf("\n");
        printf("sspcomp: Attaches calculated soundspeed values (and\n");
        printf("     optionally comparison soundspeed values based on a\n");
        printf("     comparison salinity value) to an input line that\n");
	printf("     includes depth, temp, salinity, and other possible\n");
	printf("     values as outputted from the oclfilt program.\n");
        printf("     The ssp calculation model has a limited input domain -\n");
        printf("     valid for Depth<9900 m, 0<Temp<40 C, 0<Sal<40 ppt.\n");
        printf("     If out-of-range values entered, there will be NaN\n");
        printf("     (\"not a number\") flags listed in the fields\n");
        printf("     affected by the the out-of-range value.\n");
        printf("     (last compiled: %s, %s)\n\n", __DATE__, __TIME__);
        printf("usage: sspcomp [-s <comparison_salinity> | -A [salFile] |\n");
        printf("            -S [winSalFile,sprSalFile,sumSalFile,fallSalFile]"
               " ]\n");
        printf("           [-d <depthbinsize>] [-l <labelstring>] [-t]\n");
        printf("           [-i <infilename>] [-h]\n");
	printf("     Note that no args assumes stdin & stdout.\n");
	printf("     See sspcomp.manpage for more details.\n\n");
        status=HELP_LISTING;
        break;
      case 't':  /* DON'T show title header */
        *showTitleHeader=0;
        break;
      default:
        printf("Illegal Option:  -%c\n", c);
        status=UNSPECIFIED_PROBLEM;
        break;
    }
  }
  /* catch any remaining parsing errors */
  if (argc > 0) {
    printf("There was some kind of parsing error, probably a\n");
    printf("missing dash or missing parameter value...\n");
    status=UNSPECIFIED_PROBLEM;
  }
  /* show error messsage and exit if trouble */
  if( status==UNSPECIFIED_PROBLEM ) {
    printf("For usage list, type sspcomp -h\n");
  }
  if( status!=SUCCESSFUL ) return status;



  /* assign stdin or open file depending on args */
  if( i_flag ) {
    if ((*fp_In = fopen(inFileName,"r")) == NULL) {
      printf("Unable to open file %s.\n", inFileName);
      return FAILED;
    }
  }
  else {
    *fp_In = stdin;
  }

  /* Get file data for salArray (either annual or seasonal) if specified */
  if( *compSalType==ANNUAL ) {

    if ((fpSal = fopen(salFileName[0],"r")) == NULL) {
      printf("Unable to open annual salinity file %s.\n", salFileName[0]);
      return FAILED;
    }

    /* Read the 5deg sal data into salArray */
    if( read5DegData( fpSal, salArray[0] ) == FAILED ) {
       fprintf(stderr, "parse_commandline: read5DegData failed.\n");
       return FAILED;
    }

  }
  else if( *compSalType==SEASONAL ) {

    for(i=0; i<4; i++) {

      if ((fpSal = fopen(salFileName[i],"r")) == NULL) {
        printf("Unable to open file %s.\n", salFileName[i]);
        return FAILED;
      }
  
      /* Read the 5deg sal data into salArray */
      if( read5DegData( fpSal, salArray[i] ) == FAILED ) {
         fprintf( stderr,
           "parse_commandline: read5DegData failed for salfile %d.\n", i );
         return FAILED;
      }
    }
  }

  return SUCCESSFUL;

} /* end of parse_commandline */






int read5DegData( FILE *fpIn,
   double data5Deg[MAX_SDEPTHS][MAX_LAT_INDS][MAX_LON_INDS] ) {

   int stdDepthLevel, lat, lon;

   for( stdDepthLevel=0; stdDepthLevel<MAX_SDEPTHS; stdDepthLevel++ ) {

      for( lat=0; lat<MAX_LAT_INDS; lat++ ) {

         for( lon=0; lon<MAX_LON_INDS; lon++ ) {

            /* read single value */
            if( fscanf( fpIn, "%8lf", &data5Deg[stdDepthLevel][lat][lon] ) ==
               EOF ) {
               fprintf( stderr, "read5DegData: unexpected EOF in salinfile.\n");
               return FAILED;
            }
         }
      }
   }

   return SUCCESSFUL;
}




/* get OCL stdLevel index to reference depth to salinity array */
int getStdLevelInd(double depth) {

   int i;
   int stdLevelDepth[] = { 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250,
      300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500,
      1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000,
      7500, 8000, 8500, 9000 };

   /* note this is just a quick 'n dirty approximator to find the nearest
      stdlvl depth - remember 'int' doesn't even round, it just truncates... */
   for( i=0; i<40; i++) {
      if( (int)depth==stdLevelDepth[i] ) return i;
      if( (int)depth < ( stdLevelDepth[i] + 
                       (stdLevelDepth[i+1]-stdLevelDepth[i])/2 ) ) return i;
      if( (int)depth < stdLevelDepth[i+1] ) return i+1;
   }

   fprintf(stderr,
      "getStdLevelInd: no stdlevel depth match, apparently depth > 9000m\n");
   fprintf(stderr,
      "                depth value was %lf, used stdlevel 0 for now...\n",
      depth );

}




/* get OCL 5degree latitude index for salinity array */
int getLatInd(double lat) {

/* converting lat [-90,+90] to [0,180] with 0=Spole, to correspond to
   OCL's 5degree indices */
lat+=90.;

/* sorry, this is rather silly - first I'm rounding to nearest 5, then I need
   to round to nearest int to finish conversion, but no native rounding func
   in std libs so using floor(x+0.5).
   Also, the -1 is to accomodate our 0-35 array from 1-36.  */
return (int)( floor( (lat+2.5)/5. + .5 ) - 1 );

}




/* get OCL 5degree latitude index for salinity array */
int getLonInd(double lon) {

/* converting lon (-180,+180] to [0,360) (0 remains Greenwich), to correspond
   to OCL's 5degree indices */
if( lon<0 ) lon+=360.;

/* sorry, this is rather silly - first I'm rounding to nearest 5, then I need
   to round to nearest int to finish conversion, but no native rounding func
   in std libs so using floor(x+0.5).
   Also, the -1 is to accomodate our 0-71 array from 1-72.  */
return (int)( floor( (lon+2.5)/5. + .5 ) - 1 );

}





double stdev(double *cumDiffSsp, double avg, long int N) {

  long int i;
  double sum=0.;

  for(i=0; i<N; i++)
    sum+=pow( avg-cumDiffSsp[i] , 2 );

  return sqrt( sum / N );

}


double nan() {
  double x=0;
  return sqrt(-1/x);
}
