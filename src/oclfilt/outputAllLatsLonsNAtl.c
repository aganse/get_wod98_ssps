#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ocl.h"


int main (int argc, char *argv[]) {

  int skipFlag=0, varListFlag=0, wantProfileFlag=0, dbBathyFlag=0;
  int  zeroLatLonFlag, badLatFlag=0, badLonFlag=0;
  char wmoSquare[5]="";
  long int i, stnToSkipTo=0, varList[1], numVarsOnVarList=0;
  FILE *fp_dbBathy=NULL;
  double badLat=0., badLon=0.;
  OCLStationType stnData;

  if( argc!=4 ) {
    printf("usage: outputAllLatsLons <wmo_square> <bad-lon> <bad-lat>\n");
    exit(1);
  }
  else {
    strncpy(wmoSquare,argv[1],4);
    zeroLatLonFlag=0;
    badLon=atof(argv[2]);
    badLat=atof(argv[3]);
  }


  /* loop over stations in this file */
  for (i=0; !feof(stdin); i++) {
    if( getOCLStationData( stdin, i, &stnData, wantProfileFlag,
      skipFlag, stnToSkipTo, varListFlag, varList, numVarsOnVarList,
      dbBathyFlag, fp_dbBathy, zeroLatLonFlag, wmoSquare) != SUCCESSFUL ) {
      fprintf( stderr,
         "outputAllLatsLons: failure in getOCLStationData at stn#%ld.\n", i );
      exit(1);
    }


    /* check for zero-valued lats/lons in wmo-squares where there should'nt
       be any */
    badLatFlag=badLonFlag=0;
       /* (remember can't rely on doubles being exactly equal...) */
    if( stnData.lon<0.0000001 && stnData.lon>-0.0000001 ) {  /*ie, if zero*/
       badLonFlag = !zeroLatLonOkay( wmoSquare, "lon" );
    }
    if( stnData.lat<0.0000001 && stnData.lat>-0.0000001 ) {  /*ie, if zero*/
       badLatFlag = !zeroLatLonOkay( wmoSquare, "lat" );
    }

    
    /* catch lats that are out of db bounds : >72 or <-72 */
    if( stnData.lat>72. || stnData.lat<-72. ) badLatFlag=1;


    /* substitute valid placeholder values for bad lats or lons, so that
       grdtrack will still return a line for this location - oclfilt will
       filter out the bad ones, but this way grdtrack won't completely skip
       their line */
    if(badLonFlag || badLatFlag) {
       stnData.lon=badLon;
       stnData.lat=badLat;
    }

    printf("%f  %f %ld\n", stnData.lon, stnData.lat, i);
  }

  return SUCCESSFUL;

}  /* end of main */
