/* ----------------------------------------------------------------------------
 * 
 *             Applied Physics Laboratory, University of Washington
 * 
 *          TR 9407 High Frequency Ocean Environmental Acoustic Models
 *                          Frequency Range 10-100 kHz
 * 
 * ----------------------------------------------------------------------------
 * 
 *      This code module was developed by the Applied Physics Laboratory,
 *   University of Washington, to implement algorithms supporting the models
 *   described in the APL-UW Technical Report TR 9407 "High Frequency Ocean
 *   Environmental Acoustic Models, Frequency Range 10-100 kHz".
 * 
 * ----------------------------------------------------------------------------
 * 
 * MODEL:  Volume
 * REPORT SECTION:  I
 * MODEL DEVELOPER:  C.T. Chen, F.J. Millero and Xu Li
 * PROGRAMMER:  Kristin Kulman
 * TRANSLATION FROM FORTRAN TO C:  Andy Ganse, 4/30/99
 * 
 * DESCRIPTION:  This function returns a sound speed estimate in m/s, given
 *    pressure, temperature and salinity values as inputs.  Pressure
 *    is in bars. Temperature is in degrees C.  And salinity is in
 *    parts per thousand.  This equation by Millero and Li (1994) is
 *    a recent adjustment to the Chen and Millero equation (1977).
 * 
 *    NOTE - This calculation is frequency independant.
 * 
 *    This function performs boundary checking on the inputs.  If pressure,
 *    temp, or salinity are not within the ranges shown below, the function
 *    returns an error status from which params were out of range.
 * 
 * USAGE:  status = sspcm2 (pressure, temp, salinity, &sndspd)
 * 
 * ARGUMENTS:
 *    INPUT:  pressure- (double) pressure [units = bars]
 *            temp - (double) temperature [degrees C]
 *            salinity - (double) salinity [parts per thousand]
 *   OUTPUT:  pointer to sndspd - (double) sound speed [m/s]
 *            return status - 0=good
 *                            +1=error: pressure out of range
 *                            +2=error: temperature out of range
 *                            +4=error: salinity out of range
 *            Those error values are added together and can be parsed out
 *            (with mod function) to show if there was more than one bad param.
 *            Otherwise one can just say 0 = good, >0 = out-of-range error.
 * 
 * RANGE LIMITATIONS:  0 <= PRESSURE <= 1000 bars
 *                     0 <= TEMP <= 40 deg C
 *                     0 <= SAL <= 40 parts per thousand
 * 
 * CHECK VALUE: 1745.095215 m/s for PRESSURE = 1000 bars, TEMP = 40 deg C, and
 *                                SAL = 40 ppt
 * 
 *    Note:  This check value has been provided to allow verification of the
 *           module's performance against its original form in the case of
 *           future modifications.
 * 
 * REFERENCES:
 * 
 * MODIFIED:
 * 
 * --------------------------------------------------------------------------*/

#include <math.h>

int sspcm2(double P, double T, double S, double *sndspd) {
/* note param units : Pressure (bars), Temp (degrees C), Salinity (ppt) */

  double A;       /* S**1 term */
  double A0, A1;
  double A2, A3;
  double B;       /* S**2 term */
  double B0, B1;
  double C;
  double C0, C1;
  double C2, C3;
  double CC;      /* correction term */
  double CC1;
  double CC2, CC3;
  double D;
  double SR;

  int status=0;   /* return status - 0=good, >0=bad */

  /* EQUIVALENCE (A0,B0,C0),(A1,B1,C1),(A2,C2),(A3,C3)       ???????   */


  /* Input param range validation */
     if ( (P < 0.) || (P > 1000.) ) status+=1;
     if ( (T < 0.) || (T > 40.) ) status+=2;
     if ( (S < 0.) || (S > 40.) ) status+=4;


  /* Return bad flag & exit early if input params not valid */
     if ( status /* is bad */) return status;  /* and exit */

     SR = sqrt(fabs(S));


  /* --S**2 TERM------------- */
     D = (1.727E-3) - (7.9836E-6)*P;

  /* --S**3/2 TERM----------- */
     B1 = (7.3637E-5) + (1.7945E-7)*T;
     B0 = (-1.922E-2) - (4.42E-5)*T;
     B = B0 + B1*P;

  /* --S**1 TERM------------- */
     A3 = ( (-3.389E-13)*T + (6.649E-12) )*T + (1.100E-10);
     A2 = ( ( (7.988E-12)*T - (1.6002E-10) )*T + (9.1041E-9) )*T - (3.9064E-7);
     A1 = ( ( ( (-2.0122E-10)*T + (1.0507E-8) )*T - (6.4885E-8) )*T
          - (1.2580E-5) )*T + (9.4742E-5);
     A0 = ( ( ( (-3.21E-8)*T + (2.006E-6) )*T + (7.164E-5) )*T - (1.262E-2) )*T
          + 1.389;
     A = ( (A3*P + A2)*P + A1)*P + A0;

  /* --S**0 TERM------------- */
     C3 = ( (-2.3643E-12)*T + (3.8504E-10) )*T - (9.7729E-9);
     C2 = ( ( ( (1.0405E-12)*T - (2.5335E-10) )*T + (2.5974E-8) )*T
          - (1.7107E-6) )*T + (3.1260E-5);
     C1 = ( ( ( (-6.1185E-10)*T + (1.3621E-7) )*T - (8.1788E-6) )*T
          + (6.8982E-4) )*T + 0.153563;
     C0 = ( ( ( ( (3.1464E-9)*T - (1.47800E-6) )*T + (3.3420E-4) )*T
          - (5.80852E-2) )*T + 5.03711)*T + 1402.388;

  /* --S**0 CORRECTION TERM-- */
     CC1= ( (1.4E-5)*T - (2.19E-4) )*T + 0.0029;
     CC2= ( (-2.59E-8)*T + (3.47E-7) )*T - (4.76E-6);
     CC3= (2.68E-9);
     CC = ((CC3*P+CC2)*P+CC1)*P;
     C = ((C3*P+C2)*P+C1)*P+C0-CC;

  /* --SOUND SPEED RETURN---- */
     *sndspd = C + (A+B*SR+D*S)*S;
     return 0;

}



