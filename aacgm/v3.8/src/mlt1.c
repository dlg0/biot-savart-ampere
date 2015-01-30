/* mlt1.c
   ======
   Author: R.J.Barnes
*/

/*
 Copyright 2004 The Johns Hopkins University/Applied Physics Laboratory.
 All rights reserved.
 
 This material may be used, modified, or reproduced by or for the U.S.
 Government pursuant to the license rights granted under the clauses at DFARS
 252.227-7013/7014.
 
 For any other permissions, please contact the Space Department
 Program Office at JHU/APL.
 
 This Distribution and Disclaimer Statement must be included in all copies of
 "Radar Software Toolkit - SuperDARN Toolkit" (hereinafter "the Program").
 
 The Program was developed at The Johns Hopkins University/Applied Physics
 Laboratory (JHU/APL) which is the author thereof under the "work made for
 hire" provisions of the copyright law.  
 
 JHU/APL assumes no obligation to provide support of any kind with regard to
 the Program.  This includes no obligation to provide assistance in using the
 Program or to provide updated versions of the Program.
 
 THE PROGRAM AND ITS DOCUMENTATION ARE PROVIDED AS IS AND WITHOUT ANY EXPRESS
 OR IMPLIED WARRANTIES WHATSOEVER.  ALL WARRANTIES INCLUDING, BUT NOT LIMITED
 TO, PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE ARE
 HEREBY DISCLAIMED.  YOU ASSUME THE ENTIRE RISK AND LIABILITY OF USING THE
 PROGRAM TO INCLUDE USE IN COMPLIANCE WITH ANY THIRD PARTY RIGHTS.  YOU ARE
 ADVISED TO TEST THE PROGRAM THOROUGHLY BEFORE RELYING ON IT.  IN NO EVENT
 SHALL JHU/APL BE LIABLE FOR ANY DAMAGES WHATSOEVER, INCLUDING, WITHOUT
 LIMITATION, ANY LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR
 CONSEQUENTIAL DAMAGES, ARISING OUT OF THE USE OR INABILITY TO USE THE
 PROGRAM."
 
 
 
 
 
 
*/


/*
 $Log: mlt1.c,v $
 Revision 1.6  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.5  2001/01/29 18:11:46  barnes
 Added Author Name

 Revision 1.4  1998/06/17 18:59:41  barnes
 There appeared to be a problem at certain MLT times because the order
 of the spherical harmonic expansion for the coordinate conversion was too
 high.

 Revision 1.3  1998/06/16 23:48:17  barnes
 Fixed bug that caused the interpolated mlt to be wrong.

 Revision 1.2  1998/06/05 20:46:36  barnes
 Modifications to fix problems with return codes.

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

 */

#include <stdio.h>
#include <math.h>
#include "math.h"
#include "convert_geo_coord.h"


static double sol_dec_old=0;
static double told=1e12;
static double mslon1=0;
static double mslon2=0;


double mlt1(int t0,double solar_dec,double mlon,double *mslon) {
 
  double ret_val;
  double mslat1,mslat2,slon1,slon2,height;
  int err;
 
  if ((fabs(solar_dec-sol_dec_old)>0.1) || (sol_dec_old==0)) told=1e12;
  if (fabs(mslon2-mslon1)>10) told=1e12;
    
   if ((t0>=told) && (t0<(told+600))) {
    *mslon=mslon1+(t0-told)*(mslon2-mslon1)/600.0;
   } else {
    told=t0;
    sol_dec_old=solar_dec;

    slon1 = (43200.0-t0)*15.0/3600.0;
    slon2 = (43200.0-t0-600)*15.0/3600.0;

    height = 450;

    err=convert_geo_coord(solar_dec,slon1,height,&mslat1,&mslon1,0,4);
    err=convert_geo_coord(solar_dec,slon2,height,&mslat2,&mslon2,0,4);
    *mslon=mslon1;
  }

  
  ret_val = (mlon - *mslon) /15.0 + 12.0;
  if (ret_val >=24) ret_val -=24;
  if (ret_val <0) ret_val+=24;  
  return ret_val;
} 


