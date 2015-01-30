/* mlt.c
   =====
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
 $Log: mlt.c,v $
 Revision 1.4  2004/04/26 21:18:37  barnes
 Audit - Enforced warning set to "all".

 Revision 1.3  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.2  2001/01/29 18:11:46  barnes
 Added Author Name

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "eqn_of_time.h"
#include "solar_loc.h"
#include "mlt1.h"

double mlt(int yr,int t0,double mlong, double *mslong) {

    double ret_val;
    double mean_lon;
   
    
    double apparent_time;
    double et;
    double ut;
    double dec;

    if (yr > 1900) yr-=1900;
    solar_loc(yr, t0, &mean_lon,&dec);

    et = eqn_of_time(mean_lon, yr);

    ut=(float) (t0 % (24*3600));
    apparent_time = ut + et;
   
    ret_val = mlt1(apparent_time, dec, mlong, mslong);
    return ret_val;
}






