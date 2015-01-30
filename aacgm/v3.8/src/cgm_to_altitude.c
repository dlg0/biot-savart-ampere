/* cgm_to_altitude.c
   =================
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
 $Log: cgm_to_altitude.c,v $
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
#include "math.h"


int cgm_to_altitude(double r_height_in,double r_lat_in,double *r_lat_adj) {
    double eradius = 6371.2;
    double unim=1;
    double  r1;
    double ra;
    int error=0;

   /* Compute the corresponding altitude adjusted dipole latitude. */
   /* Computing 2nd power */

    r1 = cosd(r_lat_in);
    ra = (r_height_in/ eradius+1)*(r1*r1);
    if (ra > unim) {
	ra = unim;
        error=1;
    }

    r1 = acos(sqrt(ra));
    *r_lat_adj = sgn(r1,r_lat_in)*180/PI;
    return error;
}

