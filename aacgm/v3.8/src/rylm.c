/* rylm.c
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
 $Log: rylm.c,v $
 Revision 1.7  2004/05/04 22:40:32  barnes
 Changed the header name.

 Revision 1.6  2003/08/05 20:18:08  barnes
 Changed math constant header.

 Revision 1.5  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.4  2001/01/29 18:11:47  barnes
 Added Author Name

 Revision 1.3  1999/02/12 18:07:50  barnes
 Fixed nasty bug in calculating complex products.

 Revision 1.2  1998/06/05 20:46:36  barnes
 Modifications to fix problems with return codes.

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

 */

#include <math.h>
#include "rmath.h"


int rylm(double colat,double lon,int order,
	  double *ylmval) {
   
    double d1;
    struct complex z1, z2;

    /* Local variables */
    struct complex q_fac;
    struct complex q_val;
    int l, m;
    int la,lb,lc,ld,le,lf;

    double cos_theta,sin_theta;  
    double ca,cb;
    double fac;
    double cos_lon, sin_lon;

    cos_theta = cos(colat);
    sin_theta = sin(colat);

    cos_lon = cos(lon);
    sin_lon = sin(lon);

    d1 = -sin_theta;
    z2.x = cos_lon;
    z2.y = sin_lon;

    z1.x = d1 * z2.x;
    z1.y = d1 * z2.y;

    q_fac.x = z1.x;
    q_fac.y = z1.y;

    ylmval[0] = 1;
    ylmval[2] = cos_theta;

    for (l = 1; l <= (order-1); l++) {
	la = (l - 1) * l + 1;
	lb = l * (l + 1) + 1;
	lc = (l + 1) * (l + 2) + 1;

	ca =  ((double) (l * 2 + 1)) / (l + 1);
	cb =  ((double) l) / (l + 1);

	ylmval[lc-1] = ca * cos_theta * ylmval[lb-1] - cb * ylmval[la-1];
    }

    q_val.x = q_fac.x;
    q_val.y = q_fac.y;

    ylmval[3] = q_val.x;
    ylmval[1] = -q_val.y;
    for (l = 2; l <= order; l++) {

	d1 = l*2 - 1.;
	z2.x = d1 * q_fac.x;
	z2.y = d1 * q_fac.y;
	z1.x = z2.x * q_val.x - z2.y * q_val.y;
	z1.y = z2.x * q_val.y + z2.y * q_val.x;
	q_val.x = z1.x;
	q_val.y = z1.y;

	la = l*l + (2*l) + 1;
	lb = l*l + 1;

	ylmval[la-1] = q_val.x;
	ylmval[lb-1] = -q_val.y;
    }

    for (l = 2; l <= order; l++) {

	la = l*l;
	lb = l*l - 2*(l - 1);

	lc = l*l + (2*l);
	ld = l*l + 2;

	fac = l*2 - 1;

	ylmval[lc-1] = fac * cos_theta * ylmval[la-1];
	ylmval[ld-1] = fac * cos_theta * ylmval[lb-1];
    }

    for (m = 1; m <= (order-2); m++) {

	la = (m+1)*(m+1);
        lb = (m+2)*(m+2)-1;
        lc = (m+3)*(m+3)-2;
	
	ld = la - (2*m);
	le = lb - (2*m);
	lf = lc - (2*m);

	for (l = m + 2; l <= order; l++) {
	    ca =  ((double) (2*l - 1)) / (l - m);
	    cb =  ((double) (l+m - 1)) / (l - m);

	    ylmval[lc-1] = ca * cos_theta *ylmval[lb-1] - cb *ylmval[la-1];
	    ylmval[lf-1] = ca * cos_theta *ylmval[le-1] - cb *ylmval[ld-1];

	    la = lb;
	    lb = lc;
	    lc = lb + (2*l) + 2;

	    ld = la - (2*m);
	    le = lb - (2*m);
	    lf = lc - (2*m);

	}
    }
    return 0;
} 














