/* convert_geo_coord.c
   ===================
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
 $Log: convert_geo_coord.c,v $
 Revision 1.6  2004/04/26 21:18:37  barnes
 Audit - Enforced warning set to "all".

 Revision 1.5  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.4  2001/01/29 18:11:46  barnes
 Added Author Name

 Revision 1.3  1998/06/17 18:59:41  barnes
 Changed the code so that the order of the fit can be specified.

 Revision 1.2  1998/06/17 12:51:32  barnes
 Fixed problem with sign of latitude sometimes being reversed.

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

 */

#include <stdio.h>
#include <math.h>
#include "math.h"

#include "altitude_to_cgm.h"
#include "cgm_to_altitude.h"
#include "rylm.h"

extern struct {
  double coef[121][3][5][2];
} sph_harm_model;

static double cint[121][3][2]; 

/* Table of constant values */

static double height_old[2]={-1,-1};
static double first_coeff_old=-1;

int convert_geo_coord(double lat_in,double  lon_in,
		      double height_in,double *lat_out,
		      double *lon_out,int flag,
		      int order) {

    
    int i, j, l, m, k;
    int i_err64;

    extern int rylm();
    double ylmval[121];
    double colat_temp;
    double lon_output;
    
    double lat_adj=0;
    double lat_alt=0;
    double colat_input; 
   
    double alt_var_cu=0, lon_temp=0, alt_var_sq=0, alt_var_qu=0;
    double colat_output=0, r=0, x=0, y=0, z=0;
    double alt_var=0;
    double lon_input=0;

    if (lon_in<0) lon_in+=360.0;  

    if (first_coeff_old != sph_harm_model.coef[0][0][0][0]) {
	height_old[0] = -1.0;
	height_old[1] = -1.0;
    }
    first_coeff_old= sph_harm_model.coef[0][0][0][0];

    if ((height_in < 0) || (height_in > 7200)) return -2;
    else if ((flag < 0) || (flag > 1)) return -4; 
    else if (fabs(lat_in) >90.) return -8;
    else if ((lon_in<0) || (lon_in >360)) return -16;
       
    if (height_in != height_old[flag]) {
	alt_var= height_in/7200.0;
	alt_var_sq = alt_var * alt_var;
	alt_var_cu = alt_var * alt_var_sq;
	alt_var_qu = alt_var * alt_var_cu;

	for (i=0; i<3; i++) {
	    for (j=0; j<121;j++) {
		cint[j][i][flag] =sph_harm_model.coef[j][i][0][flag]+
                sph_harm_model.coef[j][i][1][flag]*alt_var+
                sph_harm_model.coef[j][i][2][flag]*alt_var_sq+
                sph_harm_model.coef[j][i][3][flag]*alt_var_cu+
                sph_harm_model.coef[j][i][4][flag]*alt_var_qu;
	    }
	}
	height_old[flag] = height_in;
    
    }

    x = 0;
    y = 0;
    z = 0;

    lon_input =lon_in*PI/180.0;

    if (flag == 0) colat_input = (90-lat_in)*PI/180.0;
    else {
      i_err64=cgm_to_altitude(height_in, lat_in, &lat_adj);
	       
      if (i_err64 != 0) return -64;
      colat_input= (90. - lat_adj)*PI/180;

    }
    rylm(colat_input,lon_input,order,ylmval);

    for (l = 0; l <= order; l++) {
      for (m = -l; m <= l; m++) {

	    k = l * (l+1) + m+1;
	    x += cint[k-1][0][flag]*ylmval[k-1];
            y += cint[k-1][1][flag]*ylmval[k-1];
            z += cint[k-1][2][flag]*ylmval[k-1];
	}
    }
    r = sqrt(x * x + y * y + z * z);
    if ((r< 0.9) || (r > 1.1)) return -32;
   
    z /= r;
    x /= r;
    y /= r;

    if (z > 1.) colat_temp=0;
    else if (z< -1.) colat_temp =PI;
    else colat_temp= acos(z);
  
    if ((fabs(x) < 1e-8) && (fabs(y) < 1e-8)) lon_temp =0;
    else lon_temp = atan2(y,x);

    lon_output = lon_temp;

    if (flag == 0) {

	lat_alt =90 - colat_temp*180/PI;

	altitude_to_cgm(height_in, lat_alt,&lat_adj);
	colat_output = (90. - lat_adj) * PI/180;

    } else colat_output = colat_temp;


    *lat_out = (double) (90 - colat_output*180/PI);
    *lon_out  = (double) (lon_output*180/PI);

    return 0;
} 




