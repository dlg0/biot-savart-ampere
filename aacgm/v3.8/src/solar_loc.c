/* solar_loc.c
   ===========
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
 $Log: solar_loc.c,v $
 Revision 1.5  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.4  2001/01/29 18:11:47  barnes
 Added Author Name

 Revision 1.3  1999/12/15 18:57:37  barnes
 Modifications for year 2000.

 Revision 1.2  1998/06/05 21:27:38  barnes
 Fixed problem of no return value.

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

 */

#include <stdlib.h>
#include <math.h>
#include "math.h"

int solar_loc(int yr,int t1,double *mean_lon,double *dec) {

  int index,delta_yr,yr_step,i;
  double d,lambda,g,eps,L;
  double L0[12]={279.642,279.403,279.165,278.926,279.673,279.434,
		 279.196,278.957,279.704,279.465,279.226,278.982};
  double DL=0.985647;
  double G0[12]={356.892984,356.637087,356.381191,356.125295,
                 356.854999,356.599102,356.343206,356.087308,
		 356.817011,356.561113,356.31,356.05};
  double DG=0.98560028;
  double EPS0[12]={23.440722,23.440592,23.440462,23.440332,
	           23.440202,23.440072,23.439942,23.439811,
		   23.439682,23.439552,23.439422,23.439292};
  double DE=-0.00000036;
      
  d = 0;
  if (yr<1900) index = yr - 88;
  else index = yr - 1988;

  if (index<=0) delta_yr = index - 1;
  else if (index>10) delta_yr = index - 10;
  else  delta_yr = 0;
     
  if (index<=0) index = 1;
  if (index>12) index = 12;

  yr_step = sgn(1,delta_yr);
  delta_yr = fabs(delta_yr);

  for (i=1;i<=delta_yr;i++) {
    if (yr_step>0) yr=98+i-1;
    else yr=89-i;

    if (mod(yr,4)==0) d = d + 366*yr_step;
    else d = d + 365*yr_step;
  }

     
  d = d + t1/86400;
  L = L0[index-1] + DL*d;
  g = G0[index-1] + DG*d;

  while (L<0) L = L + 360;
  while (g<0) g = g + 360;
    
  L = mod(L,360.0);
  g = mod(g,360.0);

  lambda = L + 1.915*sind(g) + 0.020*sind(2*g);
  eps = EPS0[index-1] + DE*d;

  *dec = asind(sind(eps)*sind(lambda));
  *mean_lon = L;
  return 0;
}





