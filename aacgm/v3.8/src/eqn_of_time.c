/* eqn_of_time.c
   =============
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
 $Log: eqn_of_time.c,v $
 Revision 1.4  2001/06/27 20:24:39  barnes
 Added license tag.

 Revision 1.3  2001/01/29 18:11:47  barnes
 Added Author Name

 Revision 1.2  1999/12/15 18:57:37  barnes
 Modifications for year 2000.

 Revision 1.1  1998/05/12 14:04:48  barnes
 Initial revision

*/

#include <math.h>
#include "math.h"

double eqn_of_time(double mean_lon,int yr) {

  int index;
  double coef[12][7]={
     {-105.8,596.2,4.4,-12.7,-429.0,-2.1,19.3},
     {-105.9,596.2,4.4,-12.7,-429.0,-2.1,19.3},
     {-106.1,596.2,4.4,-12.7,-428.9,-2.1,19.3},
     {-106.2,596.2,4.4,-12.7,-428.9,-2.1,19.3},
     {-106.4,596.1,4.4,-12.7,-428.9,-2.1,19.3},
     {-106.5,596.1,4.4,-12.7,-428.8,-2.1,19.3},
     {-106.6,596.1,4.4,-12.7,-428.8,-2.1,19.3},
     {-106.7,596.1,4.4,-12.7,-428.7,-2.1,19.3},
     {-106.8,596.1,4.4,-12.7,-428.7,-2.1,19.3},
     {-107.0,596.1,4.4,-12.7,-428.7,-2.1,19.3},
     {-107.2,596.1,4.4,-12.7,-428.6,-2.1,19.3},
     {-107.3,596.1,4.4,-12.7,-428.6,-2.1,19.3},
  };

  if (yr<88) index = yr + 2000 - 1988;
  if ((yr>=88) && (yr<100)) index = yr - 88;
  else if ((yr>=100) && (yr<1900)) index=yr-88;
  else index=yr-1988;

   
   if (index<0) index = 1;
   if (index>12) index = 12;

    return coef[index-1][0]*sind(mean_lon)+
           coef[index-1][1]*sind(2.0*mean_lon)+
           coef[index-1][2]*sind(3.0*mean_lon)+
           coef[index-1][3]*sind(4.0*mean_lon)+
           coef[index-1][4]*cosd(mean_lon)+
           coef[index-1][5]*cosd(2.0*mean_lon)+
           coef[index-1][6]*cosd(3.0*mean_lon);
}






