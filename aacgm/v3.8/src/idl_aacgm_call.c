#if defined(WIN32)
	#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "export.h"
#include "idl_aacgm.h"
#include "mlt.h"
#include "convert_geo_coord.h"
#include "aacgm.h"
#include "convert_geo_vec.h"

/*globals here*/
static char statusBuffer[256] ;

/* function protos */
extern IDL_VPTR IDL_CDECL aacgm_mlt(int argc, IDL_VPTR argv[], char *argk);
extern IDL_VPTR IDL_CDECL aacgm_mlong(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL aacgm_help(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL aacgm_conv_coord(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL aacgm_load_coef(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL aacgm_set_path(int argc, IDL_VPTR argv[], char *argk);
extern void IDL_CDECL aacgm_conv_vec(int argc, IDL_VPTR argv[], char *argk);

static IDL_SYSFUN_DEF2 idl_aacgm_functions[] = {
	{aacgm_mlt, "AACGM_MLT", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
	{aacgm_mlong, "AACGM_MLONG", 3, 3, 0, 0},
};

static IDL_SYSFUN_DEF2 idl_aacgm_procedures[] = {
	{(IDL_FUN_RET) aacgm_help, "AACGM_HELP", 0, 0, 0, 0},
	{(IDL_FUN_RET) aacgm_conv_coord, "AACGM_CONV_COORD", 6, 6, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
	{(IDL_FUN_RET) aacgm_load_coef, "AACGM_LOAD_COEF", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
	{(IDL_FUN_RET) aacgm_set_path, "AACGM_SET_PATH", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
	{(IDL_FUN_RET) aacgm_conv_vec, "AACGM_CONV_VEC", 12, 12, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
};

/* Check which system we are one */
#if defined(WIN32) 
	#include <windows.h>
#endif

/* startup call when DLM is loaded */
int idl_aacgm_startup(void)
{
	char status_msg[256];
	char version[4];
	char year[5];

	sprintf(version,"%3.1f",IDL_AACGM_VERSION);
	sprintf(year,"%4d",IDL_AACGM_YEAR);

	strcpy(status_msg,"IDL_AACGM Version ");
	strncat(status_msg,version,3);
	strcat(status_msg,". DLM interface (c) ");
	strncat(status_msg,year,4);
	strcat(status_msg," by Haje Korth, JHU/APL.");

	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,status_msg);

    if (!IDL_SysRtnAdd(idl_aacgm_functions, TRUE, ARRLEN(idl_aacgm_functions))) {
                return IDL_FALSE;
        }
	
	if (!IDL_SysRtnAdd(idl_aacgm_procedures, FALSE, ARRLEN(idl_aacgm_procedures))) {
                return IDL_FALSE;
        }

        return(IDL_TRUE);
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

void IDL_CDECL aacgm_help(int argc, IDL_VPTR argv[], char *argk)
{
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"IDL AACGM Routines:");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_SET_PATH, path");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_LOAD_COEF, year");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_CONV_COORD, glat, glon, height, mlat, mlon, error, /TO_AACGM");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_CONV_COORD, mlat, mlon, height, glat, glon, error, /TO_GEO");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_MLT, year, t0, mlon");
	IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO,
		"AACGM_MLONG, year, t0, mlt");
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

IDL_VPTR IDL_CDECL aacgm_mlt(int argc, IDL_VPTR argv[], char *argk)
{
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		int mslong_set;
		IDL_VPTR mslong_out;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
		{"MSLONG",IDL_TYP_UNDEF,1,IDL_KW_OUT|IDL_KW_ZERO,IDL_KW_OFFSETOF(mslong_set),IDL_KW_OFFSETOF(mslong_out)},
		{NULL}
        };

	KW_RESULT kw;
	int tyr, tt0, i;
	double tmlong, tmlt_out, tmslong;
	int *pyr, *pt0;
	double *pmlong, *pmlt_out, *pmslong;
	IDL_VPTR yr, t0, mlong;
	IDL_VPTR mlt_out, mlt_out_scalar, mslong;
	IDL_MEMINT n_yr, n_t0, n_mlong;

	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,0,1,&kw);

	yr=IDL_BasicTypeConversion(1,&argv[0],IDL_TYP_LONG);
	IDL_VarGetData(yr, &n_yr,(char **) &pyr, FALSE);
	t0=IDL_BasicTypeConversion(1,&argv[1],IDL_TYP_LONG);
	IDL_VarGetData(t0, &n_t0,(char **) &pt0, FALSE);
	mlong=IDL_BasicTypeConversion(1,&argv[2],IDL_TYP_DOUBLE);
	IDL_VarGetData(mlong, &n_mlong,(char **) &pmlong, FALSE);

	if ((n_yr != n_t0) || (n_yr != n_mlong))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Array dimensions differ.");

	if (n_yr == 1) {
		pmlt_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_yr, IDL_ARR_INI_ZERO, &mlt_out);	
		pmslong=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_yr, IDL_ARR_INI_ZERO, &mslong);	
	} else {
		pmlt_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			yr->value.arr->n_dim, yr->value.arr->dim, IDL_ARR_INI_ZERO, &mlt_out);	
		pmslong=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			yr->value.arr->n_dim, yr->value.arr->dim, IDL_ARR_INI_ZERO, &mslong);	
	}

	for (i=0;i<n_yr;i++) {
		tyr=pyr[i];
		tt0=pt0[i];
		tmlong=pmlong[i];

		if ((tmlong < -180.0) || (tmlong > 360.0))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Valid longitude range is -180 < MLONG < 360.");

  		tmlt_out=mlt(tyr,tt0,tmlong,&tmslong);

		pmlt_out[i]=tmlt_out;
		pmslong[i]=tmslong;
	}

	if (kw.mslong_set) {
		if (n_yr == 1) {
			IDL_StoreScalar(kw.mslong_out, IDL_TYP_DOUBLE, (IDL_ALLTYPES *) pmslong);
		} else {
			IDL_VarCopy(mslong,kw.mslong_out);
		}
	}

	if (yr != argv[0]) IDL_Deltmp(yr);
	if (t0 != argv[1]) IDL_Deltmp(t0);
	if (mlong != argv[2]) IDL_Deltmp(mlong);
	IDL_Deltmp(mslong);

	IDL_KW_FREE;

	if (n_yr == 1) {
		mlt_out_scalar=IDL_Gettmp();
		mlt_out_scalar->type=IDL_TYP_DOUBLE;
		mlt_out_scalar->value.d=pmlt_out[0];
		IDL_Deltmp(mlt_out);
		return(mlt_out_scalar);
	} else {
		return(mlt_out);
	}
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

IDL_VPTR IDL_CDECL aacgm_mlong(int argc, IDL_VPTR argv[], char *argk)
{
	int tyr, tt0, i;
	double tmlt, tmlt0, tmslong, tmlong_out;
	int *pyr, *pt0;
	double *pmlt, *pmlong_out;
	IDL_VPTR yr, t0, mlta;
	IDL_VPTR mlong_out, mlong_out_scalar;
	IDL_MEMINT n_yr, n_t0, n_mlt;

	yr=IDL_BasicTypeConversion(1,&argv[0],IDL_TYP_LONG);
	IDL_VarGetData(yr, &n_yr,(char **) &pyr, FALSE);
	t0=IDL_BasicTypeConversion(1,&argv[1],IDL_TYP_LONG);
	IDL_VarGetData(t0, &n_t0,(char **) &pt0, FALSE);
	mlta=IDL_BasicTypeConversion(1,&argv[2],IDL_TYP_DOUBLE);
	IDL_VarGetData(mlta, &n_mlt,(char **) &pmlt, FALSE);

	if ((n_yr != n_t0) || (n_yr != n_mlt))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Array dimensions differ.");

	if (n_yr == 1) {
		pmlong_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_yr, IDL_ARR_INI_ZERO, &mlong_out);
	} else {
		pmlong_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			yr->value.arr->n_dim, yr->value.arr->dim, IDL_ARR_INI_ZERO, &mlong_out);
	}

	for (i=0;i<n_yr;i++) {
		tyr=pyr[i];
		tt0=pt0[i];
		tmlt=pmlt[i];

		if ((tmlt < 0.0) || (tmlt > 24.0))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Valid MLT range is 0 < MLT < 24.");

        tmlt0=mlt(tyr,tt0,0.0,&tmslong);
  		tmlong_out=(tmlt-tmlt0)*360.0/24.0;
		if (tmlong_out < 0) tmlong_out=tmlong_out+360.0;

		pmlong_out[i]=tmlong_out;
	}

	if (yr != argv[0]) IDL_Deltmp(yr);
	if (t0 != argv[1]) IDL_Deltmp(t0);
	if (mlta != argv[2]) IDL_Deltmp(mlta);

	if (n_yr == 1) {
		mlong_out_scalar=IDL_Gettmp();
		mlong_out_scalar->type=IDL_TYP_DOUBLE;
		mlong_out_scalar->value.d=pmlong_out[0];
		IDL_Deltmp(mlong_out);
		return(mlong_out_scalar);
	} else {
		return(mlong_out);
	}
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

void IDL_CDECL aacgm_conv_coord(int argc, IDL_VPTR argv[], char *argk)
{
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		int order;
		int to_aacgm;
		int to_geo;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
		{"ORDER",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(order)},
		{"TO_AACGM",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(to_aacgm)},
		{"TO_GEO",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(to_geo)},
		{NULL}
        };

	KW_RESULT kw;
	int flag,i,terr,order;
	int *perr;
	double tlat_in,tlon_in,theight_in;
	double tlat_out,tlon_out;
	double *plat_in,*plon_in,*pheight_in;
	double *plat_out,*plon_out;
	IDL_VPTR lat_in,lon_in,height_in;
	IDL_VPTR lat_out,lon_out,err;
	IDL_MEMINT n_lat_in,n_lon_in,n_height_in;

	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,0,1,&kw);

	lat_in=IDL_BasicTypeConversion(1,&argv[0],IDL_TYP_DOUBLE);
	IDL_VarGetData(lat_in, &n_lat_in,(char **) &plat_in, FALSE);
	lon_in=IDL_BasicTypeConversion(1,&argv[1],IDL_TYP_DOUBLE);
	IDL_VarGetData(lon_in, &n_lon_in,(char **) &plon_in, FALSE);
	height_in=IDL_BasicTypeConversion(1,&argv[2],IDL_TYP_DOUBLE);
	IDL_VarGetData(height_in, &n_height_in,(char **) &pheight_in, FALSE);

    IDL_StoreScalarZero(argv[3], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[4], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[5], IDL_TYP_LONG);

	if ((n_lat_in != n_lon_in) || (n_lat_in != n_height_in))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Array dimensions differ.");

	if (n_lat_in == 1) {
		plat_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &lat_out);	
		plon_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &lon_out);	
		perr=(int *) IDL_MakeTempVector(IDL_TYP_LONG,
			n_lat_in, IDL_ARR_INI_ZERO, &err);	
	} else {
		plat_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &lat_out);	
		plon_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &lon_out);	
		perr=(int *) IDL_MakeTempArray(IDL_TYP_LONG,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &err);	
	}

	if ((kw.to_aacgm == 0) && (kw.to_geo == 0))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Conversion direction not specified.");
	if ((kw.to_aacgm == 1) && (kw.to_geo == 1))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Only one conversion direction may be specified.");
	if (kw.to_aacgm) flag=0;
	if (kw.to_geo) flag=1;

    if (kw.order) order=kw.order;
		else order=10;

	for (i=0;i<n_lat_in;i++) {
		tlat_in=plat_in[i];
		tlon_in=plon_in[i];
		theight_in=pheight_in[i];

		if ((tlat_in < -90.0) || (tlat_in > 90.0) || 
			(tlon_in < -180.0) || (tlon_in > 360.0) ||
			(theight_in < 0.0))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Valid ranges are -90 < LAT < +90, -180 < LON < 360, Height > 0.");

		terr=convert_geo_coord(tlat_in,tlon_in,theight_in,&tlat_out,&tlon_out,flag,order);

		if (terr ==0) {
			plat_out[i]=tlat_out;
			plon_out[i]=tlon_out;
		} else{
			plat_out[i]=0;
			plon_out[i]=0;
		}
		perr[i]=terr;
	}

	if (n_lat_in == 1) {
		IDL_StoreScalar(argv[3], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) plat_out);
		IDL_StoreScalar(argv[4], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) plon_out);
		IDL_StoreScalar(argv[5], IDL_TYP_LONG, (IDL_ALLTYPES *) perr);

		IDL_Deltmp(lat_out);
		IDL_Deltmp(lon_out);
		IDL_Deltmp(err);
	} else {
		IDL_VarCopy(lat_out,argv[3]);
		IDL_VarCopy(lon_out,argv[4]);
		IDL_VarCopy(err,argv[5]);
	}

	if (lat_in != argv[0]) IDL_Deltmp(lat_in);
	if (lon_in != argv[1]) IDL_Deltmp(lon_in);
	if (height_in != argv[2]) IDL_Deltmp(height_in);

	IDL_KW_FREE;
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

void IDL_CDECL aacgm_load_coef(int argc, IDL_VPTR argv[], char *argk)
{
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		int quiet;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
		{"QUIET",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(quiet)},
		{NULL}
        };

	KW_RESULT kw;
	char *fpath;
	char fname[256];
	char yearstr[32];
	char status_msg[256];
	int year;
	int status;

	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,0,1,&kw);

    year=IDL_LongScalar(argv[0]);

    sprintf(yearstr,"%4.4d",year);  
	fpath=getenv("AACGM_PATH");
	if (fpath!=NULL) {
		strcpy(fname,getenv("AACGM_PATH"));  
		strcat(fname,"aacgm_coeffs");
		strcat(fname,yearstr);
		strcat(fname,".asc");
	} else {
		strcpy(fname,"aacgm_coeffs");
		strcat(fname,yearstr);
		strcat(fname,".asc");
	}

	status=AACGMLoadCoef(fname);

	if (status == 0) {
		if (kw.quiet == 0) {
			strcpy(status_msg,"Coefficients loaded for year ");
			strcat(status_msg,yearstr);
			strcat(status_msg,".");
			IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_INFO,
				status_msg);
		}
	} else {
		strcpy(status_msg,"Coefficient file ");
		strcat(status_msg,fname);
		strcat(status_msg," not found.");
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
			status_msg);
	}

	IDL_KW_FREE;
}

/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

void IDL_CDECL aacgm_set_path(int argc, IDL_VPTR argv[], char *argk)
{
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		int quiet;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
		{"QUIET",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(quiet)},
		{NULL}
        };

	KW_RESULT kw;
	char path[256];
	static char envvar[256];
	char status_msg[256];

	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,0,1,&kw);

	IDL_ENSURE_STRING(argv[0]);
	IDL_ENSURE_SCALAR(argv[0]);
	strcpy(path,IDL_STRING_STR(&argv[0]->value.str));
        
	strcpy(envvar,"AACGM_PATH=");
	strcat(envvar,path);
	putenv(envvar);

	if (kw.quiet == 0) {
		sprintf(status_msg,"AACGM_PATH changed to %s",path);
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_INFO,
			status_msg);
	}

	IDL_KW_FREE;
}
/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/

void IDL_CDECL aacgm_conv_vec(int argc, IDL_VPTR argv[], char *argk)
{
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		int to_aacgm;
		int to_geo;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
		{"TO_AACGM",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(to_aacgm)},
		{"TO_GEO",IDL_TYP_LONG,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(to_geo)},
		{NULL}
        };

	KW_RESULT kw;
	int i,terr;
	int *perr;
	double tlat_in,tlon_in,theight_in,tth_in,tph_in;
	double tlat_out,tlon_out,tth1_out,tph1_out,tth2_out,tph2_out;
	double *plat_in,*plon_in,*pheight_in,*pth_in,*pph_in;
	double *plat_out,*plon_out,*pth1_out,*pph1_out,*pth2_out,*pph2_out;
	IDL_VPTR lat_in,lon_in,height_in,th_in,ph_in;
	IDL_VPTR lat_out,lon_out,th1_out,ph1_out,th2_out,ph2_out,err;
	IDL_MEMINT n_lat_in,n_lon_in,n_height_in,n_th_in,n_ph_in;

	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,0,1,&kw);

	lat_in=IDL_BasicTypeConversion(1,&argv[0],IDL_TYP_DOUBLE);
	IDL_VarGetData(lat_in, &n_lat_in,(char **) &plat_in, FALSE);
	lon_in=IDL_BasicTypeConversion(1,&argv[1],IDL_TYP_DOUBLE);
	IDL_VarGetData(lon_in, &n_lon_in,(char **) &plon_in, FALSE);
	height_in=IDL_BasicTypeConversion(1,&argv[2],IDL_TYP_DOUBLE);
	IDL_VarGetData(height_in, &n_height_in,(char **) &pheight_in, FALSE);
	th_in=IDL_BasicTypeConversion(1,&argv[3],IDL_TYP_DOUBLE);
	IDL_VarGetData(th_in, &n_th_in,(char **) &pth_in, FALSE);
	ph_in=IDL_BasicTypeConversion(1,&argv[4],IDL_TYP_DOUBLE);
	IDL_VarGetData(ph_in, &n_ph_in,(char **) &pph_in, FALSE);

    IDL_StoreScalarZero(argv[5], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[6], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[7], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[8], IDL_TYP_DOUBLE);
    IDL_StoreScalarZero(argv[9], IDL_TYP_LONG);

	if ((n_lat_in != n_lon_in) || (n_lat_in != n_height_in))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Array dimensions differ.");

	if (n_lat_in == 1) {
		plat_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &lat_out);	
		plon_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &lon_out);	
		pth1_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &th1_out);	
		pph1_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &ph1_out);	
		pth2_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &th2_out);	
		pph2_out=(double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
			n_lat_in, IDL_ARR_INI_ZERO, &ph2_out);	
		perr=(int *) IDL_MakeTempVector(IDL_TYP_LONG,
			n_lat_in, IDL_ARR_INI_ZERO, &err);	
	} else {
		plat_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &lat_out);	
		plon_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &lon_out);	
		pth1_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &th1_out);	
		pph1_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &ph1_out);	
		pth2_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &th2_out);	
		pph2_out=(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &ph2_out);	
		perr=(int *) IDL_MakeTempArray(IDL_TYP_LONG,
			lat_in->value.arr->n_dim, lat_in->value.arr->dim, IDL_ARR_INI_ZERO, &err);	
	}

	if ((kw.to_aacgm == 0) && (kw.to_geo == 0))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Conversion direction not specified.");
	if ((kw.to_aacgm == 1) && (kw.to_geo == 1))
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Only one conversion direction may be specified.");

	if (kw.to_aacgm) {
		for (i=0;i<n_lat_in;i++) {
			tlat_in=plat_in[i];
			tlon_in=plon_in[i];
			theight_in=pheight_in[i];
			tth_in=pth_in[i];
			tph_in=pph_in[i];

			if ((tlat_in < -90.0) || (tlat_in > 90.0) || 
				(tlon_in < -180.0) || (tlon_in > 360.0) ||
				(theight_in < 0.0))
			IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
			"Valid ranges are -90 < LAT < +90, -180 < LON < 360, Height > 0.");

			terr=convert_geo_vec(tlat_in,tlon_in,theight_in,tth_in,tph_in,
				&tlat_out,&tlon_out,&tth1_out,&tph1_out,&tth2_out,&tph2_out);

			if (terr ==0) {
				plat_out[i]=tlat_out;
				plon_out[i]=tlon_out;
				pth1_out[i]=tth1_out;
				pph1_out[i]=tph1_out;
				pth2_out[i]=tth2_out;
				pph2_out[i]=tph2_out;
			} else {
				plat_out[i]=0;
				plon_out[i]=0;
				pth1_out[i]=0;
				pph1_out[i]=0;
				pth2_out[i]=0;
				pph2_out[i]=0;
			}
			perr[i]=terr;
		}
	}


	if (kw.to_geo) {
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,
		"Vector conversions from AACGM to GEO are not supported.");
	}

	if (n_lat_in == 1) {
		IDL_StoreScalar(argv[5], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) plat_out);
		IDL_StoreScalar(argv[6], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) plon_out);
		IDL_StoreScalar(argv[7], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) pth1_out);
		IDL_StoreScalar(argv[8], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) pph1_out);
		IDL_StoreScalar(argv[9], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) pth2_out);
		IDL_StoreScalar(argv[10], IDL_TYP_DOUBLE, (IDL_ALLTYPES *) pph2_out);
		IDL_StoreScalar(argv[11], IDL_TYP_LONG, (IDL_ALLTYPES *) perr);

		IDL_Deltmp(lat_out);
		IDL_Deltmp(lon_out);
		IDL_Deltmp(th1_out);
		IDL_Deltmp(ph1_out);
		IDL_Deltmp(th2_out);
		IDL_Deltmp(ph2_out);
		IDL_Deltmp(err);
	} else {
		IDL_VarCopy(lat_out,argv[5]);
		IDL_VarCopy(lon_out,argv[6]);
		IDL_VarCopy(th1_out,argv[7]);
		IDL_VarCopy(ph1_out,argv[8]);
		IDL_VarCopy(th2_out,argv[9]);
		IDL_VarCopy(ph2_out,argv[10]);
		IDL_VarCopy(err,argv[11]);
	}

	if (lat_in != argv[0]) IDL_Deltmp(lat_in);
	if (lon_in != argv[1]) IDL_Deltmp(lon_in);
	if (height_in != argv[2]) IDL_Deltmp(height_in);
	if (th_in != argv[3]) IDL_Deltmp(th_in);
	if (ph_in != argv[4]) IDL_Deltmp(ph_in);

	IDL_KW_FREE;
}
/*{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|*/
