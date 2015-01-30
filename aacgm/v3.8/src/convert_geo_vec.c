#include <stdio.h>
#include <math.h>
#include "convert_geo_coord.h"


void sphcar(double r, double theta, double phi, double *x, double *y, double *z)
{
	double d2r,sq;

	d2r=3.14159265/180.0;
    theta=theta*d2r;
	phi=phi*d2r;

	sq=r*sin(theta);
	*x=sq*cos(phi);
	*y=sq*sin(phi);
	*z=r*cos(theta);
}


void carsph(double x, double y, double z, double *r, double *theta, double *phi)
{
	double d2r,sq;

	d2r=3.14159265/180.0;

	sq=x*x+y*y;
	*r=sqrt(sq+z*z);
	if (sq != 0) {
		sq=sqrt(sq);
		*phi=atan2(y,x);
		*theta=atan2(sq,z);
		if (*phi < 0) *phi=*phi+6.28318531;
	} else {
		*phi=0;
		if (z < 0) *theta=3.141592654;
			else *theta=0;
	}

	*theta=*theta/d2r;
	*phi=*phi/d2r;
}


void bspcar(double theta, double phi, double br, double btheta, double bphi, 
			double *bx, double *by, double *bz)
{
	double d2r,s,c,sf,cf,be;

	d2r=3.14159265/180.0;
    theta=theta*d2r;
	phi=phi*d2r;

	s=sin(theta);
	c=cos(theta);
	sf=sin(phi);
	cf=cos(phi);
	be=br*s+btheta*c;
	*bx=be*cf-bphi*sf;
	*by=be*sf+bphi*cf;
	*bz=br*c-btheta*s;
}


void bcarsp(double x, double y, double z, double bx, double by, double bz, 
			double *br, double *btheta, double *bphi)
{
	double rho2,r,rho;	
	double cphi,sphi;
	double ct,st;

	rho2=x*x+y*y;
	r=sqrt(rho2+z*z);
	rho=sqrt(rho2);

	if (rho != 0) {
		cphi=x/rho;
		sphi=y/rho;
	} else {
		cphi=1.;
		sphi=0.;
	}

	ct=z/r;
	st=rho/r;

	*br=(x*bx+y*by+z*bz)/r;
	*btheta=(bx*cphi+by*sphi)*ct-bz*st;
	*bphi=by*cphi-bx*sphi;
}


void crossp(double v11, double v12, double v13, double v21, double v22, double v23, 
			double *v31, double *v32, double *v33)
{
	*v31=v12*v23-v13*v22;
	*v32=v13*v21-v11*v23;
	*v33=v11*v22-v12*v21;
}


int convert_geo_vec(double lat_in,double lon_in, double height_in, double th_in, double ph_in,
	double *lat_out, double *lon_out, double *th1_out, double *ph1_out, double *th2_out, double *ph2_out)
{
	int order,err;
	double dgth,dgph;
    double cglat,glon;
    double mlat,mlon,cmlat;
	double xg,yg,zg;
	double xm,ym,zm;
	double vr,vth,vph;
	double vx,vy,vz;
    double xg_th,yg_th,zg_th;
	double rg_th,cglat_th,glon_th,glat_th;
	double mlat_th,mlon_th,cmlat_th;
	double xm_th,ym_th,zm_th;
	double xg_ph,yg_ph,zg_ph;
	double rg_ph,cglat_ph,glon_ph,glat_ph;
	double mlat_ph,mlon_ph,cmlat_ph;
	double xm_ph,ym_ph,zm_ph;
	double mxyz_mg;
	double mxyz_r0,mxyz_r1,mxyz_r2;
	double mxyz_ruv0,mxyz_ruv1,mxyz_ruv2;
	double mxyz_th0,mxyz_th1,mxyz_th2;
	double mxyz_thuv0,mxyz_thuv1,mxyz_thuv2;
	double mxyz_ph0,mxyz_ph1,mxyz_ph2;
	double mxyz_phuv0,mxyz_phuv1,mxyz_phuv2;
	double mxyz_ph_gth0,mxyz_ph_gth1,mxyz_ph_gth2;
	double mrv_th,mthv_th,mphv_th;
	double mrv_ph,mthv_ph,mphv_ph;
	double mth_vec_gth,mph_vec_gth,gvec_th,gvec_ph;
	double mxyz_th_gph0,mxyz_th_gph1,mxyz_th_gph2;
	double mth_vec_gph,mph_vec_gph;

	order=10;
	dgth=1.0;
	dgph=1.0;

	cglat=90.0-lat_in;
	glon=lon_in;
	gvec_th=th_in;
	gvec_ph=ph_in;

	sphcar(height_in,cglat,glon,&xg,&yg,&zg);
    err=convert_geo_coord(lat_in,lon_in,height_in,&mlat,&mlon,0,order);
	if (mlon < 0) mlon=mlon+360.0;
	cmlat=90.0-mlat;
    sphcar(height_in,cmlat,mlon,&xm,&ym,&zm);

//  del_th in geo
	vr=0.0;          // construct a (0,theta,0) vector
	vth=dgth;
	vph=0.0;
	bspcar(cglat,glon,vr,vth,vph,&vx,&vy,&vz);
    xg_th=xg+vx;     // get geo(x,y,z) for geo(0,dth,0) (unit vec) shift
	yg_th=yg+vy;
	zg_th=zg+vz;
	carsph(xg_th,yg_th,zg_th,&rg_th,&cglat_th,&glon_th);
	glat_th=90.0-cglat_th;
    err=convert_geo_coord(glat_th,glon_th,height_in,&mlat_th,&mlon_th,0,order); // convert to aacgm mlat,mlon coords (from geo+dth)
	if (mlon_th < 0) mlon_th=mlon_th+360.0;
	cmlat_th=90.0-mlat_th;
    sphcar(height_in,cmlat_th,mlon_th,&xm_th,&ym_th,&zm_th); //	get xyz coords of aacgm+dth(geo) coords

//  del_ph in geo
    vr=0.0;          // construct a (0,0,phi) vector
	vth=0.0;
	vph=dgph;
    bspcar(cglat,glon,vr,vth,vph,&vx,&vy,&vz); // get this vector in xyz coords
    xg_ph=xg+vx;     // get geo(x,y,z) for geo(0,0,dph) (unit vector) shift
	yg_ph=yg+vy;
	zg_ph=zg+vz;
	carsph(xg_ph,yg_ph,zg_ph,&rg_ph,&cglat_ph,&glon_ph);
    glat_ph=90.0-cglat_ph;
//  convert to aacgm mlat,mlon coords (from geo+dph)
    err=convert_geo_coord(glat_ph,glon_ph,height_in,&mlat_ph,&mlon_ph,0,order);
	if (mlon_ph < 0) mlon_ph=mlon_ph+360.0;
    cmlat_ph=90.0-mlat_ph;
    sphcar(height_in,cmlat_ph,mlon_ph,&xm_ph,&ym_ph,&zm_ph); // get xyz coords of aacgm+dph(geo) coords

//	get aacgm radial unit vector
    mxyz_r0=xm;
	mxyz_r1=ym;
	mxyz_r2=zm;
    mxyz_mg=sqrt(mxyz_r0*mxyz_r0+mxyz_r1*mxyz_r1+mxyz_r2*mxyz_r2);
    mxyz_ruv0=mxyz_r0/mxyz_mg;
    mxyz_ruv1=mxyz_r1/mxyz_mg;
    mxyz_ruv2=mxyz_r2/mxyz_mg;

//  get aacgm(x,y,z) unit vec for a geo dth shift
    mxyz_th0=xm_th-xm;
	mxyz_th1=ym_th-ym;
	mxyz_th2=zm_th-zm;
    mxyz_mg=sqrt(mxyz_th0*mxyz_th0+mxyz_th1*mxyz_th1+mxyz_th2*mxyz_th2);
    mxyz_thuv0=mxyz_th0/mxyz_mg;
    mxyz_thuv1=mxyz_th1/mxyz_mg;
    mxyz_thuv2=mxyz_th2/mxyz_mg;

//  get aacgm(x,y,z) unit vec for a geo dph shift
    mxyz_ph0=xm_ph-xm;
	mxyz_ph1=ym_ph-ym;
	mxyz_ph2=zm_ph-zm;
    mxyz_mg=sqrt(mxyz_ph0*mxyz_ph0+mxyz_ph1*mxyz_ph1+mxyz_ph2*mxyz_ph2);
    mxyz_phuv0=mxyz_ph0/mxyz_mg;
    mxyz_phuv1=mxyz_ph1/mxyz_mg;
    mxyz_phuv2=mxyz_ph2/mxyz_mg;

    crossp(mxyz_ruv0,mxyz_ruv1,mxyz_ruv2,mxyz_thuv0,mxyz_thuv1,mxyz_thuv2,&mxyz_ph_gth0,&mxyz_ph_gth1,&mxyz_ph_gth2); // get dph unit vector for r x dth
//  for a geo dth shift -> aacgm_dth
    bcarsp(xm,ym,zm,mxyz_thuv0,mxyz_thuv1,mxyz_thuv2,&mrv_th,&mthv_th,&mphv_th);
//  for a geo dth shift -> aacgm_dph
    bcarsp(xm,ym,zm,mxyz_ph_gth0,mxyz_ph_gth1,mxyz_ph_gth2,&mrv_ph,&mthv_ph,&mphv_ph);
    mth_vec_gth=gvec_th*mthv_th+gvec_ph*mthv_ph;
    mph_vec_gth=gvec_th*mphv_th+gvec_ph*mphv_ph;

    crossp(mxyz_phuv0,mxyz_phuv1,mxyz_phuv2,mxyz_ruv0,mxyz_ruv1,mxyz_ruv2,&mxyz_th_gph0,&mxyz_th_gph1,&mxyz_th_gph2); // get dth unit vector for dph x r
//	for a geo dph shift -> aacgm_dth
    bcarsp(xm,ym,zm,mxyz_th_gph0,mxyz_th_gph1,mxyz_th_gph2,&mrv_th,&mthv_th,&mphv_th);
//  for a geo dph shift -> aacgm_dph
    bcarsp(xm,ym,zm,mxyz_phuv0,mxyz_phuv1,mxyz_phuv2,&mrv_ph,&mthv_ph,&mphv_ph);
    mth_vec_gph=gvec_th*mthv_th+gvec_ph*mthv_ph;
    mph_vec_gph=gvec_th*mphv_th+gvec_ph*mphv_ph;

    *lat_out=mlat;
    *lon_out=mlon;
    *th1_out=mth_vec_gth;
    *ph1_out=mph_vec_gth;
    *th2_out=mth_vec_gph;
    *ph2_out=mph_vec_gph;

	return 0;
}