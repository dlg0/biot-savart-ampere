pro iridium_to_magnetic, vecth_irid, vecph_irid, clati, loni, xavg, yavg, vecth_mag, vecph_mag, $
	CLATM=clatm, LONM=lonm, VECX=ct_x, VECY=ct_y, VECZ=ct_z

;	All inputs in radians

	npts	= n_elements(vecth_irid[*])

;	Shift vectors to Cartesian
	
	ct_x	= fltarr(npts)
	ct_y	= fltarr(npts)
	ct_z	= fltarr(npts)

	sinCLati	= sin ( cLati );
	cosCLati	= cos ( cLati );
	cosLoni	= cos ( loni );
	sinLoni	= sin ( loni );
	sinCLati_cosLoni	= sinCLati * cosLoni;
	sinCLati_sinLoni	= sinCLati * sinLoni;
	cosCLati_cosLoni	= cosCLati * cosLoni;
	cosCLati_sinLoni	= cosCLati * sinLoni;

;	rot_to_cartF	= fltArr ( 3, 3, n_elements ( cLati ) );
;	rot_to_cartF[0,0,*]	= sinCLati_cosLoni;
;	rot_to_cartF[1,0,*]	= cosCLati_cosLoni;
;	rot_to_cartF[2,0,*]	= -sinLoni;
;	rot_to_cartF[0,1,*]	= sinCLati_sinLoni;
;	rot_to_cartF[1,1,*]	= cosCLati_sinLoni;
;	rot_to_cartF[2,1,*]	= cosLoni;
;	rot_to_cartF[0,2,*]	= cosCLati;
;	rot_to_cartF[1,2,*]	= -sinCLati;
;	rot_to_cartF[2,2,*]	= 0.0;
;
;	rot_to_cartF	= reform ( rot_to_cartF, 3, 3 * n_elements ( cLati ) );
;
;	vec_iridF	= [ [ fltArr ( n_elements ( cLati ) ) ], $
;								[ vecTh_irid[*] ], $
;								[ vecPh_irid[*] ] ];
;
;	vec_cartF	= rot_to_cartF ## vec_iridF;
;
;	iiF	= indgen ( 3, n_elements ( cLati ) ) * n_elements ( cLati ) + $
;		rebin ( transpose ( indgen ( n_elements ( cLati ) ) ), 3, n_elements ( cLati ) );
;	
;	ctF	= transpose ( vec_cartF[iiF] );
		
	for jj=0,npts-1 do begin
; spherical to Cartesian coords
		rot_to_cart	= [	[	sinCLati_cosLoni[jj],	cosCLati_cosLoni[jj],	-sinLoni[jj]],$
										[	sinCLati_sinLoni[jj],	cosCLati_sinLoni[jj], cosLoni[jj]	], $
										[	cosCLati[jj],								-sinCLati[jj],		0				]]

		vec_cart	= rot_to_cart ## [[0.], [vecth_irid[jj]], [vecph_irid[jj]]]
; dB components in x,y,z coords
		ct_x[jj]	= vec_cart[0]
		ct_y[jj]	= vec_cart[1]
		ct_z[jj]	= vec_cart[2]

	endfor

;	Calculate coords in Magnetic and add pole shift
	Xcrd	= clati*!radeg*cos(loni-!pi) + Xavg[0]
	Ycrd	= clati*!radeg*sin(loni-!pi) + Yavg[0]

	clatm	= sqrt(Xcrd^2+Ycrd^2)*!dtor ; [rads]
	lonm	= (180. + atan(Ycrd,Xcrd)*!radeg)*!dtor ; [rads]

;	Shift from cartesian to magnetic
	vecth_mag		= fltarr(npts)
	vecph_mag		= fltarr(npts)

	sinCLatm	= sin ( cLatm );
	cosLonm	= cos ( lonm );
	cosCLatm	= cos ( cLatm );
	sinLonm	= sin ( lonm );
	sinCLatm_cosLonm	= sinCLatm * cosLonm;
	cosCLatm_cosLonm	= cosCLatm * cosLonm;
	sinCLatm_sinLonm	= sinCLatm * sinLonm;
	cosCLatm_sinLonm	= cosCLatm * sinLonm;

;	rot_to_sphF	= fltArr ( 3, 3, n_elements ( cLati ) );
;	rot_to_sphF[0,0,*]	= sinCLatm_cosLonm;
;	rot_to_sphF[1,0,*]	= sinCLatm_sinLonm;
;	rot_to_sphF[2,0,*]	= cosCLatm;
;	rot_to_sphF[0,1,*]	= cosCLatm_cosLonm;
;	rot_to_sphF[1,1,*]	= cosCLatm_sinLonm;
;	rot_to_sphF[2,1,*]	= -sinCLatm;
;	rot_to_sphF[0,2,*]	= -sinLonm;
;	rot_to_sphF[1,2,*]	= cosLonm;
;	rot_to_sphF[2,2,*]	= 0.0;
;
;	rot_to_sphF	= reform ( rot_to_sphF, 3, 3 * n_elements ( cLati ) );
;
;	vec_magF	= rot_to_sphF ## ctF;
;
;	vec_magF	= transpose ( vec_magF[iiF] );
	
	for jj=0,npts-1 do begin
	
		rot_to_sph	= [	[	sinCLatm_cosLonm[jj],	sinCLatm_sinLonm[jj],	cosCLatm[jj]	],$
										[	cosCLatm_cosLonm[jj],	cosCLatm_sinLonm[jj],	-sinCLatm[jj]	], $
										[	-sinLonm[jj],							cosLonm[jj],		0								]]
	
		vec_mag	= rot_to_sph ## [[ct_x[jj]], [ct_y[jj]], [ct_z[jj]]]
	
		vecth_mag[jj]	= vec_mag[1]
		vecph_mag[jj]	= vec_mag[2]
	
	endfor

;vecTh_mag	= vec_magF[*,1];
;vecPh_mag	= vec_magF[*,2];

end
