;+
;
; Calculate jPar from deltaB using SECS
; Juusola, L., et. al., Earth Planets Space, 58, 667-678, 2006
;
; [coLat] degrees
; [mlt] hours
;
;-

pro jPar_SECS, coLat_data, mlt_data,$
	bMag_data, coLat_grid, mlt_grid, $
	coLat_grid_mag, mlt_grid_mag, xAvg, yAvg, $
	NOPLOT = noPlot, JPARSECS = jPar, $
	POLELIMIT1 = poleLimit1, $
	POLELIMIT2 = poleLimit2, $
	NXY = nXY, $
	BTHSECS = bThFit_mag, $
	BPHSECS = bPhFit_mag, $
	UNITTH = unitTH, UNITPH = unitPh, $
	WEIGHT = weight, HARD = hard, $
	BPHFIT_DATA = bPhFit_data

  forward_function diffbydave
  if not keyword_set ( poleLimit1 ) then poleLimit1 = 35.0;
  if not keyword_set ( poleLimit2 ) then poleLimit2 = 20.0;
  if not keyword_set ( nXY ) then nXY = 11;

;	Set a hard boundary by adding extra zeroed data points
  if keyword_set ( hard ) then begin
    zeroCoLat	= [ fltArr ( 24 ) + hard, fltArr ( 24 ) + hard-3, fltArr ( 24 ) + hard, fltArr ( 24 ) + hard-3 ]
    zeroMlt	= [ fIndGen ( 24 ), fIndGen ( 24 )+0.5, fIndGen ( 24 ), fIndGen ( 24 ) + 0.5]

    coLat_data_Backup	= coLat_data
    mlt_data_Backup = mlt_data
    bMag_data_Backup	= bMag_data

    coLat_data	= [ coLat_data, zeroCoLat ]
    mlt_data	= [ mlt_data, zeroMlt ]
    bMag_data	= [ bMag_data, zeroMlt * 0 ]
  endif

  if keyword_set ( weight ) then begin
    unitTh_Backup = unitTh
    unitPh_Backup = unitPh
    weight_Backup = weight

    if keyword_set ( hard ) then begin
      unitTh	= [ unitTh, zeroMlt * 0 ]
      unitPh	= [ unitPh, zeroMlt * 0 + 1.0 ]
      weight	= [ weight[*], zeroMlt * 0 + 0.1 ]
    endif
  endif
;	Hard boundary done.

  u0	= 4d0 * !dpi * 10d0^(-7)
  rI	= 6356750.0 + 110d3; [m]
  r	= 6356750.0 + 800d3; [m] ( > rI )

;coLat_grid	= transpose ( rebin ( findgen ( 20 ) * 2 + 2.0, 20, 24 ) );
;mlt_grid	= rebin ( findgen ( 24 ), 24, 20 );

  maxLat = max ( coLat_data )

;	Generate SECS pole locations

  maxX	=  50 * sin ( 90.0 * !dtor - !pi )
;  nXY	= 11;21
  xPoleArray	= rebin ( ( findgen ( nXY ) - ( nXY - 1.0 ) / 2.0 ) / nXY * 2.0 * maxX, nXY, nXY )
  yPoleArray	= transpose ( xPoleArray )

  xPoleArray2	= xPoleArray + abs ( xPoleArray[1,0] - xPoleArray[0,0] ) / 2.0
  yPoleArray2	= yPoleArray + abs ( yPoleArray[0,1] - yPoleArray[0,0] ) / 2.0

  coLat_pole	= sqrt ( xPoleArray^2 + yPoleArray^2 ); [deg]
  mlt_pole	= ( 180.0 + atan ( yPoleArray, xPoleArray ) * !radeg ) / 15.0; [deg]

  coLat_pole2	= sqrt ( xPoleArray2^2 + yPoleArray2^2 ); [deg]
  mlt_pole2	= ( 180.0 + atan ( yPoleArray2, xPoleArray2 ) * !radeg ) / 15.0; [deg]

  iiPoleKeep	= where ( coLat_pole lt poleLimit1 )
  coLat_pole	= coLat_pole[iiPoleKeep]
  mlt_pole	= mlt_pole[iiPoleKeep]

  iiPoleKeep	= where ( coLat_pole2 lt poleLimit2 )
  coLat_pole2	= coLat_pole2[iiPoleKeep]
  mlt_pole2	= mlt_pole2[iiPoleKeep]

  coLat_pole	= [ coLat_pole, coLat_pole2 ]
  mlt_pole	= [ mlt_pole, mlt_pole2 ]

;	These pole locations are in magnetic coords, convert
;	to iridium coords to pass to jPar_SECS

  poleX	= coLat_pole * cos ( mlt_pole * 15 * !dtor - !pi ) - xAvg
  poleY	= coLat_pole * sin ( mlt_pole * 15 * !dtor - !pi ) - yAvg

  coLat_pole	= sqrt ( poleX^2 + poleY^2 ); [degrees]
  mlt_pole	= ( !pi + atan ( poleY, poleX ) ) * !radeg / 15.0; [degrees]

  if keyword_set ( weight ) then begin
    weight = rebin ( weight[*], n_elements ( weight ), n_elements ( coLat_pole[*] ) )
  endif

  bfnArray_data	= fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_data[*] ) )
  bfnArray_jPar	= fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_grid[*] ) )
  bfnArrayTh_grid	= fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_grid[*] ) )
  bfnArrayPh_grid	= fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_grid[*] ) )
  bfnArrayJcfTh_grid	=fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_grid[*] ) )
  bfnArrayJcfPh_grid	= fltarr ( n_elements ( coLat_pole[*] ), n_elements ( coLat_grid[*] ) )

  xShift	= coLat_pole * cos ( mlt_pole * 15.0 * !dtor - !pi ) ; These used in following loop to adj SEC pole position
  yShift	= coLat_pole * sin ( mlt_pole * 15.0 * !dtor - !pi )

; catch for where theta'=theta

  closestLat	= 0.1
  max_bph	= -u0 / ( 4.0 * !pi * r * tan ( closestLat * !dtor / 2.0 ) )
  max_jPar	= -1.0 / ( 8.0 * !pi * r^2 * sin ( closestLat * !dtor )^2 )
  max_jcf	= 1.0 / ( 4.0 * !pi * rI ) * 1.0 / tan ( closestLat * !dtor / 2.0 )

	bFnTime = 0.0

; create one cf basis fn and copy to other locations
  For ii = 0, n_elements ( coLat_pole[*] ) - 1 do begin      ; Loop over SECS poles
    xCoords	= coLat_data * cos ( mlt_data * 15.0 * !dtor - !pi ) - xShift[ii]  ; colat_data[np], mlt_data[np]
    yCoords	= coLat_data * sin ( mlt_data * 15.0 * !dtor - !pi ) - yShift[ii]

    thisCoLat	= sqrt ( xCoords^2 + yCoords^2 ) * !dtor
    thisLon	= ( 180. + atan ( yCoords , xCoords ) * !radeg ) * !dtor

    bfn_template_ph	= -u0 / ( 4.0 * !pi * r * sin ( thisCoLat>0.01*!dtor ) ) * $
      (( r - rI * cos ( thisCoLat ) ) / $
      sqrt ( r^2 - 2.0 * r * rI * cos ( thisCoLat ) + rI^2 ) - 1.0 )

	tic0 = SysTime(1)
    iridium_to_magnetic, bfn_template_ph * 0, bfn_template_ph, thisCoLat, thisLon, $
      xShift[ii], yShift[ii], $
      vThNew, vPhNew, $
      CLATM=coLatNew_rad, LONM=lonNew_rad, $
      VECX=ct_x, VECY=ct_y, VECZ=ct_z
	toc0 = SysTime(1)
	bFnTime = bFnTime + (toc0-tic0)

    if NOT keyword_set ( unitTh ) then begin
      bfnArray_data[ii,*]	= vPhNew   ; Assumed Iridium system for now

      if keyword_set ( hard ) then $
        bfnArray_data[ii,n_elements ( coLat_data[*] ) - 48 - 1 : *] = vThNew[n_elements ( coLat_data[*] ) - 48 - 1 : *]

    endif else begin
      bfnArray_data[ii,*]	= vPhNew * unitPh + vThNew * unitTh
    endelse

    xCoords_grid	= coLat_grid * cos ( mlt_grid * 15.0 * !dtor - !pi ) - xShift[ii] ; coLat_grid[24,20] ??
    yCoords_grid	= coLat_grid * sin ( mlt_grid * 15.0 * !dtor - !pi ) - yShift[ii]

    thisCoLat_grid	= sqrt ( xCoords_grid^2 + yCoords_grid^2 ) * !dtor > 0.001;
    thisLon_grid	= ( 180. + atan ( yCoords_grid , xCoords_grid ) * !radeg ) * !dtor

    bfn_jPar	= -1.0 / ( 8.0 * !pi * r^2 * sin ( thisCoLat_grid )^2 ); > max_jPar

    bfn_template_ph_grid	= -u0 / ( 4.0 * !pi * r * sin ( thisCoLat_grid ) ) * $
      ( ( r - rI * cos ( thisCoLat_grid ) ) / $
      sqrt ( r^2 - 2.0 * r * rI * cos ( thisCoLat_grid ) + rI^2 ) - 1.0 )

    iiNAN	= where ( bfn_template_ph_grid ne bfn_template_ph_grid, iiNANCnt )  ; check for probs (NaN)
    if iiNANCnt gt 0 then stop

	tic0 = SysTime(1)
    iridium_to_magnetic, bfn_template_ph_grid * 0, bfn_template_ph_grid, thisCoLat_grid, thisLon_grid, $
      xShift[ii], yShift[ii], $
      vThNew_grid, vPhNew_grid, $
      CLATM=coLatNew_rad, LONM=lonNew_rad, $
      VECX=ct_x, VECY=ct_y, VECZ=ct_z
	toc0 = SysTime(1)
	bFnTime = bFnTime + (toc0-tic0)

; bfnArray_jPar[ii,*]	= bfn_jPar; Assumed Iridium system for now

    bfnArrayPh_grid[ii,*]	= vPhNew_grid
    bfnArrayTh_grid[ii,*]	= vThNew_grid

; bfn_template_th_jcf	= 1.0 / ( 4.0 * !pi * rI ) * 1.0 / tan ( thisCoLat_grid / 2.0 ) < max_jcf;
; iridium_to_magnetic, bfn_template_th_jcf, bfn_template_th_jcf*0, $
	;	thisCoLat_grid, thisLon_grid, $
	;	xShift[ii], yShift[ii], $
	;	vThNew_grid, vPhNew_grid, $
	;	CLATM=coLatNew_rad, LONM=lonNew_rad, $
	;	VECX=ct_x, VECY=ct_y, VECZ=ct_z;

	;bfnArrayJcfPh_grid[ii,*]	= vPhNew_grid;
	;bfnArrayJcfTh_grid[ii,*]	= vThNew_grid;

    if keyword_set ( noPlot ) ne 1 then begin
		PlotDevice = 'win'
		if !version.os ne 'Win32' then PlotDevice = 'X'
		set_plot, PlotDevice

		device, decomposed = 0;
		loadct, 39;
		!p.background = 255;
		window, 0, xSize = 600, ySize = 300;
		loadct, 12, /silent;
		!p.multi = [2,2,1];
		plot_b_vectors, bfn_template_ph, bfn_template_ph * 0, $
			thisColat,thisLon, 40.0, 0.0, /plotAll, $
			0, 0, scaleFac = 1.5e-15, /multiPlot, /no_cont, /keep_ctbl, $
			/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5;
		tmp	= drawGrid ();
		!p.multi = [1,2,1];
		plot_b_vectors, vPhNew, vThNew, $
			coLat_data * !dtor, mlt_data * 15.0 * !dtor, 40.0, 0.0, $
			0, 0, scaleFac = 1.5e-15, /multiPlot, /no_cont, /keep_ctbl, $
			/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5, $
			/plotAll;
		tmp	= drawGrid ();
    endif
	stop
  endfor

	print, 'iridium_to_magnetic wallclock time: ', bFnTime 

  if keyword_set ( unitTh ) then begin
    bfnArray_data	= bfnArray_data / weight
  endif

  alpha	= transpose ( bfnArray_data ) ## bfnArray_data;
  if keyword_set ( unitTh ) then begin
    beta_	= transpose( transpose ( bfnArray_data ) ## bMag_data / weight[0,*] )
  endif else begin
    beta_	= transpose( transpose ( bfnArray_data ) ## bMag_data )
  endelse
  tic1 = SysTime(1)
  la_svd, alpha, w, u, v, status = svdStatus
  toc1 = SysTime(1)
  print, 'la_svd wall clock time: ', toc1-tic1

; if svdStatus ne 0 then goto, svdProblem

  kk	= where(w lt max ( w ) * 0.5e-4, kkcnt)
  print, 'Number of W elements set to zero:',kkcnt,' of ', n_elements( w[*] )
  if kkcnt gt 0 then w[kk]=0

  coeffs	= svsol ( u, w, v, beta_ )

;  jPar	= bFnArray_jPar ## coeffs;
;  jPar	= reform ( jPar, n_elements ( coLat_grid[*,0] ), $
;	   n_elements ( coLat_grid[0,*] ) ) * 1e-9 * 1e6; [uAm^-2]

  bThFit	= reform ( bFnArrayTh_grid ## coeffs, $
    n_elements ( coLat_grid[*,0] ), $
    n_elements ( coLat_grid[0,*] ) )
  bPhFit	= reform ( bFnArrayPh_grid ## coeffs, $
    n_elements ( coLat_grid[*,0] ), $
    n_elements ( coLat_grid[0,*] ) )

  bPhFit_data	= bFnArray_data ## coeffs

; Shift the fitted delta b's to aacgm

  iridium_to_magnetic, bThFit, bPhFit, coLat_grid * !dtor,$
    mlt_grid * 15.0 * !dtor, xAvg, yAvg, bThFit_mag, bPhFit_mag

  bThFit_mag	= reform ( bThFit_mag, $
    n_elements ( coLat_grid[*,0] ), $
    n_elements ( coLat_grid[0,*] ) );
  bPhFit_mag	= reform ( bPhFit_mag, $
    n_elements ( coLat_grid[*,0] ), $
    n_elements ( coLat_grid[0,*] ) );
;stop
;	Expand array to allow smooth derivatives at the boundaries

  bThFit_mag_expand	= fltarr ( n_elements ( bThFit_mag[*,0] ) + 4, $
    n_elements ( bThFit_mag[0,*] ) + 2 )
  bPhFit_mag_expand	= fltarr ( n_elements ( bThFit_mag[*,0] ) + 4, $
    n_elements ( bThFit_mag[0,*] ) + 2 )

  delTh	= abs ( coLat_grid_mag[ 0, 1 ] - coLat_grid_mag[ 0, 0 ] )
  delMlt	= abs ( mlt_grid_mag[ 1, 0 ] - mlt_grid_mag[ 0, 0 ] )

  coLat_grid_mag_expand	= transpose ( rebin ( findgen ( n_elements ( bThFit_mag[0,*] ) + 2 ), $
    n_elements ( bThFit_mag[0,*] ) + 2 , n_elements ( bPhFit_mag[*,0] ) + 4 ) ) * delTh
  mlt_grid_mag_expand	= rebin ( findgen ( n_elements ( bThFit_mag[*,0] ) + 4 ), $
    n_elements ( bThFit_mag[*,0] ) + 4, n_elements ( bPhFit_mag[0,*] ) + 2 ) * delMlt

  bThFit_mag_expand[ 2 : 2 + n_elements ( bThFit_mag[*,0] ) - 1,2 : * ]	= bThFit_mag
  bThFit_mag_expand[ 0 : 1, 2 : * ] = bThFit_mag[ n_elements ( bThFit_mag[*,0] ) - 2 : *, * ]
  bThFit_mag_expand[ n_elements ( bThFit_mag_expand[*,0] ) - 2 : *,2 : * ]	= bThFit_mag[ 0 : 1, * ]
  bThFit_mag_expand[ 2 : 2 + n_elements ( bThFit_mag[*,0] ) - 1, 0 : 1 ]	$
		= -1.0 * reverse ( shift ( bThFit_mag[ *, 0 : 1 ], n_elements ( bThFit_mag[*,0] ) / 2 ), 2 )
;	bThFit_mag_expand[ 2 + n_elements ( bThFit_mag[*,0] ) / 2 : $
;		2 + n_elements ( bThFit_mag[*,0] ) - 1, 0 : 1 ]	= $
;		reverse ( bThFit_mag_expand[ 2 + n_elements ( bThFit_mag[*,0] ) / 2 : $
;			2 + n_elements ( bThFit_mag[*,0] ) - 1, 0 : 1 ], 2 );
  bPhFit_mag_expand[ 2 : 2 + n_elements ( bPhFit_mag[*,0] ) - 1,2 : * ]	= bPhFit_mag
  bPhFit_mag_expand[ 0 : 1, 2 : * ] = bPhFit_mag[ n_elements ( bPhFit_mag[*,0] ) - 2 : *, * ]
  bPhFit_mag_expand[ n_elements ( bPhFit_mag_expand[*,0] ) - 2 : *, 2 : * ]	= bPhFit_mag[ 0 : 1, * ]
  bPhFit_mag_expand[ 2 : 2 + n_elements ( bPhFit_mag[*,0] ) - 1, 0 : 1 ]	$
		= -1.0 * reverse ( shift ( bPhFit_mag[ *, 0 : 1 ], n_elements ( bPhFit_mag[*,0] ) / 2 ), 2 )

  jPar_expand	= 1.0 $
		/ ( r * sin ( ( coLat_grid_mag_expand + min ( coLat_grid_mag )- 2 * delTh ) * !dtor ) ) $
		* ( sin ( ( coLat_grid_mag_expand + min ( coLat_grid_mag )- 2 * delTh ) * !dtor ) * $
		diffByDave ( bPhFit_mag_expand, coLat_grid_mag_expand * !dtor, 1 ) $
		+ cos ( ( coLat_grid_mag_expand + min ( coLat_grid_mag )- 2 * delTh ) * !dtor ) $
		* bPhFit_mag_expand $
		- diffByDave ( bThFit_mag_expand, mlt_grid_mag_expand * 15 * !dtor, 0 ) $
		) * 1e-9 * 1e6 / u0;

  jPar	= jPar_expand[ 2 : n_elements ( jPar_expand[*,0] ) - 3, 2 : * ]

; JcfTh	= bFnArrayJcfTh_grid ## coeffs
; JcfPh	= bFnArrayJcfPh_grid ## coeffs

  bPhFit_data	= bFnArray_data ## coeffs

;  Io	= reform ( coeffs, n_elements ( coLat_pole[*,0] ), $
;	  n_elements ( coLat_pole[0,*] ) );
;
;  jPar	= Io * 1.0e-9 / ( 4.0 * !pi * r^2 ) * 1e6; [uAm^-2]

; Plot SECS results: dB_in, poles, db_fit, FAC
  if keyword_set ( noPlot ) ne 1 then begin
;	set_plot, 'win'
    device, decomposed = 0;
    loadct, 13, file=colortb_path
;  !p.background = 0
    window, 1, xSize = 1200, ySize = 300
    !p.position = 0
    !p.multi = [4,4,1]

    plot_b_vectors, bMag_data, bMag_data*0, $
			coLat_data * !dtor,mlt_data * 15.0 * !dtor, maxLat, 0.0, /plotAll, $
			1, 0, scaleFac = 1e1, /multiPlot, /no_cont, /keep_ctbl, $
			/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5;
    tmp	= drawGrid ()

    !p.multi = [3,4,1]
    plot_b_vectors, bPhFit_data, bPhFit_data*0, $
			coLat_data * !dtor,mlt_data * 15.0 * !dtor, maxLat, 0.0, /plotAll, $
			1, 0, scaleFac = 1e1, /multiPlot, /no_cont, /keep_ctbl, $
			/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5;
    loadct, 12;
    plots, mlt_pole * 15.0, 90. - coLat_pole, psym = 1, color = 15*16-1
		tmp	= drawGrid ()

    !p.multi = [2,4,1];
    plot_b_vectors, bPhFit_mag, bThFit_mag, $
			coLat_grid_mag * !dtor,mlt_grid_mag * 15.0 * !dtor, maxLat, 0.0, /plotAll, $
			1, 0, scaleFac = 1e1, /multiPlot, /no_cont, /keep_ctbl, $
			/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5
    tmp	= drawGrid ()

;	plot_b_vectors, JcfTh, -JcfPh, $
;		coLat_grid * !dtor,mlt_grid * 15.0 * !dtor, 40.0, 0.0, /plotAll, $
;		1, 0, scaleFac = 1e8, /multiPlot, /no_cont, /keep_ctbl, $
;		/noGrid, color1 = 12 * 16 - 1, /noBorder, thick1 = 0.5;
;	tmp	= drawGrid ();

    !p.multi = [1,4,1];
    map_set, /ortho, 90, 0, /isotropic, $
      limit = [round(90.-42.0),0.,round(90.-0.),360.], $
      /noerase, /noborder, charsize=.01, clip=1, /advance

;	loadct, 12
;	plots, mlt_pole * 15.0, 90. - coLat_pole, psym = 1, $
;		color = 15*16-1;
	;plots, mlt_data * 15.0, 90. - coLat_data, psym = 4, $
	;	color = 14*16-1, symSize = 0.5;

	loadct, 13, file = colortb_path;

	facr	= 1.0
	levSpacing	= 0.1																; FAC plot contour spacing [uAm^-2]
	nLevs	= 50																			; Number of contour levels
	if ( nLevs MOD 2 ) ne 1 then nLevs = nLevs + 1	; Make nLevs odd
	fac_levs	= findgen(nlevs)*levSpacing - $				;	Array of contour levels for FAC plot
		(nLevs-1)*levSpacing/2;

	fac_linestyle	= (fac_levs lt 0.)*3;
	fac_colors	= bytscl(fac_levs, min=-facr, max=facr, top=253)+1;

	contour, [jpar,jpar[0,*]],[mlt_grid_mag,mlt_grid_mag[0,*]] * 15.0, $
		90.0 - [coLat_grid_mag, coLat_grid_mag[0,*]], /overplot,$
		levels=fac_levs, c_labels=fltarr(n_elements(fac_levs[*]))+1,$
		c_thick=fltarr(n_elements(fac_levs[*,0]))+1.,$
		c_colors=fac_colors, /fill;

	tmp = drawGrid ();

	endif

if keyword_set ( weight ) then begin
	unitTh = unitTh_Backup;
	unitPh = unitPh_Backup;
	weight = weight_Backup;
endif

if keyword_set ( hard ) then begin

	coLat_data	= coLat_data_Backup;
	mlt_data = mlt_data_Backup;
	bMag_data	= bMag_data_Backup;

endif

end
