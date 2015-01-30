; AMPERE Software - IDL 8.0 or later
;
; Main driver program for AMPERE dB data analysis to take Iridium vector dB
; fields and create a FAC and dB map on a AACGM and GEO regular grid.
; 
; SEC version trial - Oct 2014
; Calculate jPar from deltaB using SECS
; Juusola, L., et. al., Earth Planets Space, 58, 667-678, 2006
;
; Summary of processing:
;   Read in data (GEI)
;   Conv to GEO
;   do fit in GEO shifted
;   conv to AACGM
;
;  C.L. Waters and D.L. Green
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
; Code copyright: Johns Hopkins University Applied Physics Lab
;     Written under contract for Dr B.J. Anderson
;
pro ampfit_secs, sHr, eHr, south, $
	syr,smo,sday,sthr,stmn,stsc, $
	datFileName = datFileName, $
	nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
	mn_fac = mn_fac, mx_fac=mx_fac, $
	plt_png = plt_png, $
	aacgm_cap_coLat_deg = aacgm_cap_coLat_deg, $
	debug = debug, $
	sepoch = sepoch, eepoch = eepoch, $			  ; start/end epoch times
	extend = extend, $         							  ; get data across day if required
	SpacRes = SpacRes, status=status


  ;@d:\cwac\secs_olaf\amp_path_secs.pro				; specify dir paths
  @amp_path_secs.pro
  png_dir = './png'
  ;png_dir = 'd:\cwac\haje\'

  u0  = 4d0*!dpi*10d0^(-7)
  rI_m  = 6371000.0 + 110d3     ; [m]
  r_m = 6371000.0 + 780d3       ; [m] ( > rI )

  rI_km = rI_m/1000.0
  r_km = r_m/1000.0

; error trap routine - if error, then return to calling routine
  if (debug eq 0) then begin
	  status=0
    catch,error_status
    if (error_status ne 0) then begin
      status=1
		  Print,'ERR: ampfit_run90.pro => returning'
		  catch,/cancel
		  return
	  endif
	end

; Set AACGM coeff path
  ;aacgm_set_path,aacgmpath

; - - - - - - - - - - - 1. - - - - - - - - - - -
; Read AMPERE netcdf input data for full sphere
; - - - - - - - - - - - 1. - - - - - - - - - - -
  Print,'Reading AMPERE Data File....'
; Read GEI input data file, also do a geopack_reCalc
; Defines data_struc
; (a) read full day ncdf file
; (b) select data segment based on sHr->eHr (full sphere data)
; (c) Call geopack_recalc
  read_GEI_data90, sHr, eHr, $
    datFileName = datFileName, $
    data_struc = data_struc, dataGEI_orig = dataGEI_orig, $
    year = year, month = month, day = day, $
    avgYrSec = yrSecAvg, avgEpoch = avgEpoch, $
    debug = debug, $
    extend = extend, $
    status=status

; Trap for errors
  if (debug eq 0) then begin
    if (status ne 0) then return
  endif

; - - - - - - - - - - - 2. - - - - - - - - - - -
;   Convert the GEI coords and dB (x,y,z) into GEO
; - - - - - - - - - - - 2. - - - - - - - - - - -
; Converts full sphere data set - selection based on sHr->eHr
  amp_conv_GEI_GEO, dataGEI_orig, dataGEO_orig, avgEpoch

straighten = 1
if straighten then begin

; - - - - - - - - - - - 3. - - - - - - - - - - -
; Find the satellite track intersection point and calc rotation matrix
; - - - - - - - - - - - 3. - - - - - - - - - - -
; 'south' selects hemisphere of data for calc of intersection point
; returns cent_coLat_deg, cent_lon_deg coords
  amp_get_new_centre, dataGEO_orig, south, cent_coLat_deg, cent_lon_deg,XAvg,YAvg;, /diag
  print,'Centre coords (coLat, Lon):'
  print,'GEO ->', cent_coLat_deg, cent_lon_deg

; Calc rotation matrix
  amp_calc_rotmat_sec, r_km, cent_coLat_deg, cent_Lon_deg, south, rot_mat
  dataGEO_shifted = dataGEO_orig          ; Define dataGEO_Shifted as a copy of dataGEO

; - - - - - - - - - - - 4. - - - - - - - - - - -
; Rotate these GEO coords into shifted GEO coords
; - - - - - - - - - - - 4. - - - - - - - - - - -

  amp_rotate_coords, dataGEO_orig.px, dataGEO_orig.py, dataGEO_orig.pz, $
                     px_out, py_out, pz_out, $
                     r_out, th_out_deg, ph_out_deg, $  ; spherical coords of shifted coords
                     rot_mat

  dataGEO_shifted.px = px_out
  dataGEO_shifted.py = py_out
  dataGEO_shifted.pz = pz_out
  px_out=!null & py_out=!null & pz_out=!null

  dataGEO_shifted.r_km = r_out
  dataGEO_shifted.coLat_deg = th_out_deg
  dataGEO_shifted.lon_deg = ph_out_deg
  r_out=!null & th_out_deg=!null & ph_out_deg=!null

; - - - - - - - - - - - 5. - - - - - - - - - - -
; Rotate the dB vectors into shifted GEO coords
; - - - - - - - - - - - 5. - - - - - - - - - - -

  amp_rotate_vecs, dataGEO_shifted.px, dataGEO_shifted.py, dataGEO_shifted.pz, $
                   dataGEO_orig.dbx, dataGEO_orig.dby, dataGEO_orig.dbz, $
                   dbx_rot, dby_rot, dbz_rot, $     ; db(x,y,z) rotated
                   dbr_rot, dbth_rot, dbph_rot, $   ; db(r,t,p) rotated
                   rot_mat          ; rotation matrix
  dataGEO_shifted.dbx = dbx_rot
  dataGEO_shifted.dby = dby_rot
  dataGEO_shifted.dbz = dbz_rot
  dbx_rot=!null & dby_rot=!null & dbz_rot=!null

  dataGEO_shifted.br = dbr_rot
  dataGEO_shifted.bTheta = dbth_rot
  dataGEO_shifted.bPhi = dbph_rot
  dbr_rot=!null & dbth_rot=!null & dbph_rot=!null

; - - - - - - - - - - - 6. - - - - - - - - - - -
;   Select hemisphere data
; - - - - - - - - - - - 6. - - - - - - - - - - -
; Hemisphere data selection done on Shifted data, so need full sphere data up until this point
  tmp_dat = dataGEO_Shifted                ; dataGEO_Shifted struc here is full sphere data

  if (south eq 0) then iiSubSet = where ( dataGEO_Shifted.coLat_deg le 90.0)
  if (south eq 1) then iiSubSet = where ( dataGEO_Shifted.coLat_deg gt 90.0)
  dataGEO_Shifted = tmp_dat(iisubset)
  tmp_dat = !null

; - - - - - - - - - - - 7. - - - - - - - - - - -
; Sort data along each orbit track (required for ghosting)
; Check for coLat gaps in each track
; - - - - - - - - - - - 7. - - - - - - - - - - -

  trk_gap = intarr(2,6)
  amp_sort_trackdata, dataGEO_shifted, data_struc, south, SpacRes, trk_gap, Trk_order
; TRK_GAP is [2,  real track num]
; ---------------------------------------
; Check for data_gap exit condition
  For ii=0,5 do begin
    if trk_gap[1, trk_order[ii]] ne 0 then begin
      print,'## Data gap, track ',trk_order_sth[ii]
      status = 1
      return
    end
  end        ; Track Loop

;  Track order is now sequential and dataShifted structure is now sorted by:
;   (i) orbit plane [order specified by Tr_Num(6)]
;  (ii) distance from equator point where
; (iii) Lon start point is in same semi-circle in longitude (all < 180 deg)

  rev_rot_mat = transpose(rot_mat)

; ghost routine to go in here - to get SECS poles - for later

; unshift the dataGEO_Shifted data coords. Use these as they have Track_sort, zig_zag removal
  dataGEO_unSh = dataGEO_Shifted     ; Allocate MEM for dataGEO unshifted but sorted

  amp_rotate_coords, dataGEO_Shifted.px, dataGEO_Shifted.py, dataGEO_Shifted.pz, $ ; input dataGEO_Shifted, excluding ghosts
                     px_out, py_out, pz_out, $
                     r_out, th_out_deg, ph_out_deg, $  ; spherical coords of shifted coords
                     rev_rot_mat

  dataGEO_unSh.px = px_out
  dataGEO_unSh.py = py_out
  dataGEO_unSh.pz = pz_out
  px_out=!null & py_out=!null & pz_out=!null

  dataGEO_unSh.r_km = r_out
  dataGEO_unSh.coLat_deg = th_out_deg
  dataGEO_unSh.lon_deg = ph_out_deg
  r_out=!null & th_out_deg=!null & ph_out_deg=!null

; Second - unshift the dataGEO_Shifted dB vectors
  geopack_bspcar, dataGEO_Shifted.coLat_deg, dataGEO_Shifted.lon_deg, $
                  dataGEO_Shifted.bR, dataGEO_Shifted.bTheta, dataGEO_Shifted.bPhi, $
                  dbx_ird_sh, dby_ird_sh, dbz_ird_sh, /degree
; Rotate these dB vectors from shifted to unshifted uniform GEO
  amp_rotate_vecs, dataGEO_unSh.px, dataGEO_unSh.py, dataGEO_unSh.pz, $   ; new coord system coords
                   dbx_ird_sh, dby_ird_sh, dbz_ird_sh, $
                   dbx_ird, dby_ird, dbz_ird, $         ; db(x,y,z) rotated
                   dbR_ird, dbTheta_ird, dbPhi_ird, $   ; db(r,t,p) rotated
                   rev_rot_mat                          ; rotation matrix
  dataGEO_unSh.dbx = dbx_ird
  dataGEO_unSh.dby = dby_ird
  dataGEO_unSh.dbz = dbz_ird
  dbx_ird = !null & dby_ird = !null & dbz_ird = !null
  dbx_ird_sh = !null & dby_ird_sh = !null & dbz_ird_sh = !null

  dataGEO_unSh.br = dbR_ird
  dataGEO_unSh.btheta = dbTheta_ird
  dataGEO_unSh.bphi = dbphi_ird
  dbR_ird = !null & dbTheta_ird = !null & dbphi_ird = !null

; Use sorted GEO data
  dataGEO = dataGEO_unSh
  dataGEO = dataGEO_shifted
endif else begin ; _DO NOT SHIFT OR STRAIGHTEN_
	; Added to replace all the above, so comment this out if you uncomment
	; all the track straigtening stuff
	if (south eq 0) then iiSubSet = where ( dataGEO_orig.coLat_deg le 90.0)
	if (south eq 1) then iiSubSet = where ( dataGEO_orig.coLat_deg gt 90.0)
	dataGEO = dataGEO_orig[iiSubSet]
endelse

; - - - - - - - - - - - - - - - - - - - - - -
;      SECS section here
; - - - - - - - - - - - - - - - - - - - - - -

  print,'Calculating SECS...'

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use_cf = 0             ; bfn switch for curl_free or div_free bfns
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  polelimit1 = 50.     ; 50 in orig code
  polelimit2 = 50.

  coLat_data = dataGEO.coLat_deg
  mlt_data = (dataGEO.lon_deg/15.0 mod 24)

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Generate SECS pole locations
  maxX  =  50.0*sin(90.0*!dtor-!dpi )   ; ############ Check
  nXY  = 35 

; xPoleArray is dblarr[nXY,nXY]
  xPoleArray  = rebin ( ( findgen ( nXY ) - ( nXY - 1.0 ) / 2.0 ) / nXY * 2.0 * maxX, nXY, nXY ) ; array DBL[nXY, nXY]
  yPoleArray  = transpose ( xPoleArray )

  coLat_pole  = sqrt ( xPoleArray^2 + yPoleArray^2 ); [deg]
  mlt_pole  = ( 180.0 + atan ( yPoleArray, xPoleArray ) * !radeg ) / 15.0; [deg]

  iiPoleKeep  = where ( coLat_pole lt poleLimit1 )
  coLat_pole  = coLat_pole[iiPoleKeep]
  mlt_pole  = mlt_pole[iiPoleKeep]
	if south then coLat_pole = 180-coLat_pole

; Get a uniform grid
; -------------------
  minth = 2.0     ; in deg
  nLatGrid=30
;  nth = nLatGrid
  max_coLat = aacgm_cap_coLat_deg
  aacgm_cap_coLat_deg = 50.
  geoGrid_coLat_deg = transpose(rebin(findgen(nLatgrid)/(nLatgrid-1)*(aacgm_cap_coLat_deg-minth)+minth,nLatgrid,nLongrid)) ; deg
  ;if south eq 1 then aacgm_coLat_deg[*,*] = 180.0 - aacgm_coLat_deg[*,*]
  geoGrid_lon_deg = rebin(findgen(nLongrid)/(nLongrid)*24.*15.,nLongrid,nLatgrid)
  geoGrid_mlt = geoGrid_Lon_deg/15.0
  geoGrid_r_km = geoGrid_lon_deg*0 + R_km
  mltshift = 0.0
  if south then geoGrid_coLat_deg = 180-geoGrid_coLat_deg

	;; Set poles to be co-located with grid
	;coLat_pole = geoGrid_coLat_deg
	;mlt_pole = geoGrid_mlt

	;; Set poles to be co-located with data
	;coLat_pole = dataGEO.coLat_deg
	;mlt_pole = dataGEO.lon_deg/15.0


  n_el = n_elements(coLat_pole)            ; num SECS to use
  size_pole = fltArr(n_el) + 200e3 ; m

  	;; Overwrite the data with a basis funciton
	;ii = 12
	;ampfit_create_arb_bfn, r_m, coLat_pole[ii], mlt_pole[ii]*15, $
	;			dataGEO.coLat_deg*0+r_m, dataGEO.coLat_deg, dataGEO.lon_deg, size_pole[ii], $
	;			bfn_bR, bfn_bT, bfn_bP
	;dataGEO.bR = bfn_bR
	;dataGEO.bTheta = bfn_bT
	;dataGEO.bPhi = bfn_bP	

; Calc bfns at GEO obs locations for phi=0 SECS - in order to get coeffs
  n_obs = n_elements(dataGEO.coLat_deg)        ; obs locations to use

  Td_a = fltarr(n_el, 2*n_obs)      ; T matrix
  bd_a = fltarr(2*n_obs)            ; db vector (column matrix)

; There are 3 coord systems involved - all in GEO
; i) coords for poles of SECS : [coLat_pole,mlt_pole*15]
; ii) coords of the dB data : [coLat_data,mlt_data*15]
; iii) coords of the desired output grid [coLat_grid,mlt_grid*15] - used later

; order the dB(th,ph) components in the matrices here
  th_kk = indgen(n_obs)        ; idx for db_th data in matrices, put in top 1/2
  ph_kk = indgen(n_obs)+n_obs  ; idx for db_ph data in matrices, put in bottom half


;; Full biot-savart integral for each point
;for ii=0,n_el-1 do begin
;
;	ampfit_create_arb_bfn, r_m, coLat_pole[ii], mlt_pole[ii]*15, $
;			dataGEO.coLat_deg*0+r_m, dataGEO.coLat_deg, dataGEO.lon_deg, size_pole[ii], $
;			bfn_bR, bfn_bT, bfn_bP
;
;	scale = 1e4
;	plt_dat90, dataGEO.coLat_deg, dataGEO.lon_deg, $
;	    n_elements(dataGEO.coLat_deg[*]), $
;	    -bfn_bT*scale, bfn_bP*scale, south,title='Test', capSize=cp_sz
;
;    Td_a[ii,th_kk] = bfn_bT
;    Td_a[ii,ph_kk] = bfn_bP
;
;endfor

; Interpolate from biot-savart derived template
for ii=0,n_el-1 do begin 

	; Get rot_mat for rotation from 0,0 to coLat_pole[ii],mlt-pole[ii]*15
    amp_calc_rotmat_sec, r_km,coLat_pole[ii], mlt_pole[ii]*15.0, south, rot_mat

	; rotate input dB data coords into coord system with pole at coLat_pole[ii],mlt_pole[ii]*15
    amp_rotate_coords, dataGEO.px, dataGEO.py, dataGEO.pz, $
                     px_out, py_out, pz_out, $
                     r_out, th_out_deg, ph_out_deg, $    ; spherical components of the rotated coords
                     rot_mat

	; This allows us to get coLat array to use in bfn generation
    ecs_coLat_rad = (th_out_deg)*!dtor    ; load in coLats of dB data (in GEO)
    ecs_Lon_rad = ph_out_deg*!dtor
    ecs_r_m = ph_out_deg*0 + r_m

	; generate bfn at n_obs locations
	_nGrid = 100
	_ColatGrid = fIndGen(_nGrid)/(_nGrid-1)*max(ecs_coLat_rad*!radeg)	
	ampfit_create_arb_bfn, r_m, coLat_pole[ii]*0, mlt_pole[ii]*15*0, $
			_CoLatGrid*0+r_m, _CoLatGrid, _CoLatGrid*0, size_pole[ii], $
			bfn_bR, bfn_bT, bfn_bP

	if south then begin
		bfn_template = interpol(bfn_bP,_CoLatGrid,180-ecs_coLat_rad*!radeg)
	endif else begin
		bfn_template = interpol(bfn_bP,_CoLatGrid,ecs_coLat_rad*!radeg)
	endelse
	
	; Calc bx,by,bz components of bfn in shifted pole coords (pole with thet'=0)
	geopack_bspcar, $
		ecs_coLat_rad, ecs_Lon_rad, $
		bfn_template*0, bfn_template*0, bfn_template, $   
		bfn_dbx, bfn_dby, bfn_dbz                 

	; rotate these bfn dB back to original GEO components
    rev_rot_mat = transpose(rot_mat)
    amp_rotate_vecs, dataGEO.px, dataGEO.py, dataGEO.pz, $     ; input already rotated x,y,z
                   bfn_dbx, bfn_dby, bfn_dbz, $          ; input dB(x,y,z) of bfn dB from above
                   dbx_rot, dby_rot, dbz_rot, $          ; dB(x,y,z) rotated
                   dbr_rot, dbth_rot, dbph_rot, $        ; dB(r,t,p) rotated
                   rev_rot_mat                           ; rotation matrix to use
	;scale = 1e4
	;plt_dat90, dataGEO.coLat_deg, dataGEO.lon_deg, $
	;    n_elements(dataGEO.coLat_deg[*]), $
	;    -dbth_rot*scale, dbph_rot*scale, south,title='Test', capSize=cp_sz

    Td_a[ii,th_kk] = dbth_rot
    Td_a[ii,ph_kk] = dbph_rot

endfor

print,'Calling SVD process...'

bd_a[th_kk] = dataGEO.btheta  
bd_a[ph_kk] = dataGEO.bphi

la_svd, Td_a, w, u, v, status = svdStatus  
if svdStatus then stop

; Set tolerance on singular values
svd_tol = 1.0d-1        ; smaller number here adds more W elements -> achieves a finer fit

kk  = where(w lt max(w)*svd_tol, kkcnt)
print, 'Number of W elements set to zero:',kkcnt,' of ', n_elements( w[*] )
if kkcnt gt 0 then w[kk]=0

coeffs  = svsol ( u, w, v, bd_a )

bThFit_data = reform(Td_a[*,th_kk] ## coeffs[*])
bPhFit_data = reform(Td_a[*,ph_kk] ## coeffs[*])

; Plot original and fit data

	PlotDevice = PlotDevice 
	if !version.os ne 'Win32' then PlotDevice = 'X'

  set_plot, PlotDevice 
  device, decomposed = 0;
  loadct, 0
  wn = 0
  window, wn, xSize = 1050, ySize = 400, title='Orig (dot) and Fit (solid) dB data'
  !p.multi = [0,1,2]
  !p.charsize=1.0
  o_col = 100
  f_col=0
  plot,bThFit_data,linestyle=0,title='dB_theta',color=f_col,background=255
  oplot,dataGEO.btheta,linestyle=1,color=o_col
  rms_th = total( (bThFit_data-dataGEO.btheta)^2 )
  rms_th = sqrt(rms_th/n_obs)
  print,'RMS_bth = ',rms_th

  plot,bPhFit_data,linestyle=0,title='db_phi',color=f_col,background=255
  oplot,dataGEO.bphi,linestyle=1,color=o_col
  rms_ph = total( (bPhFit_data-dataGEO.bphi)^2 )
  rms_ph = sqrt(rms_ph/n_obs)
  print,'RMS_bph = ',rms_ph

  wn++

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;   Now get dB at output grid locations
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Create uniform grid in AACGM and output the assoc GEO grid
;
; Generate uniform spaced grid in AACGM
; Calc mltshift
; Convert these to GEO
;  mltref=0.0

; dont use just yet
;  ev_epoch = (sepoch+eepoch)/2.0
;  aacgm_cap_coLat_deg = 90.
;  aacgm_grid_GEO_secs, aacgm_cap_coLat_deg, south, $
;       aacgmGrid_coLat_deg = aacgmGrid_coLat_deg, $
;       aacgmGrid_lon_deg = aacgmGrid_lon_deg, $
;       geoGrid_coLat_deg = geoGrid_coLat_deg, $
;       geoGrid_lon_deg = geoGrid_lon_deg , $
;       geoGrid_x = geoGrid_x, $                              ; XYZ locations of uniform aacgm grid (unshifted)
;       geoGrid_y = geoGrid_y, $
;       geoGrid_z = geoGrid_z, $
;       nLat = nLatGrid, nLon = nLonGrid, $
;       year = year, ev_epoch=ev_epoch, $
;       mltShift = mltShift, $
;       mltref = mltref

; Get a uniform grid
; -------------------
  minth = 2.0     ; in deg
  nLatGrid=30
;  nth = nLatGrid
  max_coLat = aacgm_cap_coLat_deg
  aacgm_cap_coLat_deg = 50.
  geoGrid_coLat_deg = transpose(rebin(findgen(nLatgrid)/(nLatgrid-1)*(aacgm_cap_coLat_deg-minth)+minth,nLatgrid,nLongrid)) ; deg
  ;if south eq 1 then aacgm_coLat_deg[*,*] = 180.0 - aacgm_coLat_deg[*,*]
  geoGrid_lon_deg = rebin(findgen(nLongrid)/(nLongrid)*24.*15.,nLongrid,nLatgrid)
  geoGrid_mlt = geoGrid_Lon_deg/15.0
  geoGrid_r_km = geoGrid_lon_deg*0 + R_km
  mltshift = 0.0
  if south then geoGrid_coLat_deg = 180-geoGrid_coLat_deg

; - - - - -  to be changed later - - - - - - -

; Convert from (r,t,p) to (x,y,z) coords
  geoPack_sphCar,geoGrid_R_km,geoGrid_colat_deg, geoGrid_lon_deg, $
    geoGrid_x, geoGrid_y, geoGrid_z, /to_rect, /degree

  n_obs = n_elements(geoGrid_coLat_deg)

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; For testing
;;sel_ii = where(dataGEO.ipln eq 0)  ; select plane - testing
;np = n_elements(dataGEO.px)
;;np = 100
;sel_ii = indgen(np)            ; select 1st 100 points - testing
;n_obs = n_elements(sel_ii)
;twist_deg = 0.0                                    ; Change here for different twists on lon
;geoGrid_coLat_deg = dataGEO.coLat_deg
;geoGrid_Lon_deg = dataGEO.Lon_deg + twist_deg      ; twist the data a bit
;geoGrid_r_km = geoGrid_lon_deg*0 + R_km
;geoPack_sphCar, geoGrid_R_km,geoGrid_colat_deg*!dtor, geoGrid_lon_deg*!dtor, $
;    geoGrid_x, geoGrid_y, geoGrid_z, /to_rect
;geoGrid_x = dblarr(n_obs) & geoGrid_x[sel_ii] = dataGEO[sel_ii].px
;geoGrid_y = dblarr(n_obs) & geoGrid_y[sel_ii] = dataGEO[sel_ii].py
;geoGrid_z = dblarr(n_obs) & geoGrid_z[sel_ii] = dataGEO[sel_ii].pz
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  T_a = fltarr(n_el, 2*n_obs)      ; T matrix
  b_a = fltarr(2*n_obs)            ; db vector (column matrix)

; order the dB(th,ph) components in the matrices here
  th_kk = indgen(n_obs)        ; idx for db_th data in matrices, put in top 1/2
  ph_kk = indgen(n_obs)+n_obs  ; idx for db_ph data in matrices, put in bottom half

;!p.multi=0
;++wn
;window, wn, xsize=600,ysize=600

;; Full Biot-Savart Guassian implementation
;for ii=0,n_el-1 do begin
;
;	ampfit_create_arb_bfn, r_m, coLat_pole[ii], mlt_pole[ii]*15, $
;			geoGrid_colat_deg*0+r_m, geoGrid_colat_deg, geoGrid_lon_deg, size_pole[ii], $
;			bfn_bR, bfn_bT, bfn_bP
;
;	;scale = 1e4
;	;plt_dat90, geoGrid_colat_deg, geoGrid_lon_deg, $
;	;    n_elements(geoGrid_colat_deg[*]), $
;	;    -bfn_bT*scale, bfn_bP*scale, south,title='Test', capSize=cp_sz
;
;    T_a[ii,th_kk] = bfn_bT 
;    T_a[ii,ph_kk] = bfn_bP 
;
;endfor

!p.multi=0

; Interpolate biot-savart derived template per pole
for ii = 0, n_el-1 do begin 

	; Get rot_mat for rotation from 0,0 to coLat_pole[ii],mlt-pole[ii]*15
    amp_calc_rotmat_sec, r_km, coLat_pole[ii], mlt_pole[ii]*15.0, south, rot_mat

	; rotate input dB data coords into coord system with pole at coLat_pole[ii],mlt_pole[ii]*15
    amp_rotate_coords, geoGrid_x, geoGrid_y, geoGrid_z, $
                     px_out, py_out, pz_out, $
                     r_out, th_out_deg, ph_out_deg, $    ; spherical components of the rotated coords
                     rot_mat

	; This allows us to get coLat array to use in bfn generation
    ecs_coLat_rad = (th_out_deg[*])*!dtor    ; load in coLats of dB data (in GEO)
    ecs_Lon_rad = ph_out_deg[*]*!dtor
    ecs_r_m = ph_out_deg[*]*0 + r_m

	_nGrid = 100
	_ColatGrid = fIndGen(_nGrid)/(_nGrid-1)*max(ecs_coLat_rad*!radeg)	
	ampfit_create_arb_bfn, r_m, coLat_pole[ii]*0, mlt_pole[ii]*15*0, $
			_CoLatGrid*0+r_m, _CoLatGrid, _CoLatGrid*0, size_pole[ii], $
			bfn_bR, bfn_bT, bfn_bP

	if south then begin
		bfn_template = interpol(bfn_bP,_CoLatGrid,180-ecs_coLat_rad*!radeg)
	endif else begin
	   	bfn_template = interpol(bfn_bP,_CoLatGrid,ecs_coLat_rad*!radeg)
	endelse

	; Calc bx,by,bz components of bfn in shifted pole coords (pole with thet'=0)
	geopack_bspcar, $
		ecs_coLat_rad, ecs_Lon_rad, $
		bfn_template*0, bfn_template*0, bfn_template, $   ; only phi component in df bfn
		bfn_dbx, bfn_dby, bfn_dbz                         ; these are the bfn components in x,y,z

	; rotate these bfn dB back to original GEO components
    rev_rot_mat = transpose(rot_mat)
    amp_rotate_vecs, geoGrid_x, geoGrid_y, geoGrid_z, $     ; input already rotated x,y,z
                   bfn_dbx, bfn_dby, bfn_dbz, $          ; input dB(x,y,z) of bfn dB
                   dbx_rot, dby_rot, dbz_rot, $          ; dB(x,y,z) rotated
                   dbr_rot, dbth_rot, dbph_rot, $        ; dB(r,t,p) rotated
                   rev_rot_mat                           ; rotation matrix to use

    T_a[ii,th_kk] = dbth_rot
    T_a[ii,ph_kk] = dbph_rot

	;scale = 1e4
	;plt_dat90, geoGrid_colat_deg, geoGrid_lon_deg, $
	;    n_elements(geoGrid_colat_deg[*]), $
	;    -dbth_rot*scale, dbph_rot*scale, south,title='Test', capSize=cp_sz

	bfn_bT = dbth_rot

endfor

; test
;  bthfit_tst = reform(T_a[*,th_kk] ## coeffs[*])
;  bphfit_tst = reform(T_a[*,ph_kk] ## coeffs[*])
;  !p.multi=0
;  window, wn, xSize = 1050, ySize = 400, title='test'
;  !p.multi = [0,1,2]
;  o_col = 100 & f_col=0
;  plot,bthfit_tst,linestyle=0,title='dB_theta',color=f_col,background=255
;  oplot,dataGEO.btheta,linestyle=1,color=o_col
;  plot,bphfit_tst,linestyle=0,title='dB_phi',color=f_col,background=255
;  oplot,dataGEO.bphi,linestyle=1,color=o_col
;stop
  
  bThFit_geo = reform(T_a[*,th_kk] ## coeffs[*], $
     n_elements(geoGrid_coLat_deg[*,0] ), n_elements(geoGrid_coLat_deg[0,*]))
  bPhFit_geo = reform(T_a[*,ph_kk] ## coeffs[*], $)
     n_elements(geoGrid_coLat_deg[*,0] ), n_elements(geoGrid_coLat_deg[0,*]))

	_jPar = bfn_bT*0
	_jParCnt = bfn_bT*0
 	for ii = 0, n_el-1 do begin 
		ampfit_calc_fac_arb, r_m, coLat_pole[ii], mlt_pole[ii]*15, $
			geoGrid_colat_deg*0+r_m, geoGrid_colat_deg, geoGrid_lon_deg, size_pole[ii], coeffs[ii], $
			_jPar, _jParCnt
	endfor
	_i = where(_jParCnt gt 0)
	_jPar = reform(_jPar,n_elements(geoGrid_coLat_deg[*,0] ), n_elements(geoGrid_coLat_deg[0,*]))
	_jPar = _jPar * 1e-3 ; uA/m^2

; Expand array to allow smooth derivatives at the boundaries
  bThFit_geo_expand = fltarr(n_elements(bThFit_geo[*,0] ) + 4, n_elements(bThFit_geo[0,*] ) + 2 )
  bPhFit_geo_expand = fltarr(n_elements(bThFit_geo[*,0] ) + 4, n_elements(bThFit_geo[0,*] ) + 2 )

  delTh = abs(geoGrid_coLat_deg[0,1] - geoGrid_coLat_deg[0,0])
  delMlt  = abs(geoGrid_lon_deg[1,0]/15.0 - geoGrid_lon_deg[0,0]/15.0)

  coLat_grid_expand = transpose(rebin(findgen(n_elements(bThFit_geo[0,*] ) + 2 ), $
    n_elements(bThFit_geo[0,*] ) + 2 , n_elements(bPhFit_geo[*,0] ) + 4 ) ) * delTh
	if south then coLat_grid_expand = 180-coLat_Grid_Expand
  mlt_grid_expand = rebin(findgen(n_elements(bThFit_geo[*,0] ) + 4 ), $
    n_elements(bThFit_geo[*,0] ) + 4, n_elements(bPhFit_geo[0,*] ) + 2 ) * delMlt

  bThFit_geo_expand[ 2 : 2 + n_elements ( bThFit_geo[*,0] ) - 1,2 : * ] = bThFit_geo
  bThFit_geo_expand[ 0 : 1, 2 : * ] = bThFit_geo[ n_elements ( bThFit_geo[*,0] ) - 2 : *, * ]
  bThFit_geo_expand[ n_elements ( bThFit_geo_expand[*,0] ) - 2 : *,2 : * ]  = bThFit_geo[ 0 : 1, * ]
  bThFit_geo_expand[ 2 : 2 + n_elements ( bThFit_geo[*,0] ) - 1, 0 : 1 ]  $
      = -1.0 * reverse ( shift ( bThFit_geo[ *, 0 : 1 ], n_elements ( bThFit_geo[*,0] ) / 2 ), 2 )
    ; bThFit_mag_expand[ 2 + n_elements ( bThFit_mag[*,0] ) / 2 : $
    ;   2 + n_elements ( bThFit_mag[*,0] ) - 1, 0 : 1 ] = $
    ;   reverse ( bThFit_mag_expand[ 2 + n_elements ( bThFit_mag[*,0] ) / 2 : $
    ;     2 + n_elements ( bThFit_mag[*,0] ) - 1, 0 : 1 ], 2 );
  bPhFit_geo_expand[ 2 : 2 + n_elements ( bPhFit_geo[*,0] ) - 1,2 : * ] = bPhFit_geo
  bPhFit_geo_expand[ 0 : 1, 2 : * ] = bPhFit_geo[ n_elements ( bPhFit_geo[*,0] ) - 2 : *, * ]
  bPhFit_geo_expand[ n_elements ( bPhFit_geo_expand[*,0] ) - 2 : *, 2 : * ] = bPhFit_geo[ 0 : 1, * ]
  bPhFit_geo_expand[ 2 : 2 + n_elements ( bPhFit_geo[*,0] ) - 1, 0 : 1 ]  $
      = -1.0 * reverse ( shift ( bPhFit_geo[ *, 0 : 1 ], n_elements ( bPhFit_geo[*,0] ) / 2 ), 2 )

; Numerical curl for J//
  jPar_expand = 1.0 $
      / ( r_m * sin ( ( coLat_grid_expand + min ( geoGrid_coLat_deg )- 2 * delTh ) * !dtor ) ) $
      * ( sin((coLat_grid_expand + min(geoGrid_coLat_deg )- 2 * delTh ) * !dtor ) * $
      diffByDave ( bPhFit_geo_expand, coLat_grid_expand * !dtor, 1 ) $
      + cos ( ( coLat_grid_expand + min ( geoGrid_coLat_deg )- 2 * delTh ) * !dtor ) $
      * bPhFit_geo_expand $
      - diffByDave ( bThFit_geo_expand, mlt_grid_expand * 15 * !dtor, 0 ) $
      ) * 1e-9 * 1e6 / u0

  jPar  = jPar_expand[ 2 : n_elements(jPar_expand[*,0] ) - 3, 2 : * ]
  if south then jPar = -jPar
  
;	jPar = _jPar

; = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
; Plot input data: dataGEI_orig, dataGEO_orig, AACGM of input data

;  set_plot, 'win'
;  device, decomposed = 0;
;  loadct, 13, file=colortb_path
  loadct, 0
  if south eq 1 then begin
    iiSubSet = where ( dataGEI_Orig.coLat_deg ge 90.0, iiSubSetCnt)
  endif else begin             ; North hemisphere
    iiSubSet = where ( dataGEI_Orig.coLat_deg le 90.0, iiSubSetCnt )
  end
  
; get date/time for plot labels
  cdf_epoch, sepoch, yr,mnth,day,hr,minut,sec, /break
  date_str = STRING(yr,format = '(i4.4)') + $
             STRING(mnth,format = '(i2.2)') + $
             STRING(day,format = '(i2.2)')
  time_str = STRING(hr,format = '(i2.2)') + $
             STRING(minut,format = '(i2.2)') + $
             STRING(sec,format = '(i2.2)')

  window, wn, xSize = 1050, ySize = 400, title='AMPERE dB Data: '+date_str+'_'+time_str
  !p.multi = [0,3,1]
  !p.charsize=1.0
  cp_sz=50.

  plt_dat90,dataGEI_Orig[iiSubSet].coLat_deg, dataGEI_Orig[iiSubSet].lon_deg, $
      n_elements(iiSubSet), $
      -dataGEI_Orig[iiSubSet].bTheta, dataGEI_Orig[iiSubSet].bPhi, $
        south,title='Input dB: dataGEI_Original', capSize=cp_sz

  plt_dat90,dataGEO_Orig.coLat_deg, dataGEO_Orig.lon_deg, $
      n_elements(dataGEO_Orig.coLat_deg), $
      -dataGEO_Orig.bTheta, dataGEO_Orig.bPhi, $
      south,title='Input dB: dataGEO_Original',capSize=cp_sz

  dataGEO_ird = dataGEO     ; Allocate MEM for dataGEO_ird

; Third - Original input data to AACGM
  if south eq 1 then idx_sel=where( dataGEO_ird[*].coLat_deg gt 110.)
  if south eq 0 then idx_sel=where( dataGEO_ird[*].coLat_deg lt 70.)  ; ensure data within AACGM calc range
  geo_lat_deg = 90.0 - dataGEO_ird[idx_sel].coLat_deg
  hgt = replicate(780.0,n_elements(geo_lat_deg))
; Calc AACGM [lat,lon] coords of GEO input data plus AACGM dB's
  aacgm_conv_vec, $
    geo_lat_deg, dataGEO_ird[idx_sel].lon_deg, hgt, $
    dataGEO_ird[idx_sel].bTheta, dataGEO_ird[idx_sel].bPhi, $
    aacgm_ird_lat_deg, aacgm_ird_lon_deg, $
    dBTheta_aacgm_ird_gth, dBPhi_aacgm_ird_gth, $
    dBTheta_aacgm_ird_gph, dBPhi_aacgm_ird_gph, err, /to_aacgm

; get mlon for 0 MLT
;    cdf_epoch, epoch0, y,m,d,h,mn,sec ,/compute_epoch
  epoch0 = (sepoch+eepoch)/2.0
  yrsec=cnvtime(syr,smo,sday,sthr,stmn,stsc)
  zer_mlt = aacgm_mlt(syr,yrsec,0.0)
  ;zer_mlt = aacgm_mlt(epoch0,0.0)      ; get MLT for 0 deg of mlon

  plt_dat90, 90.0-aacgm_ird_lat_deg[*], ((aacgm_ird_lon_deg[*]/15.0+zer_mlt) mod 24)*15.0, $
;  plt_dat90, 90.0-aacgm_ird_lat_deg[*], ((aacgm_ird_lon_deg[*]/15.0) mod 24)*15.0, $
    n_elements(aacgm_ird_lat_deg[*] ), $
    -dBTheta_aacgm_ird_gth, dBPhi_aacgm_ird_gph, $
    south, title = 'Input dB: AACGM(Lat,MLT)', $
    capSize = abs(aacgm_cap_coLat_deg)

; PNG output
  if plt_png eq 1 then begin
    pngfname = png_dir + 'amp_secs_indB_'+date_str+'_'+time_str+'_N.png'
    image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
    write_PNG,pngfname,image,r,g,b
    Print,'PNG for FAC written to ', pngfname
  end
  wn++

; = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

; Plot fitted dB data
!p.multi = [0,2,2]

window, wn, xSize = 800, ySize = 800, title='Fit dB Data: '+date_str+'_'+time_str

plt_dat90,coLat_data, mlt_data*15.0, n_elements(bPhFit_data), $
    -bThFit_data, bPhFit_data, $
      south,title='Fitted dB_phi [GEO]', capSize=cp_sz

plt_dat90,geoGrid_coLat_deg, geoGrid_lon_deg, n_elements(bPhFit_geo), $
    -bThFit_geo, bPhFit_geo, $
    south,title='Fitted dB grid [GEO]', capSize=cp_sz

; Polar plot of FAC
if south then begin
	pole=-90
	capLimit=180-aacgm_cap_coLat_deg
end else begin
	pole = 90
	capLimit = aacgm_cap_coLat_deg
end

; Overlap for plotting
ThisJPar = [jPar,jPar[0,*]]
ThisCoLat = [geoGrid_coLat_deg,geoGrid_coLat_deg[0,*]]
ThisLon = [geoGrid_lon_deg,geoGrid_lon_deg[0,*]]

title = 'jPar (curl)'
if south then title = title + ' (through Earth from NH)'
plot_fac, ThisJPar, pole, capLimit, ThisCoLat, ThisLon, $ 
		mlt_shift = 0, title = title, south = south
xyOuts, mlt_pole*15, 90-coLat_pole, 'x', color=0


ThisJPar = [_jPar,_jPar[0,*]]
title = '_jPar (pure)'
if south then title = title + ' (through Earth from NH)'
plot_fac, ThisJPar, pole, capLimit, ThisCoLat, ThisLon, $ 
		mlt_shift = 0, title = title, south = south
xyOuts, mlt_pole*15, 90-coLat_pole, 'x', color=0


;; AACGM jPar
;    window, wn, xSize = 500, ySize = 400, title='j!d//!n [AACGM]: '+date_str+'_'+time_str
;    !p.multi=0
;; AACGM section
;    geoGrid_lat_deg = 90.0 - geoGrid_coLat_deg
;    hgt = replicate(780.0,n_elements(geoGrid_lat_deg))
;; Convert GEO grid to AACGM
;    aacgm_conv_coord, geoGrid_lat_deg, geoGrid_lon_deg, hgt, $
;      aacgm_lat_deg, aacgm_lon_deg, err, /to_aacgm
;; get mlon for 0 MLT
;;    cdf_epoch, epoch0, y,m,d,h,mn,sec ,/compute_epoch
;    epoch0 = (sepoch+eepoch)/2.0
;    zer_mlt = aacgm_mlt(epoch0,0.0)      ; get MLT for 0 mlon
;; rotate mlons
;    mlon_plt = ( (aacgm_lon_deg[*,*]/15.0+zer_mlt) mod 24 )*15.0
;
;;    capLimit = 30.
;    limit = [ pole - capLimit, 0, pole, 360 ] ; use Nth hemis for both Nth and Sth views
;    map_set, pole, 0, 0, /ortho, /iso, limit = limit, $
;       /noborder, title = 'jPar [AACGM]', /advance, $
;       yMargin = [1,3];, color = 0
;    contour, [jpar,jpar[0,*]],[mlon_plt,mlon_plt[0,*]], $
;        [aacgm_lat_deg,aacgm_lat_deg[0,*]], /overplot, $
;        levels=fac_levs, c_labels=fltarr(n_elements(fac_levs[*]))+1, $
;        c_thick=fltarr(n_elements(fac_levs[*,0]))+1.,$
;        c_colors=fac_colors, /fill;
;    nLats = fix((capLimit)/10)
;    lats = 90-(fIndGen(nLats)+1)*10
;    lnlbp=min(abs(lats),idx)
;    lnpos=lats(idx)+3
;    lt_lab=strtrim(string(lats,format='(i2)'),2)
;    if south then lt_lab=strtrim(string(-lats,format='(i3)'),2)
;    map_grid, label = 1, $
;      lonNames  = ['0','6','12','18',''], $
;      lons = [0,90,180,270,360], $
;      lonlab=lnpos, $
;      latLab = 45, $
;      lats = lats, $
;      latNames=lt_lab, $
;      color = 255
;      
;; PNG output
;    if plt_png eq 1 then begin
;      pngfname = png_dir + 'amp_secs_fac_'+date_str+'_'+time_str+'_N.png'
;      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
;      write_PNG,pngfname,image,r,g,b
;      Print,'PNG for FAC written to ', pngfname
;    end


  print,'Finished'

end
