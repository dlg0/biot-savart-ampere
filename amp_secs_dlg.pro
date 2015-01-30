; AMPERE Software - IDL 8.0 or later
;
; Main driver program for AMPERE dB data analysis to take Iridium vector dB
; fields and create a FAC and dB map on a AACGM and GEO regular grid.
; 
; SEC version trial - Oct 2014
; Calculate jPar from deltaB using SECS
; Juusola, L., et. al., Earth Planets Space, 58, 667-678, 2006
;
;  C.L. Waters and D.L. Green
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
pro amp_secs_dlg

	dlm_register, 'idl_geopack.dlm'
	dlm_load, 'Geopack'

  ;sav_pth = 'd:\cwac\secs_olaf\dlg_amp\'
  ;sav_Name = sav_pth+'20100824_0520_amp.sav'
  sav_Name = '20100824_0520_amp.sav'
  restore, sav_Name
; restores:
;  dataGEI_orig, dataGEO_orig, data_struc
;  avgepoch, sHr,eHr,south
;  nLatGrid,nLonGrid,aacgm_cap_coLat_deg
;

; Find the satellite track intersection point and calc rotation matrix
; 'south' selects hemisphere of data for calc of intersection point
; returns cent_coLat_deg, cent_lon_deg coords
  amp_get_new_centre, dataGEO_orig, south, cent_coLat_deg, cent_lon_deg,XAvg,YAvg;, /diag
  print,'Centre coords (coLat, Lon):'
  print,'GEO ->', cent_coLat_deg, cent_lon_deg

; Calc rotation matrix
  amp_calc_rotmat, cent_coLat_deg, cent_Lon_deg, south, rot_mat
  dataGEO_shifted = dataGEO_orig          ; Define dataGEO_Shifted as a copy of dataGEO

; Rotate these GEO coords into shifted GEO coords
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

; Rotate the dB vectors into shifted GEO coords
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

;   Select hemisphere data
; Hemisphere data selection done on Shifted data, so need full sphere data up until this point
  tmp_dat = dataGEO_Shifted                ; dataGEO_Shifted struc here is full sphere data
  if (south eq 0) then iiSubSet = where ( dataGEO_Shifted.coLat_deg le 90.0)
  if (south eq 1) then iiSubSet = where ( dataGEO_Shifted.coLat_deg gt 90.0)
  dataGEO_Shifted = tmp_dat(iisubset)
  tmp_dat = !null

; Sort data along each orbit track (required for ghosting)
; Check for coLat gaps in each track
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
;
; First - unshift the dataGEO_Shifted data coords. Use these as they have Track_sort, zig_zag removal
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
  XAvg = 0.0
  YAvg = 0.0

; Plot input data: dataGEI_orig, dataGEO_orig, AACGM of input data

	PlotDevice = 'win'
	if !version.os ne 'Win32' then PlotDevice = 'X'

  set_plot, PlotDevice 
  device, decomposed = 0;
  if south eq 1 then begin
    iiSubSet = where ( dataGEI_Orig.coLat_deg ge 90.0, iiSubSetCnt)
  endif else begin             ; North hemisphere
    iiSubSet = where ( dataGEI_Orig.coLat_deg le 90.0, iiSubSetCnt )
  end

  wn = 0
  window, wn, xSize = 850, ySize = 400, title='Original dB Data'
  !p.multi = [0,2,1]
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

  minth = 1.0     ; in deg
  nth = nLatGrid
  max_coLat = aacgm_cap_coLat_deg

  reg_colat_deg = transpose(rebin(findgen(nth)/(nth-1)* $
      (max_coLat-minth)+minth,nth,nLongrid))    ; deg
  if south eq 1 then reg_coLat_deg[*,*] = 180.0 - reg_coLat_deg[*,*]

  reg_lon_deg = rebin(findgen(nLongrid)/(nLongrid)*24.*15.,nLongrid,nth)

  regX  = reg_coLat_deg*cos(reg_lon_deg*!dtor-!dpi) - Xavg
  regY  = reg_coLat_deg*sin(reg_lon_deg*!dtor-!dpi) - Yavg

  coLat_grid  = sqrt(regX^2+regY^2)                   ; deg, reg_coLat_irid  (shifted)
  mlt_grid  = (!dpi + atan(regY,regX))*!radeg/15.0    ; MLT, reg_lon_irid



; SECS section here

  polelimit1 = 50.     ; 50 in orig code
  polelimit2 = 50.
  mlt_data = (dataGEO.lon_deg/15.0 mod 24)
  nXY  = 15
  noPlot = 0

  print,'Calculating SECS expansion...'
  jPar_SECS, dataGEO.coLat_deg, mlt_data, $
    dataGEO.bphi, coLat_grid, mlt_grid, $
    coLat_grid, mlt_grid, xAvg, yAvg, $
    NOPLOT = noPlot, JPARSECS = jPar, $
    POLELIMIT1 = poleLimit1, $
    POLELIMIT2 = poleLimit2, $
    NXY = nXY, $
    BTHSECS = bThFit_mag, $
    BPHSECS = bPhFit_mag, $
    UNITTH = unitTH, UNITPH = unitPh, $
    WEIGHT = weight, HARD = hard, $
    BPHFIT_DATA = bPhFit_data

;stop

stop

  print,'Finished'
end
