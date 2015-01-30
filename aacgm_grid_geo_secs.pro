; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Generate uniform AACGM grid for SECS basis function version
; D. L. Green and C.L. Waters
; Centre for Space Physics Research
; University of Newcastle, NSW,
; Australia
;
; Ver 1 : Dec, 2009
; Ver 2 : 0-90 deg fit - Sept 2011
; Ver 3 : Modified from aacgm_grid_GEO - Oct 2014
;
; to be finoshed - clw Nov 2014
;
pro aacgm_grid_GEO_secs, max_coLat, south, $       ; max_colat is +ve and in deg
  aacgmGrid_coLat_deg = aacgmGrid_coLat_deg, $
  aacgmGrid_lon_deg = aacgmGrid_lon_deg, $
  geoGrid_coLat_deg = geoGrid_coLat_deg, $
  geoGrid_lon_deg = geoGrid_lon_deg, $
  geoGrid_x = geoGrid_x, $
  geoGrid_y = geoGrid_y, $
  geoGrid_z = geoGrid_z, $
  nLat = nLat, nLon = nLon, $
  year = year, ev_epoch=ev_epoch, $
  mltShift = mltShift, $        ; output
  mltref=mltref                 ; usually zero, used to shift in longitude

;  @amp_fit_paths90

	; Default keyword variables
	; -------------------------

;  if (not keyword_set(nLat)) then nLat = calc_coLat
  if (not keyword_set(nLat)) then nLat = 20
  if (not keyword_set(nLon)) then nLon = 24
  if (not keyword_set(mltref)) then mltref=0.0

; Iridium altitude
; ----------------
  Re_km = 6371.0
  R_km = 780.0 + Re_km

; Generate uniform AACGM grid
; ---------------------------
  minth = 1.0     ; in deg
  aacgm_coLat_deg = transpose(rebin(findgen(nLat)/(nLat-1)*(max_coLat-minth)+minth,nLat,nLon)) ; deg
  if south eq 1 then aacgm_coLat_deg[*,*] = 180.0 - aacgm_coLat_deg[*,*]

  reg_lon_deg = rebin(findgen(nLon)/(nLon)*24.*15.,nLon,nLat)

;*************************
  yr=year
;  if yr gt 2010 then yr=2010        ; max year for aacgm coeffs
;  aacgm_set_path, aacgmpath         ; should be set in amp_fit_paths90
;  aacgm_load_coef, yr
; ************************
; Convert magnetic longitude to MLT
;  mltShift = aacgm_mlt(year, yrSec, mltref )  ; i.e get the MLT at 0 MLon
;  mltShift = aacgm_mlt(ev_epoch, mltref )  ; i.e get the MLT at 0 MLon
;;	mlon0 = aacgm_mlong(year,yrSec,mltref)    ; Calc magnetic longitude for MLTRef, usually 0 MLT
;;  aacgm_mlt_deg = fIndGen(nLon)*(360.0-360.0/nLon)/(nLon-1) - mltShift*15.0
  aacgm_mlt_deg = reg_lon_deg; - mltShift*15.0
  aacgm_lon_deg = (aacgm_mlt_deg + 360.0) mod 360
;  aacgmGrid_coLat_deg = rebin(aacgm_coLat_deg, nLat, nLon )
  aacgmGrid_coLat_deg = aacgm_coLat_deg
;  aacgmGrid_lon_deg = transpose (rebin(aacgm_lon_deg, nLon, nLat))
  aacgmGrid_lon_deg = aacgm_lon_deg
;stop
; Valid ranges of AACGM_CONV_COORD are:
;  -90 < LAT < +90
; -180 < LON < 360
;
  aacgmGrid_R_km = aacgmGrid_lon_deg * 0 + R_km

; Convert AACGM grid to GEO
; -------------------------
  aacgm_conv_coord, 90.0-aacgmGrid_coLat_deg, aacgmGrid_lon_deg, $
    aacgmGrid_R_km-Re_km, geoGrid_lat_deg, geoGrid_lon_deg, err, /to_geo

  geoGrid_coLat_deg = 90.0-geoGrid_lat_deg

; Convert to GEO XYZ locations (input theta, phi is radians)
; ------------------------------------------------------
  geoPack_sphCar,aacgmGrid_R_km,(90.0-geoGrid_lat_deg)*!dtor, geoGrid_lon_deg*!dtor, $
    geoGrid_x, geoGrid_y, geoGrid_z, /to_rect

end
