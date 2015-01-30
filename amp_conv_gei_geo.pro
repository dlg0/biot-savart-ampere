; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Rotate GEI to GEO coords as follows:
; Coords : GEI(x,y,z)->GEO(x,y,z)->GEO(r,t,p)
; dB vecs : GEI_db(r,t,p)->GEO_db(x,y,z)
;
; C.L. Waters
; Centre for Space Physics Research
; University of Newcastle, NSW,
; Australia
;
; Checked Feb 2012
;
pro amp_conv_gei_geo, gei_in, geo_out, epch
;
; Initially set geo_out to gei_in
  geo_out = gei_in

; Conv GEI(x,y,z) -> GEO(x,y,z)
  tme = replicate(epch,n_elements(gei_in.px))
  nn = fix(n_elements(tme)/10000)
  if nn eq 0 then geoPack_conv_coord, $          ; handles up to 10000 points
    gei_in[*].px, gei_in[*].py, gei_in[*].pz, $
    geo_x, geo_y, geo_z, epoch = tme, /from_gei, /to_geo $
  else begin           ; high resolution data
    geo_x=dblarr(n_elements(tme))
    geo_y=dblarr(n_elements(tme))
    geo_z=dblarr(n_elements(tme))
    for ii=long(0),n_elements(tme)-1 do begin
      in_x=gei_in[ii].px
      in_y=gei_in[ii].py
      in_z=gei_in[ii].pz
      tt=tme[ii]
      geoPack_conv_coord, in_x, in_y, in_z, out_x, out_y, out_z, epoch = tt, /from_gei, /to_geo
      geo_x[ii]=out_x
      geo_y[ii]=out_y
      geo_z[ii]=out_z
    end
  end

 	tme = !null

; Conv GEO(x,y,z) - > GEO(r,t,p)
  geopack_sphcar, $
		geo_x, geo_y, geo_z, $
   	geo_rkm, geo_coLat_rad, geo_lon_rad, /to_sphere

  geo_out.px = geo_x
  geo_out.py = geo_y
  geo_out.pz = geo_z

  geo_x = !null
  geo_y = !null
  geo_z = !null

  geo_out.R_km = geo_rkm
  geo_out.coLat_deg = geo_coLat_rad*180.0/!dpi
  geo_out.lon_deg = geo_lon_rad*180.0/!dpi

; Conv dB GEI(r,t,p) - > dB GEO(x,y,z)
  geopack_bspcar, $
    geo_coLat_rad, geo_lon_rad,  $
    gei_in.br, gei_in.bTheta, gei_in.bPhi, $
    dbx, dby, dbz

  geo_rkm = !null
  geo_coLat_rad = !null
  geo_lon_rad = !null

  geo_out.dbx = dbx
  geo_out.dby = dby
  geo_out.dbz = dbz

  dbx = !null
  dby = !null
  dbz = !null
end