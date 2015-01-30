; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Reads AMPERE data file (input data)
;
; Calls read_ampere_ncdf routine
; Initialises GEOPACK using GEOPACK_reCalc
; Calls amp_trk_lab90 to check for stray satellites
;
; Inputs:
;   sHr, eHr -> start and end times of requested data interval (decimal UT)
;   datFileName -> filename of *.ncdf input data
;   extend -> switch to read next file if time segment requires (goes over to the next day)
;
; Outputs:
;   data_struc-> definition of data structure type
;   dataGEI_orig -> full sphere data from Lars for given time segment
;   year,month,day,avgYrSec,avgEpoch, -> time information
;
;  D.L. Green & C.L. Waters
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
; Comments:
;
pro read_GEI_data90, sHr, eHr, $
	datFileName = datFileName, $
	data_struc = data_struc, $
	dataGEI_orig = dataGEI_orig, $
	year = year, month = month, day = day, $
	avgYrSec = yrSecAvg, $
	avgEpoch = avgEpoch, $
	debug = debug, $
	extend = extend, $
	status=status

  forward_function cnvTime

; Trap for errors
  if debug eq 0 then begin
	status=0
	catch,error_status
	if (error_status ne 0) then begin
		status=1
		Print,'## ERROR in read_GEI_data90 => returning'
		catch,/cancel
		return
	endif
  endif

; Create data structure
; ---------------------
  data_struc = { $
     utc: 0d0, $                ; UT time in dec hours
     isat: 0, $                 ; coded SV number (for Haje)
     iPln: 0, $                 ; orbit track number (0->5)
     qual: 0d0, $               ; data quality from Lars
     splice : 0, $              ; flag for spliced data - where missing data estimated
     typ : 0, $                 ; data type (for CLW)
     gcD : 0d0, $               ; Great circle distance from equator
     px: 0d0, py: 0d0, pz: 0d0, $
     dbx: 0d0, dby: 0d0, dbz: 0d0, $
     R_km : 0d0, coLat_deg: 0d0, lon_deg: 0d0, $
     br: 0d0, bTheta: 0d0, bPhi: 0d0 }

; NetCDF read section
  read_ampere_ncdf,datFileName,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice

  dateStr=''
  reads, strmid(file_baseName ( datFileName ),0,8), year, month, day, format='(i4,i2,i2)'

  if keyword_set(extend) then begin
    geopack_epoch,epoch,year,month,day,/compute
    epoch=epoch+86400.0d3
    geopack_epoch,epoch,extyear,extmonth,extday,/breakdown
    extdatfilename=file_basename(datfilename)
    strput,extdatfilename,string(extyear,format='(i4.4)'),0
    strput,extdatfilename,string(extmonth,format='(i2.2)'),4
    strput,extdatfilename,string(extday,format='(i2.2)'),6
    extdatfilename=file_dirname(datfilename)+path_sep()+extdatfilename

    if file_test(extdatfilename) then begin
      t_x_axis_frac_hour=x_axis_frac_hour
      t_pseudosvnum_total=pseudosvnum_total
      t_plane_number_total=plane_number_total
      t_pos_eci_total=pos_eci_total
      t_b_eci=b_eci
      t_pseudo_sv_quality=pseudo_sv_quality
      t_data_splice=data_splice

      read_ampere_ncdf,extdatfilename,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice

      x_axis_frac_hour=[t_x_axis_frac_hour,x_axis_frac_hour+24]
      pseudosvnum_total=[t_pseudosvnum_total,pseudosvnum_total]
      plane_number_total=[t_plane_number_total,plane_number_total]
      pos_eci_total=[[t_pos_eci_total],[pos_eci_total]]
      b_eci=[[t_b_eci],[b_eci]]
      pseudo_sv_quality=[t_pseudo_sv_quality,pseudo_sv_quality]
      data_splice=[t_data_splice,data_splice]

    endif else begin
      status=1
      return
    endelse
  endif

; Select subset of data
; ---------------------
; Get full sphere of data
   iiSubSet = where ( x_axis_frac_hour ge sHr and x_axis_frac_hour le eHr, iiSubSetCnt)

; Fill data structure
; -------------------

  dataGEI_orig = temporary ( replicate ( data_struc, iiSubSetCnt ) )

  dataGEI_orig.utc = x_axis_frac_hour[iiSubSet]
  dataGEI_orig.iSat = pseudoSVNum_total[iiSubSet]
  dataGEI_orig.iPln = plane_number_total[iiSubSet]
  dataGEI_orig.qual = pseudo_sv_quality[iiSubSet]
  dataGEI_orig.splice = data_splice[iiSubSet]

	; GEI XYZ position in km

  dataGEI_orig.px = (pos_ECI_total[0,*])[iiSubSet]*1d-3
  dataGEI_orig.py = (pos_ECI_total[1,*])[iiSubSet]*1d-3
  dataGEI_orig.pz = (pos_ECI_total[2,*])[iiSubSet]*1d-3

	; GEI XYZ db vectors

  dataGEI_orig.dbx = (B_ECI[0,*])[iiSubSet]
  dataGEI_orig.dby = (B_ECI[1,*])[iiSubSet]
  dataGEI_orig.dbz = (B_ECI[2,*])[iiSubSet]

  ; Get spherical coords of the GEI XYZ locations
  ; ---------------------------------------------

  geopack_sphcar, dataGEI_orig.px, dataGEI_orig.py, dataGEI_orig.pz, $
                GEI_R_km, GEI_coLat_deg, GEI_lon_deg, $
                /to_sphere, /degree

; coLat is 0->90 for Nth and 90->180 for Sth (deg)
  dataGEI_orig.R_km = GEI_R_km
  dataGEI_orig.coLat_deg = GEI_coLat_deg
  dataGEI_orig.lon_deg = GEI_lon_deg

; Free this memory
  GEI_R_km = !null
  GEI_coLat_deg = !null
  GEI_lon_deg = !null

 ; Rotate XYZ GEI db to spherical GEI
 ; ----------------------------------

  geopack_bcarsp, dataGEI_orig.px,  dataGEI_orig.py,  dataGEI_orig.pz, $
                  dataGEI_orig.dbx, dataGEI_orig.dby, dataGEI_orig.dbz, $
                  br_GEI, bTheta_GEI, bPhi_GEI

  dataGEI_orig.br = br_GEI
  dataGEI_orig.bTheta = bTheta_GEI
  dataGEI_orig.bPhi = bPhi_GEI

  ; Free up some memory
  br_GEI = !null
  bTheta_GEI = !null
  bPhi_GEI = !null

  x_axis_frac_hour = !null       ; utc
  pseudosvnum_total = !null      ; iSat
  plane_number_total = !null     ; iPln
  pos_eci_total = !null          ; px, py, pz
  b_eci = !null                  ; dbx, dby, dbz
  pseudo_sv_quality = !null      ; data quality from Lars
  data_splice = !null            ; splice data flag

; Get the average epoch and times for GEI to GEO and AACGM conversion
; Initialise GEOPACK here (using GEOPack_reCalc)
; -----------------------------------
  print,'Year/Month/Day Hour:Minute : ',strtrim(string(fix(year)),2),'/',$
                                        strtrim(string(fix(month)),2),'/',$
                                        strtrim(string(fix(day)),2),'  ',$
                                        strtrim(string(fix(shr)),2),':',strtrim(string(round((shr mod 1)*60)),2)
  avgHour = fix ( (sHr + eHr) / 2.0 )
  avgMin  = ( (sHr + eHr) / 2.0 mod 1 ) * 60
  yrSecAvg = cnvTime ( year, month, day, avgHour, avgMin, 0.0 )
; must be called before GEI->GEO conversion or if the date/time changes
  geoPack_reCalc, year, month, day, avgHour, avgMin, /date
  cdf_epoch, epoch0, year, month, day ,/compute_epoch
  epoch = dataGEI_orig.utc*3600d0*1d3 + epoch0
  avgEpoch = mean(epoch)

end
