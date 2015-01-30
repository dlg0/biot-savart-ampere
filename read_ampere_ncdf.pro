pro read_ampere_ncdf,ncdfname,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice
;
; Haje Korth
; JHU/APL
; June 2010
;
; error handling.

	;catch, theerror
	;if theerror ne 0 then begin
	;	catch, /cancel
	;	obj_destroy, fileobj
	;	print, 'Error message: ', !ERROR_STATE.MSG
	;	Print,'Error in netcdf file read'
	;	return
	;endif

;; initial definition of objects.
;	fileobj = obj_new()
;
;; open the source file in read-only mode.
;	fileobj = obj_new('ncdf_file', ncdfname)
;	if obj_valid(fileobj) eq 0 then message, 'Invalid file object returned from ncdf_file init.'
;
;; load variables
;	x_axis_frac_hour = fileobj -> getvardata('time')
;	pseudosvnum_total = fileobj -> getvardata('pseudo_sv_num')
;	plane_number_total = fileobj -> getvardata('plane_num')
;	pos_eci_total = fileobj -> getvardata('pos_eci')
;	b_eci = fileobj -> getvardata('b_eci')
;	pseudo_sv_quality = fileobj -> getvardata('pseudo_sv_quality')
;	data_splice = fileobj -> getvardata('data_splice')
;
;; destroy the file object.
;	obj_destroy, fileobj

	n = ncdf_open(ncdfname)
	ncdf_varget, n, 'time', x_axis_frac_hour
	ncdf_varget, n, 'pseudo_sv_num', pseudosvnum_total
	ncdf_varget, n, 'plane_num', plane_number_total
	ncdf_varget, n, 'pos_eci', pos_eci_total
	ncdf_varget, n, 'b_eci', b_eci
	ncdf_varget, n, 'pseudo_sv_quality', pseudo_sv_quality
	ncdf_varget, n, 'data_splice', data_splice


;; Test data from Haje
;  fname='d:\cwac\haje\20100216Amp_invert.sav'
;  restore,fname
;  pseudo_sv_quality = fltarr(n_elements(x_axis_frac_hour))
;  pseudo_sv_quality (*) = 1.0

; netcdf inquire sequence - CLW
; id=ncdf_open(ncdfname)
; res=ncdf_inquire(id)
; help,res
; print,ncdf_varinq(id,ii)
;
	return
end
