pro ampfit_secs_test
tic=SysTime(1,/sec)
dlm_register, 'aacgm/v3.8/idl_aacgm.dlm'
dlm_load, 'aacgm'

south = 1
datFileName = 'data/20100824Amp_invert.ncdf'
nLatGrid = 20.0
nLonGrid = 48.0
mn_fac = 0.1
mx_fac = 1.0
plt_png = 1
aacgm_cap_colat_deg = 50.0
debug = 1
SpacRes = 0
status = 0

StHr = 5
StMn = 20 
StSc = 0
GetSc = 600.0

datestr=strmid(file_basename(datFileName),0,8)
syr=fix(strmid(datestr,0,4))
smo=fix(strmid(datestr,4,2))
sday=fix(strmid(datestr,6,2))
geopack_epoch,sepoch,syr,smo,sday,StHr,StMn,StSc,/compute
eepoch=sepoch+GetSc*1.d3
geopack_epoch,eepoch,eyr,emo,eday,enHr,enMn,enSc,/breakdown

sHr = float(StHr)+float(StMn)/60.0+StSc/3600.0
eHr = float(enHr)+float(enMn)/60.0+enSc/3600.0
extend = 0                    ; default is sday=eday
if (eday ne sday) then begin
  extend = 1
  eHr = eHr+24.0
end

ampfit_secs, sHr, eHr, south, $
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

toc=SysTime(1,/sec)
print, 'Runtime: ', toc-tic
end
