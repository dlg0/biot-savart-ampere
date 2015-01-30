; Ver 2012_02
; IDL ver 8 or later
; Widget wrapper for AMPERE data fitting code - AMP_FIT using SECS basis functions
;
; C. L.  Waters
; Centre for Space Physics
; University of Newcastle
; Australia
; 
; Oct 2014
;
; Modifications:
; 
; -------------------------------------------------------------------------------------
;
; Widget event routines
Pro doParams_event,ev
 Common WidgBlk,Base1,Base2,HrSld,MnSld,ScSld,DneBut
 Common WFileBlk,path,datFileName,amp_data_path
 Common WParBlk,StHr,StMn,StSc,GetSc,mxclat,plt_png,HemSp,south,debug

 Widget_Control,ev.id,Get_UValue=UVal
 Case UVal of
 'gfi' : Begin
          Ttle='Select Input AMPERE File'
          datFileName=Dialog_PickFile(Title=Ttle,path=amp_data_path,filter='*.ncdf')
          Widget_Control,HrSld,Set_Value=StHr
          Widget_Control,MnSld,Set_Value=StMn
          Widget_Control,ScSld,Set_Value=StSc
          Widget_Control,HrSld,Sensitive=1
          Widget_Control,MnSld,Sensitive=1
          Widget_Control,ScSld,Sensitive=1
          Widget_Control,DneBut,Sensitive=1
         end
 'shr' : Widget_Control,ev.id,Get_Value=StHr
 'smn' : Widget_Control,ev.id,Get_Value=StMn
 'ssc' : Widget_Control,ev.id,Get_Value=StSc
 'asc' : Widget_Control,ev.id,Get_Value=GetSc
 'clt' : Widget_Control,ev.id,Get_Value=mxclat
 'pfc' : Widget_Control,ev.id,Get_Value=plt_png
 'dbg' : Widget_Control,ev.id,Get_Value=debug
 'hms' : Begin
          Widget_Control,ev.id,Get_Value=south
          If south eq 0 then HemSp='N'
          If south eq 1 then HemSp='S'
         end
 'dne' : Begin
          Widget_Control,ev.top,/DESTROY
          PrmFile = path + 'amp_fit_secs.par'
          OpenW,u1,PrmFile,/Get_Lun
          PrintF,u1,StHr,StMn,StSc,GetSc
          PrintF,u1,south,mxclat,plt_png,debug
          Free_Lun,u1
         end
 end
end
;
Pro doParams
  Common WidgBlk,Base1,Base2,HrSld,MnSld,ScSld,DneBut
  Common WFileBlk,path,datFileName,amp_data_path
 Common WParBlk,StHr,StMn,StSc,GetSc,mxclat,plt_png,HemSp,south,debug

  WXPos=1 & WYPos=10 & WXSz=200                ; Widget Placement
  Base1 = WIDGET_BASE(Title='AMPERE FAC Menu',XOFFSET=WXPos,YOFFSET=WYPos,/COLUMN,XSize=WXSz)
  ADFBut = Widget_Button(base1,value='Select a New File',UValue='gfi')
  Widget_Control, ADFBut,Sensitive=1

  PrmFile = path+'amp_fit_secs.par'
  OpenR,u1,PrmFile,/Get_Lun
  StHr=8 & StMn=0 & StSc=0 & GetSc=3600
  ReadF,u1,StHr,StMn,StSc,GetSc
  south=0 & mxclat=50. & plt_png=1
  ReadF,u1,south,mxclat,plt_png,debug
  debug=fix(debug)
  Free_Lun,u1

  HrSld = Widget_Slider(base1,value=StHr,UValue='shr',Minimum=0,Maximum=23,Title='Start Hour')
  Widget_Control,HrSld,Sensitive=0
  MnSld = Widget_Slider(base1,value=StMn,UValue='smn',Minimum=0,Maximum=59,Title='Start Minute')
  Widget_Control,MnSld,Sensitive=0
  ScSld = Widget_SLider(base1,value=StSc,UValue='ssc',Minimum=0,Maximum=59,Title='Start Second')
  Widget_Control,ScSld,Sensitive=0
  AvSld = Widget_SLider(base1,value=GetSc,UValue='asc',Minimum=15,Maximum=7200,Title='Averaging Time (seconds)')
  DneBut = Widget_Button(base1,value='GO FOR IT',UValue='dne')
  Widget_Control,DneBut,Sensitive=1
  WIDGET_CONTROL, base1, /REALIZE

  Base2 = Widget_Base(Title='AMPERE FAC Menu',Group_Leader=Base1,XOffset=WXPos+WXSz+8,YOffset=WYPos,/Column,XSize=WXSz)
  pltSld= Widget_SLider(base2,value=mxclat,UValue='clt',Minimum=20,Maximum=90,Title='Cap Size [Deg]')

  fit_clat=90.0
;  SpacRes = 2.0*fit_clat/KMax      ; Resolution

  z1Arr=StrArr(2)
  For i=0,1 do z1Arr(i)=StrTrim(String(i),2)
  PFac = CW_BGroup(base2,z1Arr,UValue='pfc',Label_top='PNG FAC Output',Set_Value=plt_png,/row,/Exclusive)
;
  hArr=StrArr(2)
  hArr(0)='N' & hArr(1)='S'
  HmSB = CW_BGroup(base2,hArr,UValue='hms',Label_top='Hemisphere',Set_Value=south,/row,/Exclusive)
  HemSp=hArr(south)
;
  dbgch=StrArr(2)
  For i=0,1 do dbgch(i)=StrTrim(String(i),2)
  dbSB = CW_BGroup(base2,dbgch,UValue='dbg',Label_top='Debug',Set_Value=debug,/row,/Exclusive)
;
  WIDGET_CONTROL, base2, /REALIZE
  XMANAGER, 'doParams',base1,/NO_BLOCK
  XMANAGER, 'doParams',base2
end
;
; ------------------------------------------------------------------------------------------------------------
;
; Main driver code starts here
Pro amp_widget_secs
; cpu, tpool_nthreads=4, tpool_min_elts=1000
  cpu, tpool_nthreads=2, tpool_min_elts=1000
  Common WFileBlk,path,datFileName,amp_data_path
  Common WParBlk,StHr,StMn,StSc,GetSc,mxclat,plt_png,HemSp,south,debug

; set default values
  @d:\cwac\secs_olaf\amp_path_secs.pro
  
  if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN64') then begin
	 plotDev	= 'win'
;	path	= 'd:\cwac\hi_res\davidg\'
  endif else begin                        ; other systems (linux, darwin etc.)
	 plotDev = 'X'
;	path = '~/code/ampereFit/idl/'
  endelse
  set_plot,plotdev
  loadct,0
  device,decomposed=0
  mn_fac=0.1              ; do not plot abs(FAC) lt than this
  mx_fac=1.0

  doParams
  nLatGrid	= mxclat      ; fit grid defaults
  aacgm_cap_coLat_deg = mxclat
  nLonGrid  = 24
;  nLonGrid  = 48

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

  st_run_time = systime(/seconds)
  status = 0

  ampfit_secs, sHr, eHr, south, $
		datFileName = datFileName, $
		nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
		mn_fac = mn_fac, mx_fac=mx_fac, $
		plt_png = plt_png, $
		aacgm_cap_coLat_deg = aacgm_cap_coLat_deg, $
		debug=debug, $
		sepoch = sepoch, eepoch = eepoch, $
    extend = extend, $
    SpacRes = SpacRes, status=status

  en_run_time = systime(/seconds)
  print,'Time for run : ',en_run_time-st_run_time,' sec'

  print,'Finished'
 end
