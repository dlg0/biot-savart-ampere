; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Calc centre of Iridium satellite orbit planes:
;
; Input : data_in (data structure), both Nth and Sth hemis data
;         'south' -> if set, calc rotation using Sth hemisphere orbit plane intersection
;
; output : coLat and Lon of average track intersection point
;          diag=diag -> print diagnostic messages
;
; C.L. Waters
; Centre for Space Physics Research
; University of Newcastle, NSW,
; Australia
;
; Comments:
; Revised: Sept 2011
;
Pro amp_get_new_centre, data_in, south, CRLat, RLon, XAvg,YAvg, diag=diag
; If south = 1 then rotate to south hemisphere intersection else use Nth hemisphere
;
; Error check
  status=0
  catch,error_status
  if (error_status ne 0) then begin
    status=1
    Print,'## ERROR : amp_get_new_centre => returning...'
    catch,/cancel
    return
  endif
;
; Estimate the average track intersection location
; Step 1: Calc x,y coords (funky cylindrical) of Iridium data locations
  npl=6                     ; 6 orbit planes
  mnpnts=20                 ; Min number of points in a track in order to proceed

; Get a copy and subset of data_in
  if south eq 1 then ii_sel=where(data_in.coLat_deg gt 95.0) else $
                     ii_sel=where(data_in.coLat_deg lt 85.0)
  clat_a = data_in[ii_sel].coLat_deg                   ; co_Lat in degrees [0->180]
  sel_iPln = data_in[ii_sel].iPln
  if south eq 1 then clat_a=180.0-clat_a     ; put coLat origin for smaller values

; Calc x,y coords from co_lat long location data
  xcrd_a=clat_a*cos(data_in[ii_sel].lon_deg*!dpi/180.0)
  ycrd_a=clat_a*sin(data_in[ii_sel].lon_deg*!dpi/180.0)

; Count number of data points in each Iridium orbit track
  p_cnt_pl=intarr(npl)                    ; space to store num points per orbit plane
  For mp=0,npl-1 do begin                 ; Get Merid Planes in order, ipln runs from 0->5
    idx1=where(sel_iPln eq mp,cnt1)
      p_cnt_pl[mp]=cnt1
  end
  If keyword_set(diag) then Print,'Num. Points read from Input File=',n_elements(data_in)
  If keyword_set(diag) then Print,'Num. Points per Orbit track:',p_cnt_pl

; Check for sufficient data on all orbit tracks
  If (min(p_cnt_pl) le mnpnts) then begin
    status=1                                  ; Return if insufficient data on any track
    Print,'Insufficient data on Orbit Track'
    return
;  stop
  endif

; Each orbit track is assumed to be a parabola curve in x,y
; Calc the coeffs of the parabolic eqn of each orbit track (in x,y)
  pcoef=dblarr(NPl,3)
  t_cnt=0

  For mp=0,NPl-1 do begin
; WARNING : SVDFIT crashes if we get a straight track along the 0600-1800 line: CLW
    idx=where(sel_iPln eq mp)                          ; select out each orbit track
    res_a=SVDFIT(xcrd_a[idx],ycrd_a[idx],3,status=status)  ; parabola fit to track location data
;print,'t_cnt=',t_cnt
    pcoef(t_cnt,*)=res_a                                   ; store parabola eqn coefficients
    t_cnt++
    if status ne 0 then begin
      print,'##ERROR: ampere_get_new_centre - Track intersection solver, excluding track'
      t_cnt=t_cnt-1
    end
  end
;
; Now calculate the intersection coords of these orbit track parabolas
  n_com=total(indgen(t_cnt))
  XInt=dblarr(n_com)                  ; There are n_com orbit track intersection combinations
  YInt=XInt & IntCnt=0
; Loop through every combination of the 6 planes to find intersections

  For aa=0,t_cnt-1 do begin
    For bb=1,t_cnt-1 do begin
      If (bb gt aa) then begin
        adiff=pcoef(aa,2)-pcoef(bb,2)
        bdiff=pcoef(aa,1)-pcoef(bb,1)
        cdiff=pcoef(aa,0)-pcoef(bb,0)
        sqtrm=bdiff^2-4.*adiff*cdiff
        If (sqtrm ge 0.) then begin           ; trap for -ve in sqrt of quadratic eqn
          x1=0.5*(-bdiff+sqrt(bdiff^2-4.*adiff*cdiff))/adiff
          x2=0.5*(-bdiff-sqrt(bdiff^2-4.*adiff*cdiff))/adiff
          y1=pcoef(aa,2)*x1*x1+pcoef(aa,1)*x1+pcoef(aa,0)
          y2=pcoef(aa,2)*x2*x2+pcoef(aa,1)*x2+pcoef(aa,0)
          RDst1=sqrt(x1*x1+y1*y1)
          RDst2=sqrt(x2*x2+y2*y2)
          If (RDst2 gt RDst1) then begin   ; Find the correct quadratic solution
            XInt(IntCnt)=x1
            YInt(IntCnt)=y1
          end
          If (RDst1 gt RDst2) then begin
            XInt(IntCnt)=x2
            YInt(IntCnt)=y2
          end
          IntCnt++
        end
      end     ; if bb gt aa -> MPlane combinations
    end      ; bb MPlane Loop
  end       ; aa MPlane Loop

; Calc equally weighted, average location of track intersections
  XAvg=mean(XInt(0:IntCnt-1))
  YAvg=mean(YInt(0:IntCnt-1))

; ********** Turn off shift
; XAvg=0.0 & YAvg=0.0
; ************************

; Release some memory
  ii_sel=!null
  clat_a=!null
  sel_iPln=!null
  xcrd_a=!null
  ycrd_a=!null

  If keyword_set(diag) then begin
    Print,'In amp_get_new_centre: Coordinate Shift is [x,y]: ',XAvg,YAvg
    Print,'Intersect X Values: ',XInt
    Print,'Intersect Y Values: ',YInt
    Print,'PCoeffs for each Track [x^2, x, c] : '
    For aa=0,npl-1 do print,transpose(pcoef(aa,*))
  end
;
  CRLat=sqrt(XAvg*XAvg+YAvg*YAvg)             ; Find Lat,Lon corresponding to XAvg,YAvg
  RLon=atan(YAvg,XAvg)*180.0/!dpi

  If south then CRLat=180.0-CRLat             ; Report CRLat in correct Hemisphere
  If keyword_set(diag) then Print,'Shift in X,Y=',XAvg,YAvg
  If keyword_set(diag) then Print,'Shift in Lat,Lon [deg]=',CRLat,RLon

end
