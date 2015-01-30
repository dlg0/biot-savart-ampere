; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Sort Iridium satellite track data and get sequential longitude order
; Input:
;  track intersection shifted, North or South hemisphere data structure
;
; Output:
;   replaces data_in with track sorted version
;   Trk_order -> array[6] of track number order
;   stLon_a -> array[6] of start lons for each track, which must be sequential for ghost routine to work
;   gcDist_a -> Great circle distance array
;
;  C.L. Waters
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
;  Comments:
;   Checked Sept 2011 [CLW]
;   Added SpacRes check for data gaps [CLW, May 2012]
;   Added trk_gap - record of gaps for each [hemis,trk]  [CLW, May 2012]
;   Added track sort refinement section and output gcDist array for all tracks [Sept 2012]
;
Pro amp_sort_trackdata, data_in, data_struc, south, SpacRes, trk_gap, Trk_order

; - - - debug switch - - - - 
  dbg = 0
;  - - - - - - - - - - - - -
  sortII = indGen ( n_elements(data_in.iPln) ) ; Num points in hemisphere of data
  cnt = long(0)
  stLon_a=fltarr(6)             ; Array to store start Lons
  data_in.typ = 0               ; Initialise data type with 'normal' tag

; Loop over 6 tracks of Hemis data (not full sphere)
  For trackNo = 0, 5 do begin
    iiThisTrack = where ((data_in.iPln eq trackNo), thisTrackCnt ) ; data for 1 track
    ii_s = sort(data_in[iiThisTrack].coLat_deg)  ; sort by coLat_deg
    if (south eq 1) then begin   ; Sth hemis
      strt_clat = data_in[iiThisTrack[ii_s[0]]].coLat_deg
      strt_lon  = data_in[iiThisTrack[ii_s[0]]].Lon_deg
    end else begin               ; Nth hemis
      strt_clat = data_in[iiThisTrack[ii_s[thisTrackCnt-1]]].coLat_deg
      strt_lon  = data_in[iiThisTrack[ii_s[thisTrackCnt-1]]].Lon_deg
    end      

; Calc distances from equator start point to all points along this track, in a hemisphere
    gcDistance=gc_dist( 90.-strt_clat, 90.-data_in[iiThisTrack].coLat_deg, $
                            strt_lon,      data_in[iiThisTrack].Lon_deg)
; Sort by great circle distances
    iiThisTrackSorted = sort ( gcDistance )
; get which lon is smallest of the strt:end of track
    lon_diff = data_in[iiThisTrack[iiThisTrackSorted[0]]].lon_deg - $
               data_in[iiThisTrack[iiThisTrackSorted[thisTrackCnt-1]]].lon_deg
; reverse the idx of iiThisTrackSorted if necessary
    if (lon_diff gt 0.0) then iiThisTrackSorted = reverse(iiThisTrackSorted)
    stLon_a[TrackNo] = data_in[iiThisTrack[iiThisTrackSorted[0]]].lon_deg
    sortII[cnt:cnt+thisTrackCnt-1] = iiThisTrack[iiThisTrackSorted]
    cnt += thisTrackCnt
  endfor     ; end of Loop Track_num

  trk_order=sort(stLon_a)               ; sorted by Longitude
  data_in = data_in[sortII]             ; replace with sorted structure
  sortII = !null
  
; Refine track sort data. Sometimes gcDist order gets messed at far ends of a track
; Populate gcD value in the data structure - after lon calc and piecewise sort
  n_segs = 5              ; break the track into n_segs
  For trackNo = 0, 5 do begin
    iiThisTrack = where((data_in.iPln eq trackNo), thisTrackCnt ) ; semi-sorted data for 1 track
; Get a copy of struc for this track
    tmp_struc = temporary ( replicate ( data_struc, thisTrackCnt ) )
    tmp_struc = data_in[iiThisTrack]
; Calc number of points per smaller segment
    npp_seg = fix(thisTrackCnt/n_segs)       ; num points per segment
    DistSeg_a = dblarr(npp_seg)              ; dist data for each segment
    n_left = thisTrackCnt - fix(n_segs*npp_seg)  ; left over points after n_segs*npp_seg points
;print,trackNo,npp_seg,n_left
    lat_arr = 90.0-data_in[iiThisTrack[0:npp_seg-1]].coLat_deg   ; get lats of 1st segment
    if (south eq 0 ) then mnv = min(lat_arr,iimn)
    if (south eq 1 ) then mnv = max(lat_arr,iimn)

    For nn = 0, n_segs-1 do begin  ; Loop over num of segments of the track
      sp_i = nn*npp_seg
      sp = sp_i
      strt_lat_deg = 90.0-data_in[iiThisTrack[iimn]].coLat_deg  ; start lat,lon of track
      strt_lon_deg = data_in[iiThisTrack[iimn]].Lon_deg
      if nn gt 0 then begin
        sp = sp_i-1
        strt_lat_deg = 90.0-data_in[iiThisTrack[sp]].coLat_deg  ; start lat,lon of this segment
        strt_lon_deg = data_in[iiThisTrack[sp]].Lon_deg
      end
;print,trackno,npp_seg,nn,sp_i,sp
      For np = 0, npp_seg-1 do begin
        gcD = gc_dist( strt_lat_deg, 90.0-data_in[iiThisTrack[sp_i+np]].coLat_deg, $
                       strt_lon_deg,      data_in[iiThisTrack[sp_i+np]].Lon_deg)
;        data_in[iiThisTrack[sp_i+np]].gcD = data_in[iiThisTrack[sp]].gcD + gcD
        tmp_struc[sp_i+np].gcD = tmp_struc[sp].gcD + gcD
        DistSeg_a[np] = gcD      ; Get a copy to sort
;  print,sp_i+np,gcD,tmp_struc[sp_i+np].gcD
      end
      ii_s = sort(DistSeg_a)      ; Sort the segment here
      For mm=0,npp_seg-1 do begin
;        data_in[iiThisTrack[sp_i+mm]] = data_in[iiThisTrack[sp_i+ii_s[mm]]]  ; store gcDistances for this track
        data_in[iiThisTrack[sp_i+mm]] = tmp_struc[sp_i+ii_s[mm]]  ; store gcDistances for this track
;print,sp_i+mm,data_in[iiThisTrack[sp_i+mm]].gcD
      end
    end                          ; of Loop over num segs

; Left over points at end of track sequence
    if n_left gt 0 then begin    ; deal with left over points
      DistSeg_a = dblarr(n_left)
      sp_i = n_segs*npp_seg
      sp = sp_i-1
      strt_lat_deg = 90.0-data_in[iiThisTrack[sp]].coLat_deg
      strt_lon_deg = data_in[iiThisTrack[sp]].Lon_deg
;print,trackNo,n_left,sp_i,sp,strt_lat_deg,strt_lon_deg
      For np = 0, n_left-1 do begin
        gcD = gc_dist( strt_lat_deg, 90.0-data_in[iiThisTrack[sp_i+np]].coLat_deg, $
                       strt_lon_deg,      data_in[iiThisTrack[sp_i+np]].Lon_deg)
;        data_in[iiThisTrack[sp_i+np]].gcD = data_in[iiThisTrack[sp]].gcD + gcD
        tmp_struc[sp_i+np].gcD = tmp_struc[sp].gcD + gcD
        DistSeg_a[np] = gcD      ; Get a copy to sort
;print,sp_i+np,90.0-data_in[iiThisTrack[sp_i+np]].coLat_deg,$
;      data_in[iiThisTrack[sp_i+np]].Lon_deg,gcd
      end
      ii_s = sort(DistSeg_a)
      For mm = 0, n_left-1 do begin
        data_in[iiThisTrack[sp_i+mm]] = tmp_struc[sp_i+ii_s[mm]]  ; store sorted structure
      end
    end    ; If n_left > 0
  end      ; Loop over tracks, trackNo

; Loop over Tracks to check for gaps
  For trackNo = 0, 5 do begin
; 
; DEBUG DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
if dbg eq 1 then begin
  ii_dbg=where(data_in.iPln eq trackNo)
  clat_a = data_in[ii_dbg].coLat_deg
  lon_a  = data_in[ii_dbg].Lon_deg
  ang=lon_a*!pi/180.0
  if (south eq 1) then clat_a=180.-clat_a     ; Put in Nth hemis for XY coord calc
  Xs=clat_a*cos(ang)
  Ys=clat_a*sin(ang)
; Define the test plot winodw
; Plot the cocentric circles ready for oplotting the ghost points
  device, decomposed = 0
  !p.multi = 0
  window, 15, xSize = 680, ySize = 720, $
    title='Track Num:'+strtrim(string(trackNo),2)+'  South = '+strtrim(string(south),2)
  cLatMin = 90.0             ; change lat extent here
  nlatlab = cLatMin / 10
  x = dblArr(nlatlab,180)
  y = dblArr(nlatlab,180)
  for i=0,nlatlab-1 do begin
    for j=0,179 do begin
      x(i,j)=float(i+1)*10.*cos(float(j)/90.*!pi)
      y(i,j)=float(i+1)*10.*sin(float(j)/90.*!pi)
    endfor
  endfor
  plot,x(nlatlab-1,*),y(nlatlab-1,*), /isotropic, $
    xmargin=[1,1], ymargin=[1,1], xstyle=5,ystyle=5, $
    background=255, color=0, linestyle=1, title = title
  for i=0,nlatlab-2 do oPlot,x(i,*),y(i,*), linestyle=1, color=0
  loadct, 1,/silent
  xyouts,Xs,Ys,strtrim(string(data_in[ii_dbg].iPln),2),color=0
stop
end
; DEBUG DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

    ii_sel=where(data_in.iPln eq trackNo)
;    deg_gap_arr = abs(data_in[ii_sel].coLat_deg - shift(data_in[ii_sel].coLat_deg,1))
;    miss = where(deg_gap_arr gt SpacRes)  ; coLat gap triggered for > Spatial Resolution
;    if (miss(0) gt -1) then begin
;      trk_gap[south,trackNo] = 1          ; set gap flag
;      for gg=0,n_elements(miss)-1 do begin
;        s_i=miss(gg)-1
;        if (s_i le 0) then s_i=1
;        e_i = miss(gg)+1
;        if (e_i ge (n_elements(ii_Sel)-1)) then e_i = n_elements(ii_sel)-1
;        print,'data gap (coLats)=',data_in[ii_sel[s_i:e_i]].coLat_deg
;      end
;    end

  end

end