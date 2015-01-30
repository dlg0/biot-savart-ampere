;+
; NAME:
;	plot_b_vectors
; PURPOSE:
;	Plot a vector field on a spherical cap
; EXPLANATION:
;	
; CALLING SEQUENCE:
;
; INPUT:
;
; OPTIONAL KEYWORD INPUT:
;	color_lines	Set this to get superdarn type plots.
; OUTPUT:
;
; COMMON BLOCKS:
;       None
; EXAMPLES:
;
; MODIFICATION HISTORY:
; 	Written by David L. Green (dgreen@studentmail.newcastle.edu.au),
;	Physics Department, University of Newcastle, Australia	
; 	Apr 2004
;
;;	Note:		Positive bns values are plotted equatorwar, it the theta direction
;;			is positive equatorward.
;-
;


pro plot_b_vectors, bew, bns, theta, phi, maxlat, minlat, $
      winno, overlay, bew_irid, bns_irid, theta_irid, phi_irid, AUTOSCALE=autoscale, $
      SCALEFAC=scalefac, PLOTALL = plotall, OUT2PS = out2ps, MULTIPLOT=multiplot, $
      WIDGET=widget, CONTINENTS=continents, STATLOCS=statlocs, NOBORDER=noborder, $
	 COLOR_LINES=color_lines, THINSDARN=thinsdarn, TITLE=title, ESMAG=esmag, BSMAG=bsmag, $
	 GREY=grey, COLOR2=color2, ZBUFFER=zbuffer, NO_CONT=no_cont, CONT_TIME=cont_time, iSMAG=ismag,$
	 COLOR1=color1, THICK1=thick1, HEADS1=heads1, VECLEN_BG=veclen_bg, NODATA=nodata, PLOT_YN=plot_yn, $
	 ALTERNATE_CTBL=alternate_ctbl, KEEP_CTBL=keep_ctbl, NOGRID=noGrid

titleSize = 10

;   @daves_paths
   @d:\cwac\hist_irid\davidg\idl\jpar_ver2\irid2fac_paths.pro
	
	IF KEYWORD_SET(COLOR1) THEN BEGIN
			IF N_ELEMENTS(color1) EQ 1 THEN arrow_color=FLTARR(N_ELEMENTS(bew))+color1 ELSE arrow_color=color1
	ENDIF ELSE arrow_color=fltarr(n_elements(bew))+4
	
	if keyword_set(thinsdarn) then begin
    bew_bak	= bew
		bns_bak	= bns
		theta_bak	= theta
		phi_bak	= phi
		
		theta	= transpose(rebin(dindgen(18)*1.75+.25,18,24))*!dtor
		phi		= rebin(dindgen(24),24, 18)*15.*!dtor
	
		fth	= (theta[0,*]-min(theta_bak[0,*])) / (max(theta_bak[0,*])-min(theta_bak[0,*]))*(n_elements(theta_bak[0,*])-1)
		fph	= (phi[*,0]-min(phi_bak[*,0])) / (max(phi_bak[*,0])-min(phi_bak[*,0]))*(n_elements(theta_bak[*,0])-1)

		bew	= interpolate(bew_bak, fph, fth, /grid)
		bns	= interpolate(bns_bak, fph, fth, /grid)
		if keyword_set(color1) then begin
			if n_elements(color1) gt 1 then arrow_color	= interpolate(arrow_color, fph, fth, /grid)
		endif
	endif
	 
	if overlay eq 1 then begin
		nels	= n_elements(theta_irid[*,0])*n_elements(theta_irid[0,*])
		theta_row_irid	= reform(theta_irid,nels)
		phi_row_irid	= reform(phi_irid,nels)
		bew_row_irid	= reform(bew_irid,nels)
		bns_row_irid	= reform(bns_irid,nels)
	endif

;	Create provision to plot only every second track
;	within a certain theta radius of the pole to
;	avoid the arrows becoming conjested.

	crowd_lat = 7.*!pi/180
	plot_yn = fltarr(n_elements(theta[*,0]), n_elements(theta[0,*]))
	for j=0,n_elements(theta[0,*])-1 do begin
		;for k=0,n_elements(theta[*,0])-1 do $
		;begin
		;	plot_yn(k,j) =	0
		;	if theta(k,j) lt crowd_lat then $
		;	begin
		;		if iseven(k) eq 1 then plot_yn(k,j) = 1
		;		;;if k eq 0 then plot_yn[k,j] = 0
		;	endif
		;endfor
		k=0
		while k lt n_elements(theta[*,0]) do begin
			if theta[k,j] lt crowd_lat then begin
				plot_yn[k,j]=1
			endif
			k	= k+2.
		endwhile
	endfor

	if keyword_set(plotall) then plot_yn[*,*] = 0

	nel	= n_elements(theta[*])
	theta_row	= reform(theta,nel)
	phi_row	= reform(phi,nel)
	bew_row	= reform(bew,nel)
	bns_row	= reform(bns,nel)
	
;	Convert coords to degrees and from colat to lat

	theta_row = 90. - theta_row*180/!pi
	phi_row = phi_row*180/!pi
	
	if overlay eq 1 then begin
		theta_row_irid = 90. - theta_row_irid*180/!pi
		phi_row_irid = phi_row_irid*180/!pi
	endif

	if keyword_set(out2ps) ne 1  and keyword_set(zbuffer) ne 1 then begin
		if strcmp(!version.os,'linux') eq 1 then set_plot, 'X' else set_plot, 'win'
		device, decomposed = 0
	endif

	loadct, 0, file = colortb_path, /silent
   
	if keyword_set(out2ps) ne 1 and keyword_set(zbuffer) ne 1 then begin
		if keyword_set(multiplot) then begin
			if !d.window ne winno then if keyword_set(widget) then wset, winno else $
				window, winno
		endif else begin
			if !d.window ne winno then if keyword_set(widget) then wset, winno else $
			window, winno
		endelse
	endif

	if keyword_set(multiplot) then begin
		if keyword_set(noborder) then begin
			if keyword_set(title) then $
				MAP_SET, /ortho, 90, 0, /isotropic,  TITLE = title, $
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noerase, /noborder, charsize=titleSize, clip=1 $
			else $
				MAP_SET, /ortho, 90, 0, /isotropic, $
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noerase, /noborder, charsize=titleSize, clip=1
		endif else begin
			if keyword_set(title) then $
				MAP_SET, /ortho, 90, 0, /isotropic, TITLE = title,$
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noerase, charsize=titleSize, clip=1 $
			else $
				MAP_SET, /ortho, 90, 0, /isotropic,  $
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noerase, charsize=titleSize, clip=1
		endelse
	endif else begin
		if keyword_set(noborder) then begin
			if keyword_set(title) then $
				MAP_SET, /ortho, 90, 0, /isotropic, TITLE = title,$
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noborder, charsize=titleSize, clip=1 $
			else $
				MAP_SET, /ortho, 90, 0, /isotropic, $
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, /noborder, charsize=titleSize, clip=1
		endif else begin
			if keyword_set(title) then $
				MAP_SET, /ortho, 90, 0, /isotropic,  TITLE = title,$
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, charsize=titleSize, clip=1 $
			else $
				MAP_SET, /ortho, 90, 0, /isotropic, $
				limit = [round(90.-maxlat),0.,round(90.-minlat),360.], color = 1, charsize=titleSize, clip=1
		endelse
	endelse

;	Create the blue ocean effect
	
	loadct, 0, file=colortb_path, /silent	
	if keyword_set(no_cont) ne 1 then polyfill, findgen(360),fltarr(360)+(90.-fix(maxlat)),color=8

;	Get Continent boundaries and plot polygons
	if keyword_set(no_cont) ne 1 then begin
		get_world_map_data, cont_time,aacgm_wmap_lat,aacgm_wmap_lon,$
				wmap_idx,geog_lat=geographic_wmap_lat,geog_lon=geographic_wmap_lon	
		;;get_world_map_data, [2001,11,01,04,00,00],aacgm_wmap_lat,aacgm_wmap_lon,$
		;;		wmap_idx,geog_lat=geographic_wmap_lat,geog_lon=geographic_wmap_lon	
		iipos	= 0
		for pfpf=0,n_elements(wmap_idx[*,0])-1 do begin
			polyfill, aacgm_wmap_lon[iipos:iipos+wmap_idx[pfpf]-1],$
							aacgm_wmap_lat[iipos:iipos+wmap_idx[pfpf]-1], color=15
			plots, aacgm_wmap_lon[iipos:iipos+wmap_idx[pfpf]-1],$
							aacgm_wmap_lat[iipos:iipos+wmap_idx[pfpf]-1], color=10

			iipos	= iipos+wmap_idx[pfpf]
		endfor
	endif

	;map_continents, color = 6, /fill, /hires
	;map_grid, latdel = 10, londel = 90, color = 1, /label
	
	mltlabels		= strarr(5)
	mltlabels[0]	= '';;'00'
	mltlabels[1]	= '';;'06'
	mltlabels[2]	= '';;'12'
	mltlabels[3]	= '';;'18'
	mltlabels[4]	= ''
	latlabels	= [strtrim(string(80-indgen(fix(maxlat)/10)*10),2)]
	latgvals	= 90-(indgen(fix(maxlat)/10)+1)*10
	lonlabvals		= 90.-maxlat
	;;latlabels[0]	= ''
	;;latlabels[1]	= '80'
	;;latlabels[2]	= '70'
	;;latlabels[3]	= '60'
	;;latlabels[4]	= '50'


	if keyword_set(out2ps) then begin
	
		if keyword_set(noGrid) ne 1 then begin
		map_grid, lats=latgvals, lons=[0,90,180,270,360], /label, $
 			lonnames=mltlabels, latlab=45, charsize=.7,lonlab=lonlabvals, latalign=0.5, $
			lonalign=0.5,latnames=latlabels, charthick=2.,glinestyle=0, color=9, glinethick=1. 
		map_grid, lats=latgvals, lons=[0,90,180,270,360], /label, $
 			lonnames=mltlabels, latlab=45, charsize=.7,lonlab=lonlabvals, latalign=0.5, $
			lonalign=0.5,latnames=latlabels, charthick=2.,glinestyle=0, color=1, glinethick=2.,$
			/no_grid
		endif
	endif else begin
		if keyword_set(noGrid) ne 1 then begin
		map_grid, lats=latgvals, lons=[0,90,180,270,360], /label, color=1, $
 			lonnames=mltlabels, latlab=45, charsize=2.,lonlab=lonlabvals, latalign=0.5, $
			lonalign=0.5,latnames=latlabels, charthick=1.
		endif
	endelse

	if keyword_set(noGrid) ne 1 then begin

		xyouts, 0., (90.-fix(maxlat)/10*10)*.99, '00', color=1, align=0.5, charthick=2., charsize=.7
		xyouts, 89., (90.-fix(maxlat)/10*10)*1., '06', color=1, align=0.5, charthick=2., charsize=.7
		xyouts, 180., (90.-fix(maxlat)/10*10)*1.01, '12', color=1, align=0.5, charthick=2., charsize=.7
		xyouts, 271., (90.-fix(maxlat)/10*10)*1., '18', color=1, align=0.5, charthick=2., charsize=.7

	endif	
	;;stop
	if keyword_set(continents) then map_continents, color=1, /hires
;	Set the max length of the arrows
	if keyword_set(autoscale) then $
	begin
		nan_check_ew	= where(bew ne bew,nan_check_ew_c)
		nan_check_ns	= where(bns ne bns,nan_check_ns_c)
		if nan_check_ew_c gt 0 or nan_check_ns_c gt 0 then $
		begin
			bew_good	= bew[where(bew eq bew)]
			bns_good	= bns[where(bns eq bns)]

			maxval_fit = [max(abs(bew_good)),max(abs(bns_good))]
			if overlay eq 1 then $
			begin
				maxval_irid = [max(abs(bew_irid)),max(abs(bns_irid))]
				maxval = [maxval_fit, maxval_irid]
			endif else $
			begin
				maxval = [maxval_fit, 0.]
			endelse
			maxvl = max(maxval)
			maxarrowsize = 80.
			scalef = maxvl/maxarrowsize
		endif else $
		begin
	
			maxval_fit = [max(abs(bew)),max(abs(bns))]
			if overlay eq 1 then $
			begin
				maxval_irid = [max(abs(bew_irid)),max(abs(bns_irid))]
				maxval = [maxval_fit, maxval_irid]
			endif else $
			begin
				maxval = [maxval_fit, 0.]
			endelse
			maxvl = max(maxval)
			maxarrowsize = 50.
			scalef = maxvl/maxarrowsize
			;print, scalef
		endelse
	endif else $
	begin
		if keyword_set(scalefac) then scalef = scalefac else scalef=15.
	endelse
	
	if keyword_set(out2ps) then $
	begin
		psfac=18.	;	Set this number equal to the xsize/ysize of the device keywords for the ps device (for square)
				;	18 for making movies!!
		bns_row	= temporary(bns_row)*psfac
		bew_row	= temporary(bew_row)*psfac
		if overlay eq 1 then $
		begin
			bns_row_irid	= temporary(bns_row_irid)*psfac
			bew_row_irid	= temporary(bew_row_irid)*psfac
		endif
	endif


;	Color Lines Section

	if keyword_set(color_lines) then begin
		bmag	= sqrt(bns_row^2 + bew_row^2)
		;bmag_color	= bytscl(bmag,min=0.,max=.04,top=253);+1.
		;	if keyword_set(minrange) then iridium_color, array, tv_array, minin=minrange, maxin=maxrange $
		;iridium_color, bmag, bmag_color, minin=0., maxin=scalef*100.
		;print, bmag
		;bmag_color	= bmag_color > 1		;	This is here for Haje :-P
		;bmag_color	= bmag_color < 220
		bmag_color	= bmag/scalef <180
	;	bmag_color	= (((bmag-20)>0)/35)<200
 
		;print, bmag_color
		loadct,3;,file=colortb_path, /silent
		;stop
	endif
;}

;	How long is 0.05V/m in pixels?
	
	if keyword_set(out2ps) then Escalelen	= 0.05/scalef*psfac else $
		Escalelen	= 0.05/scalef
	xaya		= convert_coord(45.,80.,/data, /to_device)
	phib		= 90.-45.
	delx		= Escalelen*sin(phib*!pi/180)
	dely		= Escalelen*cos(phib*!pi/180)
	xbyb		=[xaya[0]+delx,xaya[1]+dely]
	ESmag		= sqrt((xbyb[0]-xaya[0])^2+(xbyb[1]-xaya[1])^2)

;	How long is 500nT (veclen_bg) in pixels?
	
	if keyword_set(veclen_bg) then veclen_bg=veclen_bg else veclen_bg = 500
	
	if keyword_set(out2ps) then bscalelen	= veclen_bg/scalef*psfac else $
		bscalelen	= veclen_bg/scalef
	xaya		= convert_coord(45.,80.,/data, /to_device)
	phib		= 90.-45.
	delx		= bscalelen*sin(phib*!pi/180)
	dely		= bscalelen*cos(phib*!pi/180)
	xbyb		=[xaya[0]+delx,xaya[1]+dely]
	BSmag		= sqrt((xbyb[0]-xaya[0])^2+(xbyb[1]-xaya[1])^2)	
	
;	How long is 400mAm^-1 in pixels?
	
	if keyword_set(out2ps) then iscalelen	= 400./scalef*psfac else $
	iscalelen	= 400./scalef
	xaya		= convert_coord(45.,80.,/data, /to_device)
	phib		= 90.-45.
	delx		= iscalelen*sin(phib*!pi/180)
	dely		= iscalelen*cos(phib*!pi/180)
	xbyb		=[xaya[0]+delx,xaya[1]+dely]
	iSmag		= sqrt((xbyb[0]-xaya[0])^2+(xbyb[1]-xaya[1])^2)	


;{;;	Set the arrow color
	
		if keyword_set(color1) then arrow_color	= color1 $
				else arrow_color	= 4
		if n_elements(arrow_color[*]) ne nel then arrow_color=fltarr(nel)+arrow_color
;}

;{;;	Set the arrow thickness and head size
	
		if keyword_set(thick1) then begin
			arrow_thick	= thick1
		endif else begin
			arrow_thick	= 4
		endelse

		if keyword_set(heads1) then begin
			arrow_head	= heads1
			if n_elements(arrow_head[*]) ne nel then arrow_head = fltarr(nel)+arrow_head
		endif else begin
			if keyword_set(out2ps) then begin
				arrow_head	= fltarr(nel)+80
			endif else begin
				arrow_head = fltarr(nel)+4
			endelse
		endelse
;}

if keyword_set(keep_ctbl) ne 1 then loadct, 0	
if KEYWORD_SET(alternate_ctbl) then loadct, alternate_ctbl;, bottom=1, ncolors=254-50
if keyword_set(nodata) ne 1 then begin
	for k=0,nel-1 do $
	begin
		if theta_row(k) ge 90-maxlat and theta_row(k) le 90-minlat and theta_row(k) lt 90. and plot_yn(k) eq 0 then $
		begin
			NS			=	bns_row(k)/scalef
			ns = -ns
			EW			=	bew_row(k)/scalef
		
			x1y1		=	convert_coord([phi_row(k),theta_row(k)], /data , /to_device)
			dely1		=	NS*cos(phi_row(k)*!pi/180)
			delx1		=	NS*sin(phi_row(k)*!pi/180)
			phib		=	90-phi_row(k)
			dely2		=	EW*cos(phib*!pi/180)
			delx2		=	EW*sin(phib*!pi/180)
			delx		=	-delx1 + delx2
			dely		=	dely1	+	dely2
			x2y2		=	[x1y1(0)+delx,x1y1(1)+dely]
			;print, sqrt((x2y2[0]-x1y1[0])^2+(x2y2[1]-x1y1[1])^2), bns_row[k], bew_row[k]
			th2phi2	=	convert_coord([x2y2(0),x2y2(1)], /device, /to_data)
			
			if keyword_set(out2ps) then $
			begin
				if keyword_set(color_lines) then $
				begin
					arrow,phi_row(k),theta_row(k),th2phi2(0),th2phi2(1), /data, color = bmag_color[k], thick = 1, hsize = 0
				endif else $
				begin
					if keyword_set(grey) then $
						arrow,phi_row(k),theta_row(k),th2phi2(0),th2phi2(1), /data, $
							color = 1, thick = .1, hsize = arrow_head[k] $
					else $
						arrow,phi_row(k),theta_row(k),th2phi2(0),th2phi2(1), /data, $
							color = arrow_color[k], thick = arrow_thick, hsize = arrow_head[k]
					endelse
			endif else $
			begin
				if keyword_set(color_lines) then $
				begin
					arrow,phi_row(k),theta_row(k),th2phi2(0),th2phi2(1), /data,$
							color = bmag_color[k], thick = arrow_thick, hsize = arrow_head[k]
				endif else $
				begin
					arrow,phi_row(k),theta_row(k),th2phi2(0),th2phi2(1), /data, $
							color = arrow_color, thick = arrow_thick, hsize = arrow_head[k]
				endelse
			endelse

		endif
	endfor
	
	if keyword_set(color_lines) then $
	begin
		phi_row_sym	= phi_row[where(theta_row gt 90.-maxlat)]
		theta_row_sym	= theta_row[where(theta_row gt 90.-maxlat)]
		bmag_color	= bmag_color[where(theta_row gt 90.-maxlat)]
		if keyword_set(out2ps) then usersymbol,'circle',size=0.2,/fill else usersymbol,'circle',size=0.6,/fill
		plots,phi_row_sym,theta_row_sym,color=bmag_color,psym=8
	endif
	
;	For making movies the preferred settings are hsize=50, thick =1
	
if overlay eq 1 then begin
	for k=0,nels-1 do begin
		if theta_row_irid(k) ne 0 and bew_row_irid(k) ne 0 and $
		theta_row_irid(k) gt 90-maxlat and theta_row_irid(k) lt 90-minlat then begin
			NS		=	bns_row_irid(k)/scalef
			ns 		= -ns
			EW		=	bew_row_irid(k)/scalef
		
			x1y1		=	convert_coord([phi_row_irid(k),theta_row_irid(k)], /data , /to_device)
			dely1		=	NS*cos(phi_row_irid(k)*!pi/180)
			delx1		=	NS*sin(phi_row_irid(k)*!pi/180)
			phib		=	90-phi_row_irid(k)
			dely2		=	EW*cos(phib*!pi/180)
			delx2		=	EW*sin(phib*!pi/180)
			delx		=	-delx1 + delx2
			dely		=	dely1	+	dely2
			x2y2		=	[x1y1(0)+delx,x1y1(1)+dely]
			th2phi2	=	convert_coord([x2y2(0),x2y2(1)], /device, /to_data)
			if keyword_set(out2ps) then $
			begin
				if keyword_set(color2) then $
					arrow,phi_row_irid(k),theta_row_irid(k),th2phi2(0),th2phi2(1), /data, color = color2, thick = 1., hsize = 100 $
				else $
					arrow,phi_row_irid(k),theta_row_irid(k),th2phi2(0),th2phi2(1), /data, color = 2, thick = 2., hsize = 100
			endif else $
			begin				
				if keyword_set(color2) then $
					arrow,phi_row_irid(k),theta_row_irid(k),th2phi2(0),th2phi2(1), /data, color = color2, thick = 1, hsize = 4 $
				else $
					arrow,phi_row_irid(k),theta_row_irid(k),th2phi2(0),th2phi2(1), /data, color = 2, thick = 1, hsize = 4
			endelse
		endif
	endfor
endif
endif
	if keyword_set(thinsdarn) then begin
		bew	= bew_bak
		bns	= bns_bak
		theta	= theta_bak
		phi	= phi_bak
	endif
end
