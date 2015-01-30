function drawGrid, MAXLAT = maxLat, $
	MAPCHARSIZE = mapCharSize;
	
	loadct, 12, /silent
	
	if keyword_set ( maxLat ) then begin
		maxLat	= maxLat; 
		maxLatPlot	= maxLat + 2.0;
	endif else begin
		maxLat = 40.0;
		maxLatPlot	= 42.0;
	endelse

	latlabels	= [strtrim(string(80-indgen(fix(maxLat)/10)*10),2)]
	latgvals	= 80-indgen(fix(maxLatPlot)/10+1)*10
	lonlabvals		= 90.-maxLatPlot
	mltlabels		= strarr(5)
	if NOT keyword_set ( mapCharSize ) then mapCharSize	= 1.0;

	map_grid, lats=latgvals, lons=[0,90,180,270,360], /label, $
		lonnames=mltlabels, latlab=45, charsize=.01,lonlab=lonlabvals, latalign=0.5, $
		lonalign=0.5,latnames=latlabels, charthick=1.,glinestyle=0, color=14*16-1, glinethick=1. 
	map_grid, lats=latgvals, lons=[0,90,180,270,360], /label, $
 		lonnames=mltlabels, latlab=45, charsize=mapCharSize,lonlab=lonlabvals, latalign=0.5, $
		lonalign=0.5,latnames=latlabels, charthick=1.,glinestyle=0, color=14*16-1, glinethick=2.,$
		/no_grid
	xyouts, 0., (90.-fix(maxLatPlot)/10*10)*.99, '00', $
		color=14*16-1, align=0.5, charthick=1., charsize=mapCharSize
	xyouts, 89., (90.-fix(maxLatPlot)/10*10)*1., '06', $
		color=14*16-1, align=0.5, charthick=1., charsize=mapCharSize
	xyouts, 180., (90.-fix(maxLatPlot)/10*10)*1.01, '12', $
		color=14*16-1, align=0.5, charthick=1., charsize=mapCharSize
	xyouts, 271., (90.-fix(maxLatPlot)/10*10)*1., '18', $
		color=14*16-1, align=0.5, charthick=1., charsize=mapCharsize

	return, 1;

end


