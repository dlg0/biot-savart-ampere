pro plot_fac, jPar, pole, cap_coLat_deg, coLat_deg, lon_deg, $
		mlt_shift = mlt_shift, $
		title = title, $
		south = south

	if not keyword_set(mlt_shift) then mlt_shift = 0

	loadct, 0, /silent
	!p.background = 255

	if south then begin
			limit = [ pole, 0, pole + (180-cap_coLat_deg), 360 ]
	endif else begin
			limit = [ pole - cap_coLat_deg, 0, pole, 360 ]
	endelse

	; Make is so the FAC map is viewed from the northern hemispher
	; in all cases, i.e., for the SH you look _through_ the Earth.
	; This is consistent with the dB maps.

	reverse = 0
	if south then reverse = 2

   	map_set, pole, 0, 0, /ortho, /iso, $
			limit = limit,  /noborder, title = title, $
			/advance, yMargin = [0,3], color = 0, reverse=reverse

	facScale = 1.0
	lat = 90-coLat_deg
	lon = ((lon_deg/15.0+mlt_shift) mod 24 )*15.0

	; Just grid both north and south hemispheres
	; ------------------------------------------

	loadct, 0, /silent

	if cap_coLat_deg gt 90 then begin
		nLats = fix((180-cap_coLat_deg)/10)
		lats = 180-(fIndGen(nLats)+1)*10
		lats = -90+180-lats
	endif else begin
		nLats = fix((cap_coLat_deg)/10)
		lats = 90-(fIndGen(nLats)+1)*10
	endelse

	map_grid, label = 1, $
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360], $
			lonlab=10, $
			latLab = 45, $
			lats = lats, $
			color = 0

	; Contour positive and negative seperately.
	; This prevents zero from looking like one color.
	; Only a problem for filled contours.

	; Positive contours
	; -----------------

	nLevs = 20
	jLevels = (fIndGen(nLevs)+1)/nLevs * facScale
	colors  = (255-bytScl ( jLevels, top = 253 ) )

	loadct, 3, /silent

   	contour, JPar, Lon, Lat, /over, levels = jLevels, c_colors = colors, /cell_fill

	; Negative contours
	; -----------------

	loadct, 1, /silent
	ThisJpar = -jPar
   	contour, ThisJPar, Lon, Lat, $
			/over, levels = jLevels, $
           	c_colors = colors,  /cell_fill


	; Just grid both north and south hemispheres
	; ------------------------------------------

	loadct, 0, /silent

	if cap_coLat_deg gt 90 then begin
		nLats = fix((180-cap_coLat_deg)/10)
		lats = 180-(fIndGen(nLats)+1)*10
		lats = -90+180-lats
	endif else begin
		nLats = fix((cap_coLat_deg)/10)
		lats = 90-(fIndGen(nLats)+1)*10
	endelse

	lonLab = min(lats)+10
	if south then lonLab = max(lats)-10
	map_grid, label = 1, $
			;lonNames	= ['0/24','6','12','18',''], $
			lonNames	= ['0','90','180','270',''], $
			lons = [0,90,180,270,360], $
			lonlab = lonLab, $
			latLab = 45, $
			lats = lats, $
			color = 0
end
