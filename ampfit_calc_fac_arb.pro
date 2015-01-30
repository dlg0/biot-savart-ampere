; - Spherical coordinate system is [R,T,P]
; - Inputs are in deg
; - T coordinate is coLat

pro ampfit_calc_fac_arb, _pole_R, _pole_T, _pole_P, _grid_R, _grid_T, _grid_P, _wire_R_m, coeff, $
		jPar, jParCnt 

	geopack_sphcar_08, _pole_R, _pole_T, _pole_P, pX, pY, pZ, /degree, /to_rect	
	geopack_sphcar_08, _grid_R, _grid_T, _grid_P, gX, gY, gZ, /degree, /to_rect	

	; Get the perpendicular distance from the wire to the grid points

	pMag = sqrt(pX^2 + pY^2 + pZ^2)	
	gMag = sqrt(gX^2 + gY^2 + gZ^2)	

	; Dot product of the pole and grid position vectors to get the angle between them

	p_dot_g = pX * gX + pY * gY + pZ * gZ
	cosTh = p_dot_g / (pMag*gMag)
	th_rad = acos(cosTh)

	; Perpendicular distance (r) from wire to grid point

	r = _pole_R * sin(th_rad) 

	J = 1.0
	width = _wire_R_m ; gaussian variance	
	jPar = jPar + coeff* J * exp(-(r)^2/(2*width^2))
	++jParCnt

end
