; - Spherical coordinate system is [R,T,P]
; - Inputs are in deg
; - T coordinate is coLat

pro ampfit_create_arb_bfn, _pole_R, _pole_T, _pole_P, _grid_R, _grid_T, _grid_P, _wire_R_m, $
		bR, bT, bP 

	;_pole_T = 0.01 

	u0 = 1.2566370614e-6
	J = 1.0

	nPts = n_elements(_grid_R)

	geopack_sphcar_08, _pole_R, _pole_T, _pole_P, pX, pY, pZ, /degree, /to_rect	
	geopack_sphcar_08, _grid_R, _grid_T, _grid_P, gX, gY, gZ, /degree, /to_rect	

	p = [pX,pY,pZ]
	pU = p/mag(p)

	; Get two directions x' & y' that define a plane perp to the
	; radial position vector of the pole

	_a = cross(pU,pU+[1.0,0,0])
	_aU = _a/mag(_a)

	_x = cross(pU,_aU)
	_xU = _x/mag(_x)

	_y = cross(pU,_xU)
	_yU = _y/mag(_y)

	; Create a regtangular mesh over which to perform the dS Biot-Savart integral

	width = _wire_R_m ; gaussian variance	

	nX = 21
	nY = 21
	xPrimeGrid = 8*width*(FIndGen(nX)/(nX-1)-0.5)
	yPrimeGrid = 8*width*(FIndGen(nY)/(nY-1)-0.5)
	dX = xPrimeGrid[1]-xPrimeGrid[0] 
	dY = yPrimeGrid[1]-yPrimeGrid[0]
	dS = dX*dY

	xP = (rebin(xPrimeGrid,nX,nY))[*]
	yP = (transpose(rebin(yPrimeGrid,nY,nX)))[*]

	_x = fltArr(nX*nY)
	_y = fltArr(nX*nY)
	_z = fltArr(nX*nY)

	bx = fltArr(nPts)
	by = fltArr(nPts)
	bz = fltArr(nPts)

	; Generate 2-D plane in R3 over which to define J guassian
	for _i = 0,nX*nY-1 do begin

			; dU is +x
			; cU is +y		

			_x[_i] = pX + xP[_i]*_xU[0] 
			_y[_i] = pY + xP[_i]*_xU[1]
			_z[_i] = pZ + xP[_i]*_xU[2]
		
			_x[_i] = _x[_i] + yP[_i]*_yU[0] 
			_y[_i] = _y[_i] + yP[_i]*_yU[1]
			_z[_i] = _z[_i] + yP[_i]*_yU[2]

	endfor

	_rPMag = sqrt(xP^2+yP^2)

	for _g = 0,nPts-1 do begin

		_rMag = (sqrt((gX[_g]-_x)^2+(gY[_g]-_y)^2+(gZ[_g]-_z)^2))>(sqrt(dx^2+dy^2))
	
		_rxU = (gX[_g]-_x)/_rMag	
		_ryU = (gY[_g]-_y)/_rMag	
		_rzU = (gZ[_g]-_z)/_rMag	

		jMag = J * exp(-(_rPMag)^2/(2*width^2))

		jX = pU[0] * jMag
		jY = pU[1] * jMag
		jZ = pU[2] * jMag

		; Sum over all the current sources

		bX[_g] = u0/(2*!pi)*dS*total((jY*_rzU-jZ*_ryU)/_rMag)
		bY[_g] = u0/(2*!pi)*dS*total(-(jX*_rzU-jZ*_rxU)/_rMag)
		bZ[_g] = u0/(2*!pi)*dS*total((jX*_ryU-jY*_rxU)/_rMag)

	endfor

	; Convert b vector from XYZ back to RTP

	geopack_bcarsp_08, gX,gY,gZ, bX,bY,bZ, bR,bT,bP

	iiBad = where(bR ne bR,iiBadCnt)
   	if iiBadCnt gt 0 then stop
	iiBad = where(bT ne bT,iiBadCnt)
   	if iiBadCnt gt 0 then stop
	iiBad = where(bP ne bP,iiBadCnt)
   	if iiBadCnt gt 0 then stop

end
