; - Spherical coordinate system is [R,T,P]
; - Inputs are in deg
; - T coordinate is coLat

pro ampfit_create_bfn, _pole_R, _pole_T, _pole_P, _grid_R, _grid_T, _grid_P, _wire_R_m, $
		bR, bT, bP 

	;_pole_T = 0.01 

	u0 = 1.2566370614e-6
	I = 3e15

	geopack_sphcar_08, _pole_R, _pole_T, _pole_P, pX, pY, pZ, /degree, /to_rect	
	geopack_sphcar_08, _grid_R, _grid_T, _grid_P, gX, gY, gZ, /degree, /to_rect	

	; Get the perpendicular distance from the wire to the grid points

	pMag = sqrt(pX^2 + pY^2 + pZ^2)	
	gMag = sqrt(gX^2 + gY^2 + gZ^2)	

	; Dot product of the pole and grid position vectors to get the angle between them

	p_dot_g = pX * gX + pY * gY + pZ * gZ
	cosTh = p_dot_g / (pMag*gMag)
	th_rad = acos((cosTh<1)>(-1))

	; Perpendicular distance (r) from wire to grid point

	r = _pole_R * sin(th_rad) 

	iiInsideWire = where(r le _wire_R_m, iiInsideWireCnt)

	bMag = u0 * I / (2*!pi*r)
	;bWireEdge = u0*I/(2*!pi*_wire_R_m)
	bMag[iiInsideWire] = (u0*I*r[iiInsideWire]/(2*!pi*_wire_R_m^2));>bWireEdge*0.01

	; b unit vector [XYZ]

	gxU = gX / gMag  
	gyU = gY / gMag
	gzU = gZ / gMag
	
	pxU = pX / pMag  
	pyU = pY / pMag
	pzU = pZ / pMag

	; cU = pU x gU

	cxU = +(pyU*gzU - pzU*gyU)
	cyU = -(pxU*gzU - pzU*gxU)
   	czU = +(pxU*gyU - pyU*gxU)	

	; dU = cU x pU

	dxU = +(cyU*pzU - czU*pyU)
	dyU = -(cxU*pzU - czU*pxU)
   	dzU = +(cxU*pyU - cyU*pxU)	

	; bU = -(pU x dU)

	bxU = -(+(pyU*dzU - pzU*dyU))
	byU = -(-(pxU*dzU - pzU*dxU))
   	bzU = -(+(pxU*dyU - pyU*dxU))	

	; normalize length

	bUMag = sqrt(bxU^2+byU^2+bzU^2)

	bxU = bxU / bUMag 
	byU = byU / bUMag 
	bzU = bzU / bUMag 

	iiZero = where(bxU ne bxU,iiZeroCnt) 
	if iiZeroCnt gt 0 then begin
		bxU[iiZero] = 0 
		byU[iiZero] = 0 
		bzU[iiZero] = 0 
	endif
	
	; b = bU * bMag

	bx = bxU * bMag
	by = byU * bMag
	bz = bzU * bMag

	; Convert b vector from XYZ back to RTP

	geopack_bcarsp_08, gX,gY,gZ, bX,bY,bZ, bR,bT,bP
	;bT = -bT

	iiBad = where(bR ne bR,iiBadCnt)
   	if iiBadCnt gt 0 then stop
	iiBad = where(bT ne bT,iiBadCnt)
   	if iiBadCnt gt 0 then stop
	iiBad = where(bP ne bP,iiBadCnt)
   	if iiBadCnt gt 0 then stop



end
