;	Program to compute the partial derivative
;	of any 1 or 2 dimensional array using the
;	"deriv.pro" IDL routine.
;
;	Inputs
;		array		:		Array to be differentiated
;		dir			:		Direction in which to differentiat,
;								ie, x(0) or y(1)
;		coords	:		The x/y values to be diff against, ie
;								dB/dx or dB/dy
;
;	Outputs
;		darray	:		Differentiated array the same size
;								as the original array.
;
; David Green
;	Newcastle Uni, Physics Dept.
; Callaghan, Australia
;
; May 03
;

function diffbydave, array, coords, dir, INTERP=interp

	if keyword_set(interp) then $
	begin
		;array2 = array
		;coords2 = coords
		arraytmp	= dblarr(n_elements(array[*,0])*3,n_elements(array[0,*])*2)
		
		arraytmp[0:n_elements(array[*,0])-1,n_elements(array[0,*]):2*n_elements(array[0,*])-1]	= array
		arraytmp[n_elements(array[*,0]):2*n_elements(array[*,0])-1,n_elements(array[0,*]):2*n_elements(array[0,*])-1]	= array
		arraytmp[2*n_elements(array[*,0]):3*n_elements(array[*,0])-1,n_elements(array[0,*]):2*n_elements(array[0,*])-1]	= array
		
		arraytmp[0:n_elements(array[*,0])-1,0:n_elements(array[0,*])-1]	= $
			reverse([array[n_elements(array[*,0])/2:n_elements(array[*,0])-1,*],array[0:n_elements(array[*,0])/2-1,*]],2)
		arraytmp[n_elements(array[*,0]):2*n_elements(array[*,0])-1,0:n_elements(array[0,*])-1]	= $
			reverse([array[n_elements(array[*,0])/2:n_elements(array[*,0])-1,*],array[0:n_elements(array[*,0])/2-1,*]],2)
		arraytmp[2*n_elements(array[*,0]):3*n_elements(array[*,0])-1,0:n_elements(array[0,*])-1]	= $
			reverse([array[n_elements(array[*,0])/2:n_elements(array[*,0])-1,*],array[0:n_elements(array[*,0])/2-1,*]],2)
		
		coords2	= dblarr(n_elements(array[*,0])*3,n_elements(array[0,*])*2)
		array2	= arraytmp
		
		if dir eq 0 then $
		begin
			delph		= coords[1,0]-coords[0,0] 
			for i=0,n_elements(coords2[0,*])-1 do $
			begin
				coords2[*,i]	= indgen(n_elements(coords2[*,0]))*delph
			endfor
		endif else $
		begin
			delth		= coords[0,1]-coords[0,0] 
			for i=0,n_elements(coords2[*,0])-1 do $
			begin
				coords2[i,*]	= indgen(n_elements(coords2[0,*]))*delth
			endfor
		endelse
		;stop
	endif else $
	begin
		array2 = array
		coords2 = coords
	endelse

	m = n_elements(array2[*,0])
	n = n_elements(array2[0,*])


	darray = dblarr(m,n)

;	For partial differentiation in the x direction
	if dir eq 0 then $
	begin
		for i=0,n-1 do $
		begin
			darray(*,i) = deriv(coords2(*,i),array2(*,i))
		endfor
	endif else $
;	For partial differentiation in the y direction
	begin
		for i=0,m-1 do $
		begin
			darray(i,*) = deriv(coords2(i,*),array2(i,*))
		endfor
	endelse
	darray = darray;congrid(darray, n_elements(array[*,0]), n_elements(array[0,*]), /interp, /minus_one)

	if keyword_set(interp) then $
	begin
		darray	= darray[n_elements(array[*,0]):2*n_elements(array[*,0])-1,n_elements(array[0,*]):2*n_elements(array[0,*])-1]
	endif
;stop	
  return, darray
end
