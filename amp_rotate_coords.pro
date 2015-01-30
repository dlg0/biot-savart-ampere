; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Given coords, process through the rotation matrix
;
;  C.L. Waters and D.L. Green
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
;  Comments:
;  Checked at JHU/APL Sept 2010
;
pro amp_rotate_coords, pi_x, pi_y, pi_z, $   ; input coords
                       po_x, po_y, po_z, $   ; rotated X,Y,Z coords
                       r_out, th_out_deg, ph_out_deg, $  ; spherical components of shifted coords
                       r_mat          ; rotation matrix
;
  po_x = pi_x
  po_y = pi_y
  po_z = pi_z

; Rotate input (x,y,z) location data
  For ii=Long(0),n_elements(pi_x)-1 do begin

; rotate coord locations
    xyz_dat = [pi_x[ii], pi_y[ii], pi_z[ii]]
    xyz_shift = matrix_multiply ( r_mat, xyz_dat, /Atrans )

    po_x[ii] = xyz_shift[0]
    po_y[ii] = xyz_shift[1]
    po_z[ii] = xyz_shift[2]

  end
; Get rotated spherical coord locations
  geopack_sphcar, po_x, po_y, po_z, $
    r_out, th_out_deg, ph_out_deg, /to_sphere, /degree   ; find shifted (r,thet,phi)

end
