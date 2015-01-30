; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Given vectors, process through the rotation matrix
;
;  C.L. Waters and D.L. Green
;  Centre for Space Physics Research
;  University of Newcastle
;  New South Wales, Australia
;
;  Comments:
;  Checked at JHU/APL Sept 2010
;
pro amp_rotate_vecs, pi_x, pi_y, pi_z, $              ; coords in rotated coords (use amp_rotate_coords)
                     dbx_in, dby_in, dbz_in, $        ; Input db (x,y,z) components in old coords
                     dbx_out, dby_out, dbz_out, $     ; db(x,yz) rotated
                     dbr_out, dbth_out, dbph_out, $   ; db(r,t,p) rotated
                      r_mat          ; rotation matrix

; Initialise output arrays
  dbx_out=dbx_in & dby_out=dby_in & dbz_out=dbz_in
  dbr_out=dbx_in & dbth_out=dby_in & dbph_out=dbz_in

; Rotate input vecs
  For ii=Long(0),n_elements(dbx_in)-1 do begin
    vec_dat = [dbx_in[ii], dby_in[ii], dbz_in[ii]]
    vec_shift = matrix_multiply ( r_mat, vec_dat, /Atrans )

    dbx_out[ii] = vec_shift[0]
    dby_out[ii] = vec_shift[1]
    dbz_out[ii] = vec_shift[2]

; Get spherical dB components, calc at shifted (x,y,z) coords
    geopack_bcarsp, pi_x[ii],  pi_y[ii],  pi_z[ii], $
                  vec_shift[0], vec_shift[1], vec_shift[2], $
                  b_r, b_th, b_ph
    dbr_out[ii] = b_r
    dbth_out[ii] = b_th
    dbph_out[ii] = b_ph

  end

  b_r=!null
  b_th=!null
  b_ph=!null

end
