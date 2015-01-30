; Ver : 2012_02
; AMPERE Software - IDL 8.0 or later
;
; Procedure to calc rotation matrix given new Iridium satellite track centre (amp_get_new_centre.pro)
; Uses the Rodrigues rotation formula; R=P+(I-P)cos(thet)+Qsin(thet)
;
; C.L. Waters
; Centre for Space Pysics Research
; University of Newcastle, NSW,
; Australia
;
; Ver 1 : Sept, 2011
;
; Comments:
; Fix error if CRLat comes in as 0.0 - CLW, Nov 2014
;
Pro amp_calc_rotmat, CRLat, RLon, south, rot_mat
; If south = 1 then rotate to south hemisphere intersection else use Nth hemisphere
;
  tol = 0.01    ; tol in deg
; Error check
  status=0
  catch,error_status
  if (error_status ne 0) then begin
    status=1
    PRINT, 'Error message: ', !ERROR_STATE.MSG
    catch,/cancel
    return
  endif
;
  if (abs(CRLat) lt tol) then begin
    rot_mat=[[1.0, 0.0, 0.0],$                ; matrices for Rodrigues rotation formula; R=P+(I-P)cos(thet)+Qsin(thet)
       [0.0, 1.0, 0.0],$
       [0.0, 0.0, 1.0]]
  end else begin
    
; Get rotation axis (perp to plane defined by z and shifted point and the origin)
    r0=6371.0+780.0 ; km
    geopack_sphcar,r0,CRLat,RLon,x_sh,y_sh,z_sh,/to_rect,/degree    ; XYZ coord of new pole in old system (in correct Hemis)
    sh_mag=sqrt(x_sh^2+y_sh^2+z_sh^2)                               ; magnitude of shifted pole vector
    If south then r0=-r0                                            ; point r vec along Z axis in correct direction for Hemis

; cross product (x,y,z) so for south->(0,0,-r) is the pole
    cross_p,[x_sh,y_sh,z_sh],[0.0,0.0,r0],v_rot                     ; rotation is about this axis, v_rot should always have a zero z component
    c_arg=(r0*z_sh)/(sh_mag*abs(r0))
    rot_ang=acos(c_arg)                                             ; rotation angle, calc from dot product

    If keyword_set(diag) then Print,'Rotation angle off pole : ',rot_ang*180.0/!dpi,' deg'
    v_rot_mag=sqrt(v_rot(0)^2+v_rot(1)^2+v_rot(2)^2)                   ; magnitude of rotation vector
    urx=v_rot(0)/v_rot_mag                                             ; unit rotation axis vector (perp to plane containing Z and shifted origin)
    ury=v_rot(1)/v_rot_mag
    urz=v_rot(2)/v_rot_mag

; Matrices for rotation
    p_a=[[urx*urx, urx*ury, urx*urz],$                ; matrices for Rodrigues rotation formula; R=P+(I-P)cos(thet)+Qsin(thet)
       [urx*ury, ury*ury, ury*urz],$
       [urx*urz, ury*urz, urz*urz]]
    i_a=[[1.0, 0.0, 0.0],$                            ; identity matrix
       [0.0, 1.0, 0.0],$
       [0.0, 0.0, 1.0]]
    q_a=[[0.0, -urz, ury],$
       [urz, 0.0, -urx],$
       [-ury, urx, 0.0]]
    rot_mat = p_a + (i_a-p_a)*cos(rot_ang) + q_a*sin(rot_ang)  ; IDL orders by [column,row]
  end
end
