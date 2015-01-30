function cross,v1,v2

 	vout=[[v1(1)*v2(2)-v1(2)*v2(1)],$
 	      [v1(2)*v2(0)-v1(0)*v2(2)],$
 	      [v1(0)*v2(1)-v1(1)*v2(0)]]
 	return, transpose(vout)

end


