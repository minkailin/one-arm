function pltphi, x, y
common consts,pi
x=double(x)
y=double(y)
if x gt 0. and y ge 0. then begin
result= atan(y/x)
endif else if x lt 0. and y ge 0. then begin
result=atan(y/x)+pi
endif else if x lt 0. and y lt 0. then begin
result=atan(y/x)+pi
endif else if x gt 0. and y lt 0. then begin
result=atan(y/x)+2.*pi
endif else if x eq 0 and y gt 0. then begin
result=pi
endif else if x eq 0 and y lt 0. then begin
result=3.*pi/2.
endif
return, double(result)
end
