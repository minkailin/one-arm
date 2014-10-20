function pltphi, x, y

x=double(x)
y=double(y)
if x gt 0. and y ge 0. then begin
result= atan(y/x)
endif else if x lt 0. and y ge 0. then begin
result=atan(y/x)+!dpi
endif else if x lt 0. and y lt 0. then begin
result=atan(y/x)+!dpi
endif else if x gt 0. and y lt 0. then begin
result=atan(y/x)+2.*!dpi
endif else if x eq 0 and y gt 0. then begin
result=!dpi
endif else if x eq 0 and y lt 0. then begin
result=3.*!dpi/2.
endif
return, double(result)
end
