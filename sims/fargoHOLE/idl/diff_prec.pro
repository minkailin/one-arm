function prec, a, amin=amin, amax=amax
prec_min = 0.1
;rate = prec_min*(1.0 - a/amin);*(1.0-a/amax)
r0 = mean([amin,amax])
dr = (amax - amin)/2.0
rate = -prec_min*exp(-0.5*(a-r0)^2/dr^2)
return,rate
end 

pro diff_prec, ecc=ecc, tend=tend

nt    = 1000
tmax  = 2d0*!dpi*tend
taxis = dindgen(nt)*tmax/(nt-1.0)

amin = 1.0
amax = 2.0 
na   = 20
aaxis = amin + dindgen(na)*(amax - amin)/(na-1.0)

eaxis = dblarr(na)
for i=0, na-1 do begin
eaxis(i) = ecc*randomu(seed)
endfor 

nphi = 128
theta = 2d0*!dpi*dindgen(nphi)/(nphi-1.0)


for k=0, nt-1 do begin
ks=string(k,format='(I03)')

omp = prec(aaxis(0), amin=amin, amax=amax)

ecc = eaxis(0)
rad = aaxis(0)*(1.0 - ecc^2)/(1.0 - ecc*cos(theta - omp*taxis(k)))

set_plot, 'x'
;device, filename='diff_prec'+ks+'.ps' $
;           ,/color, bits_per_pixel=8,xsize=12, ysize=12

    plot, /polar, rad, theta,  $
    ymargin=[2,2],xmargin=[2,2], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, $
     xrange=[-1,1]*3, yrange=[-1,1]*3, xstyle=1, ystyle=1,/isotropic;,thick=4

for i=1, na - 1 do begin   
omp = prec(aaxis(i), amin=amin, amax=amax)
ecc = eaxis(i)
rad = aaxis(i)*(1.0 - ecc^2)/(1.0 - ecc*cos(theta - omp*taxis(k)))
oplot,/polar,rad,theta;,thick=4
endfor
;device,/close
;print, 'k'
endfor

end


