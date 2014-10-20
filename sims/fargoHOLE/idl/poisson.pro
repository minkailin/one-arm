pro poisson

r0 = 0.4d0
rmax = 2.5d0
B = 0.3*0.05

nr   = 256
nphi = 768

sigma = dblarr(nphi,3*nr)
sr = dblarr(nphi, 3*nr)
kr = dblarr(nphi, 3*nr)

dtheta = 2d0*!dpi/nphi
du = alog(rmax/r0)/nr

uaxis = du/2 -alog(rmax/r0) + dindgen(3*nr)*du
azi = dtheta*dindgen(nphi)

sigma(*,nr:2*nr - 1) = 2d-3

for i = 0, nphi -1 do begin
    for j = 0, 3*nr - 1 do begin
        theta =  azi(i)
        u     =  uaxis(j)
        sr(i,j) = sigma(i,j)*exp(u/2d0)
        kr(i,j) = 1d0 + B^2 - exp(-u)*cos(theta)
        kr(i,j)/= (2d0*(cosh(u) - cos(theta)) + B*B*exp(u))^(3d0/2d0)
    endfor
endfor

srfft = fft(sr, -1, /double)
krfft = fft(kr, -1, /double)

grfft = srfft*krfft*dtheta*du*3*nr*nphi

gr = fft(grfft, 1, /double)

for i=0, nphi -1 do begin
    gr(i,*) *= -exp(-uaxis/2d0)
endfor

gr += sigma*du*dtheta/B

;gr1d = dblarr(nr)

;for j=0, nr-1 do begin
;    gr1d(j) = mean(gr(*,j))
;endfor

;gr1d = alog10(abs(gr(0,0:nr-1)))

set_plot, 'ps'
device, filename='poissontest.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot, r0*exp(uaxis(nr:2*nr-1)), alog10(abs(gr(0,nr:2*nr-1))), xmargin=[6,2],ymargin=[3,2] $
  ,charsize=1.5, thick=4 $
  ,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
  ,xrange=xrange,yrange=yrange
device,/close

end



