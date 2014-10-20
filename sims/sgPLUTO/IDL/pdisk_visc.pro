function bigQ, psi, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi

arg1 = (psi - psitrans)/dpsi 
arg2 = (psi + psitrans)/dpsi
res = 1.0 + tanh(arg1) + 1.0 - tanh(arg2)

return, 1.0 + 0.5*(ampnu - 1.0)*res

end

pro pdisk_visc, nu0=nu0, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi, smallh=smallh, zmax=zmax, label=label, $
    legend=legend, xrange=xrange, yrange=yrange

nx=128

xaxis = zmax*dindgen(nx)/(nx-1.0)
xx    = tan(xaxis)/smallh 
visc  = dblarr(nx)

for i=0, nx-1 do begin
visc(i) = nu0*bigQ(xaxis(i), ampnu=ampnu(0), psitrans=psitrans(0), dpsi=dpsi(0))
endfor

fname = 'pdisk_visc.ps'
  set_plot, 'ps'
  device, filename=fname $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, xx, visc,xmargin=[7.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange,xtitle=textoidl('tan\psi/h'),$
        linestyle = 0, ytitle = textoidl('Q(\psi)'), xtickinterval=0.5, xminor=5, ytickinterval=ytickinterval,/ylog
for k=1, n_elements(psitrans)-1 do begin
for i=0, nx-1 do begin
visc(i) = nu0*bigQ(xaxis(i), ampnu=ampnu(k), psitrans=psitrans(k), dpsi=dpsi(k))
endfor
oplot, xx, visc, thick=4, linestyle=k
endfor
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for k=0, n_elements(label)-1 do begin
      oplot, [x0,x1], [y0,y0]*10.^(-dy*k), thick=4, linestyle=k
      xyouts, x1, y0*10.^(-dy*k),textoidl(label(k)),charsize=1.5
   endfor
endif
device,/close

end
