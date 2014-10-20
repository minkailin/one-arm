pro pepeos, rp=rp, hp=hp, smallh=smallh, q=q, xrange=xrange, soft=soft

eosn=3.5
hp_over_smallh = hp/smallh
epsilon = rp*(q/3d0)^(1d0/3d0)
nr=256
nz=64

array_cylind = dblarr(nr,nz)
array_sph    = dblarr(nr,nz)

xaxis = xrange(0)  +(xrange(1)-xrange(0))*dindgen(nr)/(nr-1d0)
yaxis=  2d0*smallh*(xrange(1))*dindgen(nz)/(nz-1d0) 

for i= 0, nr-1 do begin
    for j=0, nz-1 do begin
        bigR = xaxis(i)
        z    = yaxis(j)
        sphr = sqrt(bigR^2 + z^2)
        iso_eos = smallh/sqrt(bigR)

        distp = sqrt((bigR - rp)^2 + epsilon^2 + z^2)
        omegas_sq = 1d0/bigR^3
        omegap_sq = q/distp^3

        array_cylind(i,j) = smallh*bigR*hp*distp
        array_cylind(i,j)*= sqrt(omegas_sq + omegap_sq) 
        array_cylind(i,j)/= ((smallh*bigR)^eosn + (hp*distp)^eosn)^(1d0/eosn)
        array_cylind(i,j)/= iso_eos

        distp = sqrt((bigR - rp)^2 + z^2 + epsilon^2)
        omegas_sq = 1d0/sphr^3
        omegap_sq = q/distp^3

        array_sph(i,j) = smallh*sphr*hp*distp
        array_sph(i,j)*= sqrt(omegas_sq + omegap_sq) 
        array_sph(i,j)/= ((smallh*sphr)^eosn + (hp*distp)^eosn)^(1d0/eosn)
        array_sph(i,j)/= iso_eos
    endfor
endfor


loadct, 5, bottom=0
set_plot, 'ps'
device, filename='pepeos_cylind.ps' $
  ,/color, bits_per_pixel=8,xsize=18, ysize=6

plotrange=[min(array_cylind),max(array_cylind)]
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)

contour, array_cylind,xaxis,yaxis,/fill $
  ,ymargin=[3,3],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval,xtitle='R', ytitle=textoidl('z'), levels=levels
colorbar, position=[0.895, 0.18, 0.93, 0.82],/vertical,/right,format='(f5.2)', range=plotrange
device,/close


set_plot, 'ps'
device, filename='pepeos_sph.ps' $
  ,/color, bits_per_pixel=8,xsize=18, ysize=6

plotrange=[min(array_sph),max(array_sph)]
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)

contour, array_sph,xaxis,yaxis,/fill $
  ,ymargin=[3,3],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval,xtitle='R', ytitle=textoidl('z'), levels=levels
colorbar, position=[0.895, 0.18, 0.93, 0.82],/vertical,/right,format='(f5.2)', range=plotrange
device,/close
end
