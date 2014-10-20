function bigQ, psi, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi

arg1 = (psi - psitrans)/dpsi 
arg2 = (psi + psitrans)/dpsi
res = 1.0 + tanh(arg1) + 1.0 - tanh(arg2)

return, 1.0 + 0.5*(ampnu - 1.0)*res

end

function csq, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

;strictly isothermal
  soundspeed_sq = h*h
  soundspeed_sq*= bigG*mstar/g_unitLength

  return, soundspeed_sq*(g_unitLength/bigR)^q
end

function bump, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  bwidth = 0.05*g_unitLength

  return, 1.0 + (amp - 1.0)*exp(-0.5*((bigR - g_unitLength)/bwidth)^2.0)
end

function surface_density, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  sig0  = 1.0
  sig_s =sig0*(bigR/g_unitLength)^(-p)

  return, sig_s*bump(bigR)
end

function omega_k, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q
  
  omega_sq = bigG*mstar/bigR^3.0
  return, sqrt(omega_sq)
end

function bigH, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  return, sqrt(csq(bigR))/omega_k(bigR)
end

function density3D, bigR, z
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  phi = -bigG*mstar/sqrt(bigR^2 + z^2)
  phi0 = -bigG*mstar/bigR

  vertical = -(phi - phi0)
  vertical/= csq(bigR)

  sigma = surface_density(bigR)
  sigma/= sqrt(2.0*!dpi)*bigH(bigR)

  return, sigma*exp(vertical)
end

function dlogbump, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  bwidth = 0.05*g_unitLength
  
  dbump   = (bump(bigR) - 1.0)
  dbump  *=-(bigR - g_unitLength)/bwidth^2.0
  dbump  /= bump(bigR)
  
  return, dbump
end

function d2logbump, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  bwidth = 0.05*g_unitLength

  DR2 = bwidth^2.0
  B   = bump(bigR)

  d2bump = dlogbump(bigR)/B
  d2bump*= -(bigR - g_unitLength)/DR2
  d2bump += -(1.0 - 1.0/B)/DR2

  return, d2bump
end

function dlog_rho0, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q
  
  bwidth = 0.05*g_unitLength

  dsigma_s = -p/bigR

  dbump   = (bump(bigR) - 1.0)                     
  dbump  *=-(bigR - g_unitLength)/bwidth^2.0 
  dbump  /= bump(bigR)                             
  
  domega_k = -1.5/bigR 
  dbigH   =  - domega_k
  
  return, dsigma_s + dbump - dbigH
end

function azivel, bigR, z
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q
  
  cs2     = csq(bigR)
  dphi0   = bigR*omega_k(bigR)^2.0 
  
  vsq     = cs2*dlog_rho0(bigR) + $
            dphi0
  vsq    *= bigR

  return, sqrt(vsq)
end

function dlogomega_iso, bigR
  common consts, bigG, mstar
  common params, g_unitlength, p , h, amp, q

  cs2   = csq(g_unitLength)
  omega = azivel(bigR, 0.0)/bigR

  d2rho   = (1.5 + p)/bigR^2.0 + d2logbump(bigR)
  d2phi   = -2.0*bigG*mstar/bigR^3 

  result = cs2*d2rho + d2phi
  result/= bigR*omega*omega
  result-= 1.0/bigR
  result/= 2.0

  return, result*bigR
end



pro pdisk_visc2d, nu0=nu0, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi, smallh=smallh, $
                  smallp=smallp, r0=r0, zmax=zmax, bamp=bamp, rtrans=rtrans, znorm=znorm,yrange=yrange, smallq=smallq
;LOCALLY ISOTHERMAL
common consts, bigG, mstar
common params, g_unitlength, p, h, amp, q

bigG         = 1.0
mstar        = 1.0
g_unitlength = r0
p            = smallp 
h            = smallh 
nh           = zmax*smallh
amp          = bamp
q            = smallq


rmin = 0.4
rmax = 2.5
xrange = [rmin,rmax]
nrad = 128
raxis = rmin + (rmax-rmin)*dindgen(nrad)/(nrad-1.0)
temp = min(abs(raxis-r0), rgrid)

zmin = 0.0
zmax = rmax*nh
if not keyword_set(yrange) then yrange = [zmin,zmax]
nz    = 128
zaxis = zmin + (zmax-zmin)*dindgen(nz)/(nz-1.0)

upper_bound = nh*raxis

visc2d = dblarr(nrad, nz)
mdot   = dblarr(nrad, nz)
mask   = intarr(nrad, nz)


H0 = bigH(r0)

for j=0, nz-1 do begin
   z = zaxis(j)
  
   for i=0, nrad-1 do begin
      rad = raxis(i)

      if (z lt upper_bound(i)) then begin

         if not keyword_set(znorm) then begin
            psi = atan2(z, rad)
         endif else psi = z/H0
         
         BQ = bigQ(psi, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi)
         if keyword_set(rtrans) then begin
            fr = 0.5*(1.0 + tanh((rad - r0)/(rtrans*H0)))
            BQ *= fr
            BQ += (1.0-fr)*ampnu
         endif 
         
         if not keyword_set(znorm) then begin
            visc2d(i,j) = nu0*BQ*surface_density(r0) ;*dlogomega_iso(r0)
            visc2d(i,j)/= surface_density(rad)       ;*dlogomega_iso(rad)
         endif else begin 
            visc2d(i,j) = nu0*BQ*density3D(r0, z)/bump(r0) ;*dlogomega_iso(r0)
            visc2d(i,j)/= density3d(rad, z)/bump(rad)      ;*dlogomega_iso(rad)
         endelse
         
         mdot(i,j) = density3d(r0, z)*bigQ(psi, ampnu=ampnu, psitrans=psitrans, dpsi=dpsi)
         mdot(i,j)/= density3d(r0, 0.)*bigQ(0., ampnu=ampnu, psitrans=psitrans, dpsi=dpsi)
         
         mask(i,j) = 1
      endif else begin
         mask(i,j) = 0
         visc2d(i,j) = 1.0
         mdot(i,j)   = 1.0
      endelse
      
   endfor
endfor
visc2d = alog10(visc2d)

disk = where(mask eq 1)
notdisk=where(mask eq 0)

visc2d(notdisk) = max(visc2d(disk))
mdot(notdisk) = max(mdot(disk))

;if not keyword_set(znorm) then begin
plotrange0=[min(visc2d),max(visc2d)]
;endif else  begin
;plotrange0 =[min(visc2d), min(  [ max( visc2d(0,0:nz-1) ) , max(visc2d(nrad-1,0:nz-1)) ]  )  ]
;endelse

levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0

loadct, 5, /silent
title  = textoidl('log[\nu/r_0^2\Omega_K(r_0)]')
set_plot,'ps'
device, filename=strcompress('pdisk_visc2d.ps',/remove_all) $
        , bits_per_pixel=8, xsize=18, ysize=9,/color

contour,visc2d, raxis/r0, zaxis/(r0*smallh), levels=levels, xstyle=1, /fill $
        ,xmargin=[7.0,10.0], ymargin=[4,2], $
        xtickinterval=xtickinterval,xrange=xrange/r0, yrange=yrange/(r0*smallh), title=title, xtitle=textoidl('R/r_0'), $ 
        ytitle = textoidl('z/H(r_0)'), ytickinterval=1, ystyle=1;,charsize=1.2
colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange0,format='(f5.2)';,charsize=1.2

oplot, raxis, upper_bound/(r0*smallh), thick=4, linestyle=0
xyouts, 2, 5, textoidl('\alpha\sim10^{-2}'), charthick=4,charsize=1.5
xyouts, 2, 1, textoidl('\alpha\sim10^{-4}'), charthick=4,charsize=1.5, color=255
;oplot, raxis(0:rgrid+5), upper_bound(0:rgrid+5)/(r0*smallh), thick=4, linestyle=0, color=0
;oplot, raxis(rgrid+5:nrad-1), upper_bound(rgrid+5:nrad-1)/(r0*smallh), thick=4, linestyle=0, color=255
device,/close

;mdot = alog10(mdot)

plotrange0 =[min(mdot), max(mdot)]
;; plotrange0=[0,8.23]
levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0

title  = textoidl('Mdot(R,z)/Mdot(r_0,0)')
set_plot,'ps'
device, filename=strcompress('pdisk_visc2d_mdot.ps',/remove_all) $
        , bits_per_pixel=8, xsize=18, ysize=9,/color

contour,mdot, raxis/r0, zaxis/(r0*smallh), levels=levels, xstyle=1, /fill $
        ,xmargin=[7.0,10.0], ymargin=[4,2], $
        xtickinterval=xtickinterval,xrange=xrange/r0, yrange=yrange/(r0*smallh), title=title, xtitle=textoidl('R/r_0'), $ 
        ytitle = textoidl('z/H(r_0)'), ytickinterval=1;,charsize=1.2
colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange0,format='(f5.2)' ;,charsize=1.2

;oplot, raxis(0:rgrid), upper_bound(0:rgrid)/(r0*smallh), thick=4, linestyle=0, color=31
;oplot, raxis(rgrid+1:nrad-1), upper_bound(rgrid+1:nrad-1)/(r0*smallh), thick=4, linestyle=0, color=15

oplot, raxis, upper_bound/(r0*smallh), thick=4, linestyle=0


device,/close

end
