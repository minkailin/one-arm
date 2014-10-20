pro pdisk_nonaxi_evol, start=start, finish=finish, xrange=xrange, yrange=yrange, $
        label=label, legend=legend, r0=r0, max=max, basic=basic, $
        ytickinterval=ytickinterval, rrange=rrange, xtickinterval=xtickinterval, azimodes=azimodes 

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes;  ** Chombo data structure **
                                      ;  ** loaded when HDF5LOAD is called **

COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
                             bx1, bx2, bx3, $
                             Ax1, Ax2, Ax3, $
                             bx1s, bx2s, bx3s, pot, $
                         ; ----------------------------------------
                             v1, v2, v3, $   ; Kept for backward
                             b1, b2, b3, $   ; compatibility with 
                             A1, A2, A3, $   ; PLUTO 3 data name
                             b1s, b2s, b3s, $ ;
                             pr,            $ ;
                         ; -----------------------------------------
                  prs, psi_glm, fneut, fion, fh2, $
                  tmp, tr1, tr2, tr3, tr4,$
                  fhI, fheI, fheII,$
                  fcI, fcII, fcIII, fcIV, fcV,$
                  fnI, fnII, fnIII, fnIV, fnV,$
                  foI, foII, foIII, foIV, foV,$
                  fneI, fneII, fneIII, fneIV, fneV,$
                  fsI, fsII, fsIII, fsIV, fsV, vars, vname
COMMON PLUTO_RUN,  t, dt, nlast, first_call
!p.font = 0

if not keyword_set(finish) then finish=start 
if not keyword_set(r0) then r0 = 1d0 
if not keyword_set(basic) then basic=start

nmodes = n_elements(azimodes)

pload, basic

rad   = x1
theta = x2
iphi  = x3 

nrad   = nx1
ntheta = nx2
nphi   = nx3

p0  = 2d0*!dpi*r0^1.5

bigR_3D     = dblarr(nrad,ntheta,nphi)
bigR_axisym = dblarr(nrad,ntheta)
dVol_3D     = dblarr(nrad,ntheta,nphi)
dVol_axisym = dblarr(nrad,ntheta)


time     = dblarr(finish-start+1)
mode_amp = dblarr(nmodes, finish-start+1)
mode_ang = dblarr(nmodes, finish-start+1)
tot_ang  = dblarr(finish-start+1)

     if keyword_set(rrange) then begin
     temp = min(abs(rrange(0) - rad), r1)
     temp = min(abs(rrange(1) - rad), r2)
     endif else begin
     r1 = 0
     r2 = nx1-1
     endelse

;set up geometric arrays 
   for kk=0, nx3-1 do begin
   for jj=0, nx2-1 do begin
   for ii=0, nx1-1 do begin
   bigR_3D(ii,jj,kk)    = x1(ii)*sin(x2(jj))
   dVol_3D(ii,jj,kk)    = 2d0*x1(ii)^2*sin(x2(jj))*dx1(ii)*dx2(jj)*dx3(kk) ; factor of 2 for lower disk 
   endfor
   endfor
   endfor

   for jj=0, nx2-1 do begin
   for ii=0, nx1-1 do begin
   bigR_axisym(ii,jj)    = x1(ii)*sin(x2(jj))
   dVol_axisym(ii,jj)    = 2d0*x1(ii)^2*sin(x2(jj))*dx1(ii)*dx2(jj)*2d0*!dpi ; 2pi from phi-integration, factor of 2 in front for lower disk 
   endfor
   endfor
   
;normalization for mode amplitude and ang mom amplitude
   vphi_fft = fft(vx3, -1, dimension=3,/double)
   rho_fft = fft(rho, -1, dimension=3,/double)

   czero  =total(real_part(rho_fft(*,*,0))*dVol_axisym) 
;   angmom0=total(real_part(rho_fft(*,*,0)*vphi_fft(*,*,0))*bigR_axisym*dVol_axisym) 


temp = dblarr(nphi/2+1)
for m=0, nphi/2 do begin
data2d = real_part( rho_fft(r1:r2,*,m)*conj(vphi_fft(r1:r2,*,m)) )*bigR_axisym(r1:r2,*)
temp(m) = total(data2d*dVol_axisym(r1:r2,*))
if (m ne 0) then temp(m) *= 2d0
endfor
angmom0 = total(temp)





for n=start, finish do begin
   ks   = string(n,format='(I03)')

   pload, n, /silent
   vphi = vx3   

   time(n-start) = t(n)/p0 
    
   vphi_fft = fft(vphi, -1, dimension=3,/double)
   rho_fft = fft(rho, -1, dimension=3,/double)

   for mm=0, nmodes-1 do begin
      m = azimodes(mm)
;mode amplitudes
      if not keyword_set(max) then begin 
      re  = total(real_part(rho_fft(r1:r2,*,m))*dVol_axisym(r1:r2,*))
      im  = total(imaginary(rho_fft(r1:r2,*,m))*dVol_axisym(r1:r2,*)) 
      amp = dcomplex(re,im)

      mode_amp(mm, n-start) = alog10(abs(amp)/czero + 1d-16)
      endif else begin
      mode_amp(mm, n-start) = max(alog10( abs(rho_fft(r1:r2,*,m)/rho0(r1:r2,*,0)) ) )       
      endelse 

;angular momentum components  
      jtot = total(real_part(rho_fft(r1:r2,*,m)*conj(vphi_fft(r1:r2,*,m)))*bigR_axisym(r1:r2,*)*dVol_axisym(r1:r2,*))
      if (m ne 0) then jtot *= 2d0 ;extra factor of 2 for non-axisymmetric modes 
      mode_ang(mm, n-start) = jtot

    endfor

;total ang mom and total mass 
   mass  = total(real_part(rho_fft(r1:r2,*,0))*dVol_axisym(r1:r2,*))
   angmom= total(rho(r1:r2,*,*)*vx3(r1:r2,*,*)*bigR_3D(r1:r2,*,*)*dVol_3D(r1:r2,*,*))

   tot_ang(n-start) = angmom

print, n, tot_ang(n-start), total(mode_ang(*,n-start)) ,mass, format='(I03,x,3(e22.15,x))'
endfor
  if not keyword_set(max) then begin
  ytitle=textoidl('log_{10}(C_m/C_{0,t=0})')
  endif else begin
  ytitle=textoidl('max(log_{10}|\rho_m/\rho_{0,t=0}|)')
  endelse

  set_plot, 'ps'
  device, filename='nonaxi_evol.ps' $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time, mode_amp(0,*), xmargin=[6.5,1.5], ymargin=[3.5,0.5] , ytitle=ytitle, ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(mode_amp),max(mode_amp)], xstyle=1
  for j=1, nmodes-1 do begin
  oplot, time, mode_amp(j,*), thick=4, linestyle=j
  endfor
   if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
        for j=0, n_elements(label)-1 do begin
    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
    xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
    endfor
    endif
    device,/close


  tot = mode_ang(0,*)
  for i=0, finish-start do tot(i) = total(mode_ang(*,i))

  for j=0, nmodes-1 do begin
  m = azimodes(j)
  if(m lt 1) then mode_ang(j,*) = mode_ang(j,*) - angmom0  
;   mode_ang(j,*) = mode_ang(j,*) -  mode_ang(j,0)
  endfor
  mode_ang /=angmom0

  if not keyword_set(yrange) then begin
  yrange1=[min(mode_ang),max(mode_ang)]
  endif else yrange1=yrange

  set_plot, 'ps'
  device, filename='nonaxi_evol_ang.ps' $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time, mode_ang(0,*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\DeltaJ_m/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=yrange1, xstyle=1
  for j=1, nmodes-1 do begin
  oplot, time, mode_ang(j,*), thick=4, linestyle=j
  endfor
;  oplot, time, tot_ang/angmom0-1d0, thick=1, linestyle=0
  oplot, time, tot/angmom0-1d0, thick=1, linestyle=0
    if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
        for j=0, n_elements(label)-1 do begin
    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
    xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
    endfor
    endif
   device,/close

  tot_ang /= angmom0
  tot_ang -= 1d0
  set_plot, 'ps'
  device, filename='nonaxi_evol_totj.ps' $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time(*), tot_ang(*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\Delta J_{tot}/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(tot_ang),max(tot_ang)], xstyle=1
  device,/close
end
