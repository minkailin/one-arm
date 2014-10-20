pro pdisk_rz, loc=loc, start=start, finish=finish, log=log, nopert=nopert, xrange=xrange, yrange=yrange, $
        ct=ct, azislice=azislice, plotrange=plotrange, r0=r0, phi0=phi0, $
        ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, $
        red=red, arrct=arrct, arrcol=arrcol, hsize=hsize, length=length, novect=novect, basic=basic

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes;  ** Chombo data structure **
                                      ;  ** loaded when HDF5LOAD is called **

COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
                             bx1, bx2, bx3, $
                             Ax1, Ax2, Ax3, $
                             bx1s, bx2s, bx3s,$
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

if not keyword_set(ct) then ct =5   
if not keyword_set(finish) then finish=start 
if not keyword_set(zslice) then zslice= 0.99
if not keyword_set(r0) then r0=1.0
if not keyword_set(basic) then basic = 0

pload, basic
d0 = rho

if keyword_set(basic) then begin
   for j=0, nx2-1 do begin
      for i=0, nx1-1 do begin
         d0(i,j,*) = mean(d0(i,j,*))
      endfor
   endfor
endif


data     = dblarr(nx1, nx2)
data_xy  = dblarr(nx1, nx3)
vx       = dblarr(nx1, nx2)
vy       = dblarr(nx1, nx2)

;normalize axis, time
if keyword_set(mp) then begin
   rh    = (mp/3.0)^(1.0/3.0)*r0
   rad  = (x1-r0)/rh
   xtitle=textoidl('(r-r_0)/r_h')
endif else begin
   rad   = x1/r0
   xtitle=textoidl('r/r_0')
endelse

z1     = 1.0/(smallh*tan(x2))
ytitle = textoidl('tan\psi/h')
if not keyword_set(yrange) then yrange=[min(z1), max(z1)]

azi   = x3/(2d0*!dpi)
time  = t/(2d0*!dpi*r0^1.5)

if keyword_set(azislice) then begin
   nphi = nx3*azislice
   azi_string = strcompress(string(azi(nphi), format='(f4.2)'),/remove_all)
endif

if not keyword_set(red) then red=[nx1,nx2]

loadct,ct,/silent
for i=start, finish do begin
   pload, i,/silent
   
   d    = rho
   vr   = vx1
   vt   = vx2
   
   if keyword_set(nonaxi) then begin
   for jj=0, nx2-1 do begin
    for ii=0, nx1-1 do begin
        avg = mean(d(ii,jj,*))
        d0(ii,jj,*) = avg
    endfor
   endfor
   endif

   if not keyword_set(nopert) then begin
      d /= d0
      if not keyword_set(log) then   d -= 1.0
   endif
   if keyword_set(log) then d=alog10(d)
 
   if (phi0 lt 0.0) then begin
        data_xy(*,*) = d(*,nx2-1,*)
        temp = max(data_xy, maxloc)
        ind  = array_indices(data_xy, maxloc)
        nphi = ind(1)
        azi_string = strcompress(string(azi(nphi), format='(f4.2)'),/remove_all)
   endif else begin
        nphi = phi0*nx3
        azi_string = strcompress(string(azi(nphi), format='(f4.2)'),/remove_all)
   endelse
   data(*,*) = d(*,*,nphi)

   for jj=0, nx2-1 do begin
   for ii=0, nx1-1 do begin

           if keyword_set(mp) then begin
              vx(ii,jj) = vr(ii,jj,nphi)/rh   ;divide by rh if doing scaled radius
           endif else begin
              vx(ii,jj) = vr(ii,jj,nphi)/r0
           endelse

           vy(ii,jj) =  -vt(ii,jj,nphi)/(smallh*sin(x2(jj))^2)
           vy(ii,jj)/= x1(ii)     ;polar angular speed 
   endfor
   endfor

   if not keyword_set(plotrange) then begin
      temp = min(abs(xrange(0) - rad),r1)
      temp = min(abs(xrange(1) - rad),r2)

      temp = min(abs(yrange(0) - z1),y2)
      temp = min(abs(yrange(1) - z1),y1)

      plotrange0 = [min(data(r1:r2,y1:y2)),max(data(r1:r2,y1:y2))]
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0


   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1+ textoidl(', \phi/2\pi=')+azi_string
   set_plot,'ps'
   device, filename=strcompress('pdisk_rz_'+name2+'.ps',/remove_all) $
           , bits_per_pixel=8, xsize=18, ysize=9,/color

   contour,data,rad,z1,levels=levels, xstyle=1, /fill $
           , xtitle=xtitle, ytitle=ytitle,charsize=1.0,xmargin=[7.0,10.0], ymargin=[4,2], $
           xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title;, xminor=5
   colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange0,format='(f5.2)'
 
   if not keyword_set(novect) then begin
   vx_small = congrid(vx, red(0), red(1))
   vy_small = congrid(vy, red(0), red(1))
   rad_small = congrid(rad, red(0))
   zaxis = congrid(z1, red(1))

   temp = min(abs(rad_small - xrange(0)), r1)
   temp = min(abs(rad_small - xrange(1)), r2)
   temp = min(abs(zaxis - yrange(0)), t2)
   temp = min(abs(zaxis - yrange(1)), t1)

   r1 +=1
   r2 -=1

   if keyword_set(arrct) then loadct, arrct,/silent
      velovect2, vx_small(r1:r2,t1:t2), vy_small(r1:r2,t1:t2), rad_small(r1:r2), zaxis(t1:t2), color=arrcol,/overplot $
                 , length=length, hsize=hsize, thick=2
   if keyword_set(arrct) then loadct, ct,/silent
   endif

   device,/close
   endfor
end
