pro pdisk_vr, loc=loc, start=start, finish=finish, xrange=xrange, yrange=yrange, nopert=nopert $
        , ct=ct, zslice=zslice, plotrange=plotrange, r0=r0, aziline=aziline,phi0=phi0, $
        ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, basic=basic

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
d0 = vx1

if keyword_set(basic) then begin
for j=0, nx2-1 do begin
for i=0, nx1-1 do begin
d0(i,j,*) = mean(d0(i,j,*))
endfor
endfor
endif


data     = dblarr(nx1, nx3)
array    = dblarr(nx1, nx3)
dataplot = dblarr(nx1, nx3)

;normalize axis, time
if keyword_set(mp) then begin
rh    = (mp/3.0)^(1.0/3.0)*r0
rad  = (x1-r0)/rh
xtitle=textoidl('(r-r_0)/r_h')
endif else begin
rad   = x1/r0
xtitle=textoidl('r/r_0')
endelse
if keyword_set(phi0) then begin
ytitle = textoidl('(\phi-\phi_0)/\pi')
azi1 = x3/!dpi - 1.0
if not keyword_set(yrange) then yrange=[-1,1]
endif else begin
ytitle = textoidl('\phi/2\pi')
azi1 = x3/(2d0*!dpi)
if not keyword_set(yrange) then yrange=[0,1]
endelse

azi   = x3
time  = t/(2d0*!dpi*r0^1.5)

nzslice = nx2*zslice
height_string = strcompress(string(1.0/(smallh*tan(x2(nzslice))), format='(f3.1)'),/remove_all)


loadct,ct,/silent
for i=start, finish do begin
   pload, i,/silent
   
   d = vx1

   if keyword_set(nonaxi) then begin
   for jj=0, nx2-1 do begin
    for ii=0, nx1-1 do begin
        avg = mean(d(ii,jj,*))
        d0(ii,jj,*) = avg
    endfor
   endfor
   endif

   if not keyword_set(nopert) then begin
      d -= d0
   endif
   data(*,*) = d(*,nzslice,*)

   if keyword_set(phi0) then begin

        phiplt = phi0*2d0*!dpi

     if(phiplt gt !dpi) then begin
       temp2=min(abs(azi-(phiplt-!dpi)),grid2)
       dataplot(0:nx1-1,0:nx3-1-grid2)   = data(0:nx1-1,grid2:nx3-1)
       dataplot(0:nx1-1,nx3-grid2:nx3-1) = data(0:nx1-1, 0:grid2-1)
     endif
    if(phiplt lt !dpi) then begin
       temp2=min(abs(azi-(phiplt+!dpi)),grid2)
       dataplot(0:nx1-1,nx3-grid2:nx3-1) = data(0:nx1-1,0:grid2-1)
       dataplot(0:nx1-1,0:nx3-1-grid2) = data(0:nx1-1, grid2:nx3-1)
    endif
    if(phiplt eq !dpi) then dataplot = data
   endif else dataplot = data

   if not keyword_set(plotrange) then begin
      temp = min(abs(xrange(0) - rad),r1)
      temp = min(abs(xrange(1) - rad),r2)

      temp = min(abs(yrange(0) - azi1),y1)
      temp = min(abs(yrange(1) - azi1),y2)

      plotrange0 = [min(data(r1:r2,y1:y2)),max(data(r1:r2,y1:y2))]
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0


   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1+ textoidl(', (h.tan\theta)^{-1}=')+height_string
   set_plot,'ps'
   device, filename=strcompress('pdisk_vr_'+name2+'.ps',/remove_all) $
           ,/color, bits_per_pixel=8, xsize=12, ysize=14
   contour,dataplot,rad,azi1,/fill,levels=levels, xstyle=1, $
           xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
           xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title, xminor=5
   colorbar, position=[0.851, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'

   if keyword_set(aziline) then begin
   oplot, [min(rad), max(rad)], [aziline,aziline], thick=4, linestyle=1
   endif   

   device,/close
   endfor
end
