pro pdisk_pot, loc=loc, start=start, finish=finish, log=log, nopert=nopert, xrange=xrange, yrange=yrange, $
        hold=hold, label=label,float=float, ct=ct, smallq=smallq, zslice=zslice, plotrange=plotrange, r0=r0, aziline=aziline,phi0=phi0, $
        ytickinterval=ytickinterval, colxrange=colxrange, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, basic=basic, cart=cart 

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

if not keyword_set(smallq) then smallq = 1.0 
if not keyword_set(ct) then ct =5   
if not keyword_set(finish) then finish=start 
if not keyword_set(zslice) then zslice= 0.99
if not keyword_set(r0) then r0=1.0
if not keyword_set(basic) then basic = 0

if not keyword_set(float) then begin
   pload, basic
endif else pload, basic, /float
d0 = pot

if not keyword_set(xrange) then xrange=[min(x1),max(x1)]

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

azi   = x3
time  = t/(2d0*!dpi*r0^1.5)

if keyword_set(mp) then begin
   rh    = (mp/3.0)^(1.0/3.0)*r0
   rad  = (x1-r0)/rh
   xtitle=textoidl('(r-r_p)/r_h')
endif else begin
   rad   = x1/r0
   xtitle=textoidl('r/r_0')


   if keyword_set(cart) then begin
;assume uniform spacing
      dr = rad(1) - rad(0)
      dazi = x3(1) - x1(0)

      dnx1 = 512
      xaxis = -rad(nx1-1) + 2*rad(nx1-1)*dindgen(dnx1)/(dnx1-1.0)
      yaxis = xaxis
      dataxy = dblarr(dnx1, dnx1)
   endif
   
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

nzslice = nx2*zslice
height_string = strcompress(string(1.0/(smallh*tan(x2(nzslice))), format='(f3.1)'),/remove_all)

loadct,ct,/silent
for i=start, finish do begin
   if not keyword_set(float) then begin   
      pload, i,/silent
   endif else pload, i, /silent,/float  
   
   d = pot
   
   if keyword_set(nonaxi) then begin
      for jj=0, nx2-1 do begin
         for ii=0, nx1-1 do begin
            avg = mean(d(ii,jj,*))
            d0(ii,jj,*) = avg
         endfor
      endfor
   endif
   
   if not keyword_set(nopert) then begin
      d /= abs(d0) ;potential is negative, so normalize by magntude, dpot/|pot|
      if not keyword_set(log) then   d += 1.0
   endif
   
   if keyword_set(log) then d=alog10(abs(d))
   data(*,*) = d(*,nzslice,*)
   
   if keyword_set(phi0) then begin
      
      if(phi0 gt 0.0) then begin
         phiplt = phi0*2d0*!dpi
      endif else begin
         array(*,*) = d(*,nx2-1,*)
         temp = max(array, maxloc)
         ind  = array_indices(array, maxloc)
         phiplt = azi(ind(1)) + !dpi
         print, 'density center is', phiplt/(2.0*!dpi)
      endelse
      
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
      
      if keyword_set(colxrange) then begin
         temp = min(abs(colxrange(0) - rad),cr1)
         temp = min(abs(colxrange(1) - rad),cr2)
      endif else begin
         cr1 = r1
         cr2 = r2
      endelse
      
      temp = min(abs(yrange(0) - azi1),y1)
      temp = min(abs(yrange(1) - azi1),y2)
      
      plotrange0 = [min(dataplot(r1:r2,y1:y2)),max(dataplot(cr1:cr2,y1:y2))]
   dataplot(0,*) = plotrange0(1)
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0
   
   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1+ textoidl(', tan\psi=')+height_string+'h'

   if not keyword_set(cart) then begin   
   set_plot,'ps'
   device, filename=strcompress('pdisk_pot_'+name2+'.ps',/remove_all) $
           ,/color, bits_per_pixel=8, xsize=12, ysize=14
   contour,dataplot,rad,azi1,/fill,levels=levels, xstyle=1, $
           xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
           xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title
   if keyword_set(label) then begin
      xyouts, xrange(0)+(xrange(1)-xrange(0))/20., yrange(0)+(yrange(1)-yrange(0))/10., textoidl(label), charsize=2,charthick=8, color=255 
   endif
   colorbar, position=[0.865, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'
;   oplot,[0,0],[1,1]/2.,psym=7,symsize=1.5,color=!D.Table_size*1.1
   if keyword_set(aziline) then begin
      oplot, [min(rad), max(rad)], [aziline,aziline], thick=4, linestyle=1
   endif   
   device,/close
   endif else begin
      for jj=0, dnx1-1 do begin
         y = yaxis(jj)
         for ii=0, dnx1-1 do begin
            x = xaxis(ii)
            
            r_t = sqrt(x^2 + y^2)
            azi_t = pltphi(x, y)
            
            if( (r_t ge xrange(0)) and (r_t le rad(nx1-1)) )then begin
               
               temp = min(abs(rad - r_t),   x0)
               temp = min(abs(azi - azi_t),  y0)
               
               ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
               jp = y0 + (azi_t - azi(y0))/(dazi/2.0) + 0.0
               
               dataxy(ii,jj) = bilinear(dataplot, ip, jp)
            endif else begin
               dataxy(ii,jj) = 1d10
            endelse
         endfor
      endfor
      
      set_plot, 'ps'
      device, filename=strcompress('pdiskxy_pot_'+name2+'.ps',/remove_all) $
              ,/color, bits_per_pixel=8,xsize=14, ysize=12
      contour, dataxy, xaxis, yaxis, /isotropic,/fill,levels=levels,title=title $
               ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval
      colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange0,format='(f5.2)'
;      xyouts, -2.4, -2.4, 'PLUTO 3D', charsize=1.5,charthick=8, color=0
      ;; if keyword_set(aziline) then begin
      ;;     num = n_elements(aziline)
      ;;     for i = 0, num-1 do begin
      ;;         angle(*) = 2d0*!dpi*aziline(i)
      ;;         oplot, rad, angle, thick=2, color=255,/polar
      ;;     endfor
      ;; endif
      
      ;;  if keyword_set(name) then begin
      ;;     xyouts, 0, -20./r0, textoidl(name),charsize=1.5, charthick=6, color=255, alignment=0.5
      ;; endif
      device,/close
      ;stop
   endelse

endfor
end
