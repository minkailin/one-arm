function get_vorticity, vrad, vtheta, vphi, rad, theta, phi
nrad = n_elements(rad)
ntheta = n_elements(theta)
nphi = n_elements(phi)

dvphi_sintheta = vphi
dvtheta = vtheta

;get omega_r*cos(theta)
for i=0, nrad-1 do begin
    for k=0, nphi-1 do begin
        dvphi_sintheta(i,*,k) = deriv(theta(*), vphi(i,*,k)*sin(theta(*)))
        dvphi_sintheta(i,*,k) /= rad(i)*tan(theta(*))
    endfor

    for j=0, ntheta-1 do begin
        dvtheta(i,j,*) = deriv(phi(*), vtheta(i,j,*))
        dvtheta(i,j,*) /= rad(i)*tan(theta(j))
    endfor
endfor
omegar = dvphi_sintheta - dvtheta

;get omega_theta*sin(theta) !factors involving sin(theta) in the expression for the curl is absorbed with pre-factor sin(theta) here
dvrad = vrad
dr_vphi = vphi

for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
          dvrad(i,j,*) = deriv(phi(*), vrad(i,j,*))/rad(i)
        endfor

    for k=0, nphi-1 do begin
        dr_vphi(*,j,k) = deriv(rad(*), rad(*)*vphi(*,j,k))
        dr_vphi(*,j,k) *= sin(theta(j))/rad(*)
    endfor
endfor
omegatheta = dvrad - dr_vphi

omegaz = omegar - omegatheta ; geometric factors already incoroporated. 

return, omegaz
end


pro pdisk_vort, loc=loc, start=start, finish=finish, xrange=xrange, yrange=yrange, $
        ct=ct, zslice=zslice, plotrange=plotrange, r0=r0, range=range, aziline=aziline,phi0=phi0, $
        ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, vsize=vsize, $
        nobox=nobox, iline=iline, boxbound=boxbound

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
if not keyword_set(range) then range=[0.8,1.2]
if not keyword_set(vsize) then vsize=0.5; vsize*Rossby number at core is the boundary value of Ro, taken for measure aspect ratios

pload, 0 

data     = dblarr(nx1, nx3)
rho_2d   = dblarr(nx1, nx2)
vortz_2d = dblarr(nx1, nx2)
dataplot = dblarr(nx1, nx3)
array    = dblarr(nx1, nx3)

haf_nx3  = nx3/2

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

minrossby = dblarr(2, finish-start+1)

for i=start, finish do begin
   pload, i,/silent
   
   vortz = get_vorticity(vx1, vx2, vx3, x1, x2, x3)

   for jj=0, nx2-1 do begin
    for ii=0, nx1-1 do begin
        avg = mean(vx3(ii,jj,*))
        vx3(ii,jj,*) = avg/(x1(ii)*sin(x2(jj)))

        vortz_2d(ii,jj) = mean(vortz(ii,jj,*)) 
        rho_2d(ii,jj)   = mean(rho(ii,jj,*))
    endfor
   endfor   

   if keyword_set(nonaxi) then begin
   for k=0, nx3-1 do vortz(*,*,k) -= vortz_2d(*,*)
   endif

   vortz /= 2d0*vx3 

   data(*,*) = vortz(*,nzslice,*)

   if keyword_set(phi0) then begin

      if(phi0 gt 0.0) then begin
        phiplt = phi0*2d0*!dpi
      endif else begin
        ;; array(*,*) = data(*,*) ;vortz(*,nx2-1,*)
        ;; temp=min(array, minloc)
for ii=0,nx1-1 do array(ii,*) = rho(ii,nx2-1,*)/rho_2d(ii,nx2-1) - 1.0 
temp = max(array, minloc)
        ind  = array_indices(array, minloc)
        rvort  = x1(ind(0))
        phiplt = azi(ind(1))
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

      temp = min(abs(yrange(0) - azi1),y1)
      temp = min(abs(yrange(1) - azi1),y2)

      plotrange0 = [min(dataplot(r1:r2,y1:y2)),max(dataplot(r1:r2,y1:y2))]
   endif else begin
      plotrange0 = plotrange
   endelse
   minrossby(0,i-start) = time(i)
   minrossby(1,i-start) = plotrange0(0)
  
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0

   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1+ textoidl(', tan\psi=')+height_string+'h'
   set_plot,'ps'
   device, filename=strcompress('pdisk_vort_'+name2+'.ps',/remove_all) $
           ,/color, bits_per_pixel=8, xsize=12, ysize=14
   contour,dataplot,rad,azi1,/fill,levels=levels, xstyle=1, $
           xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
           xtickinterval=xtickinterval,ytickinterval=ytickinterval, xrange=xrange, yrange=yrange, title=title
   colorbar, position=[0.852, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'


   if not keyword_set(nobox) then begin
      boundary = vsize*plotrange0(0)

      levels = boundary
      contour,dataplot,rad,azi1,levels=levels, /overplot, color=255, thick=4, path_xy=xy, path_info=info, $
              /path_data_coords
      
      if not keyword_set(iline) then begin
        temp = max(info.N*(1-info.high_low), ii)
      endif else ii = iline

;       for ii=0, n_elements(info)-1 do begin

         if not keyword_set(boxbound) then begin
         S = [INDGEN(info(ii).N), 0]
         nels = n_elements(INFO(ii).OFFSET + S) 
         
         xcoord = xy(0,INFO(ii).OFFSET + S );xy(0, nels-1: n_elements(xy(0,*)) -1 )
         ycoord = xy(1,INFO(ii).OFFSET + S );xy(1, nels-1: n_elements(xy(1,*)) -1 )
         
         PLOTS, xy(*,INFO(ii).OFFSET + S ), /data, color=255, thick=4
       
         rminus = min(xcoord) 
         rplus  = max(xcoord) 
         
         phiminus = min(ycoord) 
         phiplus = max(ycoord)
         endif else begin
         rminus = boxbound(0)        
         rplus  = boxbound(1)
         phiminus = boxbound(2)
         phiplus  = boxbound(3)
         endelse  



         oplot, [rminus,rplus], [phiplus,phiplus], thick=2, linestyle=1, color=255 
         oplot, [rminus,rplus], [phiminus,phiminus], thick=2, linestyle=1, color=255
         
         oplot, [rminus,rminus], [phiminus,phiplus], thick=2, linestyle=1, color=255
         oplot, [rplus,rplus], [phiminus,phiplus], thick=2, linestyle=1, color=255
         
         rvort = mean([rminus,rplus])
         chi = rvort*(phiplus - phiminus)*!dpi/(rplus - rminus)
         print, 'aspect-ratio is', chi
         
;        xyouts, rplus, phiminus, textoidl('\chi='+string(chi, format='(f4.1)')),charsize=1.5,charthick=4
;         xyouts, rminus, phiplus, textoidl('\chi='+string(chi, format='(f4.1)')),charsize=1.5,charthick=4,align=1
          xyouts, rplus, phiplus, textoidl('\chi='+string(chi, format='(f4.1)')),charsize=1.5,charthick=4


;        endfor       
 


      
    endif
   
;     stop 
   
   device,/close
endfor

openw,1,'minrossby.dat'
for i=0, finish-start do begin
printf,1,minrossby(0,i),minrossby(1,i),format='(2(e22.15,x))'
endfor
close,1

end

