pro pdisk_vmach, start=start, finish=finish, log=log, xrange=xrange, yrange=yrange, $
                ct=ct, plotrange=plotrange, r0=r0, ytickinterval=ytickinterval, $
                smallh=smallh, smallq=smallq, mp=mp, xtickinterval=xtickinterval, mode=mode, zslice=zslice 
  
  COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
     dx1,dx2,dx3, geometry, xpos, ypos,$
     AMRLevel, AMRBoxes         ;  ** Chombo data structure **
                                ;  ** loaded when HDF5LOAD is called **
  
  COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
     bx1, bx2, bx3, $
     Ax1, Ax2, Ax3, $
     bx1s, bx2s, bx3s,$
                                ; ----------------------------------------
     v1, v2, v3, $              ; Kept for backward
     b1, b2, b3, $              ; compatibility with 
     A1, A2, A3, $              ; PLUTO 3 data name
     b1s, b2s, b3s, $           ;
     pr,            $           ;
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
  if not keyword_set(r0) then r0=1.0
  
  pload, 0
  
  data = dblarr(nx1, nx2)
  vz = dblarr(nx1,nx2,nx3)
  vz_over_vtot = dblarr(nx1, nx2, nx3)
  energy=dblarr(3,finish-start+1)
  
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
  if keyword_set(zslice) then begin
     temp = min(abs(zslice - z1),zgrid)
  endif


  azi   = x3/(2d0*!dpi)
  time  = t/(2d0*!dpi*r0^1.5)
  
  temp = min(abs(xrange(0) - rad),r1)
  temp = min(abs(xrange(1) - rad),r2)
  
  temp = min(abs(yrange(0) - z1),y2)
  temp = min(abs(yrange(1) - z1),y1)


  loadct,ct,/silent
  for i=start, finish do begin
     pload, i,/silent
          
     for jj=y1, y2 do begin
        for ii=r1, r2 do begin

           bigR = x1(ii)*sin(x2(jj))
           cs = smallh*r0^(-0.5)*(r0/bigR)^(smallq/2.0)

           vz(ii,jj,*) = vx1(ii,jj,*)*cos(x2(jj)) -vx2(ii,jj,*)*sin(x2(jj))
           data(ii,jj) = mean(abs(vz(ii,jj,*)/cs))

           vz_over_vtot(ii,jj,*) = abs(vz(ii,jj,*))/sqrt(vx1(ii,jj,*)^2+vx2(ii,jj,*)^2+vx3(ii,jj,*)^2)

        endfor
     endfor
     
     if keyword_set(log) then data = alog10(data)
     
     if keyword_set(zslice) then begin
                                ;average vmach
        energy_mid = mean(data(r1:r2,zgrid+1:nx2-1))
                                ;avearge vmach
        energy_atm = mean(data(r1:r2, 0:zgrid))
        vznorm = max(vz_over_vtot(r1:r2, 0:nx2-1))
        energy(0,i-start) = time(i)
        energy(1,i-start) = energy_mid
        energy(2,i-start) = energy_atm
        print, 'time, Mz mid, Mz atm, avg Mzmid, avg mz atm, max vz', time(i), energy_mid, energy_atm, mean(energy(1,0:i-start)), mean(energy(2,0:i-start)), vznorm
        
     endif
     
  ;;    if not keyword_set(plotrange) then begin
  ;;       plotrange0 = [min(data(r1:r2,y1:y2)),max(data(r1:r2,y1:y2))]
  ;;    endif else begin
  ;;       plotrange0 = plotrange
  ;;    endelse
  ;;    levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0
     
  ;;    name2 = string(i, format='(I03)')
  ;;    title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
  ;;    title  = title1+string(mode,format='(I1)')
  ;;    set_plot,'ps'
  ;;    device, filename=strcompress('pdisk_vmach_'+name2+'.ps',/remove_all) $
  ;;            , bits_per_pixel=8, xsize=18, ysize=9,/color
     
  ;;    contour,data,rad,z1,levels=levels, xstyle=1, /fill $
  ;;            , xtitle=xtitle, ytitle=ytitle,charsize=1.0,xmargin=[7.0,10.0], ymargin=[4,2], $
  ;;            xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title ;, xminor=5
  ;;    colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange0,format='(f5.2)'
     
  ;;    device,/close
  endfor
  
  file =strcompress('pdisk_vmach.dat',/remove_all)
  openw,1,file
  for i=0, finish-start do begin
     printf,1,energy(0:2,i),format='(3(e22.15,x))'
  endfor
  close,1
end
