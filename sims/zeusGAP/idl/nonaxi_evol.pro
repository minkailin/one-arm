function hdf, n, name
  
fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data
end

pro nonaxi_evol, loc=loc, start=start, finish=finish $
                 ,xtickinterval=xtickinterval, ytickinterval=ytickinterval $
                 ,azimodes=azimodes, xrange=xrange, yrange=yrange, rrange=rrange $
                 ,legend=legend, label=label, basic=basic, rdisk=rdisk 
  
  !p.font = 0
  
  nmodes = n_elements(azimodes)
  if not keyword_set(finish) then finish = start
 
  if not keyword_set(basic) then basic=start

;get the basic info
  location =strcompress(loc,/remove_all)
  filename = strcompress('hdfaa.'+string(basic,format='(I03)'),/remove_all)
  fileloc  = filepath(filename,root_dir='.',subdir=[location])
  
  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  
  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
  
  phi = 2d0*!dpi*dindgen(nphi)/(nphi - 1d0) 

;assume uniform log spacing in radius
  dlogr = alog(rdisk(1)/rdisk(0))/nrad
  dr    = rad*dlogr 
;assume uniform thet and phi spacing
  dtheta = theta(1) - theta(0)
  dphi   = phi(1) - phi(0) 

  bigR_3D     = dblarr(nrad,ntheta,nphi)
  bigR_axisym = dblarr(nrad,ntheta)
  dVol_3D     = dblarr(nrad,ntheta,nphi)
  dVol_axisym = dblarr(nrad,ntheta)

  ;set up geometric arrays 
  for kk=0, nphi-1 do begin
     for jj=0, ntheta-1 do begin
        for ii=0, nrad-1 do begin
           bigR_3D(ii,jj,kk)    = rad(ii)*sin(theta(jj))
           dVol_3D(ii,jj,kk)    = 2d0*rad(ii)^2*sin(theta(jj))*dr(ii)*dtheta*dphi ; factor of 2 for lower disk 
        endfor
     endfor
  endfor
  
  for jj=0, ntheta-1 do begin
     for ii=0, nrad-1 do begin
        bigR_axisym(ii,jj)    = rad(ii)*sin(theta(jj))
        dVol_axisym(ii,jj)    = 2d0*rad(ii)^2*sin(theta(jj))*dr(ii)*dtheta*2d0*!dpi ; 2pi from phi-integration, factor of 2 in front for lower disk 
     endfor
  endfor
  
  time     = dblarr(finish-start+1)
  mode_amp = dblarr(nmodes, finish-start+1)
  mode_ang = dblarr(nmodes, finish-start+1)
  tot_ang  = dblarr(finish-start+1)
  
; planet info
  nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
  planetinfo = dblarr(7,nlines)
  openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
  readf,1,planetinfo
  close,1
  torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)
  
  if keyword_set(rrange) then begin
     temp = min(abs(rrange(0) - rad), r1)
     temp = min(abs(rrange(1) - rad), r2)
  endif else begin
     r1 = 0
     r2 = nrad-1
  endelse
  
;normalization for mode amplitude and ang mom amplitude 
  rho    = hdf(19, fileloc)
  vphi   = hdf(11, fileloc)

  vphi_fft = fft(vphi, -1, dimension=3,/double)
  rho_fft = fft(rho, -1, dimension=3,/double)

  czero  =total(real_part(rho_fft(*,*,0))*dVol_axisym)
;  angmom0=total(real_part(rho_fft(*,*,0)*vphi_fft(*,*,0))*bigR_axisym*dVol_axisym)

  temp = dblarr(nphi/2+1)
  for m=0, nphi/2 do begin
  data2d = real_part( rho_fft(*,*,m)*conj(vphi_fft(*,*,m)) )*bigR_axisym
  temp(m) = total(data2d*dVol_axisym)
  if (m ne 0) then temp(m) *= 2d0
  endfor
  angmom0 = total(temp)


 
  for n=start, finish do begin
     ks   = string(n,format='(I03)')
     filename = strcompress('hdfaa.'+ks,/remove_all)
     fileloc  = filepath(filename,root_dir='.',subdir=[location])
     
     time(n-start) = planetinfo(0,n)/torb  
     
     vphi     = hdf(11, fileloc)
     rho      = hdf(19, fileloc)
     
     vphi_fft = fft(vphi, -1, dimension=3,/double) 
     rho_fft = fft(rho, -1, dimension=3,/double)
     
     for mm=0, nmodes-1 do begin
        m = azimodes(mm) 
;mode amplitudes       
        if not keyword_set(max) then begin
           re  = total(real_part(rho_fft(*,*,m))*dVol_axisym(*,*))
           im  = total(imaginary(rho_fft(*,*,m))*dVol_axisym(*,*))
           amp = dcomplex(re,im)
           
           mode_amp(mm, n-start) = alog10(abs(amp)/czero + 1d-16)
        endif else begin
           mode_amp(mm, n-start) = max(alog10( abs(rho_fft(r1:r2,*,m)/rho0(r1:r2,*,0)) ) )
        endelse
        
;angular momentum components  
        jtot = total(real_part(rho_fft(*,*,m)*conj(vphi_fft(*,*,m)))*bigR_axisym(*,*)*dVol_axisym(*,*))
        if (m ne 0) then jtot *= 2d0 ;extra factor of 2 for non-axisymmetric modes 
        mode_ang(mm, n-start) = jtot
        
     endfor 
     
;total ang mom and total mass 
     mass  = total(real_part(rho_fft(*,*,0))*dVol_axisym(*,*))
     angmom= total(rho*vphi*bigR_3D*dVol_3D)
     
     tot_ang(n-start) = angmom
      
    print, n, tot_ang(n-start),total(mode_ang(*,n-start)) ,mass, format='(I03,x,3(e22.15,x))'
endfor

  if not keyword_set(max) then begin
  ytitle=textoidl('log_{10}(C_m/C_{0,t=0})')
  endif else begin
  ytitle=textoidl('max(log_{10}|\rho_m/\rho_{0,t=0}|)')
  endelse


  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol.ps',root_dir='.',subdir=[location]) $
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
  mode_ang(j,*) = mode_ang(j,*) - mode_ang(j,0)
  endfor
  mode_ang /= angmom0

  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol_ang.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time, mode_ang(0,*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\DeltaJ_m/J_{0,t=0}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(mode_ang),max(mode_ang)], xstyle=1
  for j=1, nmodes-1 do begin
  oplot, time, mode_ang(j,*), thick=4, linestyle=j
  endfor
  oplot, time, tot/angmom0-1., thick=1, linestyle=0
;   oplot, time, tot_ang/angmom0-1d0,thick=1,linestyle=0
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
  device, filename=filepath('nonaxi_evol_totj.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time(*), tot_ang(*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\Delta J_{tot}/J_{tot,0}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(tot_ang),max(tot_ang)], xstyle=1
  device,/close


end
