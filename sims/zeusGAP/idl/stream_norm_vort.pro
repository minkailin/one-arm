function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro stream_norm_vort, loc=loc, start=start, finish=finish, azislice=azislice $
           ,xtickinterval=xtickinterval $ 
           , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, basic=basic $
            , plotrange0=plotrange0, ct=ct, r0=r0

; plot as theta-radius plane
; but modified theta co-ordinate (call psi)
; psi = 1/(h*tan(theta))
; h is aspect ratio at r=r0 (ref. rad)
; the scalar to be plotted is delta W/rho, normalized by its value at
; r=r0, z=0
  
  if not keyword_set(length) then length= 1.0
  if not keyword_set(ct) then ct=5

;get the basic info
  location =strcompress(loc,/remove_all)
  if not keyword_set(basic) then begin
     filename = strcompress('hdfaa.000',/remove_all)
  endif else begin
     name = string(basic,format='(I03)')
     filename = strcompress('hdfaa.'+name,/remove_all)
  endelse
  fileloc  = filepath(filename,root_dir='.',subdir=[location])
  
  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  
  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
    
  vort = dblarr(nrad, ntheta)
  vort3d = dblarr(nrad, ntheta, nphi)
  dvphi_dtheta = dblarr(nrad, ntheta, nphi)
  dvtheta_dphi = dblarr(nrad, ntheta, nphi)
  dvrad_dphi     = dblarr(nrad, ntheta, nphi)
  dvphi_dr     = dblarr(nrad, ntheta, nphi)


  vrad  = hdf(3, fileloc)
     vtheta= hdf(7, fileloc)
     vphi  = hdf(11, fileloc)

   azi1 = azislice*nphi

   for k=0, nphi-1 do begin
       for i=0, nrad-1 do begin
       dvphi_dtheta(i,*,k) = deriv(theta, vphi(i,*,k))
       endfor
   endfor
 

     for j=0, ntheta-1 do begin
       for i=0, nrad-1 do begin
       dvtheta_dphi(i,j,*) = deriv(phi, vtheta(i,j,*))
       dvrad_dphi(i,j,*)   = deriv(phi, vrad(i,j,*))
       endfor
     endfor

      for k=0, nphi-1 do begin
       for j=0, ntheta-1 do begin
       dvphi_dr(*,j,k) = deriv(rad, vphi(*,j,k))
       endfor
endfor

   for k=0, nphi-1 do begin
    for j=0, ntheta-1 do begin
     sint = sin(theta(j))
     cost = cos(theta(j))
     for i=0, nrad-1 do begin

      vort3d(i,j,k) = sint*dvphi_dtheta(i,j,k) - cost*vphi(i,j,k) - dvtheta_dphi(i,j,k)
      vort3d(i,j,k)*= cost/(rad(i)*sint)

      vort3d(i,j,k)-=( dvrad_dphi(i,j,k)/sint - vphi(i,j,k) - rad(i)*dvphi_dr(i,j,k) )*sint/rad(i)

     endfor
    endfor
endfor


  
  if not keyword_set(finish) then finish = start
  
; planet info
  nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
  planetinfo = dblarr(7,nlines)
  openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
  readf,1,planetinfo
  close,1
  torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)
  
  for n=start, finish do begin
     ks   = string(n,format='(I03)')
     filename = strcompress('hdfaa.'+ks,/remove_all)
     fileloc  = filepath(filename,root_dir='.',subdir=[location])
    

     vrad  = hdf(3, fileloc)
     vtheta= hdf(7, fileloc)
     vphi  = hdf(11, fileloc)
     

       for i=0, nrad-1 do begin
       dvphi_dtheta(i,*,azi1) = deriv(theta, vphi(i,*,azi1))
       endfor

     for j=0, ntheta-1 do begin
       for i=0, nrad-1 do begin
       dvtheta_dphi(i,j,*) = deriv(phi, vtheta(i,j,*))
       dvrad_dphi(i,j,*)   = deriv(phi, vrad(i,j,*)) 
       endfor
     endfor

       for j=0, ntheta-1 do begin
       dvphi_dr(*,j,azi1) = deriv(rad, vphi(*,j,azi1)) 
       endfor


    for j=0, ntheta-1 do begin
     sint = sin(theta(j))
     cost = cos(theta(j)) 
     for i=0, nrad-1 do begin
 
      vort(i,j) = sint*dvphi_dtheta(i,j,azi1) - cost*vphi(i,j,azi1) - dvtheta_dphi(i,j,azi1)
      vort(i,j)*= cost/(rad(i)*sint)
   
      vort(i,j)-=( dvrad_dphi(i,j,azi1)/sint - vphi(i,j,azi1) - rad(i)*dvphi_dr(i,j,azi1) )*sint/rad(i)


      vort(i,j)-=mean(vort3d(i,j,*))
     endfor
    endfor


     z1 = 1d0/(tan(theta)*smallh)
     rplot= rad/r0
     if not keyword_set(plotrange0) then begin
        temp = min(abs(xrange(0) - rplot(*)), r1)
        temp = min(abs(xrange(1) - rplot(*)), r2)
        
        temp = min(abs(yrange(0) - z1), t2)
        temp = min(abs(yrange(1) - z1), t1)
        
        plotrange=[min(vort(r1:r2,t1:t2)),max(vort(r1:r2,t1:t2))]
        
     endif else plotrange=plotrange0
     levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
     
     time=string(planetinfo(0,n)/torb,format='(F7.2)')
     azislicestring = string(azislice,format='(F5.3)')
     name1 = azislicestring 
     
     loadct, ct,/silent
     
     fname = strcompress('stream_norm_vort_'+ks+'.ps',/remove_all)
     set_plot, 'ps'
     device, filename=filepath(fname,root_dir='.',subdir=[location]) $
             ,/color, bits_per_pixel=8,ysize=4.5, xsize=8,xoffset=0,yoffset=0,/inches
     contour, vort, rplot, z1,title=textoidl('\delta\omega_z') $
              ,xmargin=[6,6],ymargin=[4,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
              , ytickinterval=ytickinterval,xtitle=textoidl('r_{sph}/r_0'), ytitle=textoidl('[h_{iso}tan\theta]^{-1}') $
              , charsize = 1.5,/fill, levels=levels
     colorbar, position=[0.91, 0.2, 0.94, 0.9],/vertical,/right,range=plotrange,format='(f5.2)'
     device,/close

    
  endfor

end
