function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro stream_mach, loc=loc, start=start, finish=finish, azislice=azislice $
           ,xtickinterval=xtickinterval $
           , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, plotrange0=plotrange0, ct=ct $
            , mp=mp, gmma=gmma $
            , psi=psi, azimode=azimode
  
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
  
  machz = dblarr(nrad, ntheta, nphi)
  machz_2d= dblarr(nrad, ntheta)
  angle = dblarr(nrad, ntheta)
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
    
     den  = hdf(19, fileloc)
     pres = (gmma-1d0)*hdf(23, fileloc)
     cs = sqrt(gmma*pres/den)

     vrad  = hdf(3, fileloc)
     vtheta= hdf(7, fileloc)
     
     plx=planetinfo(1,n)
     ply=planetinfo(2,n)
     plz=planetinfo(3,n)
     plrad = sqrt(plx^2 + ply^2 + plz^2)
     if keyword_set(mp) then begin
        rh = plrad*(mp/3d0)^(1d0/3d)
        rplot = (rad - plrad)/rh
        xtitle = '(r-r_p)/r_h'
        temp = min(abs(rplot), rgrid)
     endif else begin
        rplot = rad/plrad
        xtitle = 'r/r_0'
        temp = min(abs(rplot - 1d0), rgrid)
     endelse 
     
    
     azi1 = azislice*nphi
     if keyword_set(psi) then begin
        azi1 = psi*!dpi/azimode + phi(azi1) 
        if (azi1 lt 0d0) then azi1 += 2d0*!dpi
        azi1*=nphi/(2d0*!dpi)  
     endif     
         

     for k=0, nphi-1 do begin
        for j=0, ntheta-1 do begin
           sint = sin(theta(j))
           cost = cos(theta(j))
           for i=0, nrad-1 do begin                          
              
              vr = vrad(i,j,k)
              vt = vtheta(i,j,k)
              
              vely =   vr*cost - sint*vt 
              velx =   vr*sint + vt*cost
              
              machz(i,j,k) = vely/cs(i,j, k)
           endfor
        endfor
     endfor
        
     mach_2d = machz(*,*,azi1)

     for j=0, ntheta-1 do begin
        sint = sin(theta(j))
        cost = cos(theta(j))
        for i=0, nrad-1 do begin                          
           
           vr = vrad(i,j,azi1)
           vt = vtheta(i,j,azi1)
           
           vely =   vr*cost - sint*vt 
           velx =   vr*sint + vt*cost
           
           angle(i,j) = sqrt(vely^2/(vely^2+velx^2))
        endfor
     endfor

     z1 = 1d0/(tan(theta)*smallh)
     if not keyword_set(plotrange0) then begin
        temp = min(abs(xrange(0) - rplot(*)), r1)
        temp = min(abs(xrange(1) - rplot(*)), r2)
        
        temp = min(abs(yrange(0) - z1), t2)
        temp = min(abs(yrange(1) - z1), t1)
        
        plotrange=[min(mach_2d(r1:r2,t1:t2)),max(mach_2d(r1:r2,t1:t2))]
        
     endif else plotrange=plotrange0
     levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
     

     time=string(planetinfo(0,n)/torb,format='(F7.2)')  
     loadct, ct,/silent  
     set_plot,'ps'
     fname = strcompress('stream_mach_'+ks+'.ps', /remove_all)
     device, filename=filepath(fname,root_dir='.',subdir=[location]) $
             ,/color, bits_per_pixel=8,ysize=4.5, xsize=8,xoffset=0,yoffset=0,/inches
     contour, mach_2d, rplot, z1 ,title=textoidl('M_z') $
              ,xmargin=[6,6],ymargin=[4,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
              , ytickinterval=ytickinterval,xtitle=textoidl('r_{sph}/r_0'), ytitle=textoidl('[h_{iso}tan\theta]^{-1}') $
              , charsize = 1.5,/fill, levels=levels
     colorbar, position=[0.91, 0.2, 0.94, 0.9],/vertical,/right,range=plotrange,format='(f5.2)'
     device,/close
     print, 'done '+ks

     
     print, 'avg theta', mean(angle(r1:r2, t1:t2))

     ;density-weighted vertical mach number in the chosen 2D plane
     
     mean_machz_2d = mean(den(r1:r2,t1:t2,azi1)*abs(mach_2d(r1:r2,t1:t2)))/mean(den(r1:r2,t1:t2,azi1))

     print, 'mean Mz_2D=', mean_machz_2d

     ;density-weighted vertical mach number in 3D shell
     
     mean_machz_3d = mean(den(r1:r2,t1:t2,*)*abs(machz(r1:r2,t1:t2,*)))/mean(den(r1:r2,t1:t2,*))

     print, 'mean Mz_3D=', mean_machz_3d

  endfor
stop
end
