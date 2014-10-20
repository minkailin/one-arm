function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

function bigH, bigR, rref, hsmall, qsmall

height = hsmall*rref*(rref/bigR)^(0.5d0*(qsmall-3d0))

return, height
end

pro stream, loc=loc, start=start, finish=finish, azislice=azislice $
            ,xtickinterval=xtickinterval,length=length,red=red $
            , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
            , yrange=yrange ,  hsize=hsize, en=en, name=name  $
            , log=log, plotrange0=plotrange0, ct=ct, arrcol=arrcol, arrct=arrct $
            , nopert=nopert,noflux=noflux, mp=mp, basic=basic, smallq=smallq, novect=novect, nonaxi=nonaxi
  
  if not keyword_set(length) then length= 1.0
  if not keyword_set(ct) then ct=5
  if not keyword_set(red) then begin
     red = [nrad, ntheta]
  endif

;get the basic info
  location =strcompress(loc,/remove_all)
  if not keyword_set(basic) then begin
     filename = strcompress('hdfaa.000',/remove_all)
  endif else begin
     nme = string(basic,format='(I03)')
     filename = strcompress('hdfaa.'+nme,/remove_all)
  endelse
  fileloc  = filepath(filename,root_dir='.',subdir=[location])
  if keyword_set(nonaxi) then filename = strcompress('hdfaa.'+string(start,format='(I03)'),/remove_all) 
  
  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  z1 = 1d0/(smallh*tan(theta))
  
  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
  azi1 = azislice*nphi

  data0   = hdf(19, fileloc)                       ;initial density field
  if keyword_set(en) then data0 = hdf(23, fileloc) ;initial energy density
  
  vx   = dblarr(nrad, ntheta)
  vy   = dblarr(nrad, ntheta)
  data_axisymmetric = dblarr(nrad, ntheta)
  den2d = dblarr(nrad, ntheta)
  
 
  if keyword_set(basic) then begin
     for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
           data_axisymmetric(i,j) = mean(data0(i,j,*))
        endfor
     endfor
     for k=0, nphi-1 do begin
        data0(*,*,k) = data_axisymmetric(*,*)
     endfor
  endif
  
  if not keyword_set(finish) then finish = start

; planet info (or bump)
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
     den  = hdf(19, fileloc)
     if keyword_set(en) then den = hdf(23, fileloc)
     
     if not keyword_set(noflux) then begin
        vrad   *= den
        vtheta *= den
     endif
     
     plx=planetinfo(1,n)
     ply=planetinfo(2,n)
     plz=planetinfo(3,n)
     plrad = sqrt(plx^2 + ply^2 + plz^2)
     if keyword_set(mp) then begin
        rh = plrad*(mp/3d0)^(1d0/3d)
        rplot = (rad - plrad)/rh
        xtitle = '(r-r_p)/r_h'
     endif else begin
        rplot = rad/plrad
        xtitle = 'r/r_0'
     endelse 
     
;meridional flow
     for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
           vr = vrad(i,j,azi1)
           vt = vtheta(i,j,azi1)
           
           if keyword_set(mp) then begin
              vx(i,j) = vr/rh   ;divide by rh if doing scaled radius
           endif else begin
              vx(i,j) = vr/plrad
           endelse

           vy(i,j) =  -vt/(smallh*sin(theta(j))^2)
           vy(i,j)/= rad(i)     ;polar angular speed 
        endfor
     endfor
     
   if not keyword_set(nopert) then begin
      if not keyword_set(log) then den = den - data0
      den/=data0
   endif
   
   den2d(0:nrad-1,0:ntheta-1) = den(0:nrad-1,0:ntheta-1,azi1)  
   if keyword_set(log) then den2d=alog10(den2d)
   
   if not keyword_set(plotrange0) then begin
      temp = min(abs(xrange(0) - rplot(*)), r1)
      temp = min(abs(xrange(1) - rplot(*)), r2)
      temp = min(abs(yrange(0) - z1), t2)
      temp = min(abs(yrange(1) - z1), t1)
      plotrange=[min(den2d(r1:r2,t1:t2)),max(den2d(r1:r2,t1:t2))]
   endif else plotrange=plotrange0
   levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
   
   
   time=string(planetinfo(0,n)/torb,format='(F7.2)')
   azislicestring = string(azislice,format='(F5.3)')
   name1 = azislicestring 

   loadct, ct,/silent
   set_plot, 'ps'
   device, filename=filepath(strcompress('stream_'+ks+'.ps',/remove_all) $
                             ,root_dir='.',subdir=[location]), bits_per_pixel=8,xsize=18, ysize=9,/color
   contour, den2d, rplot, z1, /fill $
            ,title=time+textoidl('P_0, ')+textoidl('\phi/2\pi='+name1), levels=levels $ 
            ,ymargin=[4,2],xmargin=[7,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
            , ytickinterval=ytickinterval,xtitle=textoidl(xtitle), ytitle=textoidl('[h.tan\theta]^{-1}')
   colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange,format='(f5.2)'

   vx_small = congrid(vx, red(0), red(1))
   vy_small = congrid(vy, red(0), red(1))
   rad_small = congrid(rplot, red(0))
   zaxis = congrid(z1, red(1))
   
   temp = min(abs(rad_small - xrange(0)), r1)
   temp = min(abs(rad_small - xrange(1)), r2)
   temp = min(abs(zaxis - yrange(0)), t2)
   temp = min(abs(zaxis - yrange(1)), t1)

   r1 +=1
   r2 -=1

   ;loadct,2
   if not keyword_set(novect) then begin 
      if keyword_set(arrct) then loadct, arrct,/silent
      velovect2, vx_small(r1:r2,t1:t2), vy_small(r1:r2,t1:t2), rad_small(r1:r2), zaxis(t1:t2), color=arrcol,/overplot $
                 , length=length,/isotropic, hsize=hsize
   endif 

   if keyword_set(name) then begin
        xyouts, xrange(1)-0.2, 0.2, textoidl(name),charsize=2, charthick=6, color=255, align=1
   endif

   device,/close
  
   print, 'done '+ks
endfor
end
