function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro stream_xy2, loc=loc, start=start, finish=finish, zslice=zslice $
                ,log=log, plotrange0=plotrange0 $
                ,xtickinterval=xtickinterval, ytickinterval=ytickinterval,nopert=nopert $
                , basic=basic, xrange=xrange, yrange=yrange, mp=mp $
                , length=length, hsize=hsize, red=red $
                , ct=ct, smallh=smallh, phi0=phi0, arrcol=arrcol, arrct=arrct, nonaxi=nonaxi
  
  if not keyword_set(ct) then ct=5
  type='dens'

  loadct,ct,/silent
  
;get the basic info
  location =strcompress(loc,/remove_all)
  if not keyword_set(basic) then begin
     filename = strcompress('hdfaa.000',/remove_all)
  endif else begin
     name = string(basic,format='(I03)')
     filename = strcompress('hdfaa.'+name,/remove_all)
  endelse
  if keyword_set(nonaxi) then filename = strcompress('hdfaa.'+string(start,format='(I03)'),/remove_all) 
  fileloc  = filepath(filename,root_dir='.',subdir=[location])

  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  phinew =(phi/!dpi - 1.0)


  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
  
; assume uniform theta and phi spacing
  
  dtheta = theta(1) - theta(0)
  dphi = phi(1) - phi(0)
  
  data_axisymmetric = dblarr(nrad, ntheta)
  data2d = dblarr(nrad, nphi)
  vx     = dblarr(nrad, nphi)
  vy     = dblarr(nrad, nphi)


  dataplot=dblarr(nrad, nphi)
  vxplot  = dblarr(nrad, nphi)
  vyplot  = dblarr(nrad, nphi)

  d0   = hdf(19, fileloc)
  v0=  hdf(11, fileloc)

  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(d0(i,j,*))
     endfor
  endfor
  for k=0, nphi-1 do begin
     d0(*,*,k) = data_axisymmetric(*,*)
  endfor
      
  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(v0(i,j,*))
     endfor
  endfor
  for k=0, nphi-1 do begin
     v0(*,*,k) = data_axisymmetric(*,*)
  endfor

if not keyword_set(finish) then finish = start

; planet info
nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
planetinfo = dblarr(7,nlines)
openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
readf,1,planetinfo
close,1
torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)

nzslice = fix(zslice*ntheta)
height  = 1d0/(smallh*tan(theta(nzslice)))
print, 'height=', height

for n=start, finish do begin
    ks   = string(n,format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])

    den  = hdf(19, fileloc)
    vrad = hdf(3, fileloc) 
    vphi = hdf(11, fileloc) 
    vphi-= v0
    
    if not keyword_set(nopert) then begin
       den /= d0
       if keyword_set(log) then begin
          den = alog10(den)
       endif else den -= 1.0
    endif
   
    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2) 
    if not keyword_set(phi0) then begin
       phiplt = pltphi(plx,ply)
    endif else begin
       phiplt = phi0*2d0*!dpi
       temp = min(abs(phiplt-phi), azi1)
       temp = min(abs(rad-plrad), r1)
    endelse
    
    if keyword_set(mp) then begin
       rhill = plrad*(mp/3d0)^(1d0/3d0)
       rplot = (rad - plrad)/rhill
       xtitle=textoidl('(r-r_p)/r_h')
    endif else begin
       rplot = rad/plrad
       xtitle=textoidl('r/r_0')
    endelse
    
    data2d(0:nrad-1,0:nphi-1) = den(0:nrad-1,nzslice,0:nphi-1)
    vx(0:nrad-1,0:nphi-1)     = vrad(0:nrad-1,nzslice,0:nphi-1)
    vy(0:nrad-1,0:nphi-1)     = vphi(0:nrad-1,nzslice,0:nphi-1)
    
    if keyword_set(mp) then begin
       vx /= rhill
    endif else begin
       vx /= plrad
    endelse
  
    for i=0, nrad-1 do begin
       vy(i,*) /= rad(i)        ;this is domega
       vy(i,*) *= !dpi
    endfor   
    
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    if(phiplt gt !dpi) then begin
       temp2=min(abs(phi-(phiplt-!dpi)),grid2)
       
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1,grid2:nphi-1)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1, 0:grid2-1)
       
       vxplot(0:nrad-1,0:nphi-1-grid2)    = vx(0:nrad-1,grid2:nphi-1)
       vxplot(0:nrad-1,nphi-grid2:nphi-1) = vx(0:nrad-1, 0:grid2-1)
       
       vyplot(0:nrad-1,0:nphi-1-grid2)     = vy(0:nrad-1,grid2:nphi-1)
       vyplot(0:nrad-1,nphi-grid2:nphi-1) = vy(0:nrad-1, 0:grid2-1)
       
    endif
    if(phiplt lt !dpi) then begin
       temp2=min(abs(phi-(phiplt+!dpi)),grid2)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1,0:grid2-1)
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1, grid2:nphi-1)
       
       vxplot(0:nrad-1,nphi-grid2:nphi-1) = vx(0:nrad-1,0:grid2-1)
       vxplot(0:nrad-1,0:nphi-1-grid2) = vx(0:nrad-1, grid2:nphi-1)
       
       vyplot(0:nrad-1,nphi-grid2:nphi-1) = vy(0:nrad-1,0:grid2-1)
       vyplot(0:nrad-1,0:nphi-1-grid2) = vy(0:nrad-1, grid2:nphi-1)
    endif    
    if(phiplt eq !dpi) then begin
       dataplot = data2d
       vxplot = vx
       vyplot = vy
    endif
    
    title  = time+textoidl('P_0')
    ytitle = textoidl('(\phi-\phi_0)/\pi')
    if keyword_set(mp) then begin
       ;; title += textoidl(', r_p='+string(plrad,format='(f5.2)'))
       ytitle = textoidl('(\phi-\phi_p)/\pi')
    endif


    heightstring = string(height,format='(F4.1)')
    title += textoidl(',[h.tan\theta]^{-1}='+heightstring)

    if not keyword_set(plotrange0) then begin
       temp = min(abs(rplot - xrange(0)),r1)
       temp = min(abs(rplot - xrange(1)),r2)
       temp = min(abs(phinew - yrange(0)), t1)
       temp = min(abs(phinew - yrange(1)), t2)
       plotrange=[min(dataplot(r1:r2,t1:t2)),max(dataplot(r1:r2,t1:t2))]
    endif else plotrange=plotrange0
    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
    
    set_plot, 'ps'
    device, filename=filepath(strcompress('streamxy2_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
            ,/color, bits_per_pixel=8,xsize=4.5, ysize=8, xoffset=0, yoffset=0,/inches
    contour,dataplot,rplot, (phi/!dpi - 1.0),/fill,levels=levels,title=title, xstyle=2, $
            xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[6,4], $
            xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5, ytickinterval=ytickinterval
    colorbar, position=[0.84, 0.18, 0.89, 0.88],/vertical,/right,range=plotrange,format='(f5.2)'
    
    vxplot = congrid(vxplot, red(0), red(1))
    vyplot = congrid(vyplot, red(0), red(1))
    rplot  = congrid(rplot, red(0))
    phinew    = congrid(phinew, red(1))
 
    temp = min(abs(rplot - xrange(0)), r1)
    temp = min(abs(rplot - xrange(1)), r2)
    temp = min(abs(phinew - yrange(0)), t1)
    temp = min(abs(phinew - yrange(1)), t2)

    r1 +=1
    r2 -=1
    t1 +=1
    t2 -=1

    if keyword_set(arrct) then loadct, arrct,/silent
    velovect, vxplot(r1:r2,t1:t2), vyplot(r1:r2,t1:t2), rplot(r1:r2), phinew(t1:t2), color=arrcol,/overplot $
              , length=length,/isotropic, thick=2
    
    device,/close
    print, 'done '+ks
    
 endfor
stop
end
