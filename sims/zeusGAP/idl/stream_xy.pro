function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro stream_xy, loc=loc, start=start, finish=finish, zslice=zslice $
               , type=type, log=log, plotrange0=plotrange0 $
               ,xtickinterval=xtickinterval, ytickinterval=ytickinterval,nopert=nopert $
               ,hole=hole, aziline=aziline, basic=basic, xrange=xrange, yrange=yrange, mp=mp $
               ,name=name, length=length, hsize=hsize, arrcolor=arrcolor, red=red $
               ,gmma=gmma, ct=ct, azimode=azimode, smallh=smallh, phi0=phi0, acol=acol, act=act
  
  if not keyword_set(ct) then ct=5
  if not keyword_set(azimode) then begin
     mmmode=1.0
  endif else mmode=azimode
  loadct,ct,/silent
  
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
  p0   = (gmma-1d0)*hdf(23, fileloc)
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
        data_axisymmetric(i,j) = mean(p0(i,j,*))
     endfor
  endfor
  for k=0, nphi-1 do begin
     p0(*,*,k) = data_axisymmetric(*,*)
  endfor
  
  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(v0(i,j,*))
     endfor
  endfor
  for k=0, nphi-1 do begin
     v0(*,*,k) = data_axisymmetric(*,*)
  endfor

  csq0 = gmma*p0/d0

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
    press= (gmma-1d0)* hdf(23, fileloc) 
    vrad = hdf(3, fileloc) 
    vphi = hdf(11, fileloc) 

    bigQ = csq0*(den - d0)/d0
    bigW = (press - p0)/d0
    vphi-= v0
   
    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    phiplt = pltphi(plx,ply)
    plrad = sqrt(plx^2 + ply^2 + plz^2)
   
        if not keyword_set(phi0) then begin
       phiplt = pltphi(plx,ply)
    endif else begin
       phiplt = phi0*2d0*!dpi
                                ;normalize pert to same as lin calc.
       temp = min(abs(phiplt-phi), azi1)
       temp = min(abs(rad-plrad), r1)
      
       wnorm = abs(bigW(r1, ntheta-1, azi1))
       bigQ /= wnorm
       bigW /= wnorm
       vrad   /= wnorm
       vphi   /= wnorm 
    endelse
     


    if keyword_set(mp) then begin
       rhill = plrad*(mp/3d0)^(1d0/3d0)
       rplot = (rad - plrad)/rhill
       xtitle=textoidl('(r-r_p)/r_h')
    endif else begin
       rplot = rad/plrad
       xtitle=textoidl('r_{sph}/r_0')
    endelse
    
    if keyword_set(zslice) then begin 
       nzslice = fix(zslice*ntheta)
       height  = 1d0/(smallh*tan(theta(nzslice)))
       print, 'height=', height
       if(type eq 'Q') then data2d(0:nrad-1,0:nphi-1) = bigQ(0:nrad-1,nzslice,0:nphi-1)
       if(type eq 'W') then data2d(0:nrad-1,0:nphi-1) = bigW(0:nrad-1,nzslice,0:nphi-1)
       vx(0:nrad-1,0:nphi-1)     = vrad(0:nrad-1,nzslice,0:nphi-1)
       vy(0:nrad-1,0:nphi-1)     = vphi(0:nrad-1,nzslice,0:nphi-1)
    endif
 
    if keyword_set(mp) then begin
       vx /= rhill
    endif else begin
       vx /= plrad
    endelse
    
    for i=0, nrad-1 do begin
       vy(i,*) /= rad(i)        ;this is domega
       vy(i,*) *= mmode/!dpi
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
      
    ;; title  = time+textoidl('P_0')
     ytitle = textoidl('m(\phi-\phi_0)/\pi')
     if keyword_set(mp) then begin
    ;; title += textoidl(', r_p='+string(plrad,format='(f5.2)'))
     ytitle = textoidl('(\phi-\phi_p)/\pi')
     endif


    heightstring = string(height,format='(F4.1)')
    title = textoidl('[h_{iso}tan\theta]^{-1}='+heightstring)
    title = title + ' (nonlinear)'

    if not keyword_set(plotrange0) then begin
       temp = min(abs(rplot - xrange(0)),r1)
       temp = min(abs(rplot - xrange(1)),r2)
        phinew = mmode*(phi/!dpi - 1.0)
       temp = min(abs(phinew - yrange(0)), t1)
       temp = min(abs(phinew - yrange(1)), t2)
       plotrange=[min(dataplot(r1:r2,t1:t2)),max(dataplot(r1:r2,t1:t2))]
    endif else plotrange=plotrange0
    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
   
    

    set_plot, 'ps'
    device, filename=filepath(strcompress('streamxy_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
            ,/color, bits_per_pixel=8,xsize=4.5, ysize=8, xoffset=0, yoffset=0,/inches
    contour,dataplot,rplot,mmode*(phi/!dpi - 1.0),/fill,levels=levels,title=title, xstyle=2, $
            xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[6,4], $
            xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5, ytickinterval=ytickinterval
    colorbar, position=[0.84, 0.18, 0.89, 0.88],/vertical,/right,range=plotrange,format='(f5.2)'
    
    vxplot = congrid(vxplot, red(0), red(1))
    vyplot = congrid(vyplot, red(0), red(1))
    rplot  = congrid(rplot, red(0))
    phi    = congrid(phi, red(1))
 
   
    phinew = mmode*(phi/!dpi - 1.0)

    temp = min(abs(rplot - xrange(0)), r1)
    temp = min(abs(rplot - xrange(1)), r2)

    if not keyword_set(yrange) then yrange=[-1,1]

    temp = min(abs(phinew - yrange(0)), t1)
    temp = min(abs(phinew - yrange(1)), t2)

    if keyword_set(act) then loadct, act,/silent
    t1+=1
    t2-=1
    velovect, vxplot(r1:r2,t1:t2), vyplot(r1:r2,t1:t2), rplot(r1:r2), phinew(t1:t2), color=acol,/overplot $
               , length=length,/isotropic, thick=2

;    oplot,[0,0],[0.0,0.0]/!dpi,psym=7,symsize=1.5,color=!D.Table_size*1.1

    ;; if keyword_set(aziline) then begin
    ;;  for j=0, n_elements(aziline)-1 do begin
    ;;   oplot,[-1,1]*1d2,[aziline(j),aziline(j)], color=!D.Table_size, thick=4
    ;;   print, 'azislice=', (aziline(j)*!dpi + phiplt)/(2d0*!dpi)
    ;;  endfor
    ;; endif  

   ;; if keyword_set(name) then begin
   ;;      xyouts, xrange(1), -0.95, textoidl(name),charsize=1.5, charthick=6, color=255, align=1
   ;;  endif



  device,/close
    print, 'done '+ks

endfor

end
