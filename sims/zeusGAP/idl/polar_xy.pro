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


pro polar_xy, loc=loc, start=start, finish=finish, zslice=zslice $
           ,vertavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, ytickinterval=ytickinterval, nopert=nopert, tuniv=tuniv $
           ,hole=hole, aziline=aziline, basic=basic, xrange=xrange, yrange=yrange, mp=mp $
           ,name=name, smallq=smallq, smallh=smallh, ct=ct, nonaxi=nonaxi, phi0=phi0, cart=cart, mmode=mmode

  if not keyword_set(finish) then finish = start
  if not keyword_set(ct) then ct =5 
  loadct,ct,/silent


  !p.font = 0

;get the basic info
  location =strcompress(loc,/remove_all)
  
  if not keyword_set(basic) then begin
     filename = strcompress('hdfaa.000',/remove_all)
  endif else begin
     nme = string(basic,format='(I03)')
     filename = strcompress('hdfaa.'+nme,/remove_all)
  endelse
  if keyword_set(nonaxi) then filename = strcompress('hdfaa.'+string(start,format='(I03)'),/remove_all)
  
  fileloc  = filepath(filename,root_dir='.',subdir=[location])
  
  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  
  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
 
  phi = 2d0*!dpi*dindgen(nphi)/(nphi-1d0)
 
; assume uniform theta and phi spacing
  
  dtheta = theta(1) - theta(0)
  dphi = phi(1) - phi(0)

  if keyword_set(mp) then begin
  azi = phi/!dpi - 1d0
  endif else azi = phi/(2d0*!dpi)

   if keyword_set(cart) then begin
      dnx1 = 2*nrad
      dx   = 2*rad(nrad-1)/(dnx1-1.0)
      xaxis = -rad(nrad-1) + dx*dindgen(dnx1)
      yaxis = xaxis
      dataxy = dblarr(dnx1, dnx1)
   endif

  
  data_axisymmetric = dblarr(nrad, ntheta)
  data2d = dblarr(nrad, nphi)
  dataplot=dblarr(nrad, nphi)

  
  case type of
     'vrad':   data0   = hdf(3, fileloc)
     'vtheta': data0   = hdf(7, fileloc)
     'vphi':   data0   = hdf(11, fileloc)
     'pot':    data0   = hdf(15, fileloc)
     'dens':   data0   = hdf(19, fileloc)
     'en':     data0   = hdf(23, fileloc)
  endcase
  
  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(data0(i,j,*))
     endfor
  endfor
  for k=0, nphi-1 do begin
     data0(*,*,k) = data_axisymmetric(*,*)
  endfor

  nzslice = fix(zslice*ntheta) 
  data02d = dblarr(nrad,nphi)
  data02d(0:nrad-1,0:nphi-1) = data0(0:nrad-1,nzslice,0:nphi-1)
  data0_fft =  fft(data02d, -1, dimension=2,/double) 

; planet info
  nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
  planetinfo = dblarr(7,nlines)
  openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
  readf,1,planetinfo
  close,1
  torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)
  
  height  = 1d0/(smallh*tan(theta(nzslice)))
  heightstring = string(height,format='(F4.1)')
  
  
  
  for n=start, finish do begin
     ks   = string(n,format='(I03)')
     filename = strcompress('hdfaa.'+ks,/remove_all)
     fileloc  = filepath(filename,root_dir='.',subdir=[location])
    
     case type of
        'vrad':   data   = hdf(3, fileloc)
        'vtheta': data   = hdf(7, fileloc)
        'vphi':   data   = hdf(11, fileloc)
        'pot':    data   = hdf(15, fileloc)
        'dens':   data   = hdf(19, fileloc)
        'en':   data   = hdf(23, fileloc)
     endcase

     if keyword_set(nonaxi) then begin
     for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(data(i,j,*))
     endfor
     endfor
     for k=0, nphi-1 do begin
     data0(*,*,k) = data_axisymmetric(*,*)
     endfor
     endif
  


     plx=planetinfo(1,n)
     ply=planetinfo(2,n)
     plz=planetinfo(3,n)
     phiplt = pltphi(plx,ply)
     plrad = sqrt(plx^2 + ply^2 + plz^2)
     

     if keyword_set(mp) then begin
        rhill = plrad*(mp/3d0)^(1d0/3d0)
        rplot = (rad - plrad)/rhill
        xtitle=textoidl('(r-r_p)/r_h')
     endif else begin
        rplot = rad/plrad
        xtitle=textoidl('r/r_0')
     endelse
     
     if not keyword_set(phi0) then begin
        phiplt = pltphi(plx,ply)
     endif else begin
        phiplt = phi0*2d0*!dpi
     endelse
     
     zmax = plrad*cos(theta(nzslice)) ;height at rp
     zmax/=  bigH(plrad, plrad, smallh, smallq)
     print, 'zmax/H at plrad =', zmax
     
     if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
     endif
     
     data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
  

    if keyword_set(mmode) then begin ;replace density field by a particular fourier component
    data_fft = fft(data2d, -1, dimension=2,/double)
    for j=0,nphi-1 do begin
    for i=0, nrad-1 do begin
    data2d(i,j) = min([mmode+1,2.0])*real_part( data_fft(i,mmode)*( cos(mmode*phi(j)) + dcomplex(0,1d0)*sin(mmode*phi(j)) ) )
    data2d(i,j)/= abs(data0_fft(i,0))
    endfor
    endfor
    endif

    if keyword_set(log) then data2d=alog10(data2d)

    if not keyword_set(plotrange0) then begin
       temp = min(abs(rplot - xrange(0)),r1)
       temp = min(abs(rplot - xrange(1)),r2)
       plotrange=[min(data2d(r1:r2,*)),max(data2d(r1:r2,*))]
    endif else plotrange=plotrange0
    
    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F6.1)')
    
    title  = time+textoidl('P_0')

    if keyword_set(mp) then begin
       ;title += textoidl(', r_p='+string(plrad,format='(f5.2)'))
       if(phiplt gt !dpi) then begin
       temp2=min(abs(phi-(phiplt-!dpi)),grid2)
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1,grid2:nphi-1)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1, 0:grid2-1)
       endif
       if(phiplt lt !dpi) then begin
       temp2=min(abs(phi-(phiplt+!dpi)),grid2)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1,0:grid2-1)
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1, grid2:nphi-1)
       endif
       if(phiplt eq !dpi) then dataplot = data2d
       ytitle = textoidl('(\phi-\phi_p)/\pi')
    endif else begin
       dataplot = data2d
       ytitle = textoidl('\phi/2\pi')
    endelse
    
    title = title+strcompress(textoidl(',(h.tan\theta)^{-1}='+heightstring),/remove_all)

    if not keyword_set(cart) then begin
    set_plot, 'ps'
    device, filename=filepath(strcompress('polarxy_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
            ,/color, bits_per_pixel=8,xsize=12, ysize=14
    contour,dataplot,rplot,azi ,/fill,levels=levels,title=title, xstyle=1, $
            xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
            xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5, ytickinterval=ytickinterval
;    xyouts, -2, -0.96, 'ZEUS 3D', charsize=2,charthick=8, color=255 
    colorbar, position=[0.865, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
    oplot,[0,0],[0.0,0.0]/!dpi,psym=7,symsize=1.5,color=!D.Table_size*1.1
    if keyword_set(aziline) then begin
       for j=0, n_elements(aziline)-1 do begin
          oplot,[-1,1]*1d2,[aziline(j),aziline(j)], color=!D.Table_size, thick=4, linestyle=1
          print, 'azislice=', (aziline(j)*!dpi + phiplt)/(2d0*!dpi)
       endfor
    endif  
    
    if keyword_set(name) then begin
       xyouts, xrange(1), -0.95, textoidl(name),charsize=2, charthick=6, color=255, align=1
    endif
    endif else begin

;    for jj=0, dnx1-1 do begin
;         y = yaxis(jj)
;         for ii=0, dnx1-1 do begin
;            x = xaxis(ii)
;            r_t = sqrt(x^2 + y^2)
;            azi_t = pltphi(x, y)
;            if( (r_t ge xrange(0)) and (r_t le rad(nrad-1)) )then begin
;
;               temp = min(abs(rad - r_t),   x0)
;               temp = min(abs(phi - azi_t),  y0)
;
;                if x0 eq nrad-1 then begin
;                dr   = rad(x0) - rad(x0-1)
;                endif else if x0 eq 0 then begin
;                dr   = rad(x0+1) - rad(x0)
;                endif else begin
;                dr   = abs(rad(x0-1) - rad(x0+1))/2d0
;                endelse
;
;               ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
;               jp = y0 + (azi_t - phi(y0))/(dphi/2.0) + 0.0
;
;               dataxy(ii,jj) = bilinear(dataplot, ip, jp)
;            endif else begin
;               dataxy(ii,jj) = 1d10
;            endelse
;         endfor
;      endfor
;  

      dataxy = transpose(dataplot)
      temp = min(abs(rad - xrange(0)),x1)
      for i=0, x1 do dataxy(*,i) = 10d0*abs(plotrange(1))


          set_plot, 'ps'
      device, filename=filepath(strcompress('polarxy2_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location]) $
              ,/color, bits_per_pixel=8,xsize=14, ysize=12
;      contour, dataxy, xaxis, yaxis, /isotropic,/fill,levels=levels,title=title $
;               ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, xrange=[-1,1]*xrange(1), yrange=[-1,1]*xrange(1), xstyle=1, ystyle=1
      polar_contour, dataxy, phi, rad, /isotropic,/fill,levels=levels,title=title $
               ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, xrange=[-1,1]*xrange(1), yrange=[-1,1]*xrange(1), xstyle=1, ystyle=1,/dither

    colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange,format='(f5.2)'
    endelse 
   device,/close
     

    print, 'done '+ks
    
 endfor
end
