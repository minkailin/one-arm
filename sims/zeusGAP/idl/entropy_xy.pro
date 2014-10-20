function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro entropy_xy, loc=loc, start=start, finish=finish, zslice=zslice $
           ,vertavg=vertavg, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, ytickinterval=ytickinterval, nopert=nopert, tuniv=tuniv $
           ,hole=hole, aziline=aziline, basic=basic, xrange=xrange, yrange=yrange, mp=mp $
           ,name=name, gmma=gmma, ct=ct, azimode=azimode, smallh=smallh, phi0=phi0
if not keyword_set(azimode) then begin
mmode=1.0
endif else mmode = azimode
if not keyword_set(ct) then ct=5
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


;construct array to hold cell edge radius

;re = dblarr(nrad+1)
;re(0) = rdomain(0)
;re(nrad) = rdomain(1)
;for i=0, nrad-2 do begin
;    re(i+1) = 2d0*rad(i) - re(i)
;endfor


; assume uniform theta and phi spacing

dtheta = theta(1) - theta(0)
dphi = phi(1) - phi(0)

data_axisymmetric = dblarr(nrad, ntheta)
data2d = dblarr(nrad, nphi)
data1d = dblarr(nrad)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]
dataplot=dblarr(nrad, nphi)

den0  = hdf(19, fileloc)
eng0  = hdf(23, fileloc)


    for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(den0(i,j,*))
        endfor
    endfor
    for k=0, nphi-1 do begin
        den0(*,*,k) = data_axisymmetric(*,*)
    endfor


     for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(eng0(i,j,*))
        endfor
    endfor
    for k=0, nphi-1 do begin
        eng0(*,*,k) = data_axisymmetric(*,*)
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

    den = hdf(19, fileloc)
    eng = hdf(23, fileloc)

    data = (eng/eng0)*(den/den0)^(-gmma) - 1d0


    if keyword_set(smallh) then begin
    csq = gmma*(gmma-1d0)*eng/den
    for j=0, ntheta-1 do begin
    for i=0, nrad -1 do begin
    cylindrad = rad(i)*sin(theta(j))
    csq(i,j,*) /= smallh^2/cylindrad
    endfor
    endfor
    print, 'max csq is', max(csq)
    print, 'min csq is', min(csq)
    endif


    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2)

    if not keyword_set(phi0) then begin
        phiplt = pltphi(plx,ply)
    endif else begin
       phiplt = phi0*2d0*!dpi
       ;normalize entropy pert to same as lin calc.
       temp = min(abs(phiplt-phi), azi1)
       temp = min(abs(rad-plrad), r1)

       w0 = (gmma-1d0)*(eng-eng0)/den0
       wnorm = abs(w0(r1, ntheta-1, azi1))
       data *= (gmma-1d0)*eng0
       data /= den0*wnorm
    endelse

    if keyword_set(mp) then begin
    rhill = plrad*(mp/3d0)^(1d0/3d0)
    rplot = (rad - plrad)/rhill
    xtitle=textoidl('(r-r_p)/r_h')
    endif else begin
    rplot = rad/plrad
    xtitle=textoidl('r/r_0')
    endelse


 ;   if not keyword_set(nopert) then begin
 ;       data = data - data0
        ;data/=data0
;    endif
    
    if keyword_set(zslice) then begin 
        nzslice = fix(zslice*ntheta)
        data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
        for i = 0, nrad - 1 do data1d(i) = mean(data2d(i,*))
    endif



    if not keyword_set(plotrange0) then begin
         temp = min(abs(rplot - xrange(0)),r1)
    temp = min(abs(rplot - xrange(1)),r2)
    plotrange=[min(data2d(r1:r2,*)),max(data2d(r1:r2,*))]
    endif else plotrange=plotrange0


    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    
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

    title  = time+textoidl('P_0')
    ytitle = textoidl('(\phi-\phi_0)/\pi')
    if keyword_set(mp) then begin
    title += textoidl(', r_p='+string(plrad,format='(f5.2)'))
    ytitle = textoidl('(\phi-\phi_p)/\pi')
    endif

    if keyword_set(azimode) then begin 
        ytitle = 'm' + ytitle
        ytitle = strcompress(ytitle,/remove_all)
    endif

    set_plot, 'ps'
    device, filename=filepath(strcompress('entropyxy_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
      ,/color, bits_per_pixel=8,xsize=12, ysize=14
    contour,dataplot,rplot, mmode*(phi/!dpi - 1.0),/fill,levels=levels,title=title, xstyle=2, $
      xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
      xtickinterval=xtickinterval, ytickinterval=ytickinterval, xrange=xrange, yrange=yrange, xminor=5

;   contour,dataplot,rplot,phi/!dpi,/fill,levels=levels,title=time+textoidl('P_0, ')+ $ 
;      textoidl('r_p='+string(plrad,format='(f5.2)')), xstyle=2, $
;      xtitle=textoidl('r'), ytitle=textoidl('\phi/\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
;      xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=4

    colorbar, position=[0.851, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
    oplot,[0,0],[0.0,0.0]/!dpi,psym=7,symsize=1.5,color=!D.Table_size*1.1
    if keyword_set(aziline) then begin
     for j=0, n_elements(aziline)-1 do begin
      oplot,[-1,1]*1d2,[aziline(j),aziline(j)], color=!D.Table_size, thick=4
      print, 'azislice=', (aziline(j)*!dpi + phiplt)/(2d0*!dpi)
     endfor
    endif  

   if keyword_set(name) then begin
        xyouts, xrange(1), -0.95, textoidl(name),charsize=1.5, charthick=6, color=255, align=1
    endif



  device,/close
    print, 'done '+ks

endfor

end
