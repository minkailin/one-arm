function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro sigma, loc=loc, start=start, finish=finish $
           ,log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert $
           ,hole=hole, basic=basic

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
data1d = dblarr(nrad)

data2d = dblarr(nrad, nphi)
data2d_basic = dblarr(nrad, nphi)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]
dataplot=dblarr(nrad+2, nphi)

data0   = hdf(19, fileloc)

;if keyword_set(basic) then begin
    for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(data0(i,j,*))
        endfor
    endfor
    for i=0, nrad-1 do begin
        r = rad(i)
        data1d(i) = total(data_axisymmetric(i,0:ntheta-1)*r*dtheta)
        data2d_basic(i,0:nphi-1) = data1d(i)
    endfor
;endif

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

    data   = hdf(19, fileloc)

    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2)

    for k=0, nphi-1 do begin
        for i=0, nrad-1 do begin
            r = rad(i)
            data2d(i,k) = total(data(i,0:ntheta-1,k)*r*dtheta)*2d0
        endfor
    endfor

    if not keyword_set(nopert) then begin
        data2d = data2d - data2d_basic
        data2d /= data2d_basic + 1d-15
    endif
    
    if keyword_set(log) then data2d=alog10(data2d)

    if not keyword_set(plotrange0) then begin
        plotrange=[min(data2d),max(data2d)]
    endif else plotrange=plotrange0

    if keyword_set(hole) then begin
        for j=0, nphi-1 do begin
            for i=0, nrad-1 do begin
                if(rad(i) le hole) then data2d(i,j) = plotrange(1)*10
            endfor
        endfor 
    endif

    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F7.2)')    
    loadct, 5, bottom=0
    set_plot, 'ps'
    device, filename=filepath(strcompress('sigma_'+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=14, ysize=12 
    polar_contour,transpose(data2d),phi,rad,/isotropic,/fill,levels=levels,title=time+' orbits' $
      ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=xtickinterval
    colorbar, position=[0.85, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
    device,/close
    
    print, 'done '+ks
endfor
end
