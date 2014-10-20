function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro polar, loc=loc, start=start, finish=finish, zslice=zslice $
           ,vertavg=vertavg, type=type, log=log, plotrange0=plotrange0, rmax=rmax $
           ,xtickinterval=xtickinterval, nopert=nopert, tuniv=tuniv, xrange=xrange $
           ,hole=hole, basic=basic, smallh=smallh, fudge=fudge, aziline=aziline, ct=ct $
           , name=name, notime=notime, scale=scale, r0=r0, nonaxi=nonaxi
if not keyword_set(fudge) then fudge = 0.0
if not keyword_set(ct) then ct = 5
if not keyword_set(scale) then scale=1.0
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

angle= dblarr(nphi)

dtheta = theta(1) - theta(0)
dphi = phi(1) - phi(0)

data_axisymmetric = dblarr(nrad, ntheta)
data2d = dblarr(nrad, nphi)
data1d = dblarr(nrad)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]
dataplot=dblarr(nrad+2, nphi)

case type of
    'vrad':   data0   = hdf(3, fileloc)
    'vtheta': data0   = hdf(7, fileloc)
    'vphi':   data0   = hdf(11, fileloc)
    'pot':    data0   = abs(hdf(15, fileloc))
    'dens':   data0   = hdf(19, fileloc)
    'en':   data0   = hdf(23, fileloc)
endcase

    for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(data0(i,j,*))
        endfor
    endfor
    for k=0, nphi-1 do begin
        data0(*,*,k) = data_axisymmetric(*,*)
    endfor

;stop

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

    case type of
        'vrad':   data   = hdf(3, fileloc)
        'vtheta': data   = hdf(7, fileloc)
        'vphi':   data   = hdf(11, fileloc)
        'pot':    data   = abs(hdf(15, fileloc))
        'dens':   data   = hdf(19, fileloc)
        'en':    data   = hdf(23, fileloc)
    endcase


    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2)

;    rhill = plrad*(mp/3d0)^(1d0/3d0)
;    print, 'hill radius is', rhill

;    mass = 0.0
;    hillmass = 0.0
;    for k=0, nz-1 do begin
;        for j=0, ntheta-1 do begin
;            for i=0, nrad-1 do begin
;                dr= re(i+1) - re(i)
;                dm = data(i,j,k)*rad(i)^2*sin(theta(j))*dr*dtheta*dphi

;                xx = rad(i)*sin(theta(j))*cos(phi(k))
;                yy = rad(i)*sin(theta(j))*sin(phi(k))
;                zz = rad(i)*cos(theta(j))

;               distp = (xx - plx)^2 + (yy - ply)^2 + (zz - plz)^2
;                distp = sqrt(distp)

;                if(rad(i)*sin(theta(j)) gt hole) then mass += dm
;                if(distp lt rhill) then hillmass += dm

;            endfor
;        endfor
;    endfor

;    print, 'disc mass', mass
;    print, 'mass inside Rh', hillmass



    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif
    
    if keyword_set(zslice) then begin
        nzslice = zslice*ntheta
        height = 1d0/(smallh*tan(theta(nzslice)))
        data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)*scale
        for i = 0, nrad - 1 do data1d(i) = mean(data2d(i,*))
    endif
    
    if keyword_set(vertavg) then begin
        for i=0, nrad-1 do begin
            for k = 0, nphi-1 do begin
                data2d(i,k) = mean(data(i,*,k))
            endfor
            data1d(i) = mean(data2d(i,*))
        endfor
    endif
    
    if keyword_set(log) then data2d=alog10(data2d)

;data2(*,0:in-1)=0.
;data2(*,in:nrad+in-1)=data
    if not keyword_set(plotrange0) then begin
        temp = min(abs(rad-hole),r1)
        if keyword_set(rmax) then begin 
        temp = min(abs(xrange(1)*r0 - rad),r2)
        endif else begin
         r2 = nrad-1
       endelse
        plotrange=[min(data2d(r1:r2,*)),max(data2d(r1:r2,*))]
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
    heightstring = string(height,format='(F5.2)')
  
    title = time+textoidl('P_0,') + textoidl('(h.tan\theta)^{-1}='+heightstring)    
    if keyword_set(notime) then title = textoidl('(h.tan\theta)^{-1}='+heightstring)


;plx=info(1,k)
;ply=info(2,k)
;plrad=sqrt(plx*plx+ply*ply)
phip=pltphi(plx,ply)*(1d0)
;temp=min(abs(azi-phi),grid)
;azi1(0:nsec-1)=azi(0:nsec-1)-phi;-azi(grid)
    loadct, ct, bottom=0
    set_plot, 'ps'
;    dataplot(2:nrad+1,0:nphi-1) = data2d(0:nrad-1, 0:nphi-1)
;    dataplot(0:1, 0:nphi-1) = 10.0

    device, filename=filepath(strcompress('polar_'+type+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=14, ysize=12 
;temp=min(abs(radnew-10.0), grid);
;data2(*,grid:nrad-1) = abs(max(data2))*2.0
    polar_contour,transpose(data2d),phi,rad/r0,/isotropic,/fill,levels=levels,title=title $
      ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=xrange, xtickinterval=xtickinterval $
      , ytickinterval=xtickinterval, charsize=1.
    colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[plrad,plrad]/r0,[0.,0.],psym=7,symsize=1.25,color=512,thick=4


    if keyword_set(aziline) then begin
        num = n_elements(aziline)
        for i = 0, num-1 do begin
            angle(*) = 2d0*!dpi*aziline(i)
            oplot, rad, angle, thick=2, color=255,/polar
        endfor
    endif

     if keyword_set(name) then begin
        xyouts, 0, -20./r0, textoidl(name),charsize=1.5, charthick=6, color=255, alignment=0.5
    endif


device,/close
print, 'done '+ks
endfor
end
