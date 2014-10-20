function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro polarRZ, loc=loc, start=start, finish=finish, azislice=azislice $
           ,aziavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert $
           , hole=hole, smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, basic=basic, absolute=absolute, mp=mp, scale=scale


if not keyword_set(tuniv) then tuniv=1.0
if not keyword_set(scale) then scale=1.0

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


data2d = dblarr(nrad, ntheta)
data1d = dblarr(nrad)
data_axisymmetric = dblarr(nrad, ntheta)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]

case type of
    'vrad':   data0   = hdf(3, fileloc)
    'vtheta': data0   = hdf(7, fileloc)
    'vphi':   data0   = hdf(11, fileloc)
    'pot':    data0   = hdf(15, fileloc)
    'dens':   data0   = hdf(19, fileloc)

    'vz' : begin
        vrad  = hdf(3, fileloc)
        vtheta= hdf(7, fileloc)
        data0 = dblarr(nrad, ntheta, nphi)

        for k=0, nphi-1 do begin
            for i=0, nrad-1 do begin
                data0(i,*,k) = vrad(i,*,k)*cos(theta(*)) - vtheta(i,*,k)*sin(theta(*))
            endfor
        endfor

    end

 endcase

z1 = 1d0/(smallh*tan(theta))


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


    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2)
    if keyword_set(mp) then begin	
    rhill = (mp/3d0)^(1d0/3d0)
endif	else rhill = 0d0

    case type of
        'vrad':   data   = hdf(3, fileloc)
        'vtheta': data   = hdf(7, fileloc)
        'vphi':   data   = hdf(11, fileloc)
        'pot':    data   = hdf(15, fileloc)
        'dens':   data   = hdf(19, fileloc)

        'vz' : begin
        vrad  = hdf(3, fileloc)
        vtheta= hdf(7, fileloc)
        data = dblarr(nrad, ntheta, nphi)

        for k=0, nphi-1 do begin
            for i=0, nrad-1 do begin
                data(i,*,k) = vrad(i,*,k)*cos(theta(*)) - vtheta(i,*,k)*sin(theta(*))
            endfor
        endfor

    end


    endcase

    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif
    
    if keyword_set(azislice) then begin 
       if (azislice lt 0d0) then azislice = (pltphi(plx,ply) + azislice*rhill)/(2d0*!dpi) 
        nazislice = azislice*nphi
        data2d(0:nrad-1,0:ntheta-1) = data(0:nrad-1,0:ntheta-1,nazislice)
    endif
    
    if keyword_set(aziavg) then begin
        for i=0, nrad-1 do begin
            for j = 0, ntheta-1 do begin
                data2d(i,j) = mean(data(i,j,*))
            endfor
        endfor
    endif
    
    if keyword_set(absolute) then data2d = abs(data2d)
    if keyword_set(log) then data2d=alog10(data2d)
    

    data2d *= scale

    if not keyword_set(plotrange0) then begin
        temp = min(abs(xrange(0) - rad(*)), r1)
	  temp = min(abs(xrange(1) - rad(*)), r2)

          temp = min(abs(yrange(0) - z1), t2)
          temp = min(abs(yrange(1) - z1), t1)

        plotrange=[min(data2d(r1:r2,t1:t2)),max(data2d(r1:r2,t1:t2))]

    endif else plotrange=plotrange0


    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    azislicestring = string(azislice,format='(F5.3)')

    loadct, 5, bottom=0
    set_plot, 'ps'

    device, filename=filepath(strcompress('polarRZ_'+type+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=18, ysize=6


    contour, data2d, rad, z1, /fill,levels=levels $
      ,title=time+' orbits, '+textoidl('\phi/2\pi='+azislicestring) $ 
      ,ymargin=[3,3],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=ytickinterval,xtitle='R', ytitle=textoidl('z/H'); ,/isotropic

    colorbar, position=[0.895, 0.18, 0.93, 0.82],/vertical,/right,range=plotrange,format='(f5.2)'

;oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
device,/close

print, 'done '+ks
endfor


end
