function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro polarPZ, loc=loc, start=start, finish=finish, rslice=rslice $
           ,aziavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert $
           , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, basic=basic, absolute=absolute, scale=scale
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

data2d = dblarr(nphi, ntheta)
dataplot = dblarr(nphi, ntheta)
data_axisymmetric = dblarr(nrad, ntheta)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]

case type of
    'vrad':   data0   = hdf(3, fileloc)*scale
    'vtheta': data0   = hdf(7, fileloc)*scale
    'vphi':   data0   = hdf(11, fileloc)*scale
    'pot':    data0   = hdf(15, fileloc)*scale
    'dens':   data0   = hdf(19, fileloc)*scale
endcase

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
;what radius to slice?

temp     = min(abs(rad-rslice), grid)
radslice = rad(grid)

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
    phiplt = pltphi(plx,ply)
    plrad = sqrt(plx^2 + ply^2 + plz^2)
    
    case type of
        'vrad':   data   = hdf(3, fileloc)*scale
        'vtheta': data   = hdf(7, fileloc)*scale
        'vphi':   data   = hdf(11, fileloc)*scale
        'pot':    data   = hdf(15, fileloc)*scale
        'dens':   data   = hdf(19, fileloc)*scale
    endcase

    if not keyword_set(nopert) then begin
        data = data - data0
        data/=data0
    endif
    
    data2d(0:nphi-1,0:ntheta-1) = transpose(data(grid,0:ntheta-1,0:nphi-1))
    
    if keyword_set(absolute) then data2d = abs(data2d)
    if keyword_set(log) then data2d=alog10(data2d)
    
    if not keyword_set(plotrange0) then begin
        plotrange=[min(data2d),max(data2d)]
    endif else plotrange=plotrange0



    if(phiplt gt !dpi) then begin
        temp2=min(abs(phi-(phiplt-!dpi)),grid2)
        dataplot(0:nphi-1-grid2,0:ntheta-1) = data2d(grid2:nphi-1,0:ntheta-1)
        dataplot(nphi-grid2:nphi-1,0:ntheta-1) = data2d(0:grid2-1,0:ntheta-1)
    endif
    if(phiplt lt !dpi) then begin
        temp2=min(abs(phi-(phiplt+!dpi)),grid2)
        dataplot(nphi-grid2:nphi-1,0:ntheta-1) = data2d(0:grid2-1,0:ntheta-1)
        dataplot(0:nphi-1-grid2,0:ntheta-1) = data2d(grid2:nphi-1,0:ntheta-1)
    endif    
    if(phiplt eq !dpi) then dataplot = data2d


    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    rslicestring = string(radslice,format='(F4.2)')

    loadct, 5, bottom=0
    set_plot, 'ps'

    device, filename=filepath(strcompress('polarPZ_'+type+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=18, ysize=6

    contour, dataplot, phi/!dpi - 1.0, cos(theta)/smallh,/fill,levels=levels $
      ,title=time+' orbits, '+textoidl('r='+rslicestring) $
      ,ymargin=[4,2],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=xtickinterval,xtitle=textoidl('(\phi - \phi_p)/\pi'), ytitle=textoidl('z/H')

;    polar_contour, transpose(data2d), !dpi/2d0 - theta,rad,/fill,levels=levels,title=time+' orbits, '+textoidl('\phi/2\pi='+azislicestring) $
;      ,ymargin=[3,3],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
;      , ytickinterval=ytickinterval,xtitle='R', ytitle=textoidl('z')
    colorbar, position=[0.895, 0.24, 0.93, 0.89],/vertical,/right,range=plotrange,format='(f5.2)'

;oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
device,/close

print, 'done '+ks
endfor
end
