function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro polarxy, loc=loc, start=start, finish=finish, zslice=zslice $
           ,vertavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert, xrange=xrange, yrange=yrange $
           , hole=hole, mp=mp, basic=basic, surfden=surfden, smallq=smallq

if not keyword_set(smallq) then smallq = 1.0

;get the basic info
location =strcompress(loc,/remove_all)

if not keyword_set(basic) then begin
    filename = strcompress('hdfaa.000',/remove_all)
endif else begin
    name = string(basic,format='(I03)')
    filename = strcompress('hdfaa.'+name,/remove_all)
endelse
fileloc  = filepath(filename,root_dir='.',subdir=[location])

x  = hdf(2, fileloc)
y  = hdf(1, fileloc)
z  = hdf(0, fileloc)


dx = x(1) - x(0)
dy = y(1) - y(0)
dz = z(1) - z(0)

print, "warning: assumes UNIFORM spacing"

nx  = n_elements(x)
ny  = n_elements(y)
nz  = n_elements(z)

data2d = dblarr(nx, ny)


data_polar=dblarr(nx, 4*nx, nz)
data1d = dblarr(nx, nz)

case type of
    'vx':   data0   = hdf(3, fileloc)
    'vy':   data0   = hdf(7, fileloc)
    'vz':   data0   = hdf(11, fileloc)
    'pot':    data0   = hdf(15, fileloc)
    'dens':   data0   = hdf(19, fileloc)
endcase

rad = ((x(nx-1) - x(0))/2d0)*dindgen(nx)/(nx-1d0)
dr = rad(1) - rad(0)
azi = 2d0*!dpi*dindgen(4*nx)/(4*nx-1d0)

if keyword_set(surfden) then begin
    sigma0 = dblarr(nx, ny)
    for j=0, ny-1 do begin
        for i=0, nx-1 do begin
            sigma0(i,j) = total(data0(i,j,0:nz-1))*dz
        endfor
    endfor
endif

;construct axisymmetric basic state to measure deviations from

; if keyword_set(basic) then begin
    
;     for k=0, nz-1 do begin
;         data2d(*,*) = data0(*,*,k)
;         for j=0, 4*nx-1 do begin
;             for i=0, nx-1 do begin
                
;                 xtarg = rad(i)*cos(azi(j))
;                 ytarg = rad(i)*sin(azi(j))
                
;                 temp = min(abs(x - xtarg), grid)
;                 gridx = grid + (x(grid) - xtarg)/dx
                
;                 temp = min(abs(y - ytarg), grid)
;                 gridy = grid + (y(grid)- ytarg)/dy
                
;                 data_polar(i,j,k) = interpolate(data2d, gridx, gridy)
;             endfor
;         endfor
;     endfor
    
;     for k=0, nz-1 do begin
;         for i=0, nx-1 do begin
;             data1d(i,k) = mean(data_polar(i,*,k))
;         endfor
;     endfor    
; endif


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
        'vx':   data   = hdf(3, fileloc)
        'vy':   data   = hdf(7, fileloc)
        'vz':   data   = hdf(11, fileloc)
        'pot':    data   = hdf(15, fileloc)
        'dens':   data   = hdf(19, fileloc)
    endcase
    

;    plx=planetinfo(1,n)
;    ply=planetinfo(2,n)
;    plz=planetinfo(3,n) 
;    plrad = sqrt(plx^2 + ply^2 + plz^2)
;    rhill = plrad*(mp/3d0)^(1d0/3d0)
;    print, 'hill radius is', rhill


;    mass = 0.0
;    hillmass = 0.0
;    for k=0, nz-1 do begin
;        for j=0, ny-1 do begin
;            for i=0, nx-1 do begin
;                dm = data(i,j,k)*dx*dy*dz
 
;                r = sqrt(x(i)^2 + y(j)^2 + z(k)^2)
;                distp = (x(i) - plx)^2 + (y(j) - ply)^2 + (z(k) - plz)^2
;                distp = sqrt(distp)
;                if(r gt hole) then mass += dm
;                if(distp lt rhill) then hillmass += dm
;
;            endfor
;        endfor
;    endfor

;    print, 'disc mass', mass
;    print, 'mass inside Rh', hillmass

    if keyword_set(surfden) then begin
        for j=0, ny-1 do begin
            for i=0, nx-1 do begin
                data2d(i,j) = total(data(i,j,0:nz-1))*dz
            endfor
        endfor
        if not keyword_set(nopert) then begin
            if not keyword_set(log) then data2d -= sigma0
            data2d /= sigma0 + 1d-16
        endif
        goto, skip
    endif
    
    if not keyword_set(nopert) then begin
;        if not keyword_set(basic) then begin
            if not keyword_set(log) then data = data - data0
            data/=data0
;        endif ; else begin

;             for k=0, nz-1 do begin

;                 for j=0, ny-1 do begin
;                     for i=0, nx-1 do begin
;                         rtarg = sqrt(x(i)^2 + y(j)^2)           
;                         temp = min(abs(rad - rtarg), grid)
;                         gridr = grid + (rtarg - rad(grid))/dr
                        
;                         bg = interpolate(data1d, gridr, k)
;                         data(i,j,k) = (data(i,j,k) - bg)/bg 
;                     endfor
;                 endfor

;             endfor
;         endelse

    endif
    
    if keyword_set(zslice) then begin
;        nzslice = (ntheta/2)*(zslice + 1d0) - 1d0
;        if(nzslice lt 0) then nzslice = 0    
        nzslice = zslice*nz
        data2d(0:nx-1,0:ny-1) = data(0:nx-1,0:ny-1,nzslice)
    endif
    
    if keyword_set(vertavg) then begin
        for i=0, nx-1 do begin
            for j = 0, ny-1 do begin
                data2d(i,j) = mean(data(i,j,*))
            endfor
        endfor
    endif

skip:
    
    if keyword_set(log) then data2d=alog10(data2d + 1d-16)

;data2(*,0:in-1)=0.
;data2(*,in:nrad+in-1)=data
    if not keyword_set(plotrange0) then begin
        if not keyword_set(xrange) then begin
            plotrange=[min(data2d),max(data2d)]
        endif else begin
            temp = min(abs(x-xrange(0)),x1)
            temp = min(abs(x-xrange(1)),x2)
            temp = min(abs(y-yrange(0)),y1)
            temp = min(abs(y-yrange(1)),y2)
            plotrange=[min(data2d(x1:x2, y1:y2)),max(data2d(x1:x2,y1:y2))]
        endelse
    endif else plotrange=plotrange0

    if keyword_set(hole) then begin
        for j=0, ny-1 do begin
            for i=0, nx-1 do begin
                radius = sqrt(x(i)^2 + y(j)^2)
                if(radius le hole) then data2d(i,j) = plotrange(1)*10
            endfor
        endfor
    endif

    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
;time=string(24,format='(F7.2)')
    
    loadct, 5, bottom=0
    set_plot, 'ps'
    device, filename=filepath(strcompress('polar_'+type+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=14, ysize=12 
    contour,data2d,x,y,/isotropic,/fill,levels=levels,title=time+' orbits' $
      ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=xtickinterval
    colorbar, position=[0.86, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
device,/close

print, 'done '+ks
endfor

stop

end
