function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro profile, loc=loc, start=start, finish=finish, zslice=zslice $
             ,vertavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
             ,xtickinterval=xtickinterval, ytitle=ytitle, scale=scale $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             ,nopert=nopert, basic=basic

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

data2d = dblarr(nrad, nphi)
data_axisymmetric = dblarr(nrad, ntheta)
data1d = dblarr(nrad)

case type of
    'vrad':   data0   = hdf(3, fileloc)
    'vtheta': data0   = hdf(7, fileloc)
    'vphi':   data0   = hdf(11, fileloc)
    'pot':    data0   = hdf(15, fileloc) 
    'dens':   data0   = hdf(19, fileloc)
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
if not keyword_set(ytitle) then ytitle = ''

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
    endcase

    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif
    
    if keyword_set(zslice) then begin
;        nzslice = (ntheta/2)*(zslice + 1d0) - 1d0
;        if(nzslice lt 0) then nzslice = 0    
        nzslice = fix(zslice*ntheta)
        if(nzslice eq ntheta) then nzslice -= 1
        data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
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
    if keyword_set(scale) then data1d *= scale
    if keyword_set(log) then data1d = alog10(data1d)

    set_plot, 'ps'
    device, filename=filepath(strcompress('profile_'+type+string(ks,format='(I03)')+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    if(type eq 'vphi') then data1d *= sqrt(rad)
 
;print, deriv(alog(rad),alog(data1d))
   if(type eq 'pot') then data1d = alog10(abs(deriv(rad,data1d))) 
    
    plot, rad, data1d,xmargin=[8,2],ymargin=[3,2], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle='r', ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle)
    device,/close
    
    print, 'done '+ks
endfor
end
