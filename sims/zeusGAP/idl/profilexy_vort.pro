function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro profilexy_vort, loc=loc, start=start, finish=finish, zslice=zslice $
             ,vertavg=vertavg, type=type, log=log, plotrange0=plotrange0 $
             ,xtickinterval=xtickinterval, ytitle=ytitle, scale=scale $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             ,nopert=nopert

;get the basic info
location =strcompress(loc,/remove_all)
filename = strcompress('hdfaa.000',/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

x  = hdf(2, fileloc)
y  = hdf(1, fileloc)
z  = hdf(0, fileloc)

dx = x(1) - x(0)
dy = y(1) - y(0)

nx  = n_elements(x)
ny  = n_elements(y)
nz  = n_elements(z)

dvxdy = dblarr(nx, ny, nz)
dvydx = dblarr(nx, ny, nz)

den  = dblarr(nx, ny, nz)
vort = dblarr(nx, ny, nz)

data2d = dblarr(nx,ny)
data_polar = dblarr(nx, 4*nx)
data1d = dblarr(nx)

rad = ((x(nx-1) - x(0))/2d0)*dindgen(nx)/(nx-1d0)
azi = 2d0*!dpi*dindgen(4*nx)/(4*nx-1d0)

if not keyword_set(finish) then finish = start
if not keyword_set(ytitle) then ytitle = ''

for n=start, finish do begin
    ks   = string(n,format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])
    
    vx = hdf(3, fileloc)
    vy = hdf(7, fileloc)
    
    for k=0, nz-1 do begin
        for j=0, ny-1 do begin
            dvydx(*,j,k) = deriv(x(*), vy(*,j,k))
        endfor
        for i=0, nx-1 do begin
            dvxdy(i,*,k) = deriv(y(*), vy(i,*,k))
        endfor
    endfor
    vort = dvydx - dvxdy
    den   = hdf(19, fileloc)
    vort /= den
    data = vort

    if keyword_set(zslice) then begin
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
;    if keyword_set(scale) then data1d *= scale


    for i=0, nx-1 do begin
        for j=0, 4*nx-1 do begin
            xtarg = rad(i)*cos(azi(j))
            ytarg = rad(i)*sin(azi(j))

            temp = min(abs(x - xtarg), grid)
            gridx = grid + (x(grid) - xtarg)/dx

            temp = min(abs(y - ytarg), grid)
            gridy = grid + (y(grid)- ytarg)/dy
            
            data_polar(i,j) = interpolate(data2d, gridx, gridy)
        endfor
        data1d(i) = mean(data_polar(i,*))
    endfor

    set_plot, 'ps'
    device, filename=filepath(strcompress('profile_vort_'+string(ks,format='(I03)')+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches    
    
    plot, rad, data1d,xmargin=[8,2],ymargin=[3,2], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle='r', ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle)
    device,/close
    
    print, 'done '+ks
endfor
end
