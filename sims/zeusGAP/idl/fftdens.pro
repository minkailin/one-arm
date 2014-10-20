function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro fftdens, loc=loc, mmax=mmax, start=start, finish=finish, zslice=zslice, range=range, azimodes=azimodes, max=max, yrange=yrange 

common consts, pi
pi = !dpi
!p.font=0
if not keyword_set(azimodes) then begin
if not keyword_set(mmax) then mmax = 4
azimodes = dindgen(mmax) + 1
endif 

nmodes = n_elements(azimodes)

mstring = string(nmodes)
if not keyword_set(finish) then finish = start

;get the basic info
location =strcompress(loc,/remove_all)
filename = strcompress('hdfaa.000',/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])
data   = hdf(19, fileloc)   ; density

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)

data2d = dblarr(nrad, nphi)
data_fft_time = dblarr(finish-start + 1, nmodes+1)

; planet info
nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
planetinfo = dblarr(7,nlines)
openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
readf,1,planetinfo
close,1
torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)

openw,1,filepath('fft_time.dat',root_dir='.',subdir=[location])
format = strcompress('('+mstring+'(e22.15,1x))',/remove_all)

temp = min(abs(range(0)-rad),r1)
temp = min(abs(range(1)-rad),r2)

nzslice = zslice*ntheta

;normalizations
    data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)

    result0 = fft(data2d, -1, dimension=2,/double) 

    if not keyword_set(max) then begin
    normalisation = dcomplex(int_tabulated(rad(r1:r2), real_part(result0(r1:r2,0))), $
                             int_tabulated(rad(r1:r2), imaginary(result0(r1:r2,0))))
    normalisation = abs(normalisation)
    endif


for n=start, finish do begin
    ks   = string(n,format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])
    data   = hdf(19, fileloc)   ; density
    data_fft_time(n-start,0) = planetinfo(0,n)/torb
    
    data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
 
    result = fft(data2d, -1, dimension=2,/double)
 
    if keyword_set(max) then begin
    for i=0, nrad-1 do result(i,*) /= result0(i,0)
    endif 

    for m=0, nmodes-1 do begin
        mode = azimodes(m) 

        if not keyword_set(max) then begin
        amplitude = dcomplex(int_tabulated(rad(r1:r2), real_part(result(r1:r2,mode))), $
                             int_tabulated(rad(r1:r2), imaginary(result(r1:r2,mode))))
        amplitude = alog10(abs(amplitude)/normalisation)       
        endif else begin
        amplitude = max(alog10(abs(result(r1:r2,mode))))
        endelse

        data_fft_time(n-start,m+1) = amplitude
    endfor

    printf,1,data_fft_time(n-start,0:nmodes),format=format
    print, 'done '+ks
endfor
close,1

if not keyword_set(max) then begin
ytitle=textoidl('log_{10}(C_m/C_0)')
endif else begin 
ytitle=textoidl('max(log_{10}|\Sigma_m/\Sigma_0|)')
endelse 

set_plot, 'ps'
device, filename=filepath(strcompress('fftdens_.ps',/remove_all),root_dir='.',subdir=[location]) $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, data_fft_time(*,0), data_fft_time(*,1),xmargin=[8,2],ymargin=[3,1], ystyle=0  $
  ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('t/P_0'), ytickinterval=ytickinterval  $
  ,xtickinterval=xtickinterval, ytitle=textoidl('log_{10}(C_m/C_0)'), yrange=yrange

for m=2, nmodes do begin
    oplot, data_fft_time(*,0), data_fft_time(*,m),thick=4,linestyle=m-1
endfor
device,/close
end











