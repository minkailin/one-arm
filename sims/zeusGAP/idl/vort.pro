function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

function get_vorticity, vrad, vtheta, vphi, rad, theta, phi
nrad = n_elements(rad)
ntheta = n_elements(theta)
nphi = n_elements(phi)

dvphi_sintheta = vphi
dvtheta = vtheta

;get omega_r*cos(theta)
for i=0, nrad-1 do begin
    for k=0, nphi-1 do begin
	dvphi_sintheta(i,*,k) = deriv(theta(*), vphi(i,*,k)*sin(theta(*)))
        dvphi_sintheta(i,*,k) /= rad(i)*tan(theta(*))
    endfor
    
    for j=0, ntheta-1 do begin
	dvtheta(i,j,*) = deriv(phi(*), vtheta(i,j,*)) 
        dvtheta(i,j,*) /= rad(i)*tan(theta(j))
    endfor
endfor
omegar = dvphi_sintheta - dvtheta

;get omega_theta*sin(theta) !factors involving sin(theta) in the expression for the curl is absorbed with pre-factor sin(theta) here
dvrad = vrad
dr_vphi = vphi

for j=0, ntheta-1 do begin
	for i=0, nrad-1 do begin
          dvrad(i,j,*) = deriv(phi(*), vrad(i,j,*))/rad(i)
        endfor
     
    for k=0, nphi-1 do begin
        dr_vphi(*,j,k) = deriv(rad(*), rad(*)*vphi(*,j,k))
        dr_vphi(*,j,k) *= sin(theta(j))/rad(*)
    endfor
endfor
omegatheta = dvrad - dr_vphi

omegaz = omegar - omegatheta ; geometric factors already incoroporated. 

return, omegaz
end

pro vort, loc=loc, start=start, finish=finish, zslice=zslice $
           ,vertavg=vertavg, log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert, tuniv=tuniv $
           ,hole=hole, basic=basic, smallh=smallh $
          , absolute=absolute, scale=scale,inv=inv, aziline=aziline

;get the basic info
location =strcompress(loc,/remove_all)
if not keyword_set(scale) then scale=1.0
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


angle = dblarr(nphi)

; assume uniform theta and phi spacing

dtheta = theta(1) - theta(0)
dphi = phi(1) - phi(0)

data_axisymmetric = dblarr(nrad, ntheta)
data2d = dblarr(nrad, nphi)
data1d = dblarr(nrad)

radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]
dataplot=dblarr(nrad+2, nphi)

   vphi0   = hdf(11, fileloc)
   den0   = hdf(19, fileloc)

vtheta0=dblarr(nrad,ntheta,nphi)
vrad0  = dblarr(nrad,ntheta,nphi)

for k=0, nphi-1 do begin
      vphi0(*,*,k) = vphi0(*,*,10)
      den0(*,*,k) = den0(*,*,10)
endfor


data0 = get_vorticity(vrad0, vtheta0, vphi0, rad, theta, phi)
data0 /= den0

if keyword_set(inv) then data0 = 1d0/(data0)

if keyword_set(basic) then begin
    for k=0, nphi-1 do begin
        data0(*,*,k) = data0(*,*,10)
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


    vrad   = hdf(3, fileloc)
    vtheta   = hdf(7, fileloc)
    vphi   = hdf(11, fileloc)
    den   = hdf(19, fileloc)
    
    data = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)
    data /= den
    
    if keyword_set(inv) then data = 1d0/(data)

    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2)

    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif
    
    if keyword_set(zslice) then begin   
        nzslice = zslice*ntheta
        height = tan(!dpi/2d0 - theta(nzslice))/smallh
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
    
    if keyword_set(absolute) then data2d = abs(data2d)
    if keyword_set(log) then data2d=alog10(data2d)
    data2d *= scale

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
    heightstring = string(height,format='(F5.2)')
    phip=pltphi(plx,ply)
    loadct, 5, bottom=0
    set_plot, 'ps'

    device, filename=filepath(strcompress('vort_'+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
      ,/color, bits_per_pixel=8,xsize=14, ysize=12 

    polar_contour,transpose(data2d),phi,rad,/isotropic,/fill,levels=levels,title=time+' orbits, ' + textoidl('z/H='+heightstring) $
      ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=xtickinterval


    if keyword_set(aziline) then begin
        num = n_elements(aziline)
        for i = 0, num-1 do begin
            angle(*) = 2d0*!dpi*aziline(i)
            oplot, rad, angle, thick=2, color=255,/polar
        endfor
    endif


    colorbar, position=[0.85, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
device,/close

print, 'done '+ks
endfor
end
