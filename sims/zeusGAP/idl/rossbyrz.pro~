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


function get_vorticity_phi, vrad, vtheta, vphi, rad, theta, phi
nrad = n_elements(rad)
ntheta = n_elements(theta)
nphi = n_elements(phi)

dvr_dtheta = vrad
dvtheta_dr = vtheta

for k=0, nphi-1 do begin
    for j=0, ntheta-1 do begin
        dvtheta_dr(*,j,k) = deriv(rad(*), vtheta(*,j,k)) + vtheta(*,j,k)/rad(*)
    endfor

    for i=0, nrad-1 do begin
        dvr_dtheta(i,*,k) = deriv(theta(*), vrad(i,*,k))/rad(i)
    endfor
endfor

omegaphi = dvtheta_dr - dvr_dtheta

return, omegaphi

end


pro vortrz, loc=loc, start=start, finish=finish, azislice=azislice $
           , log=log, plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert, inv=inv $
           , hole=hole, smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, absolute=absolute, scale=scale, type=type $
            , linerange=linerange, mp=mp, name=name


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

vphi0   = hdf(11, fileloc)
den0   = hdf(19, fileloc)

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)

ztemp = 1d0/(tan(theta)*smallh)


data2d = dblarr(nrad, ntheta)

vrad0 = dblarr(nrad, ntheta, nphi)
vtheta0 = dblarr(nrad, ntheta, nphi)

for k=0, nphi-1 do begin
      vphi0(*,*,k) = vphi0(*,*,10)
      den0(*,*,k) = den0(*,*,10)
endfor

if(type eq 'omz') then begin 
    data0 = get_vorticity(vrad0, vtheta0, vphi0, rad, theta, phi)
    data0 /= den0
    
    if keyword_set(inv) then    data0 = 1d0/(data0 + 1d-15)
    
endif

if(type eq 'omphi') then data0 = get_vorticity_phi(vrad0, vtheta0, vphi0, rad, theta, phi)


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
    
    if(type eq 'omz') then begin
       data = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)
       data /= den
       
       if keyword_set(inv) then data = 1d0/(data + 1d-15)
       
    endif
    if(type eq 'omphi') then data = get_vorticity_phi(vrad, vtheta, vphi, rad, theta, phi)
    
    if not keyword_set(nopert) then begin
       if not keyword_set(log) then data = data - data0
       data/=data0
    endif
    
    if keyword_set(azislice) then begin
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
    
    if keyword_set(mp) then begin
      plx=planetinfo(1,n)
      ply=planetinfo(2,n)
      plz=planetinfo(3,n)
      plrad = sqrt(plx^2 + ply^2 + plz^2)
      rhill = plrad*(mp/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
      xtitle = textoidl('(R-r_p)/r_h')
   endif else begin
      xtitle = 'R'
   endelse
   

    if keyword_set(absolute) then data2d = abs(data2d)
    if keyword_set(log) then data2d=alog10(data2d)
    data2d *= scale


    if not keyword_set(plotrange0) then begin

       temp = min(abs(xrange(0) - rplot(*)), r1)
       temp = min(abs(xrange(1) - rplot(*)), r2)
       temp = min(abs(yrange(0) - ztemp ), t2)
       temp = min(abs(yrange(1) - ztemp ), t1)


        plotrange=[min(data2d(r1:r2,t1:t2)),max(data2d(r1:r2,t1:t2))]
    endif else plotrange=plotrange0

; 

    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
    
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    azislicestring = string(azislice,format='(F5.3)')
;
    loadct, 5, bottom=0
    set_plot, 'ps'
    device, filename=filepath(strcompress('vortRZ_'+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
            ,/color, bits_per_pixel=8,xsize=18, ysize=9
    
    contour, data2d,rplot, ztemp,levels=levels,/fill, title=time+' orbits, '+textoidl('\phi/2\pi='+azislicestring) $
             ,ymargin=[4,2],xmargin=[8,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
             , ytickinterval=ytickinterval,xtitle=xtitle, ytitle=textoidl('z/H');, charsize=1.5
    colorbar, position=[0.895, 0.18, 0.93, 0.9],/vertical,/right,range=plotrange,format='(f5.2)'
    
    if keyword_set(linerange) then begin
       levels2 = linerange(0)+(linerange(1)-linerange(0))*(dindgen(linerange(4))/(linerange(4) - 1))
       temp = min(abs(rplot - linerange(2)), r1)
       temp = min(abs(rplot - linerange(3)), r2)
       
       contour,data2d(r1:r2,*), rplot(r1:r2) , ztemp,levels=levels2,/overplot,c_color=255, thick=4
    endif
    
    
    if keyword_set(name) then begin
        xyouts, xrange(1)-0.2, 0.2, textoidl(name),charsize=1.5, charthick=6, color=255, align=1
    endif 
    


    device,/close




    
    print, 'done '+ks
 endfor
end
