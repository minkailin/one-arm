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

pro rossbyrz, loc=loc, start=start, finish=finish, azislice=azislice $
           , plotrange0=plotrange0 $
           ,xtickinterval=xtickinterval, nopert=nopert $
           , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
           , yrange=yrange, rotating=rotating, basic=basic $
            , linerange=linerange, mp=mp, name=name

if not keyword_set(finish) then finish = start

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
ztemp = 1.0/(smallh*tan(theta))

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)
nazi = fix(azislice*nphi)

vrad   = hdf(3, fileloc)
vtheta   = hdf(7, fileloc)
vphi   = hdf(11, fileloc)
data0 = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)

omega = dblarr(nrad, ntheta, nphi)
data2d = dblarr(nrad, ntheta)

for i=0, nrad-1 do begin
    for j=0,ntheta-1 do begin
        avg = mean(data0(i,j,*))
        data0(i,j,*) = avg
        
        om = mean(vphi(i,j,*))/(rad(i)*sin(theta(j)))
        data0(i,j,*) /= 2.0*om
    endfor
endfor

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
   data = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)

   plx=planetinfo(1,n)
   ply=planetinfo(2,n)
   plz=planetinfo(3,n)
   plrad = sqrt(plx^2 + ply^2 + plz^2) 
   if keyword_set(mp) then begin
      rhill = plrad*(mp/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
      xtitle=textoidl('(r-r_p)/r_h')
   endif else begin
      rplot = rad/plrad
      xtitle=textoidl('r/r_0')
   endelse
        
   if not keyword_set(rotating) then begin
      for i=0, nrad-1 do begin
         for j=0,ntheta-1 do begin
            om = mean(vphi(i,j,*))/(rad(i)*sin(theta(j)))
            omega(i,j,*) = om
         endfor
      endfor
      data /= 2.0*omega
   endif else begin
      temp = min(abs(rplot - rotating), r0)
      om = mean(vphi(r0,ntheta-1,*))/(rad(r0)*sin(theta(ntheta-1)))
      data /= 2.0*om
   endelse
   if not keyword_set(nopert) then data -= data0
   data2d(*,*) = data(*,*,nazi)


   if not keyword_set(plotrange0) then begin
      temp = min(abs(xrange(0) - rplot(*)), r1)
      temp = min(abs(xrange(1) - rplot(*)), r2)
      temp = min(abs(yrange(0) - ztemp ), t2)
      temp = min(abs(yrange(1) - ztemp ), t1)
      plotrange=[min(data2d(r1:r2,t1:t2)),max(data2d(r1:r2,t1:t2))]
   endif else plotrange=plotrange0
   levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
   
   time=string(planetinfo(0,n)/torb,format='(F7.2)')
   azislicestring = string(azislice,format='(F5.3)')

    loadct, 5, bottom=0
    set_plot, 'ps'
    device, filename=filepath(strcompress('rossbyrz_'+ks+'.ps',/remove_all) $
                              ,root_dir='.',subdir=[location]) $
            ,/color, bits_per_pixel=8,xsize=18, ysize=9
    
    contour, data2d,rplot, ztemp,levels=levels,/fill, title=time+' orbits, '+textoidl('\phi/2\pi='+azislicestring) $
             ,ymargin=[4,2],xmargin=[7,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
             , ytickinterval=ytickinterval,xtitle=xtitle, ytitle=textoidl('[h.tan\theta]^{-1}');, charsize=1.5
    ;colorbar, position=[0.895, 0.18, 0.93, 0.9],/vertical,/right,range=plotrange,format='(f5.2)'
    colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange,format='(f5.2)'
    
    if keyword_set(linerange) then begin
       levels2 = linerange(0)+(linerange(1)-linerange(0))*(dindgen(linerange(4))/(linerange(4) - 1))
       temp = min(abs(rplot - linerange(2)), r1)
       temp = min(abs(rplot - linerange(3)), r2)
       
       contour,data2d(r1:r2,*), rplot(r1:r2) , ztemp,levels=levels2,/overplot,c_color=255, thick=4
    endif
   
    if keyword_set(name) then begin
        xyouts, xrange(1)-0.2, 0.2, textoidl(name),charsize=2, charthick=6, color=255, align=1
    endif 
   
    device,/close
 
    print, 'done '+ks
 endfor
end
