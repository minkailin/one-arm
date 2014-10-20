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

pro rossbyxy, loc=loc, start=start, finish=finish, zslice=zslice $
               , plotrange0=plotrange0, basic=basic $
              ,xtickinterval=xtickinterval, xrange=xrange,yrange=yrange,nopert=nopert $
               , mp=mp, phi0=phi0, ytickinterval=ytickinterval $
              , aziline=aziline, nonaxi=nonaxi, smallh=smallh, ct=ct, rotating=rotating, name=name
  
if not keyword_set(scale) then scale=1.0
if not keyword_set(finish) then finish=start
if not keyword_set(ct) then ct=5
loadct, ct, /silent

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
phinew =(phi/!dpi - 1.0)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)
nzslice = fix(zslice*ntheta)
height  = 1d0/(smallh*tan(theta(nzslice)))
print, 'height=', height

vrad   = hdf(3, fileloc)
vtheta   = hdf(7, fileloc)
vphi   = hdf(11, fileloc)
data0 = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)

omega = dblarr(nrad, ntheta, nphi)
data2d = dblarr(nrad, nphi)
dataplot=dblarr(nrad, nphi)

for i=0, nrad-1 do begin
    for j=0,ntheta-1 do begin
        avg = mean(data0(i,j,*))
        data0(i,j,*) = avg
        
        om = mean(vphi(i,j,*))/(rad(i)*sin(theta(j)))
        data0(i,j,*) /= 2.0*om
    endfor
endfor


if keyword_set(mp) then begin
   ytitle=textoidl('(\phi-\phi_p)/\pi')
endif else begin
   ytitle=textoidl('(\phi-\phi_0)/\pi')
endelse

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
    
    plx=planetinfo(1,n)
    ply=planetinfo(2,n)
    plz=planetinfo(3,n)
    plrad = sqrt(plx^2 + ply^2 + plz^2) 
    if not keyword_set(phi0) then begin
       phiplt = pltphi(plx,ply)
    endif else begin
       phiplt = phi0*2d0*!dpi
       temp = min(abs(phiplt-phi), azi1)
    endelse

    if keyword_set(mp) then begin
       rhill = plrad*(mp/3d0)^(1d0/3d0)
       rplot = (rad - plrad)/rhill
       xtitle=textoidl('(r-r_p)/r_h')
    endif else begin
       rplot = rad/plrad
       xtitle=textoidl('r/r_0')
    endelse
    
    data = get_vorticity(vrad, vtheta, vphi, rad, theta, phi)
    
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
    
    data2d(*,*) = data(*,nzslice,*)
        
    if(phiplt gt !dpi) then begin
       temp2=min(abs(phi-(phiplt-!dpi)),grid2)
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1,grid2:nphi-1)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1, 0:grid2-1)
    endif
    if(phiplt lt !dpi) then begin
       temp2=min(abs(phi-(phiplt+!dpi)),grid2)
       dataplot(0:nrad-1,nphi-grid2:nphi-1) = data2d(0:nrad-1,0:grid2-1)
       dataplot(0:nrad-1,0:nphi-1-grid2) = data2d(0:nrad-1, grid2:nphi-1)
    endif
    if(phiplt eq !dpi) then dataplot = data2d
    
    if not keyword_set(plotrange0) then begin
       temp = min(abs(rplot - xrange(0)),r1)
       temp = min(abs(rplot - xrange(1)),r2)
       temp = min(abs(phinew - yrange(0)), t1)
       temp = min(abs(phinew - yrange(1)), t2)
       plotrange=[min(dataplot(r1:r2,t1:t2)),max(dataplot(r1:r2,t1:t2))]
    endif else plotrange=plotrange0
    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
    
    time=string(planetinfo(0,n)/torb,format='(F7.2)')
    heightstring = string(height,format='(F4.1)')
    title = textoidl(time+'P_0')+textoidl(',[h.tan\theta]^{-1}='+heightstring)

    set_plot, 'ps'
    device, filename=filepath(strcompress('rossbyxy_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
            ,/color, bits_per_pixel=8,xsize=12, ysize=14
    contour,dataplot,rplot,phinew,/fill,levels=levels,title=title, xstyle=1, $
            xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
            xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5
    colorbar, position=[0.851, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
    oplot,[0,0],[0.0,0.0]/!dpi,psym=7,symsize=1.5,color=!D.Table_size*1.1

    if keyword_set(aziline) then begin
       for j=0, n_elements(aziline)-1 do begin
          oplot,[-1,1]*1d2,[aziline(j),aziline(j)], color=!D.Table_size, thick=4, linestyle=1
          print, 'azislice=', (aziline(j)*!dpi + phiplt)/(2d0*!dpi)
       endfor
    endif

    if keyword_set(name) then begin
       xyouts, xrange(1), -0.95, textoidl(name),charsize=2, charthick=6, color=255, align=1
    endif

   
   device,/close

   print, 'done '+ks
endfor
end
