function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro stream_pz, loc=loc, start=start, finish=finish, rslice=rslice $
               ,xtickinterval=xtickinterval,length=length,red=red $
               , smallh=smallh, xrange=xrange, ytickinterval=ytickinterval $
               , yrange=yrange, hsize=hsize $
               , log=log, plotrange0=plotrange0, ct=ct, arrcolor=arrcolor $
               , name=name, nopert=nopert, mp=mp, dvphi=dvphi
           
if not keyword_set(length) then length= 1.0
if not keyword_set(ct) then ct=5

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

if not keyword_set(red) then begin
  red = [nphi, ntheta]
endif

vx   = dblarr(ntheta, nphi)
vy   = dblarr(ntheta, nphi)
den2d = dblarr(ntheta, nphi)

vx_shifted = transpose(vx)
vy_shifted = transpose(vy)
den2d_shifted = transpose(den2d)

data_axisymmetric = dblarr(nrad, ntheta)


data0   = hdf(19, fileloc);density field t=0
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

z1 = 1d0/(smallh*tan(theta))
azi1 = phi/!dpi - 1.0

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
    
   den   = hdf(19, fileloc)
   vrad  = hdf(3, fileloc)
   vtheta= hdf(7, fileloc)
   vphi  = hdf(11, fileloc)


   vrad   *= den
   vtheta *= den
   vphi   *=den
   
   plx=planetinfo(1,n)
   ply=planetinfo(2,n)
   plz=planetinfo(3,n)
   phiplt = pltphi(plx,ply)
      
   if keyword_set(mp) then begin
      plrad = sqrt(plx^2 + ply^2 + plz^2)
      rhill  = plrad*(mp/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
   endif else rplot = rad

   temp = min(abs(rplot - rslice), r1)
    
   for k=0, nphi-1 do begin
      for j=0, ntheta-1 do begin

         
            sint = sin(theta(j))
            cost = cos(theta(j))
         
            cylind_rad = rad(r1)*sint
            
            vr   = vrad(r1,j,k)
            vt   = vtheta(r1,j,k)
            vp   = vphi(r1,j,k)

            if keyword_set(dvphi) then vp -= den(r1,j,k)*cylind_rad^(-0.5)

; vx is azimuthal velocity  in shifted co-ordinates (psi_dot = phi_dot/pi)
; and phi_dot = vphi/cylind_rad
            vx(j,k) = vp/(cylind_rad*!dpi)
            
;vy is cylindrical vert. vel. but plot in z/H, so want to divide by H
            vy(j,k) = vr*cost - sint*vt
            vy(j,k)/= smallh*cylind_rad 
            
        endfor
    endfor

   if not keyword_set(nopert) then begin
      if not keyword_set(log) then den = den - data0
      den/=data0
   endif
   
   den2d(*,*) = den(r1,*,*)
   if keyword_set(log) then den2d=alog10(den2d)

   vx = transpose(vx)
   vy = transpose(vy)
   den2d = transpose(den2d)


   if(phiplt gt !dpi) then begin
        temp2=min(abs(phi-(phiplt-!dpi)),grid2)

        vx_shifted(0:nphi-1-grid2,0:ntheta-1) = vx(grid2:nphi-1,0:ntheta-1)
        vx_shifted(nphi-grid2:nphi-1,0:ntheta-1) = vx(0:grid2-1,0:ntheta-1)

        vy_shifted(0:nphi-1-grid2,0:ntheta-1) = vy(grid2:nphi-1,0:ntheta-1)
        vy_shifted(nphi-grid2:nphi-1,0:ntheta-1) = vy(0:grid2-1,0:ntheta-1)

        den2d_shifted(0:nphi-1-grid2,0:ntheta-1) = den2d(grid2:nphi-1,0:ntheta-1)
        den2d_shifted(nphi-grid2:nphi-1,0:ntheta-1) = den2d(0:grid2-1,0:ntheta-1)

    endif
    if(phiplt lt !dpi) then begin
        temp2=min(abs(phi-(phiplt+!dpi)),grid2)

        vx_shifted(nphi-grid2:nphi-1,0:ntheta-1) = vx(0:grid2-1,0:ntheta-1)
        vx_shifted(0:nphi-1-grid2,0:ntheta-1) = vx(grid2:nphi-1,0:ntheta-1)
        
        vy_shifted(nphi-grid2:nphi-1,0:ntheta-1) = vy(0:grid2-1,0:ntheta-1)
        vy_shifted(0:nphi-1-grid2,0:ntheta-1) = vy(grid2:nphi-1,0:ntheta-1)

        den2d_shifted(nphi-grid2:nphi-1,0:ntheta-1) = den2d(0:grid2-1,0:ntheta-1)
        den2d_shifted(0:nphi-1-grid2,0:ntheta-1) = den2d(grid2:nphi-1,0:ntheta-1)
    endif
    if(phiplt eq !dpi) then begin
       vx_shifted  = vx
       vy_shifted  = vy
       den2d_shifted = den2d
    endif

   if not keyword_set(plotrange0) then begin
      
      temp = min(abs(xrange(0) - azi1), t1)
      temp = min(abs(xrange(1) - azi1), t2)
      
      plotrange=[min(den2d_shifted(t1:t2, *)),max(den2d_shifted(t1:t2, *))]
      
   endif else plotrange=plotrange0
   levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)

   time=string(planetinfo(0,n)/torb,format='(F7.2)')
   if keyword_set(mp) then begin
      name1 = textoidl('r-r_p=')+ $
              string(rslice,format='(F4.1)') + $
              textoidl('r_h')
   endif else begin
      name1 = textoidl('r=')+ $
              string(rslice,format='(F4.1)')
   endelse

   loadct, ct

   set_plot, 'ps'
   device, filename=filepath(strcompress('stream_pz'+ks+'.ps',/remove_all) $
                            ,root_dir='.',subdir=[location]) $
     , bits_per_pixel=8,xsize=18, ysize=9,/color

   contour, den2d_shifted, azi1, z1, /fill $
            ,title=time+' orbits, '+name1, levels=levels $;,xstyle=2 $ 
      ,ymargin=[4,2],xmargin=[7,10],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
      , ytickinterval=ytickinterval,xtitle=textoidl('(\phi - \phi_p)/\pi'), ytitle=textoidl('z/H')

   colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange,format='(f5.2)'

   vx_small = congrid(vx_shifted, red(0), red(1))
   vy_small = congrid(vy_shifted, red(0), red(1))
   azi1_small = congrid(azi1, red(0))
   zaxis = congrid(z1, red(1))

   temp = min(abs(azi1_small - xrange(0)), t1)
   temp = min(abs(azi1_small - xrange(1)), t2)
;   temp = min(abs(zaxis - yrange(0)), t2)
;   temp = min(abs(zaxis - yrange(1)), t1)

;   r1 += 1
;   r2 -= 1
  
   a1 = 0
   a2 = red(1) - 1
   
   t1 += 1
   t2 -= 1

   velovect2, vx_small(t1:t2,a1:a2), vy_small(t1:t2,a1:a2), azi1_small(t1:t2), zaxis(a1:a2), $ 
              color=arrcolor,/overplot, length=length,/isotropic, hsize=hsize, xrange=xrange

    if keyword_set(name) then begin
       xyouts, xrange(1), 0.1, textoidl(name),charsize=1.5, charthick=6, color=255, alignment=1.2
    endif
 
device,/close

print, 'done '+ks
endfor
end
