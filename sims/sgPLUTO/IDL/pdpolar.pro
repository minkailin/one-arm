function xy_to_polar, input, x, nx, rad, azi 

  hafnx = nx/2
  nr   = hafnx
  nphi = hafnx
  
  output = dblarr(nr, nphi)
  
  xspace = x(1) - x(0)

  for j=0, hafnx-1 do begin
     for i=0, hafnx-1 do begin
        xpos = rad(i)*cos(azi(j))
        ypos = rad(i)*sin(azi(j))
        
        temp = min(abs(x-xpos), x1)
        xedge = x(x1) - xspace/2d0
        dx   = xpos - xedge
        ip  = x1 + dx/xspace
        
        temp = min(abs(x-ypos), y1)
        yedge = x(y1) - xspace/2d0
        dy   = ypos - yedge
        jp  = y1 + dy/xspace

        output(i,j) = BILINEAR(input, ip, jp)

     endfor
  endfor
  
  return, output
end


pro dpolar, loc=loc, start=start, finish=finish, log=log, nopert=nopert, xrange=xrange, yrange=yrange, $
        hold=hold, mp=mp, ct=ct, zslice=zslice, plotrange=plotrange, nonaxi=nonaxi, $
        ytickinterval=ytickinterval, xtickinterval=xtickinterval, r0=r0, azishift=azishift
if not keyword_set(ct) then ct =5   
if not keyword_set(finish) then finish=start 
if not keyword_set(zslice) then zslice= 0.5
if not keyword_set(r0) then  r0=1.0

COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi

location =strcompress(loc,/remove_all)

name = string(0)
file  = filepath(strcompress('rt'+name+'.vtk',/remove_all),root_dir='.',subdir=[location])
readvtk, file

d0 = d

rmin = 0.01*r0
rmax = x(nx-1) 

hafnx = nx/2
rad = rmin + (rmax - rmin)*dindgen(hafnx)/(hafnx-1d0)
azi = 2d0*!dpi*dindgen(hafnx)/hafnx

dataplot = dblarr(hafnx, hafnx)

if keyword_set(azishift) then begin
;planet info
loc_hist =strcompress(loc+'/id0',/remove_all)
name = string(0)
file  = filepath(strcompress('rt.hst',/remove_all),root_dir='.',subdir=[loc_hist])
nlines= file_lines(file) - 3
header= strarr(3)
info  = dblarr(11, nlines)
openr, 1, file
readf,1,header
readf, 1, info
close,1
endif

loadct,ct,/silent
for i=start, finish do begin
   name = string(i)
   file  = filepath(strcompress('rt'+name+'.vtk',/remove_all),root_dir='.',subdir=[location])
   readvtk, file

if keyword_set(azishift) then begin
   time = info(0, i)/(2.0*!dpi*r0^1.5)
   plx  = info(9, i)
   ply  = info(10,i)   
   plrad= sqrt(plx^2 + ply^2)
   phiplt = pltphi(plx,ply)  
   title  = string(time,format='(f5.2)')+textoidl('P_{k0}')
   if keyword_set(mp) then rhill = (mp/3.0)^(1.0/3.0)*plrad
endif

   if not keyword_set(nopert) then begin
      d /= d0 + 1d-16
      if not keyword_set(log) then   d -= 1.0
   endif
   if keyword_set(log) then d=alog10(d)

   data = d(*,*,nz*zslice)
  
   data_polar = xy_to_polar(data, x, nx, rad, azi)

   if keyword_set(nonaxi) then begin
      for j=0, hafnx-1 do begin
         data_polar(j,*) -= mean(data_polar(j,*))
      endfor
   endif


   if keyword_set(mp) then begin
   radplot = (rad - plrad)/rhill
   xtitle = textoidl('(r-r_0)/r_h')
   endif else begin
   radplot = rad/r0
   xtitle = textoidl('r/r_0')
   endelse

   if not keyword_set(plotrange) then begin
      temp = min(abs(radplot-xrange(0)), r1)
      temp = min(abs(radplot-xrange(1)), r2)
      plotrange0 = [min(data_polar(r1:r2,*)),max(data_polar(r1:r2,*))]
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0


   if keyword_set(azishift) then begin
     nrad = hafnx
     nphi = hafnx
     phi  = azi
     data2d= data_polar 
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
    ytitle = textoidl('(\phi-\phi_0)/\pi')
    azinew = azi/!dpi - 1.0
   endif else begin
      dataplot = data_polar
      ytitle = textoidl('\phi/2\pi')
      azinew = azi/(2d0*!dpi)
   endelse
   
   name2 = string(i, format='(I03)')
   set_plot,'ps'
   device, filename= filepath(strcompress('dpolar_'+name2+'.ps',/remove_all),root_dir='.',subdir=[location]) $
           ,/color, bits_per_pixel=8, xsize=12, ysize=14
   ;; contour, data_polar, rad/r0, azi/(2d0*!dpi), xrange=xrange, yrange=yrange, /fill,/isotropic  $
   ;;          ,xtickinterval=xtickinterval, ytickinterval=ytickinterval,levels=levels ,ymargin=[2,2],xmargin=[4,8]
   ;;  colorbar, position=[0.87, 0.06, 0.9, 0.84],/vertical,/right,range=plotrange0,format='(f5.2)'
   contour,dataplot,radplot, azinew,/fill,levels=levels, xstyle=2, $
           xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
           xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title, xminor=5
    colorbar, position=[0.851, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'
   device,/close
   
endfor
end
