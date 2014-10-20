pro polarxy, type=type, loc=loc, start=start, finish=finish, r0=r0 $
             , ct=ct,plotrange0=plotrange0, log=log, xrange=xrange, nopert=nopert, basic=basic $
             ,xtickinterval=xtickinterval,mp=mp,yrange=yrange, scale=scale, nonaxi=nonaxi, cart=cart, mdisk=mdisk, mmode=mmode
if not keyword_set(ct) then ct=5
if not keyword_set(finish) then finish=start
common consts, pi, nrad, time
!p.font = 0
pi=3.141592654
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
    nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[location]))
    info=dblarr(11,nlines)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
if not keyword_set(r0) then begin 
a0=info(1,0)
endif else a0 = r0 
;dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dblarr(nsec)
azi1=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

if not keyword_set(xrange) then xrange=[min(rad), max(rad)]


if keyword_set(cart) then begin
   dazi = azi(1) - azi(0)
   dlogr = alog(radtmp(nrad)/radtmp(0))/nrad

   dnx1 = 2*nrad
   xaxis = -rad(nrad-1) + 2*rad(nrad-1)*dindgen(dnx1)/(dnx1-1.0)
   yaxis = xaxis
   dataxy = dblarr(dnx1, dnx1)
endif

;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
data0=dblarr(nsec,nrad)
dataplot = dblarr(nrad,nsec)
pattern  = dblarr(2,finish-start+1)
loadct,ct, bottom=0

;if not keyword_set(nopert) then begin
    if not keyword_set(basic) then begin 
        basic = 0
    endif else basic = 0

    openr,2,filepath(strcompress('gas'+type+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
    readu,2,data0
    close,2
    
    dataT = transpose(data0)
    
    dataplot0 = dataT
;endif

data0_fft = fft(data0, -1, dimension=1,/double) 


for k=start, finish do begin
      
if (k lt 1000) then begin
   ks=string(k,format='(I03)')  
endif else ks=string(k,format='(I04)')

openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2

if keyword_set(mmode) then begin ;replace density field by a particular fourier component
data_fft = fft(data, -1, dimension=1,/double) 
for i=0,nsec-1 do begin
for j=0, nrad-1 do begin
data(i,j) = min([mmode+1,2.0])*real_part( data_fft(mmode,j)*( cos(mmode*azi(i)) + dcomplex(0,1d0)*sin(mmode*azi(i)) ) )
data(i,j)/= abs(data0_fft(0,j)) 
endfor
endfor

temp = max(data,grid) 
result = array_indices(data, grid)
rmax   = rad(result(1))
phimax = azi(result(0)) 
pattern(0,k-start) = info(7,k)
pattern(1,k-start) = phimax  
endif 




if keyword_set(mdisk) then begin
mass = 0.0
dphi = 2d0*!dpi/nsec
dlogr = alog(radtmp(nrad)/radtmp(0))/nrad
for j=0, nsec-1 do begin
for i=0, nrad-1 do begin
ds = rad(i)^2*dlogr*dphi

mass+= data(j,i)*ds

endfor
endfor
print, 'mdisk=', mass
endif

if keyword_set(scale) then data *= scale

if keyword_set(nonaxi) then begin
for i=0, nrad-1 do begin
data0(*,i) = mean(data(*,i))
endfor
dataplot0 = transpose(data0)
endif


dataT = transpose(data)
 

plx=info(1,k)
ply=info(2,k)
if not keyword_set(r0) then begin
 plrad=sqrt(plx*plx+ply*ply)
endif else plrad = r0 

rplot =  rad/plrad

dataplot = dataT

if not keyword_set(nopert) then begin
    dataplot /= dataplot0
if not keyword_set(log) then  dataplot -= 1.0
endif

if keyword_set(log) then dataplot=alog10(dataplot)

if not keyword_set(plotrange0) then begin 
    temp = min(abs(rplot - xrange(0)),r1)
    temp = min(abs(rplot - xrange(1)),r2)
    plotrange=[min(dataplot(r1:r2,*)),max(dataplot(r1:r2,*))]
endif else plotrange=plotrange0
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
levels2= (plotrange(1)/3d0)*(dindgen(2)/1.)
time=string(info(7,k)/p0,format='(F7.2)')
title = time+textoidl('P_0, ')+textoidl('r_p='+string(plrad,format='(f5.2)'))
xtitle = textoidl('r/r_0');textoidl('(r-r_p)/r_h')

if not keyword_set(cart) then begin
set_plot, 'ps'
device, filename=filepath(strcompress('polarxy_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
        ,/color, bits_per_pixel=8,xsize=12, ysize=14
contour,dataplot,rplot,azi/(2.0*pi),/fill,levels=levels,title=title, $
        xtitle=xtitle, ytitle=textoidl('\phi/2\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
        xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, xminor=5, xstyle=1
;contour,dataplot,rplot, azi/pi - 1.0 ,levels=levels2, /overplot
;xyouts, 0.55, 0.02, 'FARGO 2D', charsize=3.5,charthick=8, color=255
colorbar, position=[0.865, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[0.0,0.0],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*1.0
;oplot,[0.,0.],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.5
;xyouts,plrad,phi/pi,'X',charsize=1.5
device,/close
endif else begin

;   for jj=0, dnx1-1 do begin
;      y = yaxis(jj)
;      for ii=0, dnx1-1 do begin
;         x = xaxis(ii)
;         
;         r_t = sqrt(x^2 + y^2)
;         azi_t = pltphi(x, y)
;         
;         if( (r_t ge xrange(0)) and (r_t le rad(nrad-1)) )then begin
;            
;            temp = min(abs(rad - r_t),   x0)
;            temp = min(abs(azi - azi_t),  y0)
;            
;            dr = rad(x0)*dlogr
;            ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
;            jp = y0 + (azi_t - azi(y0))/(dazi/2.0) + 0.0
;            
;            dataxy(ii,jj) = bilinear(dataplot, ip, jp)
;         endif else begin
;            dataxy(ii,jj) = 1d10
;         endelse
;      endfor
;   endfor
  
   dataxy = transpose(dataplot)
   temp = min(abs(rad - xrange(0)),x1)
   for i=0, x1 do dataxy(*,i) = 10d0*abs(plotrange(1))
   
   set_plot, 'ps'
   device, filename=filepath(strcompress('polarxy2_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
           ,/color, bits_per_pixel=8,xsize=14, ysize=12
;   contour, dataxy, xaxis, yaxis, /isotropic,/fill,levels=levels,title=title $
;            ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, $
;            xrange=[-1,1]*xrange(1), yrange=[-1,1]*xrange(1), xstyle=1, ystyle=1

     polar_contour, dataxy, azi, rad, /dither, /isotropic,/fill,levels=levels,title=title $
    ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, $
     xrange=[-1,1]*xrange(1), yrange=[-1,1]*xrange(1), xstyle=1, ystyle=1

   colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange,format='(f5.2)'
;   for i=1, 4 do begin
;   rc = i 
;   omp = rc^(-1.5)
if keyword_set(mmode) then begin
   xp  = rmax*cos(phimax)/plrad
   yp  = rmax*sin(phimax)/plrad
   oplot, [xp,xp], [yp,yp], thick=4, psym=7
endif
;   endfor


   device,/close
                                ;stop
endelse

print, 'done '+ks



;overdensity
;rin=5.5
;rout=6.5
;temp = min(abs(rad-rin),grid1)
;temp = min(abs(rad-rout),grid2)
;subgrid = dataplot(grid1:grid2,*)
;overdense = where(subgrid gt 0d0)

;print, 'average overdensity', max(subgrid)
endfor
if ( keyword_set(mmode) and (finish - start ge 2) ) then begin
pattern_speed = deriv(pattern(0,*),pattern(1,*))
print,  mean(pattern_speed)^(-2./3.)

set_plot, 'ps'
  device, filename=filepath(strcompress('polarxy_pattern.ps',/remove_all),root_dir='.',subdir=[location]) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, pattern(0,*)/p0, pattern(1,*)/(2d0*!dpi),xmargin=[8,2],ymargin=[3.5,0.5], ytitle=textoidl('\phi_{max}/2\pi')  $
        ,charsize=1.5, thick=4, xtitle=textoidl('t/P_0') 
device,/close

endif
end
