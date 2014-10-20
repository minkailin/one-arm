; THIS PROCEDURE COMPARES THE AVERAGE RADIAL FLUID VELOCITY IN THE
; NEIGHBOURHOOD OF PLANET, AND THE PLANET MIGRATION RATE. WE MAY WANT
; TO CHANGE TO DENSITY-WEIGHTED RADIAL VELOCITY IN A LATER
; VERSION. THIS VERSION IS VERY SPECIFIC FOR THE EXPERIMENT IN MIND
; (ZERO VISC AND ALL ITS OUTPUT SETUP), MANY PARTS HARD-WIRED. NEED TO GENERALISE LATER.
pro uradot, loc=loc, segment=segment, xs=xs, ufit=ufit
common consts, pi
pi=3.141592654
;;;;;;;;;;;;;;;;;;;;;;;;
;PART I: INFO ABOUT RUN;
;;;;;;;;;;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
info=dblarr(9,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PART II: GETTING AVERAGE RADIAL VELOCITY IN VINCINITY OF PLANET;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Get semi-major axis. There are 10,000 small time steps but
;only 250 outputs for radial velocity. later we could use planet
;radius rather than semimajor axis to get rhill.
orbit=READ_ASCII(filepath('orbit0.dat',root_dir='.',subdir=[location]))
orbit=orbit.(0)
lines=n_elements(orbit(0,*))
step=fix(lines/nout)
semia=dblarr(nout+1)
for i=0, nout do semia(0:i)=orbit(2,i*step)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET ARRAY OF RADIUS VALUES, AND ARRAY OF AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rplot=radtmp(0:nrad-1)
;endif else begin 
;rplot=0.5*(radtmp(0:nrad-1)+radtmp(1:nrad))
;endelse
azifix=dindgen(nsec)*2.*pi/nsec
;Array to hold raw data
data=dblarr(nsec,nrad)
;setup array to hold average radial velocity in vincinity of
;planet. and factor such that rhill=factor*semimajor axis. standard to
;include extra factor of 2.5, as influence region is few times rhill.
avgurad=dblarr(2,nout+1)
;fill the time column
avgurad(0,0:nout)=info(7,0:nout)
factor=double((0.00028/3.)^(1./3.))
if keyword_set(xs) then begin
factor=factor*xs
endif else factor=factor*2.5
for k=0, 250 do begin
;Read raw data.
openr,2,filepath(strcompress('gas'+'vrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2
;physical size of half width. find planet radius and azimuth, and find
;range of grid points that correspond to sector
width=factor*semia(k)
plx=info(1,k)
ply=info(2,k)
plrad=sqrt(plx*plx+ply*ply)
;width=factor*plrad
arr1=abs(rplot-(plrad-width))
arr2=abs(rplot-(plrad+width))
rmns=where(arr1 eq min(arr1))
rpls=where(arr2 eq min(arr2))
phi=pltphi(plx,ply)
arr3=abs(azifix-(phi-width/plrad))
arr4=abs(azifix-(phi+width/plrad))
phimns=where(arr3 eq min(arr3))
phipls=where(arr4 eq min(arr4))
;get the average radial velocity in this region. 
avgurad(1,k)=mean(data(phimns:phipls,rmns:rpls))
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PART III: FIT CURVE TO a(t) AND GET CURVE OF ADOT(T);
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
adot=dblarr(2,lines)
adot(0,0:lines-1)=orbit(0,0:lines-1)
up=intarr(n_elements(segment)-1)
segment=double(segment*p0)
lower=0
for j=1, n_elements(segment)-1 do begin
;lower=abs(orbit(0,0:lines-1)-segment(j-1))
;lower=where(lower eq min(lower))
upper=abs(orbit(0,0:lines-1)-segment(j))
upper=where(upper eq min(upper))
xx=orbit(0,lower:upper)
yy=orbit(2,lower:upper)
quadfit=poly_fit(xx,yy,2,/double)
adot(1,lower:upper-1)=quadfit(1)+2.*quadfit(2)*adot(0,lower:upper-1)
lower=upper
up(j-1)=upper-1
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PART III-b: OPTION TO DO LINEAR FIT FOR URAD; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(ufit) then begin
lower=0
up2=intarr(n_elements(segment)-1)
for l=1, n_elements(segment)-1 do begin
upper=abs(avgurad(0,0:nout)-segment(l))
upper=where(upper eq min(upper))
xx=avgurad(0,lower:upper)
yy=avgurad(1,lower:upper)
linear=linfit(xx,yy,/double)
avgurad(1,lower:upper-1)=linear(0)+linear(1)*avgurad(0,lower:upper-1)
lower=upper
up2(l-1)=upper-1
endfor
endif
;;;;;;;;;;;;;;;;;;;
;PART IV: PLOTTING;
;;;;;;;;;;;;;;;;;;;
avgurad(0,*)=avgurad(0,*)/p0
avgurad(1,*)=avgurad(1,*)*100.
adot(0,*)=adot(0,*)/p0
adot(1,*)=adot(1,*)*100.
set_plot, 'ps'
device, filename='uradot.ps',xsize=8,ysize=6,xoffset=0,yoffset=0$
,/inches,/color,bits_per_pixel=8
loadct, 40
plot, avgurad(0,*), avgurad(1,*),xrange=[0,175],xtickinterval=25,xminor=5,xmargin=[8,2]$
,ymargin=[3,3],xtitle='Orbits',ytitle='Velocity/'+textoidl('10^{-2}'),/nodata,yrange=[-0.1,0.]
if keyword_set(ufit) then begin
oplot, avgurad(0,0:up2(0)), avgurad(1,0:up2(0))
for k=1, n_elements(up2)-1 do oplot, avgurad(0,up2(k-1)+1:up2(k)), avgurad(1,up2(k-1)+1:up2(k))
endif else oplot, avgurad(0,*), avgurad(1,*)

oplot, adot(0,0:up(0)),adot(1,0:up(0)),color=64,thick=4
for m=1, n_elements(up)-1 do begin
oplot, adot(0,up(m-1)+1:up(m)), adot(1,up(m-1)+1:up(m)),color=64,thick=4
vline, adot(0,up(m-1)+1),col=64
endfor
device,/close
end

