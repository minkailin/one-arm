; THIS PROCEDURE PLOTS SIGMA OR VTHETA OR VRAD AS WE GO AROUND
; A CIRCLE OF RADIUS RAD. DEFAULT IS TO NOT SET RAD, AND USE PLANET
; RADIUS. WE ALSO SET PLANET AZIMUTH=0 (CENTER OF PLOT). NOTE VTHETA
; AND SIGMA ARE STORED AT RADIUS OF MIDDLE OF CELLS, AND VRAD IS
; STORED AT RADIUS OF CELL INTERFACE.
pro radavg, type=type, loc=loc, start=start, finish=finish, rad=rad
common consts, pi
pi=3.141592654
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
info=dblarr(9,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
dtheta=2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET ARRAY OF RADIUS VALUES, AND ARRAY OF AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
if type eq 'vrad' then begin
rplot=radtmp(0:nrad-1)
endif else begin 
rplot=0.5*(radtmp(0:nrad-1)+radtmp(1:nrad))
endelse
azifix=dindgen(nsec)*2.*pi/nsec
;Array to hold raw data and data for plotting.
data=dblarr(nsec,nrad)
data1d=dblarr(nsec)
;;;Option to choose fixed radius, find grid point whose radius is
;;;closest
if keyword_set(rad) then begin
rad=double(rad)
array=abs(rplot-rad)
grid=where(array eq min(array))
plpos=string(rplot(grid),format='(F4.2)')
endif
xxtitle=textoidl('\phi')
case type of
'dens':yytitle=textoidl('\Sigma/10^{-4}')
'vrad':yytitle=textoidl('u_r')
'vtheta':yytitle=textoidl('u_{\phi}')
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;SETUP DONE, BEGIN PLOTTING;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for k=start, finish do begin
ks=string(k,format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2
plx=info(1,k)
ply=info(2,k)
;;;;default is to follow planet's radius. later we can modify this to
;;;;follow radius of inner/outer edge of co-orbital region. this is
;;;;roughly plus/minus planet hill radius
if not keyword_set(rad) then begin
plrad=sqrt(plx*plx+ply*ply)
array=abs(rplot-plrad)
grid=where(array eq min(array))
plpos=string(plrad,format='(F4.2)')
endif
data1d(0:nsec-1)=data(0:nsec-1,grid)
data1d=data1d*10.^4.
;plotrange=[min(data1d),max(data1d)]
plotrange=[5.,8.]
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress('rnd'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
; set planet at zero azimuth. this involves shifting azimuthal axis
; values and plotting order needs changing. if do not need to have
; planet at theta=0 then just use the last plot option (the else
; if).
phi=pltphi(plx,ply)
azi=azifix-phi
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg1=where(azi eq min(azi))
plot,azi(seg1:nsec-1),data1d(seg1:nsec-1),title=time+' orbits, r='+plpos,xmargin=[6,4],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
oplot,azi(0:seg1-1),data1d(0:seg1-1)
endif else if (0. le phi and  phi lt pi) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg2=where(azi eq max(azi))
plot,azi(seg2+1:nsec-1),data1d(seg2+1:nsec-1),title=time+' orbits, r='+plpos,xmargin=[6,4],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
oplot,azi(0:seg2),data1d(0:seg2)
endif else begin
plot,azi,data1d,title=time+' orbits, r='+plpos,xmargin=[6,4],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
endelse
;do not plot points in order of index. may be gap so manually insert
;line segment
arrow, azi(0),data1d(0),azi(nsec-1),data1d(nsec-1),hsize=0,/data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPTION TO MARK PLANET AZIMUTH;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; if keyword_set(markazi) then begin
; ;pos=0.90*(plotrange(1)-plotrange(0))
; phi=pltphi(plx,ply)
; ;xyouts,phi*(1.01),pos, 'planet', charsize=1.5, charthick=1.5
vline, 0.
; endif
device,/close
print, 'done '+ks
endfor
end

