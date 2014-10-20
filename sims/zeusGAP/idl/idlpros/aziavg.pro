; THIS VERSION IS FOR 1D PLOT OF SIGMA OR VR OR VTHETA
; AS A FUNCTION OF RADIUS. DEFAULT IS TO AVERAGE VARIABLE OVER
; AZIMUTH AT EACH RADIUS. OPTION TO PLOT VARIABLE ALONG THE AZIMUTH
; WHICH INTERSECTS THE PLANET. 
pro aziavg, type=type, loc=loc, start=start, finish=finish, plazi=plazi, opp=opp
common consts, pi, nrad, time
common array, avgsig, data1d, rad, radius
pi=3.141592654
if keyword_set(opp) then begin
opp=double(opp*pi)
endif else opp=0.
xxtitle='r'
if keyword_set(plazi) then begin
case type of
'dens':yytitle=textoidl('\Sigma/10^{-4}')
'vrad':yytitle=textoidl('u_r/10^{-4}')
'vtheta':yytitle=textoidl('u_{\phi}')
end
endif else begin
case type of
'dens':yytitle=textoidl('<\Sigma>_{\phi}/10^{-4}')
'vrad':yytitle=textoidl('<u_r>_{\phi}/10^{-4}')
'vtheta':yytitle=textoidl('<u_{\phi}>_{\phi}')
end
endelse
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nrad=fix(dims(6))
nsec=fix(dims(7))
nout=fix(dims(5))
info=dblarr(9,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR RADIUS AND AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
if type eq 'vrad' then begin
rad(0:nrad-1)=radtmp(0:nrad-1)
endif else begin
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.
endelse
azi=dindgen(nsec)*2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
avgsig=dblarr(nrad)
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE BEGIN PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
for k=start, finish do begin
ks=string(k,format='(I03)')
openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2
plx=info(1,k)
ply=info(2,k)
if keyword_set(plazi) then begin
phi=(pltphi(plx,ply)+opp) mod (2.*pi)
array=abs(azi-phi)
grid=where(array eq min(array))
avgsig(0:nrad-1)=data(grid,0:nrad-1)
endif else begin
for m=0, nrad-1 do avgsig(m)=mean(data(*,m))
endelse
avgsig=avgsig*10.^4.
plotrange=[min(avgsig),max(avgsig)]
;plotrange=[4.,9.]
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress('avg'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
plot,rad,avgsig,title=time+' orbits',xmargin=[8,4],ymargin=[3,3] $
,xtitle=xxtitle,ytitle=yytitle $
,charsize=1.5,xrange=[min(rad),max(rad)],xtickinterval=0.4,xminor=4,yrange=plotrange
;;;;;;;;;;;;;;;;;;;;
;MARK PLANET RADIUS;
;;;;;;;;;;;;;;;;;;;;
pos=0.25*(plotrange(1)-plotrange(0))
radius=sqrt(plx*plx+ply*ply)
vline, radius
device,/close
print, 'done '+ks
endfor
; if keyword_set(getur) then begin
; array2=abs(rad-radius)
; plgrid=where(array2 eq min(array2))
; print,'velocity at planet radius is = ', avgsig(plgrid)
; endif
end
