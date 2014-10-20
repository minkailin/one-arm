; THIS VERSION IS FOR 1D PLOT OF SIGMA OR VR OR VTHETA
; AS A FUNCTION OF RADIUS. DEFAULT IS TO AVERAGE VARIABLE OVER
; AZIMUTH AT EACH RADIUS. OPTION TO PLOT VARIABLE ALONG THE AZIMUTH
; WHICH INTERSECTS THE PLANET. 
pro aziavg, type=type, loc=loc, start=start, finish=finish $
, plazi=plazi, opp=opp, noplot=noplot, toomre=toomre, h=h, xrange=xrange, out=out
common consts, pi, nrad, time
common array, avgsig, rad, radius, plx, ply
pi=3.141592654
if not keyword_set(finish) then finish=start
if keyword_set(opp) then begin
opp=double(opp*pi)
endif else opp=0.
xxtitle='r'
if keyword_set(plazi) then begin
case type of
'dens':yytitle=textoidl('\Sigma\times10^{3}')
'vrad':yytitle=textoidl('u_r/10^{-4}')
'vtheta':yytitle=textoidl('u_{\phi}')
end
endif else begin
case type of
'dens':yytitle=textoidl('<\Sigma>_{\phi}\times10^{3}')
'vrad':yytitle=textoidl('<u_r>_{\phi}/10^{-4}')
'vtheta':yytitle=textoidl('<u_{\phi}>_{\phi}')
end
endelse
if not keyword_set(h) then h=0.05
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nrad=fix(dims(6))
nsec=fix(dims(7))
nout=fix(dims(5))
if keyword_set(out) then nout = out
info=dblarr(11,nout+1)
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
avgsig=avgsig
;stop
plotrange=[min(avgsig),max(avgsig)]
time=string(k*dt/p0,format='(F7.2)')
if not keyword_set(noplot) then begin
set_plot, 'ps'
device, filename=filepath(strcompress('avg'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot,rad,avgsig,title=time+' orbits',xmargin=[8,4],ymargin=[3,3] $
,xtitle=xxtitle,ytitle=yytitle $
,charsize=1.5,xrange=xrange,xminor=4,yrange=yrange, thick=4

openw,1,filepath(strcompress('avg'+type+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, nrad-1 do printf,1,rad(i), avgsig(i), format='(2(e12.6,2x))'
close,1


;for i=0, nrad-1 do print, rad(i), avgsig(i)
;;;;;;;;;;;;;;;;;;;;
;MARK PLANET RADIUS;
;;;;;;;;;;;;;;;;;;;;
pos=0.25*(plotrange(1)-plotrange(0))
radius=sqrt(plx*plx+ply*ply)
vline, radius
device,/close
endif
if keyword_set(toomre) then begin
omega = rad^(-3./2.)
cs=h*rad^(-1.0/2.0) 
q=omega*cs/(!dpi*avgsig)
set_plot, 'ps'
device, filename=filepath(strcompress('Q'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot,rad,q,title=time+' orbits',xmargin=[8,4],ymargin=[3,3] $
,xtitle=xxtitle,ytitle='Q' $
,charsize=1.5, xrange=xrange, thick=4
device,/close
endif
print, 'done '+ks
endfor
end
