; THIS VERSION IS FOR 1D PLOT OF SIGMA OR VR OR VTHETA
; AS A FUNCTION OF RADIUS. DEFAULT IS TO AVERAGE VARIABLE OVER
; AZIMUTH AT EACH RADIUS. OPTION TO PLOT VARIABLE ALONG THE AZIMUTH
; WHICH INTERSECTS THE PLANET. 
pro radavg, type=type, loc=loc, start=start, finish=finish $
            , plazi=plazi, opp=opp, noplot=noplot, toomre=toomre, h=h, xrange=xrange, out=out $
            , width = width, yrange=yrange
common consts, pi, nrad, time
common array, avgsig, rad, radius, plx, ply
pi=3.141592654
mp=3d-4
f0=(mp/3.0)^(1.0/3.0)
if not keyword_set(finish) then finish=start
case type of
'dens':yytitle=textoidl('\Sigma\times10^{5}')
'vrad':yytitle=textoidl('u_r/10^{-4}')
'vtheta':yytitle=textoidl('u_{\phi}')
end
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
data_1d = dblarr(nsec)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
avgsig=dblarr(nrad)
maxmode = dblarr(2,finish-start+1)
count = intarr(nsec)
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
    plrad = sqrt(plx*plx + ply*ply)
    rhill = f0*plrad
    phi=pltphi(plx,ply) 
    array=abs(azi-phi)
    temp=min(array, grid)
    temp=min(abs(rad-(plrad+width(0)*rhill)), r0)
    temp=min(abs(rad-(plrad+width(1)*rhill)), r1)
    for j=0, nsec-1 do data_1d(j) = mean(data(j,r0:r1))
    average = mean(data_1d)
;     dsigma = deriv(azi, data_1d)
;     ddsigma = deriv(azi, dsigma)
;     count(*)  = 0
;     for j=1, nsec-1 do begin
;         if(dsigma(j)*dsigma(j-1) lt 0d0) then begin
;             if(ddsigma(j)/data_1d(j) lt -1d-1) then begin
;                 if(data_1d(j) gt average) then count(j) = j
;             endif
;         endif
;     endfor
    
;     maxmode(0,k-start) = info(7,k)/p0
;     maxmode(1,k-start) = n_elements(where(count ne 0))

    fourier = fft(data_1d,/double)
    fourier = abs(fourier)
    temp = where(fourier(1:nsec-1) eq max(fourier(1:nsec-1)))
    maxmode(0,k-start) = info(7,k)/p0
    maxmode(1,k-start) = temp(0) + 1

    time=string(info(7,k)/p0,format='(F7.2)')
    if not keyword_set(noplot) then begin
        set_plot, 'ps'
        device, filename=filepath(strcompress('ravg'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
        plot,azi/pi, data_1d*1d5,title=time+' orbits',xmargin=[8,4],ymargin=[3,3] $
          ,xtitle=textoidl('\phi/\pi'),ytitle=yytitle $
          ,charsize=1.5,xrange=xrange,xminor=4,yrange=yrange, thick=4
;        for j =0, nsec-1 do begin
;            if(count(j) ne 0) then oplot, [azi(j),azi(j)]/pi, [data_1d(j),data_1d(j)]*1d5 $
;              , psym=2
;        endfor
;;;;;;;;;;;;;;;;;;;;
;MARK PLANET RADIUS;
;;;;;;;;;;;;;;;;;;;;
vline, phi/pi
device,/close
endif

if keyword_set(toomre) then begin
    csomega=dblarr(nsec)
    csomega = h*mean(rad(r0:r1))^(-2.0)
    q=csomega/(!dpi*data_1d)
    set_plot, 'ps'
    device, filename=filepath(strcompress('rQ'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    plot,azi/pi,q,title=time+' orbits',xmargin=[8,4],ymargin=[3,3] $
      ,ytitle='Q',xtitle=textoidl('\phi/\pi') $
      ,charsize=1.5, xrange=xrange, thick=4, yrange=[min(q),max(q)]
    device,/close
endif
print, 'done '+ks
endfor
print, maxmode

 set_plot, 'ps'
        device, filename=filepath(strcompress('maxmode_.ps',/remove_all),root_dir='.',subdir=[location])$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
        plot,maxmode(0,*), maxmode(1,*),xmargin=[8,4],ymargin=[3,3] $
          ,xtitle='t/orbits',ytitle='number of maxima',psym=2 $
          ,charsize=1.5,xrange=xrange,xminor=4,yrange=yrange, thick=4
device,/close


end
