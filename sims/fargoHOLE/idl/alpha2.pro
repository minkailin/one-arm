pro alpha2, loc=loc, start=start, finish=finish, out=out, xrange=xrange, mp=mp $
            ,yrange=yrange, legend=legend, label=label
common consts, pi

pi=!dpi
if not keyword_set(mp) then mp=3d-4
f0 = (mp/3.0)^(1.0/3.0)
h = 0.05
if not keyword_set(finish) then finish = start
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
if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info = dblarr(11,out+1)
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
rmed = (radtmp(0:nrad-1) + radtmp(1:nrad))/2.0
azi = 2d0*pi*dindgen(nsec)/nsec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vrad   = dblarr(nsec,nrad)
vtheta = dblarr(nsec,nrad)
vtheta_e = dblarr(nsec,nrad)

sigma0 =dblarr(nsec,nrad)
sigma = dblarr(nsec,nrad)
dsigma = dblarr(nsec, nrad)
sigma_e = dblarr(nsec, nrad)
rflux_2d = dblarr(nsec, nrad)
rflux_1d = dblarr(nrad)

omega = dblarr(nsec, nrad)
domega = dblarr(nsec, nrad)
shear_1d = dblarr(nrad)
shear_2d = dblarr(nsec, nrad)

vrad_1d=dblarr(nrad)

alpha_2d = dblarr(nsec, nrad)
alpha_1d = dblarr(nrad)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PREP FOR G_PHI CALCULATION ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
r0   = radtmp(0)
rmax = radtmp(nrad)  
B = 0.3*h

nr = nrad
rmed_small = congrid(rmed, nr)
uaxis = dblarr(3*nr)

dtheta = 2d0*pi/nsec
du = alog(rmax/r0)/nr

uaxis(0:nr-1) = -alog(rmax/r0) + du*dindgen(nr)
uaxis(nr:2*nr-1) = alog(rmed_small/rmed_small(0))
uaxis(2*nr: 3*nr - 1) = du/2d0 + alog(rmed_small(nr-1)/rmed_small(0)) + du*dindgen(nr)

kernel_phi = dblarr(nsec, 3*nr)
sarray_phi = dblarr(nsec, 3*nr)

gtorque_2d = dblarr(nr,nsec)
gtorque_1d = dblarr(nr)
gflux_1d = dblarr(nr)

for i = 0, nsec-1 do begin
    for j =  0, 3*nr-1 do begin        
        theta = azi(i)
        u     = uaxis(j)
        kernel_phi(i,j) = sin(theta)
        kernel_phi(i,j)/= (2d0*(cosh(u) - cos(theta)) + B*B*exp(u))^(3d0/2d0)
    endfor
endfor
kernel_phifft = fft(kernel_phi, -1,/double)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE BEGIN PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma0
close,2
for k=start, finish do begin
    ks=string(k,format='(I03)')

    openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vrad
    close,2
    vrad_old = vrad
    
    openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vtheta
    close,2
    
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2
     
                                ;reynolds flux
    vtheta_e(0:nsec-2,1:nrad-1) = (vtheta(0:nsec-2,1:nrad-1) + vtheta(1:nsec-1,1:nrad-1) + vtheta(0:nsec-2, 0:nrad-2) $ 
                                   + vtheta(1:nsec-1,0:nrad-2))/4.0
    vtheta_e(nsec-1,1:nrad-1) =  (vtheta(nsec-1,1:nrad-1) + vtheta(0,1:nrad-1) + vtheta(nsec-1, 0:nrad-2) + vtheta(0,0:nrad-2))/4.0 
    sigma_e(0:nsec-1, 1:nrad-1) = 0.5*(sigma(0:nsec-1, 0:nrad-2) + sigma(0:nsec-1, 1: nrad-1))
   
     for i = 1, nrad-1 do begin
         mean_vtheta = mean(vtheta_e(0:nsec-1,i))
         mean_vrad = mean(vrad(0:nsec-1,i))
         vtheta_e(0:nsec-1,i) -= mean_vtheta
         vrad(0:nsec-1,i)     -= mean_vrad
     endfor    
     rflux_2d(0:nsec-1, 1:nrad-1) = sigma_e(0:nsec-1, 1:nrad-1)*vtheta_e(0:nsec-1,1:nrad-1)*vrad(0:nsec-1,1:nrad-1) 

     for i =1, nrad-1 do  begin
         rflux_1d(i) = mean(rflux_2d(*,i))*rmed(i)
     endfor

                                ;gravity flux
                                ;azimuthal gravity acceleration
     sigma_small = congrid(sigma, nsec, nr)
     for j = nr, 2*nr-1 do begin
         sarray_phi(*,j) = sigma_small(*, j-nr)*exp(3d0*uaxis(j)/2d0)
     endfor
     
     sfft_phi = fft(sarray_phi, -1, /double)
     gphifft = sfft_phi*kernel_phifft*dtheta*du*3*nr*nsec
     gphi = fft(gphifft, 1, /double)
     
     for i=0, nsec -1 do begin
         gphi(i,*) *= -exp(-3d0*uaxis/2d0)
     endfor
    gtorque_2d = gphi(*,nr:2*nr-1)*sigma_small
    for j=0, nr-1 do gtorque_1d(j) = mean(gtorque_2d(*,j))*rmed_small(j)^2d0

    for i=0, nr-3 do begin       
        gflux_1d(i) = int_tabulated(rmed_small(i:nr-1), gtorque_1d(i:nr-1))/rmed_small(i)
    endfor


        
    set_plot, 'ps'
    time=string(info(7,k)/p0,format='(F7.2)')  
    ytitle = textoidl('Ang. mom. flux')

    device, filename=filepath(strcompress('alpha2_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    plot, rmed_small, gflux_1d*1d10,title=time+' orbits',xmargin=[6,2],ymargin=[3,2] $ 
      ,xtitle='r',charsize=1.5,xrange=xrange, xminor=4, thick=4, ytitle=ytitle, yrange=yrange
;    oplot, r0*exp(uaxis(nr:2*nr-1)), gflux_1d, thick=4, linestyle=1
    if keyword_set(legend) then begin
        x0=legend(0)
        x1=legend(1)
        y0=legend(2)
        dy=legend(3)
        for j=0, 2 do begin
            oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
            xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
        endfor
    endif
    device,/close



;     set_plot, 'ps'
;     gtorque_1d /= rmed_small
    
;     norm_reytorq = ceil(mean(alog10(abs(rtorque_1d(1:nrad-1)))))
;     norm_grvtorq = ceil(mean(alog10(abs(gtorque_1d))))

;     rtorque_1d /= 10d0^(norm_reytorq)
;     gtorque_1d /= 10d0^(norm_grvtorq)

;     ytitle1 = 'T_R/'+strcompress('10^{'+string(norm_reytorq)+'}',/remove_all)+', solid'
;     ytitle1 = textoidl(ytitle1)
    
;     ytitle2 = 'T_G/'+strcompress('10^{'+string(norm_grvtorq)+'}',/remove_all)+', dashed'
;     ytitle2 = textoidl(ytitle2)

;     device, filename=filepath(strcompress('alpha2_tq_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;       ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;     plot, rmed, rtorque_1d, title=time+' orbits',xmargin=[7,7],ymargin=[3,2] $
;       ,xtitle='r',charsize=1.5,xrange=xrange, xminor=4, thick=4, /nodata, ystyle=4 $
;       , xtickinterval=1
;     axis,yaxis=0,/save,ytitle=ytitle1,charsize=1.5, yrange=[min(rtorque_1d),max(rtorque_1d)]
;     oplot, rmed, rtorque_1d , thick=4, linestyle=2
;     axis,yaxis=1,/save,ytitle=ytitle2,charsize=1.5, yrange=[min(gtorque_1d),max(gtorque_1d)]
;     oplot, rmed_small, gtorque_1d, thick=4, linestyle=1
;     if keyword_set(legend) then begin
;         x0=legend(0)
;         x1=legend(1)
;         y0=legend(2)
;         dy=legend(3)
;         for j=0, 2 do begin
;             oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;             xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
;         endfor
;     endif
;     device,/close

    print, 'done '+ks
    
endfor

end






























































;    dataT = transpose(ratio)
;     dataplot = dataT

;     plx=info(1,k)
;     ply=info(2,k)
;     plrad=sqrt(plx*plx+ply*ply)
;     phi=pltphi(plx,ply)
;     plx=info(1,k)
   
;     temp=min(abs(azi-phi),grid)

;     if(phi gt !dpi) then begin
;         temp2=min(abs(azi-(phi-!dpi)),grid2)
;         dataplot(*,0:nsec-1-grid2) = dataT(*,grid2:nsec-1)
;         dataplot(*,nsec-grid2:nsec-1) = dataT(*, 0:grid2-1) 
;     endif
    
;     if(phi lt !dpi) then begin
;         temp2=min(abs(azi-(phi+!dpi)),grid2)
;         dataplot(*,nsec-grid2:nsec-1) = dataT(*,0:grid2-1)
;         dataplot(*,0:nsec-1-grid2) = dataT(*, grid2:nsec-1) 
;     endif
    
;     if(phi eq !dpi) then dataplot = dataT
        
;     if not keyword_set(plotrange0) then begin 
;         temp = min(abs(rmed_small - xrange(0)), grid1)
;         temp = min(abs(rmed_small - xrange(1)), grid2)


;         plotrange=[min(dataplot(grid1:grid2, *)),max(dataplot(grid1:grid2, *))]
;     endif else plotrange=plotrange0
;     levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
;     time=string(info(7,k)/p0,format='(F7.2)')

    
;     loadct,5
;     set_plot, 'ps'
;     device, filename=filepath(strcompress('alpha2_ratio_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;       ,/color, bits_per_pixel=8,xsize=12, ysize=14
;     contour,dataplot,rmed_small,azi/!dpi - 1.0 ,/fill,levels=levels,title=time+' orbits', $
;       xtitle=textoidl('r'), ytitle=textoidl('(\phi-\phi_p)/\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
;       xtickinterval=1,xrange=xrange, yrange=yrange
;     colorbar, position=[0.85, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
;     oplot,[plrad,plrad],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.5

; temperature = dblarr(nsec,nrad)
; temperature_e = dblarr(nsec,nrad)

; openr,2,filepath(strcompress('gasTemperature'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location])
; readu,2,temperature
; close,2

; temperature_e(0:nsec-2,1:nrad-1) = (temperature(0:nsec-2,1:nrad-1)+ temperature(0:nsec-2,0:nrad-2))/2.0
; temperature_e(nsec-1,1:nrad-1) = (temperature(nsec-1,1:nrad-1)+
; temperature(0,0:nrad-2))/2.0



;     alpha_time(0,k-start) = info(7,k)/p0
;     alpha_time(1,k-start) = mean(alpha) ; over entire disc

;     if not keyword_set(width) then begin
;         rpmxs = 0
;         rppxs = nrad-1
;     endif else begin
;         temp = min(abs(rmed - width(0)),rpmxs)
;         temp = min(abs(rmed - width(1)),rppxs)
;     endelse
;     alpha_time(2,k-start) = mean(alpha(*,rpmxs:rppxs)) ; over region near planet

;     ;alpha from 1d alpha-disc model
;     if keyword_set(visc1d) then begin
;         for i = 0, nrad-1 do begin
;             omega(0:nsec-1,i) = vtheta(0:nsec-1,i)/rmed(i) 
;             ;sigma_1d(i) = mean(sigma(*,i))
;             ;vrad_1d(i) = mean(vrad_old(*,i))
;         endfor
;         for i = 1, nrad-1 do begin
;             omega_e(0:nsec-1,i) = (omega(0:nsec-1,i) + omega(0:nsec-1,i-1))/2.0   
;             ;sigma_1d_e(i) = (sigma_1d(i) + sigma_1d(i-1))/2.0 
;         endfor

;         ;hprime = deriv(radtmp(1:nrad-1), omega_e(1:nrad-1)*radtmp(1:nrad-1)*radtmp(1:nrad-1))
;         ;g = radtmp(1:nrad-1)*sigma_1d_e(1:nrad-1)*vrad_1d(1:nrad-1)*hprime

;         ;omegaprime = deriv(rmed, omega)
;         ;alpha_visc1d(rpmxs) = (alpha_1d(rpmxs)/1d3)*h*h*sigma_1d(rpmxs)*rmed(rpmxs)^(3.5)*omegaprime(rpmxs)
;         ;for i= rpmxs+1, nrad-1 do begin
;         ;    alpha_visc1d(i) = g(i-1)*(rmed(i) - rmed(i-1)) + alpha_visc1d(i-1)
;         ;endfor

;         ;f = h*h*rmed^(3.5)*sigma_1d*omegaprime
;         ;alpha_visc1d /= f
        
;         ;print, alpha_visc1d

;         for i=0, nsec-1 do begin
;             alpha_visc2d(i, 1:nrad-1) = omega_e(i,1:nrad-1)*vrad_old(i, 1:nrad-1)/deriv(radtmp(1:nrad-1), omega_e(i, 1:nrad-1))
;         endfor
;         for i=0, nrad-1 do alpha_visc1d(i) = mean(alpha_visc2d(*,i))
        
;         ;alpha_visc1d(1:nrad-1)/=h*h*sqrt(radtmp(1:nrad-1))

;         ;print, alpha_visc1d

;     endif
