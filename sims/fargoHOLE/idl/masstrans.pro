pro masstrans, loc=loc, start=start, finish=finish, out=out, xrange=xrange, width=width $
               ,yrange=yrange

pi=3.141592654
G = 1.0
mstar = 1.0
nu=1d-5
mp = 3d-4
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
rad = dblarr(nrad+1)
rmed = dblarr(nrad)
azi    = dblarr(nsec)

openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,rad
close,1
rmed(0:nrad-1) = (rad(0:nrad-1) + rad(1:nrad))/2.0 

dlogr = alog(rad(nrad)/rad(0))/double(nrad)
dphi = 2.0*!dpi/double(nsec)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sigma =  dblarr(nsec,nrad)
vrad   = dblarr(nsec,nrad)
vrad_c = dblarr(nsec, nrad)

vtheta = dblarr(nsec, nrad)
alpha = dblarr(nsec, nrad)

sigma_1d = dblarr(nrad)
vrad_1d  = dblarr(nrad)
vtheta_1d = dblarr(nrad)
alpha_1d = dblarr(nrad)

; vtheta = dblarr(nsec,nrad)

mdot = dblarr(nsec, nrad)
mdot_1d  = dblarr(nrad)

hdot = dblarr(nsec, nrad)
hdot_1d = dblarr(nrad)

mdot_visc = dblarr(nrad)
mdot_alpha = dblarr(nrad)

massloss = dblarr(6,finish-start+1)
angmomflux = dblarr(3,finish-start+1)

openr,2,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vtheta
close,2 

vtheta0=vtheta

for l=start, finish do begin
    ks=string(l,format='(I03)')
    openr,2,filepath(strcompress('gasdens'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2 

    openr,2,filepath(strcompress('gasvrad'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vrad
    close,2 

    openr,2,filepath(strcompress('gasvtheta'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vtheta
    close,2 

for i = 0, nsec-1 do begin
    for j = 1, nrad-2 do begin
        vrad_c(i,j) = 0.5*(vrad(i,j) + vrad(i,j+1))
    endfor
endfor

for j=0, nrad-1 do begin
    sigma_1d(j)= mean(sigma(*,j))
    vrad_1d(j) = mean(vrad(*,j))
    vtheta_1d(j) = mean(vtheta(*,j))
endfor

i = nsec - 1
for j = 1, nrad-1 do begin
;    if(vrad(i,j) ge 0d0) then begin
;        urad = vrad(i,j) - vrad_1d(j)
;        uphi = 0.5*(vtheta(i, j-1) + vtheta(i+1,j-1)) - vtheta_1d(j-1)


;        mdot(i,j) = vrad(i,j)*sigma(i, j-1)*rad(j)*dphi ; mass flux        
;        hdot(i,j) = urad(i,j)*sigma(i, j-1)*uphi
        ;hdot(i,j)*=  rad(j)*dphi
;    endif else begin
        urad = vrad(i,j) - vrad_1d(j)
        uphi = 0.5*(vtheta(i, j-1) + vtheta(0,j-1)) - vtheta_1d(j-1)

        mdot(i,j) = vrad(i,j)*sigma(i, j)*rad(j)*dphi
        hdot(i,j) = urad*sigma(i, j)*uphi
        ;hdot(i,j)*=  rad(j)*dphi
;    endelse
endfor

mdisc = 0d0 
mdisc_inner = 0d0

for i = 0, nsec-2 do begin
    for j = 1, nrad-1 do begin
        dm = sigma(i,j)*rad(j)*rad(j)*dlogr
        mdisc += dm
        if(rad(j) le 3.0) then mdisc_inner += dm
;        if(vrad(i,j) ge 0d0) then begin
;            urad = vrad(i,j) - vrad_1d(j)
;            uphi = 0.5*(vtheta(i, j-1) + vtheta(i+1,j-1)) - vtheta_1d(j-1)

;            mdot(i,j) = urad*sigma(i, j-1)*rad(j)*dphi                
;            hdot(i,j) = vrad(i,j)*sigma(i, j-1)*uphi;*rmed(j-1)*rad(j)*dphi             
;        endif else begin
            urad = vrad(i,j) - vrad_1d(j)
            uphi = 0.5*(vtheta(i, j) + vtheta(i+1,j)) ;- vtheta_1d(j)

            uphi+= 0.5*(vtheta(i, j-1) + vtheta(i+1,j-1)) ;- vtheta_1d(j-1)

            uphi /= 2.0
 
            uphi -= 0.5*(vtheta_1d(j) + vtheta_1d(j-1))

            mdot(i,j) = urad*sigma(i, j)*rad(j)*dphi
            hdot(i,j) = urad*sigma(i, j)*uphi
;        endelse
    endfor
endfor

for j = 0, nrad - 1 do begin
    mdot_1d(j) = total(mdot(*,j)) ; mass across r per time
    hdot_1d(j) = mean(hdot(*,j)) ; ang. mom. across r per time
endfor

massloss(0,l-start) = info(7,l)/p0
massloss(1,l-start) = mdot_1d(1)*dt/mdisc
massloss(2,l-start) = mdot_1d(2)*dt/mdisc
massloss(3,l-start) = mdisc_inner/mdisc

;print, mdot(*,1)
temp = min(abs(width(0) - rad), grid1)
temp = min(abs(width(1) - rad), grid2)

mdot_1d /= mdisc*1d-7

massloss(4,l-start) = mdot_1d(grid1)
massloss(5,l-start) = mdot_1d(grid2)

angmomflux(0,l-start) = info(7,l)/p0
angmomflux(1,l-start) = hdot_1d(grid1)
angmomflux(2,l-start) = hdot_1d(grid2)


; set_plot, 'ps'
; time=string(info(7,l)/p0,format='(F7.2)')
; device, filename=filepath(strcompress('masstrans_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
; plot,rad(1:nrad-1), mdot_1d(1:nrad-1),title=time+' orbits',xmargin=[6,1],ymargin=[3,2],yrange=yrange $
;   ,ytitle='Fluxes',xtitle='r', charsize=1.5,xrange=xrange, thick=4, xtickinterval=1
; device,/close

openw,1,filepath(strcompress('masstrans_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
for i = 1, nrad-1 do printf,1, rad(i), mdot_1d(i), hdot_1d(i)
close,1

print, 'done '+ks
endfor

openw,1,filepath(strcompress('masstrans_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,massloss(0,i),massloss(1,i), massloss(2,i), massloss(3,i), massloss(4,i), massloss(5,i) $
, format='(6(e22.15,2x))'
close,1

openw,1,filepath(strcompress('masstrans_ang_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,angmomflux(0,i), angmomflux(1,i), angmomflux(2,i), format='(3(e22.15,2x))'
close,1
end

; mdot_visc = nu*deriv(rmed, sigma_1d*rmed^3*deriv(rmed,vtheta_1d/rmed))/deriv(rmed, rmed*vtheta_1d)

; for j= 0, nrad-1 do begin
;     mean_vtheta = mean(vtheta(0:nsec-1,j))
;     mean_vrad = mean(vrad_c(0:nsec-1,j))
;     alpha(0:nsec-1, j) = (vtheta(0:nsec-1, j) - mean_vtheta)*(vrad_c(0:nsec-1, j) - mean_vrad)*rad(j)^(3./2.) ; nu due to alpha
;     alpha_1d(j) = mean(alpha(*,j))
; endfor

;mdot_alpha = deriv(rmed,
;alpha_1d*sigma_1d*rmed^3*deriv(rmed,vtheta_1d/rmed))/deriv(rmed,
;rmed*vtheta_

;  set_plot, 'ps'

;  device, filename=filepath(strcompress('masstrans_time.ps',/remove_all),root_dir='.',subdir=[location])$
;    ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;  plot,mdot_time(0,*),mdot_time(2,*),xmargin=[6,1],ymargin=[3,2] $
;    ,ytitle=textoidl('Mdot\times10^6'),xtitle='t/orbits' $
;    ,charsize=1.5,xminor=25, thick=4, xrange=[25,200];, yrange=[0.08,0.10]
; ;  oplot,mdot_time(0,*),mdot_time(3,*), linestyle=1, thick=4
; ;  oplot,mdot_time(0,*),mdot_time(4,*), linestyle=2, thick=4
;  device,/close



; temp = min(abs(rad(1:nrad) - width(0)),rpmxs)
; temp = min(abs(rad(1:nrad) - width(1)),rppxs)

; mdot_time(0,l-start) = info(7,l)/p0
; mdot_time(1,l-start) = mean(mdot) ; over entire disc
; mdot_time(2,l-start) =   total(mdot(*,rpmxs:rppxs))         ;mean(mdot(*,rpmxs:rppxs)) ; over region near planet
; mdot_time(3,l-start) = mean(mdot_visc(rpmxs:rppxs))
; ;mdot_time(4,l-start) = mean(mdot_alpha(rpmxs:rppxs))
