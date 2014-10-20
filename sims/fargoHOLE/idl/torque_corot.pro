pro torque_corot, loc=loc, start=start, finish=finish, out=out, mp=mp

if not keyword_set(finish) then finish=start
common consts, pi
pi=!dpi
if not keyword_set(mp) then mp=2d-3
f0=(mp/3.0)^(1.0/3.0)
h=0.05

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
if not keyword_set (out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
period=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.0
dlogr = alog(radtmp(nrad)/radtmp(0))/nrad

for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
azi1 = azi/pi - 1d0
dphi = 2.0*!dpi/nsec

;;;;;;;;;;;;;
;TORQUE CALC;
;;;;;;;;;;;;;

density = dblarr(nsec,nrad)
torque2d = dblarr(nsec, nrad)
torque2d_nohill =  dblarr(nsec, nrad)

torque2d_shifted = dblarr(nsec, nrad)
torque2d_shifted_nohill = dblarr(nsec, nrad)


density1d = dblarr(nrad)
density1d_0 = dblarr(nrad)

torque1d_rad = dblarr(nrad)
torque1d_azi = dblarr(nsec)

torque1d_rad_nohill = dblarr(nrad)
torque1d_azi_nohill = dblarr(nsec)



openr,2,filepath(strcompress('gasdens0'+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,density
close,2
for i=0, nrad-1 do density1d_0(i) = mean(density(*,i))


for k=start, finish do begin
    ks=string(k,format='(I03)')
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,density
    close,2
    for i=0, nrad-1 do density1d(i) = mean(density(*,i))
    dsigma = (density1d - density1d_0)/density1d_0

    
    plx=info(1,k)
    ply=info(2,k)
    plrad = sqrt(plx*plx + ply*ply)
    phi=pltphi(plx,ply)

    temp = min(abs(rad-plrad),r0)    
    rhill = plrad*f0
    rplot = (rad - plrad)/rhill

    ;find inner gap edge
    for i = r0, 1, -1 do begin
        minus = dsigma(i)
        plus  = dsigma(i - 1)
        if(plus*minus lt 0d0) then begin
            rpmxs = i
            break
        endif
    endfor
  
    ; find outer gap edge
    for i = r0, nrad-2 do begin
        minus = dsigma(i)
        plus  = dsigma(i + 1)
        if(plus*minus lt 0d0) then begin
            rppxs = i
            break
        endif
    endfor

    
    ; fill torque array
    for i=0, nsec-1 do begin
        for j=0,  nrad-1 do begin
            dm = rad(j)*rad(j)*dphi*dlogr*density(i,j)           
            x = rad(j)*cos(azi(i))
            y = rad(j)*sin(azi(i))
            d = sqrt((x-plx)^2.0 + (y-ply)^2.0 + (0.6*rad(j)*h)^2.0)
            
            dm2 = dm/d^3.0
            fx = (x-plx)*dm2
            fy = (y-ply)*dm2
            torque = plx*fy - ply*fx                 

            torque2d(i,j) = torque

            if (d gt rhill) then torque2d_nohill(i,j) = torque
        endfor
    endfor

    ;azimuthal shift the torque array to make planet at center
    if(phi gt !dpi) then begin
        temp2=min(abs(azi-(phi-!dpi)),grid2)
        torque2d_shifted(0:nsec-1-grid2, 0:nrad-1) = torque2d(grid2:nsec-1, 0:nrad-1)
        torque2d_shifted(nsec-grid2:nsec-1, 0:nrad-1) = torque2d(0:grid2-1, 0:nrad-1 )

 ;       torque2d_shifted_nohill(0:nsec-1-grid2, 0:nrad-1) = torque2d_nohill(grid2:nsec-1, 0:nrad-1)
 ;       torque2d_shifted_nohill(nsec-grid2:nsec-1, 0:nrad-1) = torque2d_nohill(0:grid2-1, 0:nrad-1 )
    endif
    
    if(phi lt !dpi) then begin
        temp2=min(abs(azi-(phi+!dpi)),grid2)
        torque2d_shifted(nsec-grid2:nsec-1, 0:nrad-1) = torque2d(0:grid2-1, 0:nrad-1)
        torque2d_shifted(0:nsec-1-grid2, 0:nrad-1) = torque2d(grid2:nsec-1, 0:nrad-1)

        
 ;        torque2d_shifted_nohill(nsec-grid2:nsec-1, 0:nrad-1) = torque2d_nohill(0:grid2-1, 0:nrad-1)
 ;       torque2d_shifted_nohill(0:nsec-1-grid2, 0:nrad-1) = torque2d_nohill(grid2:nsec-1, 0:nrad-1)
    endif
    
    if(phi eq !dpi) then begin
        torque2d_shifted = torque2d
 ;       torque2d_shifted_nohill = torque2d_nohill
    endif



    for i=0, nrad-1 do begin
        torque1d_rad(i) = total(torque2d_shifted(*,i))/(rad(i)*dlogr)
 ;       torque1d_rad_nohill(i) = total(torque2d_shifted_nohill(*,i))
    endfor

    for j = 0, nsec-1 do begin
        torque1d_azi(j) = total(torque2d_shifted(j, rpmxs : rppxs))/dphi
 ;       torque1d_azi_nohill(j) = total(torque2d_shifted_nohill(j, rpmxs : rppxs))
    endfor

    openw,1,filepath(strcompress('torque1d_rad_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 0, nrad-1 do printf,1,rplot(i),torque1d_rad(i)
    close,1
    
 ;   openw,1,filepath(strcompress('torque1d_rad_nohill'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
 ;   for i = 0, nrad-1 do printf,1,rplot(i),torque1d_rad_nohill(i)
 ;   close,1

    openw,1,filepath(strcompress('torque1d_azi_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    for j = 0, nsec-1 do printf,1,azi1(j),torque1d_azi(j)
    close,1
    
 ;   openw,1,filepath(strcompress('torque1d_azi_nohill'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
 ;   for j = 0, nsec-1 do printf,1,azi1(j),torque1d_azi(j)
 ;   close,1
        
print, 'done '+ks

endfor

end



; openw,1,filepath(strcompress('tq_corot_time.dat',/remove_all),root_dir='.',subdir=[location])
; for i = 0, finish-start do printf,1,torque_time(0,i),torque_time(1,i),torque_time(2,i), torque_time(3,i) $
;   ,torque_time(4,i)
; close,1


; set_plot, 'ps'
; device, filename=filepath(strcompress('tq_corot_out_exl.ps',/remove_all),root_dir='.',subdir=[location])$
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
; plot,torque_time(0,*),torque_time(2,*)*1d4,xmargin=[6,2],ymargin=[3,1] $
;   ,ytitle=textoidl('Torque from outer gap'),xtitle='t/orbits' $
;   ,charsize=1.5,xminor=4, thick=4
; device,/close
