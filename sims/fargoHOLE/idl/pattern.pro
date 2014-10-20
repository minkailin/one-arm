function epicyclesq, rad, vtheta
kappasq = deriv(rad, rad*rad*vtheta*vtheta)/(rad^3.0)
return, kappasq
end

function get_vortensity, rad, vtheta, sigma
vortensity = epicyclesq(rad, vtheta)/(2.0*vtheta/rad)
return, vortensity
end

function soundspeed, rad
h=5d-2
return, h/sqrt(rad)
end


pro pattern, loc=loc, start=start, finish=finish, legend=legend, xrange=xrange, modes=modes, yrange=yrange $
             , plotmodes=plotmodes, out=out, xtickinterval=xtickinterval, job=job $
             , getpatt=getpatt, mp=mp
common consts, pi
pi=3.141592654
f0=(mp/3.0)^(1.0/3.0)
if not keyword_set(finish) then finish=start
if not keyword_set(modes) then modes = 10

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION,GET RUN INFO, READ PLANET ORBIT INFO.;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)

azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5

azi=dindgen(nsec)*2.*pi/nsec
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)

case job of
'calc': begin
;;;;;;;
;SETUP;
;;;;;;;

;Arrays to hold raw data, and center-valued velocities.
sigma  = dblarr(nsec,nrad)
vtheta = dblarr(nsec,nrad)
vrad   = dblarr(nsec,nrad)

sigma0  = dblarr(nsec,nrad)
vtheta0 = dblarr(nsec,nrad)
vrad0   = dblarr(nsec,nrad)

sigma_1d = dblarr(nrad) 
sigma_1d0= dblarr(nrad)
vtheta_1d = dblarr(nrad)

vthetac = dblarr(nsec,nrad)
vradc   = dblarr(nsec,nrad)

fft_sigma      = dcomplexarr(nsec,nrad)
fft_sigmavrad  = dcomplexarr(nsec,nrad)
fft_sigmavtheta= dcomplexarr(nsec,nrad)

am    = dcomplexarr(modes+1, nrad)
amdot = dcomplexarr(modes+1, nrad)

integrated_amp = dcomplexarr(modes+1)
integrated_ampdot = dcomplexarr(modes+1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DATA AND INTERPOLATE TO GET CELL-CENTER VELOCITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
for i = 0, nrad -1 do begin 
    sigma_1d0(i) = mean(sigma(*,i))
endfor


openw,10,filepath(strcompress('mode_amplitudes.dat',/remove_all),root_dir='.',subdir=[location])
openw,20,filepath(strcompress('mode_growthrate.dat',/remove_all),root_dir='.',subdir=[location])

for k=start, finish do begin
time = info(7,k)/p0
ks=string(k,format='(I03)')

;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,3,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vtheta
close,3
openr,4,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,4,vrad
close,4

for i = 0, nsec-2 do vthetac(i,*) = (vtheta(i,*) + vtheta(i+1,*))/2.0
vthetac(nsec-1,*) = (vtheta(nsec-1,*) + vtheta(0,*))/2.0

for i = 0, nrad -1 do begin 
    vtheta_1d(i) = mean(vthetac(*,i))
    sigma_1d(i) = mean(sigma(*,i))
endfor

for j = 0, nrad-2 do vradc(*,j) = (vrad(*,j) + vrad(*,j+1))/2.0
vradc(*,nrad-1) = vrad(*,nrad-1)

sigmavrad   = radius*sigma*vradc

for i = 0, nsec-1 do begin
    sigmavrad(i,*) = deriv(rmed(*), sigmavrad(i,*)) ;d/dr(r sigma u_r)
endfor
sigmavtheta = sigma*vthetac

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FFT w.r.t. phi of sigma, sigmavrad and sigmavtheta;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for i=0, nrad-1 do begin
    fft_sigma(*,i)       = fft(sigma(*,i),/double) 
    fft_sigmavrad(*,i)   = fft(sigmavrad(*,i),/double) 
    fft_sigmavtheta(*,i) = fft(sigmavtheta(*,i),/double)
endfor

;Extract the m=modes components

am = fft_sigma(0:modes,*)

;Get amdot 
for i = 0, nrad-1 do begin
    amdot(0:modes, i) = -(1d0/rmed(i))*fft_sigmavrad(0:modes,i) $
      - dcomplex(0d0, 1d0)*dindgen(modes+1)*fft_sigmavtheta(0:modes,i)/rmed(i)
endfor

;Integrate the modes, outer disc

;dsigma = (sigma_1d - sigma_1d0)/sigma_1d0

plx = info(1, k)
ply = info(2, k)
plrad = sqrt(plx*plx + ply*ply)
temp = min(abs(rmed - plrad), grid)

;for i = grid, nrad-2 do begin
;    minus = dsigma(i)
;    plus  = dsigma(i + 1)
;    if(plus*minus lt 0d0) then begin
;        r0 = i
;        break
;    endif
;endfor
;halfwidth  = 0.5*(rout - plrad)
;temp = min(abs(rmed - (plrad + halfwidth)), r0)

r0 = grid
r1 = nrad - 1

for i=0, modes do begin
        integrated_amp(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(am(i,r0:r1)),/double) $
                                     ,int_tabulated(rmed(r0:r1),imaginary(am(i,r0:r1)),/double))
        integrated_ampdot(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(amdot(i,r0:r1)),/double) $
                                     ,int_tabulated(rmed(r0:r1),imaginary(amdot(i,r0:r1)),/double))
endfor

;Mode amplitudes
printf,10, time, abs(integrated_amp(1:modes)/integrated_amp(0)) $ 
  ,format=strcompress('('+string(modes+1)+'(e10.4,2x))',/remove_all) 

;Mode growth rate
cm = integrated_amp/integrated_amp(0)
cmdot = integrated_ampdot/integrated_amp(0) - integrated_ampdot(0)*integrated_amp/(integrated_amp(0)^2)
printf, 20, time, real_part(cmdot(1:modes)/cm(1:modes)) $
  , format=strcompress('('+string(modes+1)+'(1e12.5,2x))',/remove_all)

if keyword_set(getpatt) then begin
;Pattern speeds
    x = real_part(am)
    y = imaginary(am)
    xdot = real_part(amdot)
    ydot = imaginary(amdot)
    pattspeed = (xdot*y - ydot*x)/(x*x + y*y)

 kappa = sqrt(epicyclesq(rmed, vtheta_1d))
 ilr = vtheta_1d/rmed + kappa/plotmodes(0)
 olr = vtheta_1d/rmed - kappa/plotmodes(0)
    set_plot, 'ps'
    device, filename=filepath(strcompress('pattern_speed'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    plot, rmed,pattspeed(plotmodes(0),*)/plotmodes(0) ,xmargin=[8,2],ymargin=[3,2] $
      ,xtitle='r',ytitle=textoidl('\Omega_{pat}'),charsize=1.5,thick=4, xrange=xrange, yrange=yrange $
      ,xtickinterval=xtickinterval, title=textoidl('m=')+string(plotmodes(0),format='(I1)')
     oplot, rmed, vtheta_1d/rmed, linestyle=2, thick=1
     oplot, rmed, ilr, linestyle=2, thick=1
     oplot, rmed, olr, linestyle=2, thick=1

     xyouts, 6, 0.1, textoidl('\Omega + \kappa/3'), charsize=1.5
     xyouts, 5.05, 0.09, textoidl('\Omega'), charsize=1.5
     xyouts, 3.55, 0.1, textoidl('\Omega - \kappa/3'), charsize=1.5
     device,/close
 endif
 
 print, 'done'+ks
endfor

close, 10
close, 20
end

'plot': begin

    temp=READ_ASCII(filepath('mode_amplitudes.dat',root_dir='.',subdir=[location]))
    temp=double(temp.(0))
    
     num = n_elements(temp(0,*))
     time_avg = dblarr(num)
     for j = 0, num-1 do time_avg(j) = mean(temp(plotmodes(0), j))
  
    set_plot, 'ps'
    device, filename=filepath('mode_amp_time.ps',root_dir='.',subdir=[location]) $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

    plot,temp(0,*), alog10(temp(plotmodes(0),*)) ,xmargin=[6,2],ymargin=[3,2] $
      ,ytitle=textoidl('log |C_2|'),xtitle='t/orbits' $
      ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
      ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval

    numcases = n_elements(plotmodes)
    for i=1, numcases-1 do begin
        oplot,temp(0,*) ,alog10(temp(plotmodes(i),*)),thick=4,linestyle=i
    endfor
 
  if keyword_set(legend) then begin
      x0=legend(0)
      x1=legend(1)
      y0=legend(2)
      dy=legend(3)
      for j=0, numcases-1 do begin
          oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
          xyouts, x1, y0-dy*j,label(j),charsize=1.5
      endfor
  endif
  device,/close
end
endcase

end

