function get_vortensity, rad, vtheta, sigma
kappasq = deriv(rad, rad*rad*vtheta*vtheta)/(rad^3.0)
vortensity = kappasq/(2.0*vtheta/rad)
return, kappasq
end

pro compare_profiles, type=type, cases=cases, start=start, out=out, legend=legend,label=label $
                      ,xrange=xrange, yrange=yrange,mass=mass,mode=mode

pi = !dpi
; mp = mass(0)
; f0=(mp/3.0)^(1.0/3.0)

cases=strcompress(cases,/remove_all)
numcases=n_elements(cases)

;GRID INFO

dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[cases(0)]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))

if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11, out(0)+1)

openr,3,filepath('planet0.dat',root_dir='.',subdir=[cases(0)])
readf,3,info
close,3

a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
time = string(info(7,start(0))/p0,format='(F7.2)')
; plx = info(1,start(0))
; ply = info(2,start(0))
; plrad = sqrt(plx*plx + ply*ply)
; rhill = f0*plrad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi   =dblarr(nsec)
radtmp=dblarr(nrad+1)
rad   =dblarr(nrad)
data  =dblarr(nsec,nrad)

sigma0=dblarr(nrad)
sigma =dblarr(nrad)
vtheta0=dblarr(nrad)
vtheta = dblarr(nrad)
 
;data_1d = dblarr(nrad)
;data0_1d = dblarr(nrad)
;data2_1d=dblarr(nrad)

openr,1,filepath('used_rad.dat',root_dir='.',subdir=[cases(0)])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

openr,2,filepath(strcompress('gasdens'+'0'+'.dat',/remove_all),root_dir='.',subdir=[cases(0)])
readu,2,data
close,2
for i=0, nrad-1 do sigma0(i) = mean(data(*,i))

openr,3,filepath(strcompress('gasdens'+string(start(0))+'.dat',/remove_all),root_dir='.',subdir=[cases(0)])
readu,3,data
close,3
for i=0, nrad-1 do sigma(i) = mean(data(*,i))

openr,3,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[cases(0)])
readu,3,data
close,3    
for j=0, nrad-1 do vtheta0(j) = mean(data(*,j))
   
openr,3,filepath(strcompress('gasvtheta'+string(start(0))+'.dat',/remove_all),root_dir='.',subdir=[cases(0)])
readu,3,data
close,3
for i=0, nrad-1 do vtheta(i) = mean(data(*,i))

; if (type eq'vtheta') then begin
;     data_1d /= rad^(-0.5)
;     data_1d -= 1.0
;     yytitle = textoidl('(<v_\phi>_\phi-v_{kep})/v_{kep}')
;yytitle = textoidl('<\Delta\eta>_\phi\times 10^{-2}')
;yytitle = textoidl("relative change in |\Omega^{'}|")
;     yytitle = textoidl('(<\Omega_\phi>_\phi-\Omega_{kep})\times10^4') 
; endif else begin
;     data_1d /= data0_1d
;     data_1d -= 1.0
     yytitle = textoidl('\Delta <\Sigma>_\phi/\Sigma_{t=0}')
; ;    yytitle = 'relative gap depth'
; endelse
  
;domega = deriv(rad,vtheta/rad)
;ddomega = deriv(rad,domega)

;dsigma = deriv(rad,sigma)

;temp = min(abs(rad - 5.0),plgrid)
;omega_p = omega(plgrid)

;vortensity0 = get_vortensity(rad, vtheta0, sigma0)
vortensity = get_vortensity(rad,vtheta,sigma)
delsq_omega = deriv(rad, rad*deriv(rad, vortensity))/rad

;kappa = get_vortensity(rad,vtheta,sigma)

set_plot, 'ps'
device, filename=strcompress('compare_profiles_'+type+string(start(0),format='(I03)')+'.ps',/remove_all) $ 
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches


plot,rad, (sigma-sigma0)/sigma0, title=time+' orbits',xmargin=[8,2],ymargin=[3,2] $
,xtitle='r',ytitle=yytitle  $
,charsize=1.5, thick=4, xrange=xrange, yrange=yrange;,/nodata

;oplot, rad, abs(dsigma/sigma), thick=2

;oplot, rad, dblarr(nrad), thick=1

;  temp=min(abs(rad- 5.0),plgrid)
;  temp=min(abs(rad- 5.58),beg)
;  temp=min(abs(rad- 8.0), fin) 
;  x = rad(beg:fin)/5.0
;  domega=deriv(rad, data_1d)/(-1.5*rad^(-2.5))
;    dfdx = deriv(x, data_1d(beg:fin)*5.0^1.5) 
;    y1 = (2.0/3.0)*int_tabulated(x,x*x*dfdx*(x-1.0)^(-4.0),/double)
;    y2 = int_tabulated(x,x^(-0.5)*(x-1.0)^(-4.0),/double) 
;    print, 'rel change at rp', data_1d(plgrid)*rad(plgrid)^(1.5)
;    print, 'rel change at re', domega(beg)
;    print, 'rel change at rinfty', domega(fin)
;    print, 'ratio of integrals', y1/y2
; print, 'dens at rp', data_1d(beg)/data_1d(fin)

for i=1, numcases-1 do begin
    dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[cases(i)]))).(0)
    nout=fix(dims(5))
    nrad=fix(dims(6))
    nsec=fix(dims(7))
    if not keyword_set(out) then begin
        info=dblarr(11,nout+1)
    endif else info=dblarr(11, out(i)+1)

    radtmp=dblarr(nrad+1)
    rad = dblarr(nrad)
    openr,1,filepath('used_rad.dat',root_dir='.',subdir=[cases(i)])
    readf,1,radtmp
    close,1
    rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

    openr,3,filepath('planet0.dat',root_dir='.',subdir=[cases(i)])
    readf,3,info
    close,3
    
;    mp = mass(i)
;     f0=(mp/3.0)^(1.0/3.0)

;     plx = info(1,start(i))
;     ply = info(2,start(i))
;     plrad = sqrt(plx*plx + ply*ply)
;     temp = min(abs(rad-plrad),plgrid)
;     rhill = f0*plrad

    data  =dblarr(nsec,nrad)
    sigma0 =dblarr(nrad)
    sigma = dblarr(nrad)
    vtheta0 =dblarr(nrad)
    vtheta = dblarr(nrad)

    openr,2,filepath(strcompress('gasdens'+'0'+'.dat',/remove_all),root_dir='.',subdir=[cases(i)])
    readu,2,data
    close,2    
    for j=0, nrad-1 do sigma0(j) = mean(data(*,j))

    openr,3,filepath(strcompress('gasdens'+string(start(i))+'.dat',/remove_all),root_dir='.',subdir=[cases(i)])
    readu,3,data
    close,3    
    for j=0, nrad-1 do sigma(j) = mean(data(*,j))

    openr,3,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[cases(i)])
    readu,3,data
    close,3    
    for j=0, nrad-1 do vtheta0(j) = mean(data(*,j))

    openr,3,filepath(strcompress('gasvtheta'+string(start(i))+'.dat',/remove_all),root_dir='.',subdir=[cases(i)])
    readu,3,data
    close,3    
    for j=0, nrad-1 do vtheta(j) = mean(data(*,j))

    vortensity0 = get_vortensity(rad, vtheta0, sigma0)
    vortensity = get_vortensity(rad,vtheta,sigma)


;    omega = vtheta/rad
;    temp = min(abs(rad - 5.0),plgrid)
;    omega_p = omega(plgrid)

;    kappa = sqrt(get_vortensity(rad,vtheta,sigma))

;    if(i eq 1) then sigma_nsg = data_1d

;    if(i eq 2) then sigma_nsg2 = data_1d

;     if (type eq 'vtheta') then begin
;          data_1d /= rad^(-0.5)
;          data_1d -= 1.0
;      endif else begin   
;          data_1d /= data0_1d
;          data_1d -= 1.0
;      endelse


;     omega_nsg = data_1d/rad

;     temp=min(abs(rad- 5.0),plgrid)
;     temp=min(abs(rad- 5.58),beg)
;     temp=min(abs(rad- 8.0), fin)
    
;     x = rad(beg:fin)/5.0
;     domega =  omega_sg - omega_nsg
;     domegadr_nsg = deriv(rad, omega_nsg)
;     ddomegadr = deriv(rad,omega_sg - omega_nsg)
;     dgdx = deriv(x, omega_nsg(beg:fin)*5.0^1.5)
;     dfdx = deriv(x, domega(beg:fin)*5.0^1.5)
;     y1 = int_tabulated(x,x^(-3.0)*dfdx*dgdx^(-2.0)*(x-1.0)^(-4.0),/double)
;     y2 = int_tabulated(x,x^(-3.0)*dgdx^(-1.0)*(x-1.0)^(-4.0),/double)
;     print, 'rel change at rp', domega(plgrid)/omega_nsg(plgrid)
;     print, 'rel change at re', ddomegadr(beg)/domegadr_nsg(beg)
;     print, 'rel change at rinfty', ddomegadr(fin)/domegadr_nsg(fin)
;     print, 'ratio of integrals', -y1/y2

;      sigma_nsg = data_1d    
;      openr,3,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[cases(i)])
;      readu,3,data
;      close,3
;      for j=0, nrad-1 do data_1d(j) = mean(data(*,j))
;      omega_nsg = data_1d/rad
;      delta_omega = omega_sg - omega_nsg
;      delta_omegap = delta_omega(plgrid)
;      delta_sigma = sigma_sg - sigma_nsg   
;      lhs = 3.0*(rad-plrad)^4.0*!dpi*1e-6*sqrt(rad)*deriv(rad,delta_sigma)/sigma_nsg
;      lhs+= ((rad-plrad)^4.0*!dpi*1e-6/sqrt(rad) - 6.0*0.836*mp^2.0*plrad^4.0)*delta_sigma/sigma_nsg
;      rhs = 6.0*0.836*mp^2.0*plrad^4.0*delta_omegap/omega_nsg(plgrid)
;      rhs+=
;      0.5*(rad-plrad)^4.0*!dpi*1e-6/sqrt(rad)*(4.0*rad^(7.0/2.0)*(deriv(rad,delta_omega)* $
;      (deriv(rad, alog(sigma_nsg)) +3.0/rad)+ deriv(rad, deriv(rad, delta_omega))))   
    

;    domega = deriv(rad, vtheta/rad)                                           
;    ddomega = deriv(rad,domega)
    
;    dsigma = deriv(rad,sigma)
               

    vortensity = get_vortensity(rad,vtheta,sigma)
    delsq_omega = deriv(rad, rad*deriv(rad, vortensity))/rad

    oplot, rad,  (sigma-sigma0)/sigma0, thick=4, linestyle=i

;    oplot, rad, abs(dsigma/sigma), thick=2, linestyle=i

;oplot, rad, dblarr(nrad), thick=1

endfor

;oplot, rad, (sigma_sg-sigma_nsg)/sigma_nsg, thick=2, linestyle=0
;oplot, rad, (sigma_nsg2 - sigma_nsg)/sigma_nsg, thick=2, linestyle=1

if keyword_set(legend) then begin
x0=legend(0)
x1=legend(1)
y0=legend(2)
dy=legend(3)
for j=0, n_elements(label)-1 do begin
    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j+1
    xyouts, x1, y0-dy*j,label(j),charsize=1.5
endfor
endif

device,/close
end



