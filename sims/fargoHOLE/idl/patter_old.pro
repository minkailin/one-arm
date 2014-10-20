pro pattern, loc=loc, start=start, finish=finish, legend=legend, xrange=xrange, modes=modes, yrange=yrange $
             , plotmodes=plotmodes, out=out, xtickinterval=xtickinterval, width=width
common consts, pi
pi=3.141592654
mp=3d-4
f0=(mp/3.0)^(1.0/3.0)
if not keyword_set(finish) then finish=start
if not keyword_set(modes) then modes = 16

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

;;;;;;;
;SETUP;
;;;;;;;
;Arrays to hold raw data, and center-valued velocities.
sigma  = dblarr(nsec,nrad)
vtheta = dblarr(nsec,nrad)
vrad   = dblarr(nsec,nrad)

vthetac = dblarr(nsec,nrad)
vradc   = dblarr(nsec,nrad)

fft_sigma      = dcomplexarr(nsec,nrad)
fft_sigmavrad  = dcomplexarr(nsec,nrad)
fft_sigmavtheta= dcomplexarr(nsec,nrad)

am    = dcomplexarr(modes+1, nrad)
amdot = dcomplexarr(modes+1, nrad)

am_conj    = dcomplexarr(modes+1, nrad)
amdot_conj = dcomplexarr(modes+1, nrad)

int_am   = dcomplexarr(modes+1)
int_amdot= dcomplexarr(modes+1)

int_amconj    = dcomplexarr(modes+1)
int_amdotconj = dcomplexarr(modes+1)



rate = dblarr(modes+2, finish-start+1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DATA AND INTERPOLATE TO GET CELL-CENTER VELOCITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for k=start, finish do begin
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

;Get d(a_1)/dt 
for i = 0, nrad-1 do begin
    amdot(0:modes, i) = -(1d0/rmed(i))*fft_sigmavrad(0:modes,i) $
      - dcomplex(0d0, 1d0)*dindgen(modes+1)*fft_sigmavtheta(0:modes,i)/rmed(i)
endfor

am_conj = conj(am)
amdot_conj = conj(amdot)

;amabs = abs(am)
;amabsdot = real_part(am*amdot_conj + amdot*am_conj)/(2d0*amabs)

plx=info(1,k)
ply=info(2,k)
plrad=sqrt(plx*plx+ply*ply)
rhill = plrad*f0
rplot = (rmed - plrad)/rhill

if not keyword_set(width) then begin
    r0=1
    r1=nrad-2
endif else begin
    beg = min(abs(width(0)-rplot),r0)
    fin = min(abs(width(1)-rplot),r1)
endelse

for i=0, modes do begin
    int_am(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(am(i,r0:r1)),/double) $
                                 ,int_tabulated(rmed(r0:r1),imaginary(am(i,r0:r1)),/double))
    int_amdot(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(amdot(i,r0:r1)),/double) $
                                 ,int_tabulated(rmed(r0:r1),imaginary(amdot(i,r0:r1)),/double))
;int_am(i) = int_tabulated(rmed(r0:r1), amabs(i,r0:r1))
;int_amdot(i) = int_tabulated(rmed(r0:r1), amabsdot(i,r0:r1))
endfor

int_amconj    = conj(int_am)
int_amdotconj = conj(int_amdot)

rate(0, k - start) = info(7,k)/p0
rate(1:modes+1, k-start) = 0.5*real_part( int_amdotconj/int_amconj + int_amdot/int_am $
                                      -int_amdotconj(0)/int_amconj(0) - int_amdot(0)/int_am(0))
print, 'done', k

endfor

openw,1,filepath(strcompress('growth_rates.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,rate(*,i),format=strcompress("("+string(modes+2)+"(e12.5,1x))",/remove_all)
close,1

;PLOTTING
if keyword_set(plotmodes) then begin
    set_plot, 'ps'
    device, filename=filepath(strcompress('pattern_.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    plot,rate(0,*),rate(plotmodes(0)+1,*),xmargin=[8,4],ymargin=[3,3] $
      ,xtitle='t/orbits',ytitle=textoidl('growth rate'),charsize=1.5,thick=4, xrange=xrange, yrange=yrange $
      ,xtickinterval=xtickinterval
    
    nelements = n_elements(plotmodes)
    
    if( nelements gt 1) then begin
        for i = 1, nelements - 1 do begin
            oplot,rate(0,*),rate(plotmodes(i)+1,*),thick=4, linestyle=i
        endfor
    endif
    
    if keyword_set(legend) then begin
        x0 = legend(0)
        x1 = legend(1)
        y0 = legend(2)
        dy = legend(3)

        for i=0, nelements-1 do begin 
            oplot, [x0,x1], [y0,y0] - i*dy, thick=4, linestyle=i
            xyouts,x1,y0-i*dy,label(i),charsize=1.5
        endfor
    endif
    device,/close
endif

end
