pro stress, loc=loc, start=start, finish=finish, out=out, xrange=xrange, mcut=mcut

pi=3.141592654
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
nrad0=fix(dims(6))
nsec0=fix(dims(7))
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

;reduce sizes
nrad = nrad0
nsec = nsec0/8

;if not keyword_set(mcut) then mcut = nsec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR RADIUS AND AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp = dblarr(nrad0+1)
azi    = fltarr(nsec)
rsource = dblarr(nrad)
rfield  = dblarr(nrad)

rsource_3d = dblarr(nsec,nrad,nrad)
rfield_3d = dblarr(nsec,nrad,nrad)
azi_3d = dblarr(nsec,nrad,nrad)


openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rsource(0:nrad-1) = rebin((radtmp(0:nrad0-1) + radtmp(1:nrad0))/2.0, nrad)
rfield = rsource

azi = 2.0*pi*dindgen(nsec)/double(nsec)


for i=0, nsec-1 do begin
    for k=0, nrad-1 do begin
        rsource_3d(i,0:nrad-1,k)=rsource(0:nrad-1)
        rfield_3d(i,k,0:nrad-1) =rfield(0:nrad-1)
    endfor
endfor

for j=0, nrad-1 do begin
    for k=0, nrad-1 do begin
        azi_3d(*,j,k) = azi(*)
    endfor
endfor

print, 'here1'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;KERNEL ARRAYS for self-gravity;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
kernel_rad = complexarr(nsec,nrad,nrad)
kernel_phi = complexarr(nsec,nrad,nrad)

kernel_int_rad = complexarr(nsec,nrad,nrad)
kernel_int_phi = complexarr(nsec,nrad,nrad)

g_rad = complexarr(nsec,nrad)
g_phi = complexarr(nsec,nrad)

; for i=0, nsec-1 do begin
;     phi = azi(i)
;     for j=0, nrad-1 do begin
;         r = rfield(j)
;         for k=0, nrad-1 do begin
;             rp = rsource(k)
;             temp = r*r + rp*rp - 2.0*r*rp*cos(phi) + (0.3*h*rp)^2.0
;             temp = temp^(3.0/2.0)
            
;             kernel_rad(i,j,k) = rp*(r - rp*cos(phi))/temp
;             kernel_phi(i,j,k) = rp*(rp*sin(phi))/temp
;         endfor
;     endfor
; endfor

r = rfield_3d
rp = rsource_3d
phi = azi_3d

temp = r*r + rp*rp - 2.0*r*rp*cos(phi) + (0.3*h*rp)^2.0
temp = temp^(3.0/2.0)
kernel_rad = rp*(r - rp*cos(phi))/temp
kernel_phi = rp*(rp*sin(phi))/temp

print, 'here2'

; for j = 0, nrad-1 do begin
;     for k = 0, nrad-1 do begin
;         kernel_rad(0:nsec-1,j,k) = fft(kernel_rad(0:nsec-1,j,k), -1, /double, /overwrite)
;         kernel_phi(0:nsec-1,j,k) = fft(kernel_phi(0:nsec-1,j,k), -1, /double, /overwrite)
;     endfor
; endfor

kernel_rad = fft(kernel_rad, -1, /double, /overwrite, dimension=1)
kernel_phi = fft(kernel_phi, -1, /double, /overwrite, dimension=1)

print, 'here3'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sigma =  dblarr(nsec0,nrad0)
sigma_reduced = dblarr(nsec,nrad)
sigma_fft = dcomplexarr(nsec,nrad)

; vrad   = dblarr(nsec,nrad)
; vtheta = dblarr(nsec,nrad)

alpha_1d  = dblarr(nrad)
alpha_time=dblarr(3,finish-start+1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET THE SELF-GRAVTY ACCELERATIONS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for l=start, finish do begin

    ks=string(l,format='(I03)')

     openr,2,filepath(strcompress('gasdens'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,sigma
     close,2 

sigma_reduced = rebin(sigma,nsec,nrad)
   
;for j=0, nrad-1 do sigma_fft(0:nsec-1,j) = fft(sigma_reduced(0:nsec-1,j), -1, /double)

sigma_fft = fft(sigma_reduced, -1, /double, dimension=1)

;for i = 0, nsec-1 do begin
    for j = 0, nrad - 1 do begin
 ;       for k = 0, nrad - 1 do begin
            kernel_int_rad(0:nsec-1,j,0:nrad-1) = kernel_rad(0:nsec-1,j,0:nrad-1)*sigma_fft(0:nsec-1,0:nrad-1)
            kernel_int_phi(0:nsec-1,j,0:nrad-1) = kernel_phi(0:nsec-1,j,0:nrad-1)*sigma_fft(0:nsec-1,0:nrad-1)
 ;       endfor
    endfor
;endfor


print, 'here4'

for i = 0, nsec - 1 do begin
    for k = 0, nrad - 1 do begin
        real_int = int_tabulated(rsource(*), real_part(kernel_int_rad(i,*,k)), /double)
        cplx_int = int_tabulated(rsource(*), imaginary(kernel_int_rad(i,*,k)), /double)
        g_rad(i,k) = dcomplex(real_int, cplx_int)

        real_int = int_tabulated(rsource(*), real_part(kernel_int_phi(i,*,k)), /double)
        cplx_int = int_tabulated(rsource(*), imaginary(kernel_int_phi(i,*,k)), /double)
        g_phi(i,k) = dcomplex(real_int, cplx_int)
    endfor
endfor

;kernel_int_rad(i,*,k) = int_tabulated(rsource_3d(i,*,k), kernel_int_rad(i,*,k) , /double)

print, 'here5'

;for k = 0, nrad-1 do begin
;    g_rad(*,k) = fft(g_rad(*,k),/double,/inverse)
;    g_phi(*,k) = fft(g_phi(*,k),/double,/inverse)
;endfor
g_rad = fft(g_rad,/double,/inverse,/overwrite,dimension=1)
g_phi = fft(g_phi,/double,/inverse,/overwrite,dimension=1)

alpha = - real_part(g_rad*g_phi)/(4.0*!dpi)

for i = 0, nrad - 1 do alpha_1d(i) = mean(alpha(*,i))*1d3/(h*rfield(i))

print, 'here6'

 plx = info(1,l)
 ply = info(2,l)
 plrad = sqrt(plx*plx + ply*ply)
 rhill = f0*plrad

 ;rfield = (rfield-plrad)/rhill
 if not keyword_set(xrange) then xrange = [rfield(0),rfield(nrad-1)]
 temp = min(abs(rfield - xrange(0)),rpmxs)
 temp = min(abs(rfield - xrange(1)),rppxs)


 alpha_time(0,l-start) = info(7,l)/p0
 alpha_time(1,l-start) = mean(alpha) ; over entire disc
 alpha_time(2,l-start) = mean(alpha(*,rpmxs:rppxs)) ; over region near planet


set_plot, 'ps'
time=string(info(7,l)/p0,format='(F7.2)')
device, filename=filepath(strcompress('stress'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot,rfield,alpha_1d ,title=time+' orbits',xmargin=[9,1],ymargin=[3,2] $
,ytitle=textoidl('<\alpha>_\phi\times10^{3}'),xtitle='r' $
,charsize=1.5,xrange=xrange, thick=4
device,/close

openw,1,filepath(strcompress('stress_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
for i = 1, nrad-1 do printf,1, rfield(i), alpha_1d(i)
close,1

print, 'done '+ks

endfor



; set_plot, 'ps'

; device, filename=filepath(strcompress('alpha_time.ps',/remove_all),root_dir='.',subdir=[location])$
; ,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
; plot,alpha_time(0,*),alpha_time(1,*),xmargin=[8,4],ymargin=[3,3] $
; ,ytitle=textoidl('<\alpha>\times10^{3}'),xtitle='t' $
; ,charsize=1.5,xminor=4, thick=4, xrange=[20,200]
; device,/close

 openw,1,filepath(strcompress('alphaSG_time.dat',/remove_all),root_dir='.',subdir=[location])
 for i = 0, finish-start do printf,1,alpha_time(0,i),alpha_time(1,i), alpha_time(2,i)
 close,1

end
