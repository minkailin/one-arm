pro stress2, loc=loc, start=start, finish=finish, out=out, xrange=xrange, mcut=mcut

pi=3.141592654
G = 1.0
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

nrad = nrad0
nsec = nsec0/8

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR RADIUS AND AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp = dblarr(nrad+1)
rsource = dblarr(nrad)
azi    = dblarr(nsec)
darea = dblarr(nsec,nrad)

openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rsource(0:nrad-1) = (radtmp(0:nrad-1) + radtmp(1:nrad))/2.0

dphi = 2.0*pi/double(nsec)
azi = dphi*dindgen(nsec)

dlogr = alog(radtmp(1)/radtmp(0))/double(nrad)

for i=0, nsec-1 do begin
    for j=0, nrad-1 do begin
        darea(i,j) = rsource(j)*rsource(j)*dlogr*dphi
    endfor
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ARRAYS for self-gravity;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pot_integrand = dblarr(nsec, nrad)
potential = dblarr(nsec, nrad)
g_rad = dblarr(nsec, nrad)
g_phi = dblarr(nsec, nrad)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sigma =  dblarr(nsec0,nrad0)
sigma_reduced = dblarr(nsec,nrad)

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

;kernel_int_rad(i,*,k) = int_tabulated(rsource_3d(i,*,k), kernel_int_rad(i,*,k) , /double)

for i = 0, nsec-1 do begin
    for j=0, nrad-1 do begin
        r  = rsource(j)
        phi= azi(i)
        ;fill the integrand
        for iprime = 0, nsec-1 do begin
            for jprime = 0, nrad-1 do begin
                rprime = rsource(jprime)
                phiprime = azi(iprime)
                pot_integrand(iprime, jprime) = -G*sigma(iprime, jprime)*darea(iprime, jprime)
                pot_integrand(iprime, jprime)/= sqrt(r*r + rprime*rprime - 2.0*r*rprime*cos(phi - phiprime) $
                                                     + (0.3*h*rprime)^2)
            endfor
        endfor
        potential(i,j) = total(pot_integrand)
    endfor
endfor

for i = 0, nsec - 1 do begin
    g_rad(i, *) = -deriv(rsource(*), potential(i,*))
endfor

for j = 0, nrad - 1 do begin
    g_phi(*,j) = -deriv(azi(*), potential(*,j))/rsource(j)
endfor

alpha = - g_rad*g_phi/(4.0*!dpi)

for i = 0, nrad - 1 do alpha_1d(i) = mean(alpha(*,i))*1d3/(h*rfield(i))

print, 'here6'

 plx = info(1,l)
 ply = info(2,l)
 plrad = sqrt(plx*plx + ply*ply)
 rhill = f0*plrad

 ;rfield = (rfield-plrad)/rhill
 if not keyword_set(xrange) then xrange = [rfield(0),rfield(nrad-1)]
 temp = min(abs(rsource - xrange(0)),rpmxs)
 temp = min(abs(rsource - xrange(1)),rppxs)


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
