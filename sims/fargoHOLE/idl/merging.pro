pro merging, loc=loc, start=start, finish=finish, width=width, out=out, modes=modes, plotmode=plotmode $
             ,xrange=xrange,xtickinterval=xtickinterval, yrange=yrange

pi=3.141592654
mp=3d-4
f0=(mp/3.0)^(1.0/3.0)
if not keyword_set(modes) then  modes = 16
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
endif else info=dblarr(11,out+1)
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

rad(0:nrad-1) = 0.5*(radtmp(0:nrad-1) + radtmp(1:nrad))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sigma   = dblarr(nsec,nrad)
fft_sigma      = dcomplexarr(nsec,nrad)
integrated_amp = dcomplexarr(modes+1)
vortex_time = dblarr(modes+1,finish-start+1); time, m=1...m=modes
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE BEGIN PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
for k=start, finish do begin
    ks=string(k,format='(I03)')

    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location])
    readu,2,sigma
    close,2

    if not keyword_set(width) then begin
        r0=1
        r1=nrad-2
    endif else begin
        beg = min(abs(width(0)-rad),r0)
        fin = min(abs(width(1)-rad),r1)
    endelse
   
    fft_sigma = fft(sigma,/double,dimension = 1)

    for i=0, modes do begin
        integrated_amp(i) = dcomplex(int_tabulated(rad(r0:r1),real_part(fft_sigma(i,r0:r1)),/double) $
                                     ,int_tabulated(rad(r0:r1),imaginary(fft_sigma(i,r0:r1)),/double))
    endfor
    
    test = abs(integrated_amp(1:modes)/integrated_amp(0))
         
    vortex_time(0,k-start) = info(7,k)/p0
    for j = 1, modes do vortex_time(j,k-start) = test(j-1)

    openw,1,filepath(strcompress('merging_mode.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 0, modes do printf,1, double(i), integrated_amp(i), format='(2(e10.4,2x))'
    close,1
    
    print, 'done '+ string(k) 
endfor

set_plot, 'ps'

; device, filename=filepath(strcompress('merging_.ps',/remove_all),root_dir='.',subdir=[location])$
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
; plot,vortex_time(0,*), alog10(vortex_time(plotmode(0),*)),xmargin=[8,4],ymargin=[3,3] $
;   ,ytitle=textoidl('log(amplitude)'),xtitle='t/orbit', yrange=yrange $
;   ,charsize=1.5,xminor=4, thick=4,xrange=xrange,xtickinterval=xtickinterval
; oplot, vortex_time(0,*), alog10(vortex_time(plotmode(1),*)), linestyle=1, thick=4
; oplot, vortex_time(0,*), alog10(vortex_time(plotmode(2),*)), linestyle=2, thick=4

device,/close

openw,1,filepath(strcompress('vortex_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,vortex_time(*,i), $
  format=strcompress("("+string(modes+1)+"(e12.5,1x))",/remove_all)
close,1

end
