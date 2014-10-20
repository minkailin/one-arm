pro polar_old, type=type, loc=loc, start=start, finish=finish
pi=3.141592654
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nrad=fix(dims(6))
nsec=fix(dims(7))
info=dblarr(9,2)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
in=4
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(in+nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
rad(0:in-1)=(radtmp(0)/in)*(dindgen(in)+0.5)
for j=in, nrad+in-1 do rad(j)=(radtmp(j-in)+radtmp(j-in+1))/2.
;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
data2=dblarr(nsec,in+nrad)
loadct,3, bottom=0
for k=start, finish do begin
ks=string(k,format='(I03)')
openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2
data=alog(data)
data2(*,0:in-1)=0.
data2(*,in:nrad+in-1)=data
levels=min(data)+(max(data)-min(data))*(dindgen(32)/32.)
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress(type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=16, ysize=12
polar_contour,data2,azi,rad,/isotropic,/fill,levels=levels,title=time+' orbits',ymargin=[2.5,2.5]
colorbar, position=[0.85, 0.1, 0.9, 0.9],/vertical,/right,range=[min(data),max(data)]
device,/close
print, 'done '+ks
endfor
end
