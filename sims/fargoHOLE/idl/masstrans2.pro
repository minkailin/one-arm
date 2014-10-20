pro masstrans2, loc=loc, start=start, finish=finish, out=out, width=width

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

dm = dblarr(nsec,nrad)
for i=0,  nsec-1 do begin
    dm(i,*)  = dlogr*dphi*rmed(*)^2
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sigma =  dblarr(nsec,nrad)
vrad   = dblarr(nsec,nrad)
vrad_c = dblarr(nsec, nrad)

mdot = dblarr(nsec, nrad)
mdot_1d  = dblarr(nrad)

gap_1d = dblarr(nrad)

dx = dblarr(nsec)
data = dblarr(nsec)

gapevol = dblarr(5,finish-start+1)


openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2 
sigma0 = sigma
init_mass = total(sigma0*dm)


for l=start, finish do begin
    ks=string(l,format='(I03)')
    openr,2,filepath(strcompress('gasdens'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2 

    openr,2,filepath(strcompress('gasvrad'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vrad
    close,2 
    
    for i = 0, nsec-1 do begin
        for j = 1, nrad-2 do begin
            vrad_c(i,j) = 0.5*(vrad(i,j) + vrad(i,j+1))
        endfor
    endfor
    

    mass = total(sigma*dm)/init_mass

;co-orbital region    
    plx=info(1,l)
    ply=info(2,l)
    plrad=sqrt(plx*plx+ply*ply)
    rhill = f0*plrad
    xs = 2.5*rhill
    temp = min(abs(rmed - (plrad - xs)), x1)
    temp = min(abs(rmed - (plrad + xs)), x2)

    temp = min(abs(rmed - plrad), plgrid)

;get change in density, 1d gap profile and gap depth
    dsigma = sigma - sigma0
    dsigma /= sigma0

    for j=0, nrad-1 do gap_1d(j) = mean(dsigma(*,j))
    gapdepth = -mean(dsigma(*,x1:x2))

;     for i=0, nsec-1 do begin
;         array = dsigma(i,*)
;         test = 1d0
;         beg = plgrid
;         while (test gt 0d0) do begin
;             test = array(beg)*array(beg-1)
;             beg -= 1
;         endwhile
;         minrad = mean(rmed(beg:beg+1))

;         test = 1d0
;         beg = plgrid
;         while (test gt 0d0) do begin
;             test = array(beg)*array(beg+1)
;             beg += 1
;         endwhile
;         maxrad = mean(rmed(beg-1:beg))

;         temp = min(abs(rmed - minrad), x1)
;         temp = min(abs(rmed - maxrad), x2)

;         dx(i) = (maxrad - minrad)/(rhill*2d0)
;         data(i) = mean(dsigma(x1:x2))
;     endfor

         test = 1d0
         beg = plgrid
         while (test gt 0d0) do begin
             test = gap_1d(beg)*gap_1d(beg-1)
             beg -= 1
         endwhile
         minrad = mean(rmed(beg:beg+1))

         test = 1d0
         beg = plgrid
         while (test gt 0d0) do begin
             test = gap_1d(beg)*gap_1d(beg+1)
             beg += 1
         endwhile
         maxrad = mean(rmed(beg-1:beg))    
         gapwidth = (maxrad- minrad)/(2d0*rhill)

;get mass transport 1d, and averaged around outer gap edge

    mdot = sigma*vrad_c/sigma0
    for j=0, nrad-1 do begin
        mdot(*,j) /= rmed(j)^(-0.5) ;normalise by keplerian speed 
        mdot_1d(j) = mean(mdot(*,j))
    endfor

    temp = min(abs(rmed -  width(0)), x1)
    temp = min(abs(rmed -  width(1)), x2)

    massflux = mean(mdot(*, nrad-1)) ;mean(mdot(*, x1:x2))

    gapevol(0,l-start) = info(7,l)/p0
    gapevol(1,l-start) =  gapdepth
    gapevol(2,l-start) =  massflux
    gapevol(3,l-start) =  gapwidth
    gapevol(4,l-start) =  mass

    openw,1,filepath(strcompress('masstrans2_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 0, nrad-1 do printf,1, rmed(i), gap_1d(i), mdot_1d(i), format='(3(e22.15,2x))'
    close,1
    
    print, 'done '+ks
endfor

openw,1,filepath(strcompress('masstrans2_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,gapevol(0,i),gapevol(1,i), gapevol(2,i),  gapevol(3,i), gapevol(4,i) $
  , format='(5(e22.15,2x))'
close,1

end

