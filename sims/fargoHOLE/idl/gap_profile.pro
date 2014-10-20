pro gap_profile, loc=loc, start=start, finish=finish , out=out, mp=mp, xs=xs

f0 = (mp/3d0)^(1d0/3d0)

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
if keyword_set(out) then nout = out
info=dblarr(11,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*!dpi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR RADIUS AND AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data  = dblarr(nsec,nrad)
avgsig =dblarr(nrad)
avgsig0=dblarr(nrad)
;;;;;;;;;;;;
;SETUP DONE;
;;;;;;;;;;;;
openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2
for i=0, nrad-1 do begin
    avgsig0(i) = mean(data(*, i))
endfor


openw,10, filepath('gap_profile.dat',root_dir='.',subdir=[location])

for k=start, finish do begin
    ks=string(k,format='(I03)')
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,data
    close,2

    for i=0, nrad-1 do begin
        avgsig(i) = mean(data(*, i))
    endfor

    dsigma = (avgsig - avgsig0)/avgsig0

    plx=info(1,k)
    ply=info(2,k)
    plrad = sqrt(plx*plx + ply*ply)
    rhill = f0*plrad
    rplot = (rad - plrad)/rhill

    temp = min(abs(rplot), rp) ; planet grid

    ;find inner gap width
    for i = rp, 1, -1 do begin
        minus = dsigma(i)
        plus  = dsigma(i - 1)
        if(plus*minus lt 0d0) then begin
            in_gap_width =( 0.5*(rad(i) + rad(i-1)) - plrad )/rhill
            rpmxs = i
            break
        endif
    endfor
    ; find inner gap depth
    in_gap_depth = mean(dsigma(rpmxs:rp))

 

    ; find outer gap width
    for i = rp, nrad-2 do begin
        minus = dsigma(i)
        plus  = dsigma(i + 1)
        if(plus*minus lt 0d0) then begin
            out_gap_width =( 0.5*(rad(i) + rad(i+1)) - plrad )/rhill
            rppxs = i
            break
        endif
    endfor
    ;find outer gap width
    out_gap_depth= mean(dsigma(rp:rppxs))


    printf,10, info(7,k)/p0, in_gap_depth, in_gap_width, out_gap_depth, out_gap_width, format='(5(e22.15,x))'

print, 'done', k
endfor
close,10
end



;    temp = min(abs(rplot - xs), rppxs)
;   temp = min(abs(rplot + xs), rpmxs)
