PRO orbit, start=start, finish=finish, ntime=ntime
;;;;;;;;;;;;;;;;;
;SURFACE DENSITY;
;;;;;;;;;;;;;;;;;
; sigma=dblarr(finish-start+1)
; openr,lun,'sigmas.dat', /get_lun
; readf, lun, sigma, format='(d5.1)'
; close,lun
;;;;;;;;;;;;;;;;;;
;SET UP TIME AXIS;
;;;;;;;;;;;;;;;;;;
;data=filepath('dims.dat',root_dir='.',subdir=['out1'])
;dims=(read_ascii(data)).(0)
;time=findgen(dims(0,6))
;time=(time+1)*(2.*3.141592654)
;ntime=n_elements(time)
;;;;;;;;;;;;;;;;;
;GET PLANET DATA;
;;;;;;;;;;;;;;;;;
orbs=dblarr(2,ntime,finish-start+1)
for i=start, finish do begin
    num=string(i)
    temp=READ_ASCII(filepath('orbit0.dat',root_dir='.',subdir=[strcompress('out'+num,/remove_all)]))
    temp=temp.(0)
    a=n_elements(temp(0,*))/ntime
    a=fix(a)
     for j=0, ntime-1 do begin
     orbs(0,j,i-start)=temp(0,j*a)
      orbs(1,j,i-start)=temp(2,j*a)
 endfor
endfor
;;;;;;;;;;
;PLOTTING;
;;;;;;;;;;
set_plot, 'ps'
device, filename='satacc.ps', xsize=8, ysize=6, xoffset=0, yoffset=0, /inches
plot, orbs(0,*,0)/(2.*3.141592654), orbs(1,*,0), xmargin=[14,6], ymargin=[6,6], xtitle='Orbits', ytitle='a', yrange=[0.4, 1.0]
for k=1, finish-start do begin
oplot, orbs(0,*,k)/(2.*3.141592654), orbs(1,*,k)
endfor
device,/close
END
