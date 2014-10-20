pro pattern_torque_compare, loc=loc, modes=modes, a0=a0, plotmode=plotmode, begavg=begavg $
                            , xrange=xrange, xtickinterval=xtickinterval, xminor=xminor $
                            , legend = legend, label=label, title=title
scale = 1d4
period = 2.0*!dpi*a0^1.5
numcases = n_elements(plotmodes)
location=strcompress(loc,/remove_all)

;mode amplitude information
lines = file_lines(filepath('mode_amplitudes.dat',root_dir='.',subdir=[location]))
mode_amps = dblarr(1 + modes, lines)
openr,1,filepath('mode_amplitudes.dat',root_dir='.',subdir=[location])
readf,1,mode_amps
close,1
amp_time = mode_amps(0,*)

;get running time averages
nmax = n_elements(amp_time)
temp = min(abs(begavg - amp_time), grid)
amplitude_avg = dblarr(modes, nmax)

for j = 1, modes  do begin
    for i=grid, nmax-1 do begin
        amplitude_avg(j-1,i) = alog10(mean(mode_amps(j, grid:i)))
    endfor
endfor


;torque
lines = file_lines(filepath('tqwk0.dat',root_dir='.',subdir=[location]))
data = dblarr(10, lines)
openr,1,filepath('tqwk0.dat',root_dir='.',subdir=[location])
readf,1,data
close,1
tq_time = data(9,*)/period
tq_inc = (data(2,*) + data(1,*))*scale ; with hill
tq_exl = (data(3,*) + data(4,*))*scale ; no hill
torque = tq_inc
torque_corot = tq_inc - tq_exl


;get running time averages
temp = min(abs(begavg - tq_time), grid)
tqnmax = n_elements(tq_time)

torque_avg = dblarr(tqnmax)
torque_corot_avg = dblarr(tqnmax)

for i=grid, tqnmax-1 do begin
    torque_avg(i) = mean(torque(grid:i))
    torque_corot_avg(i) = mean(torque_corot(grid:i))
endfor


; plot mode amp and outer torque v.s. time (running time avg)
  
set_plot, 'ps'
device, filename=filepath('mode_torque_avg.ps',root_dir='.',subdir=[location]) $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot,amp_time, amplitude_avg(plotmode(0)-1,*),xmargin=[8,8],ymargin=[3,2] $
  ,ystyle=4,xtitle='t/orbits', title=textoidl(title) $
  ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval,/nodata 
axis, /save, yaxis=0, ytitle='Mode amplitude (solid, dotted)' $ ;strcompress('m='+string(plotmode(0),format='(I1)')+'(solid) and '+ $
                                ;            'm='+string(plotmode(1),format='(I1)')+'(dotted)' + $
                                ;            'amplitude') $
,charsize=1.5, yrange=[-2,-1.5];[min(amplitude_avg(plotmode(0)-1,*)), max(amplitude_avg(plotmode(0)-1,*))]
oplot, amp_time, amplitude_avg(plotmode(0)-1,*), linestyle=0, thick=4
oplot, amp_time, amplitude_avg(plotmode(1)-1,*), linestyle=1, thick=4

axis, /save, yaxis=1, ytitle='Torque (dashed, dashed-dot)',charsize=1.5, yrange=[-0.5,0.5];[min(torque_avg), max(torque_avg)]
oplot, tq_time, torque_avg, linestyle=2, thick=4
oplot, tq_time, torque_corot_avg, linestyle=3, thick=4



if keyword_set(legend) then begin
    x0=legend(0)
    x1=legend(1)
    y0=legend(2)
    dy=legend(3)
    for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
    endfor
endif
device,/close
end

;lines = file_lines(filepath('tq_corot_time.dat',root_dir='.',subdir=[location]))
;data = dblarr(5, lines)
;openr,1,filepath('tq_corot_time.dat',root_dir='.',subdir=[location])
;readf,1,data
;close,1
;tq_time = data(0,*)
;outer_tq_inc = (data(2,*) + data(4,*))*scale ;  with hill contribution
;outer_tq_exl = (data(3,*) + data(1,*))*scale ; no hill contribution
;corot_torque = outer_tq_inc - outer_tq_exl


; get torques at the same time as mode amplitude output
;outer_tq_reduced = amplitude
;outer_tq_avg_reduced = amplitude_avg
;for i = 0, nmax-1 do begin
;    temp =  min(abs(amp_time(i) - tq_time), grid)
;    outer_tq_reduced(i) = outer_tq(grid)
;    outer_tq_avg_reduced(i) = outer_tq_avg(grid)
;endfor

; plot mode amp and outer torque v.s. time (instant)
  
; set_plot, 'ps'
; device, filename=filepath('mode_torque.ps',root_dir='.',subdir=[location]) $
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
; plot,amp_time, alog10(amplitude),xmargin=[8,8],ymargin=[3,1] $
;   ,ystyle=4,xtitle='t/orbits' $
;   ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;   ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval,/nodata 
; axis, /save, yaxis=0, ytitle=strcompress('m='+string(plotmode,format='(I1)')+' amplitude (solid)') $
;   ,charsize=1.5, yrange=alog10([min(amplitude), max(amplitude)])
; oplot, amp_time, alog10(amplitude), linestyle=0, thick=4
; axis, /save, yaxis=1, ytitle='Outer torque (dotted, dashed)',charsize=1.5, yrange=[min(outer_tq_inc), max(outer_tq_inc)]
; oplot, tq_time, outer_tq_inc, linestyle=1, thick=4
; oplot, tq_time, outer_tq_exl, linestyle=2, thick=4
; device,/close


 ;plot mode amp v.s outer torque (instant)
  
;set_plot, 'ps'
;device, filename=filepath('mode_torque_correl.ps',root_dir='.',subdir=[location]) $
;  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;plot, alog10(amplitude), outer_tq_reduced , xmargin=[8,2],ymargin=[3,1] $
;  ,ystyle=2,ytitle='Outer torque', xtitle =  strcompress('m='+string(plotmode,format='(I1)')+' amplitude') $
;  ,charsize=1.5, psym=2, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;  ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval
;device,/close

; plot mode amp v.s outer torque (running time average)
  
;set_plot, 'ps'
;device, filename=filepath('mode_torque_correl_avg.ps',root_dir='.',subdir=[location]) $
;  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;plot,   alog10(amplitude_avg), outer_tq_avg_reduced, xmargin=[8,2],ymargin=[3,1] $
;  ,ystyle=2, psym=2, ytitle='Outer torque', xtitle =  strcompress('m='+string(plotmode,format='(I1)')+' amplitude') $
;  ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;  ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval
;device,/close
