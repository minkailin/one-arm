function run_time_average, input, start=start

wind = 10

output = input
num = n_elements(input(0,*))

temp = min(abs(start - input(0,*)), t1)

for j=t1, num-1 do begin
   for i = 1, 6 do begin
      output(i,j) = mean(input(i,t1:j));mean(input(i,j-wind:min([j+wind, num-1])))
   endfor
endfor

return, output
end


pro torque3, cases=cases, r0=r0, scale=scale $
             , xtickinterval=xtickinterval, ytickinterval=ytickinterval $
             , xrange = xrange, yrange = yrange, legend = legend, label=label $
             , yyrange=yyrange, twodim=twodim, tavg=tavg, custom_margin=custom_margin

ncase = n_elements(cases)

torb = 2d0*!dpi*r0^(3d0/2d0)
if not keyword_set(scale) then scale = 1d0

location =strcompress(cases(0),/remove_all)
nlines = file_lines(filepath('dptorque.dat',root_dir='.',subdir=[location]))
tq = dblarr(8,nlines)
openr, 1, filepath('dptorque.dat',root_dir='.',subdir=[location])
readf,1,tq
close,1

time = tq(0,*)/torb
tq(0,*) = time 

if keyword_set(tavg) then tq = run_time_average(tq, start=tavg)

tq_in = tq(1,*)*scale
tq_out= tq(2,*)*scale
tq_tot= tq(3,*)*scale
tq_in_ex  = tq(4,*)*scale
tq_out_ex = tq(5,*)*scale
tq_tot_ex = tq(6,*)*scale

;loadct, 5

set_plot, 'ps'
device, filename='torque3_tqex.ps' $
        ,/color, bits_per_pixel=8,xsize=8, ysize=4,xoffset=0,yoffset=0,/inches

plot, time, tq_tot_ex, thick=4 ,ymargin=[3.5,0.5],xmargin=[8,2],xrange=xrange $
      ,yrange=yrange, xtickinterval=xtickinterval, charsize=1.5 $
      , ytickinterval=ytickinterval, ytitle='Torque (tapered)', xtitle=textoidl('t/P_0')
;if keyword_set(tavg) then oplot, time, tqavg(6,*)*scale, thick=8, linestyle=0

for j=1, ncase-1 do begin
   
   location =strcompress(cases(j),/remove_all)
   nlines = file_lines(filepath('dptorque.dat',root_dir='.',subdir=[location]))
   tq = dblarr(8,nlines)
   openr, 1, filepath('dptorque.dat',root_dir='.',subdir=[location])
   readf,1,tq
   close,1

   time = tq(0,*)/torb
   tq(0,*) = time 

   if keyword_set(tavg) then tq = run_time_average(tq, start=tavg)

   tq_in = tq(1,*)*scale
   tq_out= tq(2,*)*scale
   tq_tot= tq(3,*)*scale
   tq_in_ex  = tq(4,*)*scale
   tq_out_ex = tq(5,*)*scale
   tq_tot_ex = tq(6,*)*scale

   oplot, time, tq_tot_ex, thick=4, linestyle=j
;   if keyword_set(tavg) then oplot, time, tqavg(6,*)*scale, thick=8, linestyle=j
endfor


if keyword_set(legend) then begin
    for j=0, ncase-1 do begin
       
       x0 = legend(0)
       x1 = legend(1)
       y0 = legend(2)
       dy = legend(3)
       
        oplot, [x0, x1], [y0, y0]- j*dy, thick=4, linestyle=j
        xyouts, x1, y0-j*dy, textoidl(label(j)), charsize=1.5
    endfor
endif
device,/close



end

