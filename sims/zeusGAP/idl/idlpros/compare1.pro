;this procedure uses aziavg (with density option) and omgrad to
;compare surface density and vorticity as a function of radius. set
;avg to take azimuthal average values. otherwise will use azimuth at
;pi relative to planet
pro compare1, loc=loc, slice=slice, avg=avg
common consts, pi, nrad,time
common array, avgsig, data1d, rad, radius
if keyword_set(avg) then begin
aziavg, type='dens', loc=loc, start=slice, finish=slice
sigma=avgsig
sigma=sigma-7.
omgrad, loc=loc, start=slice, finish=slice
vorticity=data1d
vorticity=(vorticity-0.5*rad(1:nrad-2)^(-1.5))*10.
set_plot, 'ps'
device, filename='omgsig'+strcompress(string(slice,format='(I03)'),/remove_all)+'_phi.ps'$
, xsize=8, ysize=6, xoffset=0, yoffset=0, /inches, bits_per_pixel=8,/color
plot, rad, sigma, xmargin=[6,6], ymargin=[3,3], xtitle='r',/nodata $
,xtickinterval=0.4,xminor=4,ystyle=4,charsize=1.5,thick=2,title=time+' orbits'
loadct, 40
axis,yaxis=0,/save,ytitle=textoidl('(<\omega>_\phi-\omega_0)/10^{-1}')+', blue',yrange=[min(vorticity),max(vorticity)]$
,ystyle=2,charsize=1.5
oplot,rad(1:nrad-2), vorticity, color=60,thick=2
axis,yaxis=1,/save,ytitle=textoidl('(<\Sigma>_\phi-\Sigma_0)/10^{-4}')+', black',yrange=[min(sigma),max(sigma)]$
,ystyle=2,charsize=1.5
oplot,rad,sigma,thick=2
vline, radius
device,/close
endif else begin
aziavg, type='dens', loc=loc, start=slice, finish=slice, plazi=1, opp=1
sigma=avgsig
sigma=sigma-7.
omgrad, loc=loc, start=slice, finish=slice, plazi=1, opp=1
vorticity=data1d
vorticity=(vorticity-0.5*rad(1:nrad-2)^(-1.5))*10.
set_plot, 'ps'
device, filename='omgsig'+strcompress(string(slice,format='(I03)'),/remove_all)+'.ps'$
, xsize=8, ysize=6, xoffset=0, yoffset=0, /inches, bits_per_pixel=8,/color
plot, rad, sigma, xmargin=[6,6], ymargin=[3,3], xtitle='r',/nodata $
,xtickinterval=0.4,xminor=4,ystyle=4,charsize=1.5,thick=2,title=time+' orbits'
loadct, 40
axis,yaxis=0,/save,ytitle=textoidl('(\omega-\omega_0)/10^{-1}')+', blue',yrange=[min(vorticity),max(vorticity)]$
,ystyle=2,charsize=1.5
oplot,rad(1:nrad-2), vorticity, color=60,thick=2
axis,yaxis=1,/save,ytitle=textoidl('(\Sigma-\Sigma_0)/10^{-4}')+', black',yrange=[min(sigma),max(sigma)]$
,ystyle=2,charsize=1.5
oplot,rad,sigma,thick=2
vline, radius
device,/close
endelse
end
