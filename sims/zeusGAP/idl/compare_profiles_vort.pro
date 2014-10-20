function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro compare_profiles_vort, cases=cases, start=start, finish=finish, zslice=zslice $
             , type=type, log=log, plotrange0=plotrange0 $
             ,xtickinterval=xtickinterval, ytitle=ytitle, scale=scale $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             , mp=mp $
             ,legend=legend, label=label

;get the basic info

location =strcompress(cases(0),/remove_all)	
if not keyword_set(basic) then begin
    filename = strcompress('hdfaa.000',/remove_all)
endif else begin
    name = string(basic,format='(I03)')
    filename = strcompress('hdfaa.'+name,/remove_all)
endelse

fileloc  = filepath(filename,root_dir='.',subdir=[location])

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)

data2d = dblarr(nrad, nphi)
data_axisymmetric = dblarr(nrad, ntheta)

den1d = dblarr(nrad)
vphi1d= dblarr(nrad)
data1d = dblarr(nrad)

    ks   = string(start(0),format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])

    vphi    = hdf(11, fileloc)
    den   = hdf(19, fileloc)

    nzslice = fix(zslice*ntheta)
    if(nzslice eq ntheta) then nzslice -= 1
    for j=1, 2 do begin
    if(j eq 1) then data = den
    if(j eq 2) then data = vphi
    data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
    if(j eq 1) then begin
    for i = 0, nrad - 1 do den1d(i) = mean(data2d(i,*))
    endif else begin
    for i = 0, nrad - 1 do vphi1d(i) = mean(data2d(i,*))
    endelse   
    endfor

    if keyword_set(log) then data1d = alog10(data1d)

    ksq = deriv(rad, rad^2*vphi1d^2)/rad^3
    data1d = ksq/(2d0*(vphi1d/rad)*den1d)

    if keyword_set(mp) then begin
    ; planet info
      nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
      planetinfo = dblarr(7,nlines)
      openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
      readf,1,planetinfo
      close,1

      plx=planetinfo(1,start(0))
      ply=planetinfo(2,start(0))
      plz=planetinfo(3,start(0))
      plrad = sqrt(plx^2 + ply^2 + plz^2)
      rhill = plrad*(mp(0)/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
      xtitle = textoidl('(r - r_p)/r_h')
     endif else begin
      rplot = rad
      xtitle = textoidl('r')
     endelse


     if not keyword_set(ytitle) then ytitle = ''

    set_plot, 'ps'
    device, filename=strcompress('compare_profiles_vort'+string(start(0),format='(I03)')+'.ps',/remove_all) $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

    plot, rplot, data1d, xmargin=[8,1],ymargin=[3.5,0.5], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle=xtitle, ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle), linestyle=0

;overplot the other cases

    for k=1, n_elements(cases) - 1 do begin
	
     	location =strcompress(cases(k),/remove_all)

	ks   = string(start(k),format='(I03)')
        filename = strcompress('hdfaa.'+ks,/remove_all)
        fileloc  = filepath(filename,root_dir='.',subdir=[location])

       phi    = hdf(11, fileloc)
       den   = hdf(19, fileloc)

     for j=1, 2 do begin
    if(j eq 1) then data = den
    if(j eq 2) then data = vphi
    data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
    if(j eq 1) then  begin
    for i = 0, nrad - 1 do den1d(i) = mean(data2d(i,*))
    endif else begin
    for i = 0, nrad - 1 do vphi1d(i) = mean(data2d(i,*))
    endelse
    endfor

    if keyword_set(log) then data1d = alog10(data1d)

    ksq = deriv(rad, rad^2*vphi1d^2)/rad^3
    data1d = ksq/(2d0*(vphi1d/rad)*den1d)


    if keyword_set(mp) then begin
    ; planet info
      nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
      planetinfo = dblarr(7,nlines)
      openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
      readf,1,planetinfo
      close,1

      plx=planetinfo(1,start(k))
      ply=planetinfo(2,start(k))
      plz=planetinfo(3,start(k))
      plrad = sqrt(plx^2 + ply^2 + plz^2)
      rhill = plrad*(mp(k)/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
     endif else begin
      rplot = rad
     endelse

     oplot, rplot, data1d, thick=4, linestyle=k 
    
    endfor


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
;stop
end
