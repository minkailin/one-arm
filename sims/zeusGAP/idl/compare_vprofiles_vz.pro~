function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro compare_vprofiles, cases=cases, start=start, finish=finish, rslice=rslice $
             , type=type, log=log, plotrange0=plotrange0, azislice=azislice $
             ,xtickinterval=xtickinterval, ytitle=ytitle, scale=scale $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             ,nopert=nopert, basic=basic, mp=mp $
             ,legend=legend, label=label, smallh=smallh

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

data2d = dblarr(ntheta,nphi)
data_axisymmetric = dblarr(nrad, ntheta)
data1d = dblarr(ntheta)

;plot the first case

case type of
    'vrad':   data0   = hdf(3, fileloc)
    'vtheta': data0   = hdf(7, fileloc)
    'vphi':   data0   = hdf(11, fileloc)
    'pot':    data0   = hdf(15, fileloc) 
    'dens':   data0   = hdf(19, fileloc)
endcase

if keyword_set(basic) then begin
    for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(data0(i,j,*))
        endfor
    endfor
    for k=0, nphi-1 do begin
        data0(*,*,k) = data_axisymmetric(*,*)
    endfor
endif

    ks   = string(start(0),format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])

    case type of
        'vrad':   data   = hdf(3, fileloc)
        'vtheta': data   = hdf(7, fileloc)
        'vphi':   data   = hdf(11, fileloc)
        'pot':    data   = hdf(15, fileloc)
        'dens':   data   = hdf(19, fileloc)
    endcase

    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif

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
      rhill = plrad*(mp/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
     endif else begin
      rplot = rad
     endelse

    temp = min(abs(rplot-rslice(0)),nrslice1)
    temp = min(abs(rplot-rslice(1)),nrslice2)

    for i=0, ntheta-1 do begin
       for j=0, nphi-1 do begin
         data2d(i,j) = mean(data(nrslice1:nrslice2,i,j))
       endfor
    endfor

    if not keyword_set(azislice) then begin
    for i = 0, ntheta - 1 do data1d(i) = mean(data2d(i,*))
    endif else begin
    azi1=fix(azislice(0)*nphi)
    azi2=fix(azislice(1)*nphi)
    for i= 0, ntheta-1 do data1d(i) = mean(data2d(i, azi1:azi2))
    endelse

    if keyword_set(log) then data1d = alog10(data1d)

     if not keyword_set(ytitle) then ytitle = ''

    zaxis = 1d0/(smallh*tan(theta))
    set_plot, 'ps'
    device, filename=strcompress('compare_vprofiles_'+type+string(start(0),format='(I03)')+'.ps',/remove_all) $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

    plot, zaxis, data1d, xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle='z/H(R)', ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle), linestyle=0

;overplot the other cases

    for k=1, n_elements(cases) - 1 do begin
	
     	location =strcompress(cases(k),/remove_all)
	if not keyword_set(basic) then begin
    	filename = strcompress('hdfaa.000',/remove_all)
	endif else begin
    	name = string(basic,format='(I03)')
    	filename = strcompress('hdfaa.'+name,/remove_all)
	endelse

	fileloc  = filepath(filename,root_dir='.',subdir=[location])

	case type of
    	'vrad':   data0   = hdf(3, fileloc)
    	'vtheta': data0   = hdf(7, fileloc)
    	'vphi':   data0   = hdf(11, fileloc)
    	'pot':    data0   = hdf(15, fileloc)
    	'dens':   data0   = hdf(19, fileloc)
	endcase

    	  if keyword_set(basic) then begin
    	for j=0, ntheta-1 do begin
        for i=0, nrad-1 do begin
            data_axisymmetric(i,j) = mean(data0(i,j,*))
        endfor
    	endfor
    	for k=0, nphi-1 do begin
      	  data0(*,*,k) = data_axisymmetric(*,*)
    	endfor
	endif

	ks   = string(start(k),format='(I03)')
        filename = strcompress('hdfaa.'+ks,/remove_all)
        fileloc  = filepath(filename,root_dir='.',subdir=[location])

        case type of
        'vrad':   data   = hdf(3, fileloc)
        'vtheta': data   = hdf(7, fileloc)
        'vphi':   data   = hdf(11, fileloc)
        'pot':    data   = hdf(15, fileloc)
        'dens':   data   = hdf(19, fileloc)
        endcase

    if not keyword_set(nopert) then begin
        if not keyword_set(log) then data = data - data0
        data/=data0
    endif

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
      rhill = plrad*(mp/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
     endif else begin
      rplot = rad
     endelse

     temp = min(abs(rplot-rslice(k+1)),nrslice1)
     temp = min(abs(rplot-rslice(k+2)),nrslice2)

     for i=0, ntheta-1 do begin    
       for j=0, nphi-1 do begin
         data2d(i,j) = mean(data(nrslice1:nrslice2,i,j))
       endfor
     endfor


     if not keyword_set(azislice) then begin
     for i = 0, ntheta - 1 do data1d(i) = mean(data2d(i,*))
     endif else begin
     azi1=fix(azislice(k+1)*nphi)
     azi2=fix(azislice(k+2)*nphi)
     for i= 0, ntheta-1 do data1d(i) = mean(data2d(i, azi1:azi2))
     endelse





     if keyword_set(log) then data1d = alog10(data1d)

     oplot, zaxis, data1d, thick=4, linestyle=k 
    
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

end
