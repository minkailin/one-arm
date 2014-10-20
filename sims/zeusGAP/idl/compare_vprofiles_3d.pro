function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro compare_vprofiles_3d, cases=cases, start=start, rslice=rslice $
             ,xtickinterval=xtickinterval, ytitle=ytitle $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             , mp=mp,scale=scale, mmode=mmode $
             ,legend=legend, label=label, smallh=smallh 

if not keyword_set(scale) then scale=1.0

location =strcompress(cases(0),/remove_all)	
ks   = string(start(0),format='(I03)')
filename = strcompress('hdfaa.'+ks,/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)

vrad     = hdf(3, fileloc)
vtheta   = hdf(7, fileloc)
den      = hdf(19, fileloc)

vrad   *= den ;work with momentum densities
vtheta *= den

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
   xtitle = textoidl('(r-r_p)/r_h')
endif else begin
   rplot = rad
   xtitle = 'r'
endelse

temp = min(abs(rplot-rslice(0)),nrslice1)
temp = min(abs(rplot-rslice(1)),nrslice2)
newnrad = nrslice2 - nrslice1 + 1

vz_3d = dblarr(newnrad, ntheta, nphi)
vr_3d = dblarr(newnrad, ntheta, nphi)

vz_m = dcomplexarr(newnrad, ntheta)
vr_m = dcomplexarr(newnrad, ntheta)

avg_vz = dblarr(ntheta)
avg_vr = dblarr(ntheta)

for k=0, nphi-1 do begin
   for j=0, ntheta-1 do begin
      for i=nrslice1, nrslice2 do begin
      
         sint = sin(theta(j))
         cost = cos(theta(j))

         vr = vrad(i,j,k) 
         vt = vtheta(i,j,k) 

         vz_3d(i-nrslice1,j,k) = (vr*cost - sint*vt)
         vr_3d(i-nrslice1,j,k) = (vr*sint + vt*cost)

      endfor
   endfor
endfor

;FFT. pick out azimuthal mode mmode

result = fft(vz_3d, -1, dimension=3, /double)
vz_m(*,*) = result(*,*, mmode(0)) 

result = fft(vr_3d, -1, dimension=3, /double)
vr_m(*,*) = result(*,*, mmode(0)) 

;average over radius
for j=0,ntheta-1 do begin
   avg_vz(j) = mean(abs(vz_m(*,j))^2)
   avg_vr(j) = mean(abs(vr_m(*,j))^2)
endfor

dataplot = scale*sqrt( avg_vz /(avg_vz + avg_vr) )

if not keyword_set(ytitle) then ytitle = ''

zaxis = 1d0/(smallh*tan(theta))

set_plot, 'ps'
device, filename=strcompress('compare_vprofiles_3d'+string(start(0),format='(I03)')+'.ps',/remove_all) $ 
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot, zaxis, dataplot, xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle='z/H(R)', ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle), linestyle=0

;overplot the other cases

for k=1, n_elements(cases) - 1 do begin
   location =strcompress(cases(k),/remove_all)	
   ks   = string(start(k),format='(I03)')
   filename = strcompress('hdfaa.'+ks,/remove_all)
   fileloc  = filepath(filename,root_dir='.',subdir=[location])

   vrad     = hdf(3, fileloc)
   vtheta   = hdf(7, fileloc)
   den      = hdf(19, fileloc)

   vrad   *= den                ;work with momentum densities
   vtheta *= den

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
      xtitle = textoidl('(r-r_p)/r_h')
   endif else begin
      rplot = rad
      xtitle = 'r'
   endelse
 
   temp = min(abs(rplot-rslice(2*k)),nrslice1)
   temp = min(abs(rplot-rslice(2*k+1)),nrslice2)
 
   newnrad = nrslice2 - nrslice1 + 1
   
   vz_3d = dblarr(newnrad, ntheta, nphi)
   vr_3d = dblarr(newnrad, ntheta, nphi)
   
   vz_m = dcomplexarr(newnrad, ntheta)
   vr_m = dcomplexarr(newnrad, ntheta)

   
   for l=0, nphi-1 do begin
      for j=0, ntheta-1 do begin
         for i=nrslice1, nrslice2 do begin
            
            sint = sin(theta(j))
            cost = cos(theta(j))
            
            vr = vrad(i,j,l) 
            vt = vtheta(i,j,l) 
            
            vz_3d(i-nrslice1,j,l) =( vr*cost - sint*vt)
            vr_3d(i-nrslice1,j,l) =( vr*sint + vt*cost)
            
         endfor
      endfor
   endfor
   
;FFT. pick out azimuthal mode mmode

   result = fft(vz_3d, -1, dimension=3, /double)
   vz_m(*,*) = result(*,*, mmode(k)) 

   result = fft(vr_3d, -1, dimension=3, /double)
   vr_m(*,*) = result(*,*, mmode(k)) 

;average over radius
   for j=0,ntheta-1 do begin
      avg_vz(j) = mean(abs(vz_m(*,j))^2)
      avg_vr(j) = mean(abs(vr_m(*,j))^2)
   endfor
   dataplot = scale*sqrt( avg_vz /(avg_vz + avg_vr) )
   
   oplot, zaxis, dataplot, thick=4, linestyle=k     
endfor


if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
      oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
      xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
;      xyouts, x1, y0-dy*j,label(j),charsize=1.5
   endfor
endif

device,/close

end

