function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro compare_vprofiles_vz, cases=cases, start=start, finish=finish, rslice=rslice $
             , azislice=azislice $
             ,xtickinterval=xtickinterval, ytitle=ytitle $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             , mp=mp,scale=scale, mag=mag $
             ,legend=legend, label=label, smallh=smallh, vrscale=vrscale

if not keyword_set(scale) then scale=1.0

;get the basic info

location =strcompress(cases(0),/remove_all)	
filename = strcompress('hdfaa.000',/remove_all)

fileloc  = filepath(filename,root_dir='.',subdir=[location])

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)

vz_2d = dblarr(ntheta,nphi)
vz_1d = dblarr(ntheta)

vr_2d = dblarr(ntheta,nphi)
vr_1d = dblarr(ntheta)

ks   = string(start(0),format='(I03)')
filename = strcompress('hdfaa.'+ks,/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

vrad   = hdf(3, fileloc)
vtheta   = hdf(7, fileloc)

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

if not keyword_set(azislice) then begin
   azi1 = 0
   azi2 = nphi-1
endif else begin
   azi1=fix(azislice(0)*nphi)
   azi2=fix(azislice(1)*nphi)
endelse

newnphi = azi2-azi1+1
newnrad = nrslice2 - nrslice1 + 1

vr_3d = dblarr(newnrad, ntheta, newnphi)
vz_3d = dblarr(newnrad, ntheta, newnphi)

for k=azi1, azi2 do begin
   for j=0, ntheta-1 do begin
      for i=nrslice1, nrslice2 do begin

      
         sint = sin(theta(j))
         cost = cos(theta(j))

         cs = smallh*rad(i)*sint
         cs*= rad(i)^(-1.5d0)
         
         vr = vrad(i,j,k)
         vt = vtheta(i,j,k)
         
         vr/=cs
         vt/=cs

         vr_3d(i-nrslice1,j,k-azi1) = vr*sint + vt*cost
         vz_3d(i-nrslice1,j,k-azi1) = vr*cost - sint*vt
         
         


      endfor
   endfor
endfor

for j=0, ntheta - 1 do begin
 if keyword_set(mag) then begin   
   vr_1d(j) = sqrt(mean(vr_3d(*,j,*)^2))
   vz_1d(j) = sqrt(mean(vz_3d(*,j,*)^2))
  endif else begin
   vr_1d(j) = mean((vr_3d(*,j,*)))
   vz_1d(j) = mean((vz_3d(*,j,*)))
  endelse
endfor

dataplot = vz_1d*scale
if keyword_set(vrscale) then dataplot /= mean(vr_1d)

if not keyword_set(ytitle) then ytitle = ''

zaxis = 1d0/(smallh*tan(theta))
set_plot, 'ps'
device, filename=strcompress('compare_vprofiles_vz'+string(start(0),format='(I03)')+'.ps',/remove_all) $ 
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot, zaxis, dataplot, xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle='z/H(R)', ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle=textoidl(ytitle), linestyle=0

;overplot the other cases

for k=1, n_elements(cases) - 1 do begin
   
   location =strcompress(cases(k),/remove_all)	
   fileloc  = filepath(filename,root_dir='.',subdir=[location])

   ks   = string(start(k),format='(I03)')
   filename = strcompress('hdfaa.'+ks,/remove_all)
   fileloc  = filepath(filename,root_dir='.',subdir=[location])

   vrad   = hdf(3, fileloc)
   vtheta   = hdf(7, fileloc)

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

   temp = min(abs(rplot-rslice(2*k)),nrslice1)
   temp = min(abs(rplot-rslice(2*k+1)),nrslice2)

   if not keyword_set(azislice) then begin
      azi1 = 0
      azi2 = nphi-1
   endif else begin
      azi1=fix(azislice(2*k)*nphi)
      azi2=fix(azislice(2*k+1)*nphi)
   endelse
   
   newnphi = azi2-azi1+1
   newnrad = nrslice2 - nrslice1 + 1

   vr_3d = dblarr(newnrad,  ntheta, newnphi)
   vz_3d = dblarr(newnrad,  ntheta, newnphi)

   for l=azi1, azi2 do begin
      for j=0, ntheta-1 do begin
         for i=nrslice1, nrslice2 do begin
            
            sint = sin(theta(j))
            cost = cos(theta(j))
            
            cs = smallh*rad(i)*sint
            cs*= rad(i)^(-1.5d0)

            vr = vrad(i,j,l)
            vt = vtheta(i,j,l)

            vr /=cs
            vt /=cs

            vr_3d(i-nrslice1,j,l-azi1) = vr*sint + vt*cost
            vz_3d(i-nrslice1,j,l-azi1) = vr*cost - sint*vt
         endfor
      endfor
   endfor
   
   for j=0, ntheta - 1 do begin
    if keyword_set(mag) then begin   
   vr_1d(j) = sqrt(mean(vr_3d(*,j,*)^2))
   vz_1d(j) = sqrt(mean(vz_3d(*,j,*)^2))
     endif else begin
    vr_1d(j) = mean((vr_3d(*,j,*)))
    vz_1d(j) = mean((vz_3d(*,j,*)))
   endelse
   endfor
   
   dataplot = vz_1d*scale
   if keyword_set(vrscale) then dataplot /= mean(vr_1d)
 
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
   endfor
endif

device,/close

end
