function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

function aziavg, data3d, nzslice
  nr = n_elements(data3d(*,0,0))
  data1d = dblarr(nr)

  
  for i = 0, nr - 1 do begin
     data1d(i) = mean(data3d(i,nzslice,*))
  endfor
return, data1d
end

function bigH, bigR, rref, hsmall, qsmall

height = hsmall*rref*(rref/bigR)^(0.5d0*(qsmall-3d0))

return, height
end



pro compare_toomreq3d, cases=cases, start=start, zslice=zslice $
             ,xtickinterval=xtickinterval $
             ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval $
             , mp=mp, smallh=smallh, smallq=smallq, gmma=gmma  $
             ,legend=legend, label=label, r0=r0

;get the basic info

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


vphi1d = dblarr(nrad)
den1d  = dblarr(nrad)
cs1d   = dblarr(nrad)
bigH1d = dblarr(nrad)


vphi  = hdf(11, fileloc)
den   = hdf(19, fileloc)
if(gmma gt 1.0) then begin
   pres  = (gmma-1d0)*hdf(23, fileloc) 
endif else begin
   pres =dblarr(nrad, ntheta, nphi)
   for j=0, ntheta-1 do begin
   for i=0, nrad-1 do begin
   bigR = rad(i)*sin(theta(j))
   omegak = bigR^(-3d0/2d0)
   csq =( bigH(bigR, r0, smallh, smallq)*omegak )^2
   pres(i,j,*) = csq*den(i,j,*)
   endfor
   endfor
endelse

nzslice = fix(zslice*ntheta)
if(nzslice eq ntheta) then nzslice -= 1

vphi1d = aziavg(vphi, nzslice)
den1d  = aziavg(den, nzslice)
cs1d   = aziavg(sqrt(gmma*pres/den), nzslice)

ksq = deriv(rad, rad^2.*vphi1d^2.)/rad^3.
kappa = sqrt(ksq)

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

if keyword_set(mp) then begin
   rhill = plrad*(mp(0)/3d0)^(1d0/3d0)
   rplot = (rad - plrad)/rhill
   xtitle = textoidl('(r - r_p)/r_h')
endif else begin
   rplot = rad/plrad
   xtitle = textoidl('r/r_0')
endelse

bigH1d = smallh*r0*(r0/rad)^(0.5d0*(smallq-3d0))

toomreQ = cs1d*kappa/(!dpi*2d0*den1d*bigH1d)



set_plot, 'ps'
device, filename=strcompress('compare_toomreq3d_'+string(start(0),format='(I03)')+'.ps',/remove_all) $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot, rplot, toomreq, xmargin=[8,1],ymargin=[3.5,0.5], ystyle=0  $
      ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle=xtitle, ytickinterval=ytickinterval  $
      ,xtickinterval=xtickinterval, ytitle='Q', linestyle=0

;overplot the other cases

for k=1, n_elements(cases) - 1 do begin
   
   location =strcompress(cases(k),/remove_all)
   
   ks   = string(start(k),format='(I03)')
   filename = strcompress('hdfaa.'+ks,/remove_all)
   fileloc  = filepath(filename,root_dir='.',subdir=[location])

        
   vphi  = hdf(11, fileloc)
   den   = hdf(19, fileloc)
   if(gmma gt 1.0) then begin
   pres  = (gmma-1d0)*hdf(23, fileloc)
endif else begin
   for j=0, ntheta-1 do begin
   for i=0, nrad-1 do begin
   bigR = rad(i)*sin(theta(j))
   omegak = bigR^(-3d0/2d0)
   csq =( bigH(bigR, r0, smallh, smallq)*omegak )^2
   pres(i,j,*) = csq*den(i,j,*)
   endfor
   endfor
endelse



   nzslice = fix(zslice*ntheta)
   vphi1d = aziavg(vphi, nzslice)
   den1d  = aziavg(den, nzslice)
   cs1d   = aziavg(sqrt(gmma*pres/den), nzslice)

   ksq = deriv(rad, rad^2.*vphi1d^2.)/rad^3.
   kappa = sqrt(ksq)

   toomreQ = cs1d*kappa/(!dpi*2d0*den1d*bigH1d)

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
   if keyword_set(mp) then begin
      rhill = plrad*(mp(k)/3d0)^(1d0/3d0)
      rplot = (rad - plrad)/rhill
   endif else begin
      rplot = rad/plrad
   endelse
   
   oplot, rplot, toomreq, thick=4, linestyle=k 
   
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
