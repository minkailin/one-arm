function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro torque, loc=loc, start=start, finish=finish, rdomain=rdomain

mp = 3d-4


;get the basic info
location =strcompress(loc,/remove_all)
filename = strcompress('hdfaa.000',/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)


nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)


;construct array to hold cell edge radius

re = dblarr(nrad+1)
re(0) = rdomain(0)
re(nrad) = rdomain(1)
for i=0, nrad-2 do begin
    re(i+1) = 2d0*rad(i) - re(i)
endfor


; assume uniform theta and phi spacing

dtheta = theta(1) - theta(0)
dphi = phi(1) - phi(0)


if not keyword_set(finish) then finish = start


; planet info
nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
planetinfo = dblarr(7,nlines)
openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
readf,1,planetinfo
close,1
torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)



openw,1,filepath('torque_info.dat',root_dir='.',subdir=[location])


for n=start, finish do begin
    ks   = string(n,format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])

    data   = hdf(19, fileloc)
 


     plx=planetinfo(1,n)
     ply=planetinfo(2,n)
     plz=planetinfo(3,n)
     plrad = sqrt(plx^2 + ply^2 + plz^2)
     rhill = plrad*(mp/3d0)^(1d0/3d0)
;     print, 'hill radius is', rhill

     mass = 0.0
     hillmass = 0.0
 
     tqin    = 0.0
     tqin_ex = 0.0

     tqout   = 0.0	
     tqout_ex= 0.0

     totaltq = 0.0
     totaltq_ex = 0.0

     for k=0, nphi-1 do begin
         for j=0, ntheta-1 do begin
             for i=0, nrad-1 do begin
                 dr= re(i+1) - re(i)
                 dm = data(i,j,k)*rad(i)^2*sin(theta(j))*dr*dtheta*dphi

                 xx = rad(i)*sin(theta(j))*cos(phi(k))
                 yy = rad(i)*sin(theta(j))*sin(phi(k))
                 zz = rad(i)*cos(theta(j))

                 dsq = (xx - plx)^2 + (yy - ply)^2 + (zz - plz)^2
                 dsqsoft = dsq + (0.1*rhill)^2
 
                 hillcutfactor = 1d0 - exp(-0.5*dsq/(rhill^2))                 

		 fx = (xx - plx)*dm/(dsqsoft^1.5d0)
                 fy = (yy - ply)*dm/(dsqsoft^1.5d0)	
                 
                 fx_ex = fx*hillcutfactor
		 fy_ex = fy*hillcutfactor	

                 tqcell   = plx*fy - ply*fx
                 tqcell_ex= plx*fy_ex - ply*fx_ex       
                 
                 totaltq    += tqcell
                 totaltq_ex += tqcell_ex

		 if(rad(i) lt plrad) then begin
			tqin   += tqcell
                        tqin_ex+= tqcell_ex
 		 endif
        
                 if(rad(i) gt plrad) then begin
                        tqout   += tqcell
                        tqout_ex += tqcell_ex
                 endif
             endfor
         endfor
     endfor

    print, '-----------------------------'
    print, 'torque inc rhill', totaltq
    print, 'torque exc rhill', totaltq_ex


    printf,1,planetinfo(0,n),tqin,tqout,totaltq,tqin_ex,tqout_ex,totaltq_ex, format='(1x,7(e22.15,1x))'


;     print, 'disc mass', mass
;     print, 'mass inside Rh', hillmass



;    if not keyword_set(nopert) then begin
;        data = data - data0
;        data/=data0
;    endif
    
;    if keyword_set(zslice) then begin
;        nzslice = (ntheta/2)*(zslice + 1d0) - 1d0
;        if(nzslice lt 0) then nzslice = 0    
;        nzslice = zslice*ntheta
;        data2d(0:nrad-1,0:nphi-1) = data(0:nrad-1,nzslice,0:nphi-1)
;        for i = 0, nrad - 1 do data1d(i) = mean(data2d(i,*))
;    endif
    
;    if keyword_set(vertavg) then begin
;        for i=0, nrad-1 do begin
;            for k = 0, nphi-1 do begin
;                data2d(i,k) = mean(data(i,*,k))
;            endfor
;            data1d(i) = mean(data2d(i,*))
;        endfor
;    endif
    
;    if keyword_set(log) then data2d=alog10(data2d)

;    if not keyword_set(plotrange0) then begin
;        plotrange=[min(data2d),max(data2d)]
;    endif else plotrange=plotrange0

;    if keyword_set(hole) then begin
;        for j=0, nphi-1 do begin
;            for i=0, nrad-1 do begin
;                if(rad(i) le hole) then data2d(i,j) = plotrange(1)*10
;            endfor
;        endfor 
;    endif

;    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
;    time=string(planetinfo(0,n)/torb,format='(F7.2)')

;    loadct, 5, bottom=0
;    set_plot, 'ps'
;    device, filename=filepath(strcompress('polar_'+type+ks+'.ps',/remove_all) $
;                              ,root_dir='.',subdir=[location]) $
;      ,/color, bits_per_pixel=8,xsize=14, ysize=12 
;    polar_contour,transpose(data2d),phi,rad,/isotropic,/fill,levels=levels,title=time+' orbits' $
;      ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
;      , ytickinterval=xtickinterval
;    colorbar, position=[0.85, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
; device,/close

;print, 'done '+ks
endfor

close,1

end
