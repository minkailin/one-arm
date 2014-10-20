function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro torquexy, loc=loc, start=start, finish=finish  $
            , mp=mp, xrange=xrange, yrange=yrange $
            , plotrange0=plotrange0, scale=scale, size=size $
              , log=log
if not keyword_set(scale) then scale=1d0   


;get the basic info
location =strcompress(loc,/remove_all)
filename = strcompress('hdfaa.000',/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

x  = hdf(2, fileloc)
y  = hdf(1, fileloc)
z  = hdf(0, fileloc)


nx  = n_elements(x)
ny  = n_elements(y)
nz  = n_elements(z)

; assume UNIFORM SPACING

dx = x(1) - x(0)
dy = y(1) - y(0)
dz = z(1) - z(0)

; data output

tq3d = dblarr(nx, ny, nz)
tq3d_ex = dblarr(nx, ny, nz)

tq2d = dblarr(nx, ny)
tq2d_ex = dblarr(nx, ny)

if not keyword_set(finish) then finish = start

; planet info
nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
planetinfo = dblarr(7,nlines)
openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
readf,1,planetinfo
close,1
torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)

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

     for k=0, nz-1 do begin

         for j=0, ny-1 do begin
             for i=0, nx-1 do begin
  
               dm = data(i,j,k)*dx*dy*dz
                 xx = x(i)
                 yy = y(j)
                 zz = z(k)

                 dsq = (xx - plx)^2 + (yy - ply)^2 + (zz - plz)^2
                 dsqsoft = dsq + (0.1*rhill)^2
 
                 hillcutfactor = 1d0 - exp(-0.5*dsq/(rhill^2))                 

		 fx = (xx - plx)*dm/(dsqsoft^1.5d0)
                 fy = (yy - ply)*dm/(dsqsoft^1.5d0)	
                 
                 fx_ex = fx*hillcutfactor
		 fy_ex = fy*hillcutfactor	

                 tq3d(i,j,k)   = plx*fy - ply*fx
                 tq3d_ex(i,j,k)= plx*fy_ex - ply*fx_ex     
               endfor
         endfor

     endfor

    print, '-----------------------------'
    print, 'torque inc rhill', total(tq3d)
    print, 'torque exc rhill', total(tq3d_ex)
    
    for j=0, ny-1 do begin
        for i=0, nx-1 do begin
            tq2d(i,j)   =  total(tq3d(i,j,0:nz-1))
            tq2d_ex(i,j)=  total(tq3d_ex(i,j,0:nz-1))
        endfor
    endfor

    tq2d    *= scale
    tq2d_ex *= scale

    if keyword_set(size) then begin
;        filter = where(tq2d le 0d0)
;        tq2d(filter) = 1d-30
        
;        filter = where(tq2d_ex le 0d0)
;        tq2d_ex(filter) = 1d-30

        tq2d = abs(tq2d)
        tq2d_ex = abs(tq2d_ex)

        if keyword_set(log) then begin
            tq2d = alog10(tq2d)
            tq2d_ex = alog10(tq2d_ex)
        endif
    endif


    for count=0, 1 do begin
        if(count eq 0) then begin
            filename = 'torquexy_'
            data2d = tq2d
        endif
        if(count eq 1) then begin
            filename = 'torquexy_ex_'
            data2d = tq2d_ex
        endif

        if not keyword_set(plotrange0) then begin
            if not keyword_set(xrange) then begin
                plotrange=[min(data2d),max(data2d)]
            endif else begin
                temp = min(abs(x-xrange(0)),x1)
                temp = min(abs(x-xrange(1)),x2)
                temp = min(abs(y-yrange(0)),y1)
                temp = min(abs(y-yrange(1)),y2)
                plotrange=[min(data2d(x1:x2, y1:y2)),max(data2d(x1:x2,y1:y2))]
            endelse
        endif else plotrange=plotrange0
        
        levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
        time=string(planetinfo(0,n)/torb,format='(F7.2)')
        
        loadct, 5, bottom=0
        set_plot, 'ps'
        device, filename=filepath(strcompress(filename+ks+'.ps',/remove_all) $
                                  ,root_dir='.',subdir=[location]) $
          ,/color, bits_per_pixel=8,xsize=14, ysize=12
        contour,data2d,x,y,/isotropic,/fill,levels=levels,title=time+' orbits' $
          ,ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
          , ytickinterval=xtickinterval
        colorbar, position=[0.85, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
        device,/close
    endfor
print, 'done', ks
endfor

end
