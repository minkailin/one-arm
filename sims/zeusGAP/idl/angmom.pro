function hdf, n, name

fileid = hdf_sd_start(name,/read)
sds = hdf_sd_select(fileid,n)
hdf_sd_getdata,sds,data

return, data

end

pro angmom, loc=loc, start=start, finish=finish

;get the basic info
location =strcompress(loc,/remove_all)
filename = strcompress('hdfaa.000',/remove_all)
fileloc  = filepath(filename,root_dir='.',subdir=[location])

rad  = hdf(2, fileloc)
theta= hdf(1, fileloc)
phi  = hdf(0, fileloc)


dr = rad(1)- rad(0)
dtheta = theta(1) - theta(0)
dphi   = phi(1) - phi(0)

nrad  = n_elements(rad)
ntheta= n_elements(theta)
nphi  = n_elements(phi)


if not keyword_set(finish) then finish = start

for n=start, finish do begin
    ks   = string(n,format='(I03)')
    filename = strcompress('hdfaa.'+ks,/remove_all)
    fileloc  = filepath(filename,root_dir='.',subdir=[location])

    vphi   = hdf(11, fileloc)
    density   = hdf(19, fileloc)

    momentum = 0.0
    for k=0, nphi-1 do begin
        for j=0, ntheta-1 do begin
            for i=0, nrad-1 do begin
                
                bigR = rad(i)*sin(theta(j))
                dV = dr*rad(i)*dtheta*bigR*dphi
                dm = dV*density(i,j,k)
                
                little_j = bigR*vphi(i,j,k)
                momentum += little_j*dm
            endfor
        endfor
    endfor

print, 'angular momentum is = ', momentum 

endfor
end
