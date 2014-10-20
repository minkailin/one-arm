function epicyclesq, rad, vtheta
  kappasq = deriv(rad, rad*rad*vtheta*vtheta)/(rad^3.0)
  return, kappasq
end

function vortensity, rad, vtheta, sigma
  vortensity = epicyclesq(rad, vtheta)/(2.0*vtheta/rad)
  return, vortensity/sigma
end

function entropy, pressure, sigma, gmma=gmma
  return, pressure/sigma^gmma
end

function gen_vortensity, rad, vtheta, sigma, pressure, gmma=gmma
  ent = entropy(pressure, sigma, gmma=gmma)
  vorten = vortensity(rad, vtheta, sigma)
  return, vorten*ent^(-2.0/gmma)
end

function soundspeed, pressure, sigma
  return, sqrt(pressure/sigma)
end

function scaleheight, rad, pressure, sigma
  omega_k = rad^(-1.5)
  cs = soundspeed(pressure, sigma)
  return, cs/omega_k
end

function toomreQ, rad, vtheta, sigma, pressure
  cs    = soundspeed(pressure, sigma)
  kappa = sqrt(epicyclesq(rad,vtheta))
  return, cs*kappa/(!dpi*sigma)
end

function alpha, rad, vrad, vtheta, sigma, pressure
cs    = soundspeed(pressure, sigma)
bigH  = scaleheight(rad,pressure,sigma)

omega = vtheta/rad 
domega= deriv(rad,omega) 
dang= deriv(rad, omega*rad*rad)

nrad = n_elements(rad)
visc = dblarr(nrad) 

visc(0:1) = 0.0

for i=2, nrad-1 do begin
visc(i) = int_tabulated(rad(0:i), rad(0:i)*sigma(0:i)*vrad(0:i)*dang(0:i))
endfor

visc /= sigma*rad^3*domega

return, visc/(cs*bigH)
end


function get_grid, loc, start, r0 , hnorm=hnorm, smallh=smallh
  
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[loc]))).(0)
  nout=fix(dims(5))
  nrad=fix(dims(6))
  nsec=fix(dims(7))
  
  radtmp=dblarr(nrad+1)
  rad   =dblarr(nrad)
  rplot =dblarr(nrad)
  
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[loc]))
  info=dblarr(11,nlines)
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[loc])
  readf,3,info
  close,3
 
  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[loc])
  readf,1,radtmp
  close,1
  rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.
  
  if keyword_set(hnorm) then begin 
     data  =dblarr(nsec,nrad)
     sigma =dblarr(nrad)
     pressure  = dblarr(nrad)
     
     ;density field
     openr,3,filepath(strcompress('gasdens'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
     readu,3,data
     close,3
     for i=0, nrad-1 do sigma(i) = mean(data(*,i))
     
     ;pressure field 
     if(gmma gt 0.0) then begin        
        openr,3,filepath(strcompress('gasTemperature'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
        readu,3,data
        close,3
        for i=0, nrad-1 do pressure(i) = mean(data(*,i))*sigma(i)
     endif else begin
        pressure = (smallh*smallh/rad)*sigma
     endelse
     
     temp = min(abs(rad - r0), r1)
     bigH = scaleheight(rad, pressure, sigma)
     
     rplot(r1-2)= (rad(r1-2) - rad(r1))/mean(bigH(r1-2:r1))
     rplot(r1-1)= (rad(r1-1) - rad(r1))/mean(bigH(r1-1:r1))
     rplot(r1)  = 0.0 
     rplot(r1+1)= (rad(r1+1) - rad(r1))/mean(bigH(r1:r1+1))
     rplot(r1+2)= (rad(r1+2) - rad(r1))/mean(bigH(r1:r1+2))
     
     for i=0, r1-3 do begin
        rplot(i) = -int_tabulated(rad(i:r1), 1.0/bigH(i:r1))
     endfor
     for i=r1+3, nrad-1 do begin
        rplot(i) =  int_tabulated(rad(r1:i), 1.0/bigH(r1:i))
     endfor
     
  endif else begin
     rplot = rad/r0
  endelse 
  
  return, rplot 
end

function get_data, loc, start, type, gmma, smallh=smallh
  
  common share, refrad, pert, diskmass 

  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[loc]))).(0)
  nout=fix(dims(5))
  nrad=fix(dims(6))
  nsec=fix(dims(7))
    
  radtmp=dblarr(nrad+1)
  rad   =dblarr(nrad)
  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[loc])
  readf,1,radtmp
  close,1
  rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.

  data  =dblarr(nsec,nrad)
  
  sigma0=dblarr(nrad)
  sigma =dblarr(nrad)
  
  vtheta0=dblarr(nrad)
  vtheta = dblarr(nrad)
  
  vrad0    = dblarr(nrad) 
  vrad     = dblarr(nrad) 
  
  pressure0 = dblarr(nrad)
  pressure  = dblarr(nrad)

  angmom0 = dblarr(nrad)

  ;density field
  openr,2,filepath(strcompress('gasdens'+'0'+'.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,2,data
  close,2
  for i=0, nrad-1 do sigma0(i) = mean(data(*,i))

  temp = min(abs(rad - refrad), grid)
  den0 = sigma0(grid) 

  openr,3,filepath(strcompress('gasdens'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,3,data
  close,3
  for i=0, nrad-1 do sigma(i) = mean(data(*,i))
  diskmass = int_tabulated(rad, 2d0*!dpi*rad*sigma)

  ;radial vel 
  openr,3,filepath(strcompress('gasvrad'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,3,data
  close,3
  for i=0, nrad-1 do vrad(i) = mean(data(*,i))
  
  ;azi vel
  openr,3,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,3,data
  close,3    
  for i=0, nrad-1 do vtheta0(i) = mean(data(*,i))
 
  openr,3,filepath(strcompress('gasvtheta'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,3,data
  close,3
  for i=0, nrad-1 do vtheta(i) = mean(data(*,i))
 
  ;pressure field 
  if(gmma gt 0.0) then begin
     openr,3,filepath(strcompress('gasTemperature0.dat',/remove_all),root_dir='.',subdir=[loc])
     readu,3,data
     close,3
     for i=0, nrad-1 do pressure0(i) = mean(data(*,i))*sigma0(i)
     
     openr,3,filepath(strcompress('gasTemperature'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
     readu,3,data
     close,3
     for i=0, nrad-1 do pressure(i) = mean(data(*,i))*sigma(i)
  endif else begin
     pressure0 = (smallh*smallh/rad)*sigma0
     pressure = (smallh*smallh/rad)*sigma
  endelse
  
  case type of 
     'Q': begin
        data1d = smallh*rad^(-2d0)/!dpi/sigma;toomreQ(rad, vtheta, sigma, pressure)
        if(gmma gt 0.0) then begin  
           data1d *= sqrt(gmma) ; factor of sqrt(gamma) to convert cs to adiabatic sound speed 
        endif
     end
     'dens': begin
        if (pert gt 0) then begin 
        data1d = (sigma - sigma0)/sigma0
        endif else begin
        data1d = sigma/den0
        endelse
     end
     'ent': begin
        data0 =entropy(pressure0, sigma0, gmma=gmma)
        data1d = entropy(pressure, sigma, gmma=gmma)/data0 - 1.0
     end
     'vort': begin
        data0 =  vortensity(rad, vtheta0, sigma0)
        if (pert gt 0) then begin 
        data1d = vortensity(rad, vtheta, sigma)/data0 - 1.0 
        endif else begin
        data1d = vortensity(rad, vtheta, sigma)/data0(grid)
        endelse 
     end
     'gvort': begin
        data0 = gen_vortensity(rad, vtheta0, sigma0, pressure0, gmma=gmma)
        if (pert gt 0) then begin
         data1d = gen_vortensity(rad, vtheta, sigma, pressure, gmma=gmma)/data0 - 1.0
        endif else begin
         data1d = gen_vortensity(rad, vtheta, sigma, pressure, gmma=gmma)/data0(grid)
        endelse
     end
     'h': begin
        data1d = scaleheight(rad, pressure, sigma)/rad
     end
     'nr':begin;radial brunt freq
        s = entropy(pressure, sigma, gmma=gmma)
        invlp = deriv(rad, alog(pressure))/gmma
        invls = deriv(rad, alog(s))/gmma
        csq = gmma*pressure/sigma ;adiabatic sound speed
        nrsq = -csq*invlp*invls
        data1d = nrsq*rad^3 ;normalize by kep freq^2
     end
     'alpha':begin
        data1d = alpha( rad, vrad, vtheta, sigma, pressure)
     end
     'mdot':begin
      data1d = -vrad*sigma
      cs0 =  soundspeed(pressure0, sigma0)
      data1d/=  cs0*sigma0
      end
      'j':begin
       data0     = vtheta0*rad*sigma0 ;initial ang mom dens
       tot_angmom0 = int_tabulated(rad, data0*2d0*!dpi*rad)
        
       data1d     = vtheta*rad*sigma   ; current ang mom dens  
       tot_angmom = int_tabulated(rad, data1d*2d0*!dpi*rad)

       print, 'tot ang mom/init ang mom', tot_angmom/tot_angmom0
 
       data1d = (data1d - data0)/data0
       end
     'prec': begin ;precession rate
      kappa = sqrt(epicyclesq(rad, vtheta))
      omega = vtheta/rad
      data1d = (omega - kappa)/omega
      end
 
  endcase

  return, data1d
end

pro compare_profiles, type=type, cases=cases, start=start, legend=legend,label=label $
                      ,xrange=xrange, yrange=yrange, ytickinterval=ytickinterval, r0=r0, hnorm=hnorm $
                      ,xtickinterval=xtickinterval, gmma=gmma, smallh=smallh, nopert=nopert, mass=mass 
  common share, refrad, pert, diskmass 
  !p.font = 0
 
  cases=strcompress(cases,/remove_all)
  numcases=n_elements(cases)
  
  location=strcompress(cases(0),/remove_all)

  if not keyword_set(r0) then r0 = 1.0 
  refrad = r0 
  if keyword_set(nopert) then begin
	pert = -1 
  endif else pert = 1

  if keyword_set(hnorm) then begin
     rplot = get_grid(cases(0), start(0), r0, /hnorm, smallh=smallh)
     xtitle=textoidl('N_H')
  endif else begin
     rplot = get_grid(cases(0), start(0), r0, smallh=smallh)
     xtitle=textoidl('r/r_0') 
  endelse 

  data1d = get_data(cases(0), start(0), type, gmma, smallh=smallh)
  
  case type of 
     'Q': begin
        ytitle = 'Q'
     end
     'dens': begin
        if not keyword_set(nopert) then begin
        ytitle = textoidl('\Delta\Sigma/\Sigma')
        endif else begin
        str = 'r='+string(r0,format='(f3.1)')+'r_0,t=0'
        str = '\Sigma('+str+')'
        str = strcompress('\Sigma/'+str, /remove_all)
        ytitle =  textoidl(str)
        endelse

     if keyword_set(mass) then print, diskmass

     end
     'ent': begin
        ytitle = textoidl('\DeltaS/S')
     end
     'vort': begin
        if not keyword_set(nopert) then begin
        ytitle = textoidl('\Delta\eta/\eta')
        endif else begin
        str = 'r='+string(r0,format='(f3.1)')+'r_0,t=0'
        str = '\eta('+str+')'
        str = strcompress('\eta/'+str, /remove_all)
        ytitle =  textoidl(str)
        endelse
     end
     'gvort': begin
        s   = texsyms()
        eta = s.tilde+'!8'+textoidl('\eta')+'!X'
        if not keyword_set(nopert) then begin
        ytitle = textoidl('\Delta'+eta+'/'+eta)
        endif else begin
        str = 'r='+string(r0,format='(f3.1)')+'r_0,t=0'
        str = eta+'('+str+')'
        str = strcompress(eta+'/'+str, /remove_all)
        ytitle =  textoidl(str)
        endelse
     end
     'h': begin       
        ytitle = textoidl('H/r')
     end
     'nr': begin
       ytitle = textoidl('N_r^2/\Omega_k^2') 
     end
     'alpha': begin
       ytitle = textoidl('\alpha')
     end 
     'mdot': begin
       ytitle = textoidl('-v_r\Sigma/(c_s\Sigma)_{t=0}')
     end
      'j': begin
       ytitle = textoidl('\Deltaj/j')
     end
     'prec': begin
       ytitle = textoidl('(\Omega-\kappa)/\Omega')
     end
     'disp': begin
       ytitle = textoidl('1 - \nu^2')
     end 
 endcase
  
  file    = strcompress('compare_profiles_'+type+string(start(0),format='(I03)')+'.ps',/remove_all)
  filename=filepath(file,root_dir='.',subdir=[location])
 
  set_plot, 'ps'
  device, filename=filename $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, rplot, data1d,xmargin=[8,2],ymargin=[3.5,0.5], ytitle=ytitle, ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle=xtitle, ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  
  for j=1, numcases-1 do begin

     if keyword_set(rh) then begin
        rplot = get_grid(cases(j), start(j), /rh)
     endif
     if keyword_set(hnorm) then begin
        rplot = get_grid(cases(j), start(j), /hnorm)
     endif
     
     data1d = get_data(cases(j), start(j), type, gmma, smallh=smallh)
     if keyword_set(mass) then print, diskmass  
    
     oplot, rplot,  data1d, thick=4, linestyle=j
   
  endfor

;     rc = 7.16
;     omp = rc^(-1.5)
     oplot, [0,1e10],[0,0]


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

