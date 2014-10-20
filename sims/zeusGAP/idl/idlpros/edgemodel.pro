; first define the profile function, it takes in an array
function profile, array
array=double(array)
r=array(0)
F=array(1)
xi=array(2)
re=array(3)
return, 1.0+((F-1)/(F+1))*tanh(xi*(r-re))
end
;this is dln profile /d ln r, beta is minus this because sigma propto 1/profile
function dprofile, array
array=double(array)
r=array(0)
F=array(1)
xi=array(2)
re=array(3)
top=(F-1.0)*r*xi*(cosh(xi*(r-re)))^(-2.0)
return, top/((F+1.0)*profile(array))
end
;ddprofile is r*dprofile/dr
function ddprofile, array
array=double(array)
r=array(0)
F=array(1)
xi=array(2)
re=array(3)
term1=(F-1.0)*r*xi*(cosh(xi*(r-re)))^(-2.0)
term2=(-1.+F)*r*xi*(cosh(xi*(r-re)))^(-2.0)
term3=(1.+F+(-1.+F)*tanh(xi*(r-re)))*(-1.+2.*r*xi*tanh(xi*(r-re)))
term4=(1.+F+(-1.+F)*tanh(xi*(r-re)))^(2.0)
return, term1*(term2+term3)/term4
end
;main plotting program
PRO edgemodel, params=params, range=range, var=var
; this procedure plot the initial disc configuration. variables is the
; 1D array of input parameters=[nuplus,sigmaplus,F,xi,re].range is
; range of disc to plot. first define the parameters
params=double(params)
range=double(range)
nuplus=params(0)
sigmaplus=params(1)
F=params(2)
xi=params(3)
re=params(4)
rin=range(0)
rout=range(1)
;generate array of r values and arrays to hold viscosity, sigma and ur
r=rin+(dindgen(100)/100.)*(rout-rin)
viscosity=dblarr(100)
sigma=dblarr(100)
ur=dblarr(100)
beta=dblarr(100)
rdbeta=dblarr(100)
dur=dblarr(100)
;fill in arrays
;calculate prefactor for dur here for convinience
urpref=-3.0*nuplus*(F+1.0)/(4.0*F)
for i=0,99 do begin
array2=[r(i),F,xi,re]
prof=profile(array2)
viscosity(i)=nuplus*(F+1.0)*prof/(2.0*F)
sigma(i)=2.0*F*sigmaplus/((F+1.0)*prof)
ur(i)=urpref*prof/r(i)
beta(i)=-dprofile(array2)
rdbeta(i)=-ddprofile(array2)
dur(i)=urpref*(dprofile(array2)-1.0)*prof/(r(i)*r(i))
endfor
viscosity=viscosity*10.^5.
sigma=sigma*10.^4.
set_plot, 'ps'
device, filename=var+'.ps', xsize=8, ysize=6, xoffset=0, yoffset=0, /inches
 if (var eq 'nusig') then begin
plot, r, viscosity, xmargin=[7,6], ymargin=[3,3], xtitle='r',/nodata,xrange=range $
,xtickinterval=0.4,xminor=2,ystyle=4,charsize=1.5,thick=2
axis,yaxis=0,/save,ytitle=textoidl('\nu/10^{-5}')+', solid',yrange=[min(viscosity),max(viscosity)]$
,ystyle=2,charsize=1.5
oplot,r,viscosity
axis,yaxis=1,/save,ytitle=textoidl('\Sigma/10^{-4}')+', dashed',yrange=[min(sigma),max(sigma)]$
,ystyle=2,charsize=1.5
oplot,r,sigma,linestyle=2
 endif else if (var eq 'beta') then begin
 plot, r, beta, xmargin=[7.5,6], ymargin=[3,3], xtitle='r',/nodata,xrange=range $
 ,xtickinterval=0.4,xminor=2,ystyle=4,charsize=1.5,thick=2
 axis,yaxis=0,/save,ytitle=textoidl('\beta')+', solid',yrange=[min(beta),max(beta)]$
 ,ystyle=2,charsize=1.5
 oplot,r,beta
 axis,yaxis=1,/save,ytitle=textoidl('r(d\beta/dr)')+', dashed',yrange=[min(rdbeta),max(rdbeta)]$
 ,charsize=1.5
 oplot,r,rdbeta,linestyle=2
 device,/close
 endif else if (var eq 'ur') then begin
ur=ur*10.^6.0
dur=dur*10.^6.0
 plot, r, ur, xmargin=[7.5,7.5], ymargin=[3,3], xtitle='r',/nodata,xrange=range $
 ,xtickinterval=0.4,xminor=2,ystyle=4,charsize=1.5,thick=2
 axis,yaxis=0,/save,ytitle=textoidl('u_r/10^{-6}')+', solid',yrange=[min(ur),max(ur)]$
 ,charsize=1.5
 oplot,r,ur
 axis,yaxis=1,/save,ytitle=textoidl('(du_r/dr)/10^{-6}')+', dashed',yrange=[min(dur),max(dur)]$
 ,charsize=1.5
 oplot,r,dur,linestyle=2
 endif
device,/close
END

