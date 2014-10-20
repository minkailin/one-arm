;THIS IS A PROCEDURE FOR ESTIMATING POST-SHOCK INVERSE VORTENSITY
;VALUES FOR THE PLANET-DISC SCENARIO. METHOD FOR VORTICITY JUMP
;FOLLOWS KEVLAHAN (1997), WITH LOCAL ISOTHERMAL EQUATION OF STATE AND
;ADDITIONAL CORIOLIS TERM. 
;Basic state
function vy_k, x;y-velocity in shearing box. Keplerian.
common params, r0, h, y_s
result=double((r0+x)^(-0.5)-r0^(-0.5))
return, result
end
function dvy_k, x;dvy_k/dx
common params, r0, h, y_s
result=double(-0.5*(r0+x)^(-1.5))
return, result
end
function cs, x;Sound speed
common params, r0, h, y_s
result=double(h*(r0+x)^(-0.5))
return, result
end
function dcs, x;dcs/dx
common params, r0, h, y_s
result=double(-0.5*h*(r0+x)^(-1.5))
return, result
end
;Shock curve
function yshock, x; Height of shock in shearing box co-ordinates
common output, coeff
x=double(x)
for i=0, n_elements(coeff)-1 do y_s=y_s+coeff(i)*x^i
return, y_s
end
function alpha,x;Acute angle between shock and x-axis. -tan(alpha)=ydash
common output, coeff
x=double(x)
ydash=0.0
for i=1, n_elements(coeff)-1 do ydash=ydash+i*coeff(i)*x^(i-1)
angle=atan(-ydash)
return, angle
end
function ddy, x;d2y/dx2, shock curvature
common output, coeff
x=double(x)
yddash=0.0
for i=2, n_elements(coeff)-1 do yddash=yddash+i*(i-1)*coeff(i)*x^(i-2)
return, yddash
end
;Pre-shock velocity field, vy and assume vx=0
function vy, x; y-velocity after acceleration due to planet.
common disc, sigma, mp
x=double(x)
vyk=vy_k(x)
result=vyk^2.0+2.0*mp*(x^2.0+yshock()^2.0)^(-0.5)
result=sqrt(result)
if vyk le 0 then result=-result
return, result
end
function dvy, x; dvy/dx
common disc, sigma, mp
x=double(x)
a=double(alpha())
result=vy_k(x)*dvy_k(x)-mp*cos(a)*(x-tan(a)*yshock())*(x^2.0+yshock()^2.0)^(-3.0/2.0)
result=result/vy(x)
return, result
end
function dA, x; dA/dS
x=double(x)
a=double(alpha())
term1=dvy(x)*(cos(a))^2.0
term2=ddy(x)*vy(x)*sin(a)*(cos(a))^3.0
return, term1+term2
end
function BdB, x; BdB/dS
x=double(x)
a=double(alpha())
term1=-vy(x)*vy(x)*(cos(a))^4.0*sin(a)*ddy(x)
term2=vy(x)*dvy(x)*cos(a)*(sin(a))^2.0
return, term1+term2
end
function dM, x; dM/dS
a=double(alpha())
result=dA(x)/cs(x)-vy(x)*(cos(a))^2.0*dcs(x)/(cs(x)*cs(x))
return, result
end
;Vorticity stuff
function vorticity, x; Basic state vorticity distribution for Keplerian flow.
common params, r0, h, y_s
x=double(x)
return, 0.5*(r0+x)^(-1.5)
end
function domega, x; Vorticity jump at x.
common params, r0, h, y_s
common disc, sigma, mp
a=double(alpha())
aa=vy(x)*cos(a)
mach=aa/cs(x)
f1=mach^(-2.0)-1.0
f2=-2.0*cs(x)*f1
f3=aa/(cs(x)*cs(x))-1.0/aa
term1=f1*(dA(x)+2.0*r0^(-1.5)/sigma);Assumed frame speed is Keperial value at r0 (=planet speed).
term2=f2*dM(x)
term3=f3*(2.0*cs(x)*cos(a)*dcs(x)+BdB(x))
return, term1+term2+term3
end
pro jump, x_shock=x_shock, y_shock=y_shock, x_obs=x_obs
common params, r0, h, y_s
common disc, sigma, mp
common output, coeff
;Constants
r0=double(1.98)
h=double(0.05)
sigma=double(7.0*10.0^(-4.0))
mp=double(0.00028)
rhill=(mp/3.0)^(1.0/3.0)*r0
y=double(alpha())
y_s=double(y_shock*!pi*r0)
;Shock curve from data
shock_fit,loc='14',start=17,xrange=[0,4],yrange=[-.1,.1], eta=0.9, order=5,/nplot
;Convert inputs
x_shock=double(x_shock*rhill)
x_obs=double(x_obs*rhill)
;Post-shock vorticity, and observed vortensity change at x_obs
omega_new=vorticity(x_shock)+domega(x_shock)
mach_sq=(vy(x_shock)*cos(y)/cs(x_shock))^2.0
;print, mach_sq
vorten=sigma*mach_sq/omega_new
;no_shock=where(mach_sq le 1.0)
dvort=vorten-sigma/vorticity(x_obs)
shocks=where(mach_sq gt 1.0)
;dvort(no_shock)=0.0
set_plot, 'ps'
device, filename='jump_.ps', xsize=8, ysize=6, xoffset=0, yoffset=0, /inches, bits_per_pixel=8,/color
plot, y(shocks)/!pi, dvort(shocks)*1000., xtitle=textoidl('\alpha_s/\pi (shock angle)'), ytitle=textoidl('\Delta(\Sigma/\omega)\times1000') $
, title=textoidl('r_{obs}=')+string(r0+x_obs,format='(f4.2)')+textoidl(', r_{shock}=') +string(r0+x_shock,format='(f4.2)') $
,charsize=1.5
device,/close
end







































