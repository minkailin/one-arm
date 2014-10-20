function toomre, radius, density, a
h=0.05
cs=h/sqrt(radius)
omega=radius^(-1.5)
sigma=density*radius^(-a)
q=cs*omega/(!dpi*sigma)
return, q
end
