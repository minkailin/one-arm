pro resolution, r1=r1, r2=r2, nh=nh, smallh=smallh
nrad = nh*(r2 - r1)/smallh
nphi = 2d0*!dpi*nh/smallh

print, 'ideal r-phi resolution is', nrad, nphi
end
