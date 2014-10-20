pro gap, yrange=yrange, xtickinterval=xtickinterval, ytickinterval=ytickinterval

array_nsg=dblarr(4,1024)
openr,1,'gap_theory_nsg.dat'
readf,1,array_nsg
close,1

radius_nsg = array_nsg(0,*)
basic_nsg = array_nsg(1,*)
gap_nsg = array_nsg(2,*)
theory_nsg = array_nsg(3,*)

rel_gap_nsg = gap_nsg/basic_nsg - 1.0
rel_theory_nsg = theory_nsg/basic_nsg - 1.0

array_sg=dblarr(4,1024)
openr,1,'gap_theory_sg.dat'
readf,1,array_sg
close,1

radius_sg = array_sg(0,*)
basic_sg = array_sg(1,*)
gap_sg = array_sg(2,*)
theory_sg = array_sg(3,*)

rel_gap_sg = gap_sg/basic_sg - 1.0
rel_theory_sg = theory_sg/basic_sg - 1.0

set_plot, 'ps'
device, filename='gap_plot.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,radius_sg,rel_gap_sg,xmargin=[8,2],ymargin=[3,2] $
,ytitle='relative gap depth',xtitle='r' $
,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval,title='100 orbits' $
,xrange=xrange,yrange=yrange
oplot, radius_sg, rel_theory_sg, thick=2

oplot, radius_nsg, rel_gap_nsg, thick=4, linestyle=1
oplot, radius_nsg, rel_theory_nsg, thick=2, linestyle=1

oplot, radius_sg,(gap_sg-gap_nsg)/gap_nsg, thick=4, linestyle=2
oplot, radius_sg, (theory_sg - theory_nsg)/theory_nsg, thick=2, linestyle=2




device,/close
end
