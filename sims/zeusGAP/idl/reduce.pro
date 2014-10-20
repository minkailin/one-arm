pro reduce, loc=loc

location =strcompress(loc,/remove_all)
nlines = file_lines(filepath('dptorque.dat',root_dir='.',subdir=[location]))
tq = dblarr(8,nlines)
openr, 1, filepath('dptorque.dat',root_dir='.',subdir=[location])
readf,1,tq
close,1

for i=0, nlines-2 do begin
time = tq(0,i)
timep1 = tq(0,i+1)
if(timep1 le time) then print, i
endfor

end

