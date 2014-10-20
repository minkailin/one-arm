pro reduce, loc=loc

location=strcompress(loc,/remove_all)

nlines = file_lines(filepath('bigplanet0.dat',root_dir='.',subdir=[location]))

array = dblarr(11, nlines)

openr,3,filepath('bigplanet0.dat',root_dir='.',subdir=[location])
readf,3,array
close,3

for i=0, nlines-2 do begin
    time = array(7,i)
    timep1 = array(7,i+1)
    if(timep1 le time) then begin
        print, i
    endif
endfor

end
