pro reduce_tqwk, loc=loc

location=strcompress(loc,/remove_all)

nlines = file_lines(filepath('tqwk0.dat',root_dir='.',subdir=[location]))

array = dblarr(10, nlines)

openr,3,filepath('tqwk0.dat',root_dir='.',subdir=[location])
readf,3,array
close,3

for i=0, nlines-2 do begin
    time = array(9,i)
    timep1 = array(9,i+1)
    if(timep1 le time) then begin
        print, i
    endif
endfor

end
