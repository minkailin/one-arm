pro compare, files=files, frame=frame
common array, data1d
files=string(files)
num=n_elements(files)
for i=0, num-1 do begin
vort1d, start=frame, finish=frame, loc=files(i)
print, data1d(128)
endfor
end
