pro pdisk_m1_amp, mode=mode, r0=r0

if not keyword_set(mode) then mode = 1
if not keyword_set(r0)   then r0 = 1.0

p0 = 2.0*!dpi*r0^1.5

file = 'pdisk_modes_amp.dat'
lines = file_lines(file)

array=dblarr(7,lines)

openr,1,file
readf,1,array
close,1

;find max 
temp = max(array(mode,*), grid)

print, 'max m=', mode, ' amp at t=', array(0,grid)/p0

avg = mean(array(mode,grid:lines-1))
;avg = mean(array(1,5:10))
print, 'average mode amp after max is', avg
end
