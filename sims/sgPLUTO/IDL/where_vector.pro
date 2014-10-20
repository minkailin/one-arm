function exist,filename,count
;+
; NAME:
;   EXIST
; PURPOSE:
;   A very simple check to see if a file exists...
; CALLING SEQEUNCE:
;   tmp = Exist('STARS.DAT')
; INPUT:
;   FILENAME  This is the filename or search spec. that should be checked.
; OUTPUT:
;   TMP       The returned result is 1 is the file or files exist and 0 if
;               the file of files do(es) not exist.
; OPTIONAL OUTPUT:
;   COUNT     Number of occurances of the given search spec.
; EXAMPLE:
;   if exist('tmp.tmp') then print,'Yes' else print,'No'
;   if not exist('tmp.tmp') then print,'Create'
;   if exist('*.hhh',count) then print,strn(count),' Header files available'
; HISTORY:
;   27-JUL-1992 Header added to old routine  (E. Deutsch)
;
;-

a=findfile(filename,count=count)

  return,count<1

end


function where_vector,vector,array,count,nosort=nosort,$
       trim_string=trim_string,case_sens=case_sens,rest=rest,rcount=rcount

;temp = exist(vector)

;if not temp  or not temp then return,-1

count=0
np=n_elements(array)
rcount=np
rest=lindgen(np)

;-- protect inputs and modify

trim_string=keyword_set(trim_string)
case_sens=keyword_set(case_sens)

svec=vector & sarr=array

if not keyword_set(nosort) then begin
 rs=uniq([svec],sort([svec]))
 svec=svec(rs) 
endif

if datatype(vector) eq 'STR' then begin
 if trim_string then svec=strtrim(svec,2)
 if not case_sens then svec=strupcase(svec)
endif
if datatype(array) eq 'STR' then begin
 if trim_string then sarr=strtrim(sarr,2)
 if not case_sens then sarr=strupcase(sarr)
endif

state=''
nvecs=n_elements(svec)
pieces=strarr(nvecs)
v=svec & s=sarr
for i=0,nvecs-1 do begin
 index=strtrim(string(i),2)
 pieces(i)='(v('+index+') eq s)'
 if i eq 0 then pieces(i)='clook=where('+pieces(i)
 if i eq (nvecs-1) then pieces(i)=pieces(i)+',count)'
 if (nvecs eq 1) or (i eq 0) then conn='' else conn=' or '
 state=state+conn+pieces(i)
endfor

status=execute(strcompress(strtrim(state,2)))

if count gt 0 then begin
 rest(clook)=-1
 rlook=where(rest gt -1,rcount)
 if rcount gt 0 then rest=rest(rlook) else rest=-1
endif

return,clook & end
