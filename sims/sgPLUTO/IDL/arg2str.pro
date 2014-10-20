; ***************************************************************
;   string format conversion routine
; ***************************************************************
FUNCTION ARG2STR, x

  szx = SIZE(x)

  IF (szx[1] EQ 4 OR szx[1] EQ 5) THEN BEGIN  ; arguments is either double 
    sgnx = x GT 0.0 ? 1.0:-1.0                ; or float
    n    = 0
    IF (abs(x) GT 1.e-7) THEN n = floor(alog10(abs(x)))

    IF (n GT 2) THEN BEGIN
      y = x/10.0^n
      sn = strcompress(string(n,format='(i8)'),/remove_all)
      sx = string(y,format='(f18.3)')
      sx += " x 10^"+sn
    ENDIF ELSE BEGIN
      sx = string(x,format='(f18.3)')
    ENDELSE
    RETURN,strcompress(sx,/remove_all)
  ENDIF ELSE BEGIN
    RETURN, strcompress(string(x),/remove_all)
  ENDELSE

END


