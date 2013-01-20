pro primitive, filename = filename

  if not(keyword_set(filename)) then filename = 'primitive.nc'

  id = NCDF_CREATE(filename, /CLOBBER)
; Create a new NetCDF file with the filename inquire.nc.

  NCDF_CONTROL, id, /FILL	; Fill the file with default values.

 

  data = FLTARR(1,32)     

  data = {nx:128, ny:128}
  names = ['nx','ny']
     
  ;grdpdb = file_export('primitive',data,custom_name= names)
  ;filename = 'primitive.nc'
  
  handle = file_open(filename, /create)

  ;; FOR i=0, N_TAGS(data)-1 DO BEGIN
  ;;   var = data.(i)
  ;;   status = file_write(handle, names[i], var)
  ;; ENDFOR
  
; grdpdb = file_export(filename,data,custom_name= names)
 
  FOR i=0, N_TAGS(data)-1 DO BEGIN
    var = data.(i)
    status = file_write(handle, names[i], var)
  
 ENDFOR

end
