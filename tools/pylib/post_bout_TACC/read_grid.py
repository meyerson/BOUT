#modelled after boutdata.collect
def read_grid(gridfile='/home/cryosphere/Q-Machine/grids/Q3_short.nc'):
   try :
      from netCDF4 import Dataset
    #from Scientific.IO.Netcdf . .  a less 
    #  complete module than netCDF4, no Dataset support
   except ImportError:
      print "ERROR: netcdf4-python module not found"
      raise
   
   
   def read_var(gridfile, name):
        var = file.variables[name]
        return var[:]
   
   print gridfile
   f = Dataset(gridfile, "r")
   
   # print "File format    : " + f.file_format
   # print f
   # print f.variables
   return f
