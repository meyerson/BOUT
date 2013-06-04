try :
      from netCDF4 import Dataset
    #from Scientific.IO.Netcdf . .  a less 
    #  complete module than netCDF4, no Dataset support
except ImportError:
      print "ERROR: netcdf4-python module not found"
      raise


def write_cyl_grid(gridfile='standard.nc',nx = 32,ny=3,dx=1,dy=1,
                   ni0=1.0e17,Te0=10.0,Ti0 =1.0,
                   Bz0=.01,bphi0=.1,
                   rMin = 1, rMax=1.5, Zmax=2.0,
                   ni_profile_type = 1,te_profile_type = 2,
                   ti_profile_type = 1,phi_profile_type = 0,
                   expr_prof = None,mxg = 2):

      #Tempratures are in eV
      #densities in m^-3
      #phi in V?

      nx = nx+2*mxg
      

  # mesh->get(I,    "sinty");// m^-2 T^-1
  # mesh->get(Psixy, "psixy");//get Psi             
  # mesh->get(Psiaxis,"psi_axis");//axis flux      
  # mesh->get(Psibndry,"psi_bndry");//edge flux 
  
  # When iysptrx=Nr:  periodic boundary condition
  # When iysptrx=0:   sheath plate, material wall condition
      

      variable_list =['Jpar0','pressure','bxcv','Rxy','Bpxy','Bxy','hthe','I']
      
      norms = ['Lbar','Bbar']
      
      f = Dataset(gridfile, "w",format='NETCDF3_CLASSIC')
      x = f.createDimension('x',nx)
      y = f.createDimension('y',ny)
      
      nxx = f.createVariable('nx','i4')
      nyy = f.createVariable('ny','i4')
      nxx[:] = nx
      nyy[:] = ny
      
      dxx = f.createVariable('dx','f8',('x','y'))
      dyy = f.createVariable('dy','f8',('x','y'))
      dxx[:] = dx
      dyy[:] = dy

  

  f.close()