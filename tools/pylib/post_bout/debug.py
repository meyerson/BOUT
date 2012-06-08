import matplotlib.pyplot as plt
from netCDF4 import Dataset
plt.figure()
f_re = Dataset("BOUT.restart.13.nc","r")
f_re2 = Dataset("../data_eig64_.5/BOUT.restart.13.nc","r")
f = Dataset("BOUT.dmp.13.nc","r")




Ni = f.variables['Ni']
Ni.shape
plt.contour(Ni[0,:,3,:])
plt.show()


Ni_re2 = f_re2.variables['Ni']
plt.contour(Ni_re2[:,3,:])
plt.show()

phi_re = f_re.variables['phi']
plt.contour(phi_re[:,3,:])
Ni_re = f_re.variables['Ni']
plt.contour(Ni_re[:,3,:])
plt.show()


plt.show()


