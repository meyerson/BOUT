#! /usr/bin/env python
# note - these commands are only run by default in interactive mode
import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')



the_world_is_flat = 1
if the_world_is_flat:
    print "Be careful not to fall off!"

#some usefull variables
path='/home/cryosphere/BOUT/examples/Q3/data_hd'
var='Ni'

print "Trying to import GTK..."
import gobject
widget = "gtk"

try:
    import matplotlib
    if widget == "gtk":
        matplotlib.use('GTKAgg')
    else:
        matplotlib.use('WXAgg') # do this before importing pylab
        
    import numpy as np
    import matplotlib.pyplot as plt
        
except ImportError:
    print "ERROR: Showdata needs numpy, matplotlib and gobject modules"
    raise

#lets collect the data
print var 
print path
print sys.path

import boutdata
import boututils

#from boutdata import *
#from boututils import *


ni = boutdata.collect(var,path=path)

data = ni[:,:,5,:]
boututils.showdata(data)
size = data.shape
ndims = len(size)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Title')
cmap = None
m = plt.imshow(data[0,:,:], interpolation='bilinear', cmap=cmap, animated=True,aspect='auto')
plt.show()



#print ni.ndims

#help(ni)

#import os
#if os.path.isfile(os.environ['PYTHONSTARTUP']):
#    execfile(os.environ['PYTHONSTARTUP'])


