##################################################
#            BOUT++ data package
#
# Routines for examining simulation results for BOUT++
#
##################################################

print "Loading BOUT++ post processing routines"

# Load routines from separate files
import sys

try:
    sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    
    import matplotlib
    import gobject
    import numpy as np
except ImportError:
    print "can't find the modules I need, you fail"
    sys.exit() #no point in going on
      
      
#import some bout specific modules
try:
    import boutdata
    import boututils
except:
    print "can't find bout related modules, you fail"
         
#import some home-brewed modules
         
         
#create some aliases


try:
    from read_grid import read_grid
except:
    print "Sorry, no read_grid"

try:
    from read_inp import parse_inp, read_inp, read_log,metadata
except:
    print "Sorry no parse_inp"

try:
    from post_bout import show, save, read
except:
    print "Sorry, no show"

try:
    from basic_info import basic_info, fft_info
except:
    print "Sorry, no basic_info"

try:
    from corral import corral,LinRes
except:
    print "No corral"
