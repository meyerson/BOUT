def setup:
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
         import matplotlib.pyplot as plt
         
         from matplotlib.backends.backend_pdf import PdfPages
        
