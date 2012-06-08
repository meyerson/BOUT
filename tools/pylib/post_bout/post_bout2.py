import sqlite3 as sql
import sys
import pickle
import os
import json

try:
    sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
    import matplotlib
    import gobject
    import numpy as np
    from ordereddict import OrderedDict
except ImportError:
    print "can't find the modules I need, you fail"
    sys.exit() #no point in going on
            
try:
    import boutdata
    import boututils
except:
    print "can't find bout related modules, you fail"
    sys.exit() #no point in going on, shoot yourself 
 
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except:
    print "can' find matplotlib"

from read_inp import parse_inp, read_inp, read_log,metadata
from read_grid import read_grid
from basic_info import basic_info, fft_info
from corral import corral


def save(path='/home/cryosphere/BOUT/examples/Q3/data_short'): 
    print 'in post_bout.save'
    boutpath = path

    print path
    print sys.path
  
    meta = metadata(path=path)
   
    print 'ok got meta'
    output = OrderedDict()
    data = OrderedDict()

    con = sql.connect('test.db')
   
    all_modes = []
    print meta['evolved']['v'] #evolved fields for the given run

    for i,active in enumerate(meta['evolved']['v']): #loop over the fields
        print path, active
        data[active] = boutdata.collect(active,path=path)
        modes,ave = basic_info(data[active],meta)
        
