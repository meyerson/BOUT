#! /usr/bin/python
import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
sys.path.append('/usr/local/pylib')

import matplotlib
matplotlib.use('Agg')
from read_inp import metadata
import sys
import os

path=sys.argv[1]
key=sys.argv[2]

cache='/tmp/hlmk/'+key
if os.path.exists(cache):
    os.rmdir(cache)
else:
    os.makedirs(cache)

meta = metadata(path=path)

from boutdata import collect
from boututils import savemovie
gamma = collect("Gammax",yind=[5,5],path=path)
ni = collect("Ni",yind=[5,5],path=path)

#one movie per cpu

savemovie(gamma[:,:,0,:],ni[:,:,0,:],moviename='movie_'+key+'.avi',
          meta=meta,cache=cache+"/",overcontour=True)

os.rmdir(cache)

