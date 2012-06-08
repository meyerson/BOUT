#! /usr/bin/env python
# note - these commands are only run by default in interactive mode

import sys
sys.path.append('/home/cryosphere/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout_TACC')
sys.path.append('/usr/local/pylib')

import post_bout_TACC as post_bout
from ListDict import ListDictKey, ListDictFilt
from read_inp import parse_inp, read_inp, read_log
from basic_info import weighted_avg_and_std
from read_cxx import read_cxx, findlowpass



import os
import numpy as np
import pickle
import subprocess  


def corral(cached=True,refresh=False,debug=False,IConly=1,logname='status.log'):

   log  = read_log(logname=logname) 
   #done = log['done']
   runs = log['runs'] #a list of all directories, we need this,
   # only need 'runs' if the simulation is done 
  
   current = log['current'] #always return the last data_dir    
   
   print log
   print 'current:', current
   
   if refresh==True:
      
      for i,path in enumerate(runs):
         print i,path
         a = post_bout.save(path=path,IConly=IConly) #re post-process a run
         
         
      
   elif cached == False: #if all the ind. simulation pkl files are in place skip this part
  
      a = post_bout.save(path=current) #save to current dir
      #here is really where you shoudl write to status.log
      #write_log('status.log',
      cached = True
  
   #if done:
     
   alldata = []
   alldb =[]
   print 'last_one: '
   for i,val in enumerate(runs):
      print val
      db = post_bout.read(path=val)
      #alldata.append(array)
      alldb.append(db)
         
         #build the end database
         
         
      #remove the read in pickle

   #return alldb
   def islist(input):
      return isinstance(input,list)

   alldb = filter(islist,alldb)

   alldb = sum(alldb,[])
   #alldata = np.array(alldata)
   
   nt = []
   for mode in alldb:
      nt.append(len(mode['amp']))
      nt = [max(nt)]
   
   nt = nt[0]
   t = range(nt)             
   i = 0
    
   #return alldb
   result = LinRes(alldb)
   
  
      #s = subset(alldata,alldb,'dz',[1.0])
   if debug:
      return result
   else:
      result.show2()
      return 0

class LinRes(object):
   def __init__(self,alldb):
      #self.raw = data
      self.db = alldb

      #self.modekeys = data[0]['fields']['Ni']['modes'][0].keys()
      print len(alldb)
      
      self.meta = np.array(ListDictKey(alldb,'meta'))[0]
      
      self.keys = alldb[0].keys()
      #self.avekeys = data[0]['fields']['Ni']['ave'].keys()
      
      #self.nrun = len(alldb) #number of runs
      
      self.path= np.array(ListDictKey(alldb,'path'))
      self.cxx =[]
      self.maxN =[]

      [self.cxx.append(read_cxx(path=elem,boutcxx='2fluid.cxx.ref')) for elem in self.path]
      [self.maxN.append(findlowpass(elem)) for elem in self.cxx]
      

      #self.maxN = findlowpass(self.cxx) #low pass filt from .cxx

      self.nx = np.array(ListDictKey(alldb,'nx'))[0]
      self.ny = np.array(ListDictKey(alldb,'ny'))[0]
      self.nz = np.array(ListDictKey(alldb,'nz'))[0]

      #self.nt = int(data[0]['meta']['NOUT']['v']+1)
      self.Rxy = np.array(ListDictKey(alldb,'Rxy'))
      self.Rxynorm = np.array(ListDictKey(alldb,'Rxynorm'))
      self.nt = np.array(ListDictKey(alldb,'nt'))

     
      self.dt = np.array(ListDictKey(alldb,'dt')) 
      self.nfields = np.array(ListDictKey(alldb,'nfields')) 

      self.field = np.array(ListDictKey(alldb,'field'))
      
      self.k = np.array(ListDictKey(alldb,'k'))
      self.mn = np.array(ListDictKey(alldb,'mn'))

      #return ListDictKey(alldb,'phase')

      #self.phase = np.array(ListDictKey(alldb,'phase'))
      self.phase = ListDictKey(alldb,'phase')

      # self.amp= np.array(ListDictKey(alldb,'amp'))
      # self.amp_n=np.array(ListDictKey(alldb,'amp_n'))
      # self.dc= []
      # self.freq = np.array(ListDictKey(alldb,'k'))
      # self.gamma = np.array(ListDictKey(alldb,'gamma'))
      
      self.amp= ListDictKey(alldb,'amp')
      self.amp_n = ListDictKey(alldb,'amp_n')
      self.dc= []
      self.freq = np.array(ListDictKey(alldb,'k'))
      self.gamma = np.array(ListDictKey(alldb,'gamma'))
     
      self.freq  = ListDictKey(alldb,'freq')
      
      self.IC = np.array(ListDictKey(alldb,'IC'))
      self.dz = np.array(ListDictKey(alldb,'dz'))
      self.meta['dz'] = np.array(list(set(self.dz).union()))

      self.nmodes = self.dz.size
    
      self.MN = np.array(ListDictKey(alldb,'MN'))
      #self.MN = np.float32(self.mn) 
      #self.MN[:,1] = self.mn[:,1]/self.dz
      self.nrun =  len(set(self.path).union())
      self.L = np.array(ListDictKey(alldb,'L'))
      #self.C_s = 
      self.modeid = np.array(ListDictKey(alldb,'modeid'))

      self.model()
      
      #self.model = model(self.k,self.L)
   
   def _amp(self,tind,xind):
      #first select modes that actually have valid (tind,xind)
      #indecies
      #s = subset(self.db,'modeid',modelist)
      return np.array([self.amp[i][tind,xind] 
                       for i in range(self.nmodes)])
   

               
 
   
   def model(self,field='Ni',plot=False):
      
      #enrich the object
      allk = self.k[:,1,self.nx/2] #one location for now
   
      
      self.M = []
      self.eigsys = []
      self.gammaA = []
      self.omegaA = []
      self.eigvec = [] 
      self.gammamax = []
      self.omegamax = []

      #allk = np.arange(0.1,100.0,.1)
     # allk=  np.sort(list(set(allk).union()))

      
      for i,k in enumerate(allk):
         print i
         #M =np.matrix(np.random.rand(3,3),dtype=complex)
         M = np.zeros([2,2],dtype=complex)
         M[0,0] = 0
         M[0,1] = k/(self.L[i,self.nx/2,self.ny/2])
         M[1,0] = (2*np.pi/self.meta['lpar'][self.nx/2])**2 * self.meta['sig_par'][0]*complex(0,k**-2)
         M[1,1]= -(2*np.pi/self.meta['lpar'][self.nx/2])**2 * self.meta['sig_par'][0]*complex(0,k**-2)
         
         eigsys= np.linalg.eig(M)  
         gamma = (eigsys)[0].imag
         omega =(eigsys)[0].real
         eigvec = eigsys[1]
         
         self.M.append(M)
         self.eigsys.append(eigsys)
         self.gammaA.append(gamma)
         self.gammamax.append(max(gamma))
         where = gamma == gamma.max()
         self.omegamax.append(omega[where])
         self.eigvec.append(eigvec)
         self.omegaA.append(omega)

class subset(LinRes):
   def __init__(self,alldb,key,valuelist,model=False):
      selection = ListDictFilt(alldb,key,valuelist)
      if len(selection) !=0:
         LinRes.__init__(self,selection)
         self.skey = key
         if model==True:
            self.model()
      else:
         LinRes.__init__(self,alldb)
         if model==True:
            self.model()

         

