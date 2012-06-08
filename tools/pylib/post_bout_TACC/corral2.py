#! /usr/bin/env python
# note - these commands are only run by default in interactive mode

import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')

import post_bout
from read_inp import parse_inp, read_inp, read_log
import os
import numpy as np
import pickle

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
  
#called after after every run


class corral(object):

   def __init__(self,cached=False,refresh=False):
      self.cached = cached
      self.refresh = refresh
      self.print = 0
      self.status = read_log()
   def add(self,path=self.status['current']):
      a = post_bout.save(path=path) #save to current dir
   def refresh(self):
      self.refresh = True
      for i,path in enumerate(self.status['runs']):
         a = post_bout.save(path=path) #re post-process a run
   def info(self):
      print self.status
      print 'current:', self.status['current']
   def publish(self):
      alldata = []
      
      for i,val in enumerate(self.status['runs']):
         print val
         alldata.append(post_bout.read(path=val))
     
      alldata = np.array(alldata)
      print "alldata.shape: ", alldata.shape
      print  alldata[0]['Ni']['ave']['amp'].shape
      data = alldata[0]['Ni']['ave']['amp'] #Nt long array
      pp = PdfPages('output.pdf')
      plt.figure()
      cs = plt.plot(data)
      plt.title('amplitude ')
      plt.savefig(pp, format='pdf')

      gamma = alldata[2]['Ni']['modes'][0]['gamma'] #single value
      phase = alldata[2]['Ni']['modes'][0]['phase'] #Nt x Nx array 

      plt.figure()
      cs = plt.contour(phase)
      plt.clabel(cs, inline=1, fontsize=10)
      plt.title('phase')
      plt.savefig(pp, format='pdf')

      plt.close() 
      pp.close()
   
      explain = alldata[0]['meta']
      print explain
      
      allpickle = open('allpickled.pkl','wb')
      pickle.dump(alldata,allpickle)
      allpickle.close()

#determine the status of the run and if this is the last run of a batch 
#produce some pdfs

#we can look in $PWD, to look at BOUT.inp and ls output and decide 
#if we have enough run the script that will produce the final output

#we have the submission script produce a file that signal when to 
#run the final script rather than parse it with python and then figure
# it out by hand



