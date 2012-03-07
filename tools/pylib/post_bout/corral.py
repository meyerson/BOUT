#! /usr/bin/env python
# note - these commands are only run by default in interactive mode

import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')

import post_bout
from read_inp import parse_inp, read_inp, read_log
from basic_info import weighted_avg_and_std
import os
import numpy as np
import pickle
import json

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
  
#called after after every run

   
#post_bout.corral() will post-process the current run of a simulation and corral if done=True, this is the typically usage from bout

#post_bout.corral(cache=True) will assume that all the .pkl are in place and corral them into a pdf

#post_bout.corral(refresh=True) will rerun the post processign script for each simulation found in status.log and corral them all together into a pdf

def corral(cached=False,refresh=False):

   log  = read_log() 
   done = log['done']
   runs = log['runs'] #a list of all directories, we need this,
   # only need 'runs' if the simulation is done 
   current = log['current'] #always return the last data_dir    
   
   print log
   print 'current:', current
   
   if refresh==True:
      
      for i,path in enumerate(runs):
         print i,path
         a = post_bout.save(path=path) #re post-process a run
         
      
   elif cached == False: #if all the ind. simulation pkl files are in place skip this part
  
      a = post_bout.save(path=current) #save to current dir
      cached = True
  
   if done:
      alldata = []
      print 'last_one: '
      for i,val in enumerate(runs):
         print val
         alldata.append(post_bout.read(path=val))
      #remove the read in pickle

      alldata = np.array(alldata)
     
      result = LinRes(alldata)
      result.modes()
      a = result.show()

      return 0
   


class LinRes(object):
   def __init__(self,data):
      self.raw = data
      
      self.modekeys = data[0]['fields']['Ni']['modes'][0].keys()
      self.avekeys = data[0]['fields']['Ni']['ave'].keys()
      
      self.nrun = len(data) #number of runs
      
      self.nx = data[0]['meta']['Ni0']['v'].shape[0]
      self.ny = data[0]['meta']['Ni0']['v'].shape[1]
      self.nz = data[0]['meta']['MZ']['v']-1
      self.nt = data[0]['meta']['NOUT']['v']+1

      self.dt = [data[i]['meta']['dt']['v'] for i in range(self.nrun)]

      self.nmodes = [len(data[i]['fields']['Ni']['modes']) for i in range(self.nrun)]
      self.nfields = len(data[0]['fields'])

      self.k = []# np.zeros([nmodes,nx,2])
      self.nm =[]
      self.phase = []
      self.amp= []
      self.amp_n=[]
      self.dc= []
      self.freq =[]
      self.gamma = []

      self.key = ['run','field','harmonic','x','(t)']
   
   def modes(self): #
      runxmode = 0
      
      #collect key values across all runs and modes
      
      for i,elem in enumerate(self.raw): #runs, i isan index, elem is the value
         self.k.append(np.zeros([self.nfields,self.nmodes[i],self.nx,2]))
         self.nm.append(np.zeros([self.nfields,self.nmodes[i],self.nx,2]))
         
         self.phase.append(np.zeros([self.nfields,self.nmodes[i],self.nx,self.nt]))
         self.amp.append(np.zeros([self.nfields,self.nmodes[i],self.nx,self.nt]))
         self.amp_n.append(np.zeros([self.nfields,self.nmodes[i],self.nx,self.nt]))
         self.freq.append(np.zeros([self.nfields,self.nmodes[i],self.nx,2])) 
         self.gamma.append(np.zeros([self.nfields,self.nmodes[i],self.nx,2]))

       
         for f,fname in  enumerate(elem['fields']): #field,fname is a key
            for m,mode in enumerate(elem['fields'][fname]['modes']): #modes - mode is a dict

               print i,f,m
               
               self.k[i][f,m,:,:] = (np.array(mode['k'][1]).transpose())
               self.nm[i][f,m,:,:] =  (np.array(mode['k'][0]).transpose())

               self.phase[i][f,m,:,:] =  (np.array(mode['phase']).transpose())
               self.amp[i][f,m,:,:] =  (np.array(mode['amp']).transpose())
               self.amp_n[i][f,m,:,:] =  (np.array(mode['amp_n']).transpose())

               self.freq[i][f,m,:,:] = (np.array(mode['freq']).transpose())
               self.gamma[i][f,m,:,:] = (np.array(mode['gamma']).transpose())
               
               print type(self.k[i][f,m,:,:]),(self.k[i][f,m,:,:]).shape
      
         print 'len(self.k[i]): ',(self.k[i]).shape

      self.k = np.array(self.k,dtype = object) #in general every run will have a diff number of modes
      # self.nm = np.array(self.nm,dtype = object)
      # self.phase = np.array(self.phase,dtype = object)
      # self.amp = np.array(self.amp,dtype = object)
      # self.amp_n = np.array(self.amp_n,dtype = object)
      # self.dc = np.array(self.dc,dtype = object)
      # self.freq = np.array(self.freq,dtype = object)
      # self.gamma = np.array(self.gamma,dtype = object)

   def show(self):
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      pp = PdfPages('output.pdf')
      
      #runs,field, mods,modes are the first 3 indicies, x t follow
      
      # plt.figure()
      fig = Figure()
      canvas = FigureCanvas(fig)
      ax =fig.add_subplot(221)
      ax.scatter(self.k[0][0,:,:,1].flatten(),self.amp[0][0,:,:,0].flatten())
      ax.set_title(' Ni spectrum at t=0, all x run 0',fontsize=10)
      ax.set_xlabel(u'k \u03c1',fontsize=10)
      ax =fig.add_subplot(2,2,2)
      ax.scatter(self.k[0][0,0,:,1].flatten(),self.amp[0][0,0,:,0].flatten())
      ax.set_title(' Ni spectrum at t=0, all x run 0')
      ax.set_xlabel(u'k \u03c1')
      
      fig.savefig(pp, format='pdf')

      plt.figure()
      plt.scatter(self.nm[0][0,:,:,1].flatten(),self.amp[0][0,:,:,0].flatten())
      plt.title(' Ni spectrum at t=0, all x')
      plt.xlabel('n')
      plt.savefig(pp, format='pdf')

      plt.figure()
      for j in range(self.nrun):
         plt.scatter(self.nm[j][0,:,:,0].flatten(),self.nm[j][0,:,:,1].flatten(),c=colors[j])
         
      plt.title(' Ni spectrum at t=0, all x')
      plt.ylabel('m')
      plt.xlabel('n')
      plt.savefig(pp, format='pdf')

      #all runs
      
      plt.figure()
      
      #print self.amp.shape,self.k.shape,self.amp[7,0,0,:,0].flatten()

   
      for i in range(self.nrun):    
         print i
         plt.scatter(self.k[i][0,:,:,1].flatten(),
                     self.amp[i][0,:,:,0].flatten(),c=colors[i])

      plt.yscale('log')
      plt.title(' Ni spectrum ,IC- All runs')
      plt.xlabel(u'k \u03c1')
      plt.savefig(pp, format='pdf')

      plt.figure()
      for i in range(self.nrun):
         plt.scatter(self.k[i][0,:,self.nx/2,1].flatten(),
                     self.amp[i][0,:,self.nx/2,0].flatten(),c=colors[i])

      plt.xscale('log')
      plt.yscale('log')
      #plt.ylim(self.amp[0][:,:,:,:].min(),5*self.amp[0][:,:,:,:].max())
      plt.title(' Ni spectrum ,IC- All runs at x = nx/2')
      plt.xlabel(u'k \u03c1')
      plt.savefig(pp, format='pdf')
           
      fig = Figure()
      plt.figure()
      canvas = FigureCanvas(fig)
      #for i in range(self.nfields):
      
      Nplots = sum(self.nmodes)

      for j in range(self.nrun):
         for i in range(self.nmodes[j]):
            ax =fig.add_subplot(int(Nplots)/3,3,j+1)
            ax.scatter(range(self.nt),self.amp[j][0,i,int(self.nx/2),:].flatten(),c=colors[j])
         #ax.scatter(range(10),self.amp[j][0,0,int(self.nx/2),0:10].flatten(),c=colors[j])
            plt.scatter(range(self.nt),self.amp[j][0,i,int(self.nx/2),:].flatten(),c=colors[j])

         #print self.amp[j][0,0,int(self.nx/2),:]
            ax.set_yscale('log')
            ax.set_title(' Ni amp, all runs',fontsize=10)
            ax.set_xlabel('t')
            ax.set_ylim(self.amp[j][0,:,int(self.nx/2),:].min(),5*self.amp[j][0,:,int(self.nx/2),:].max())
      plt.yscale('log')
      plt.title(' Ni amp, all runs',fontsize=10)
      plt.xlabel('t')
      plt.ylim(self.amp[0][:,:,:,:].min(),5*self.amp[0][:,:,:,:].max())


      fig.savefig(pp, format='pdf')
      plt.savefig(pp, format='pdf')
      # gamma = np.diff(self.amp[j,0,0,int(self.nx/2),:].flatten())/
      #                self.amp[j,0,0,int(self.nx/2),:].flatten())
      

      plt.figure()
      for j in range(self.nrun):
         for i in range(self.nmodes[j]):
            #print len(np.diff(self.amp[j][0,0,int(self.nx/2),:].flatten()))
            #print self.nt
            plt.plot(np.diff(self.amp[j][0,i,int(self.nx/2),:].flatten())/self.amp[j][0,i,int(self.nx/2),0:self.nt-1].flatten(),c=colors[j])

      plt.title(' Ni gamma estimate all runs')
      plt.xlabel('t')
      plt.savefig(pp, format='pdf')

      

      
      gamma_est =[]
      for j in range(self.nrun):
         for i in range(self.nmodes[j]):
            gamma_est.append((np.diff(self.amp[j][0,i,int(self.nx/2),:].flatten())/
                              self.amp[j][0,i,int(self.nx/2),0:self.nt-1].flatten())/
                             np.array(self.dt[j]))
         
         
      
      gamma_est2 = np.array(gamma_est)

      print np.array(gamma_est)[:,-1].shape,self.dt
      gamma_est = np.array(gamma_est)[:,-1]#/np.array(self.dt)
      
      
#(np.diff(self.amp[:,0,0,int(self.nx/2),:],axis=1)/self.amp[:,0,0,int(self.nx/2),0:self.nt-1])
      
      print gamma_est2.shape
      print gamma_est2[0].shape

      gamma_w =  np.diff(gamma_est2,axis=1)
      
      np.abs(gamma_w).argmin(1) #index of the minimum for any given run
      

      
      #gamma_est2 = weighted_avg_and_std(gamma_est2,self.amp_n[:,0,0,int(self.nx/2),0:self.nt-1])
      #gamma_est2 = np.average(gamma_est2,1)/np.array(self.dt) #shoudl be nruns long
      
     

         #print len(np.diff(self.amp[j][0,0,int(self.nx/2),:].flatten()))
         #print self.nt
      
      #print len(gamma_est2[0]),self.k[:,0,0,int(self.nx/2),1].shap
      plt.figure()
      for j in range(self.nrun):
         plt.scatter(range(self.nt-2),gamma_w[j,:],c=colors[j])
      #plt.xscale('log')
      plt.title(' Ni weight')
      plt.xlabel('t')
      plt.savefig(pp, format='pdf')
                  
      # plt.figure()
      # plt.scatter(self.k[:,0,0,int(self.nx/2),1],gamma_est2)
      # #plt.xscale('log')
      # plt.title(' Ni gamma 2 estimate static')
      # plt.xlabel('t')
      # plt.savefig(pp, format='pdf')

      # plt.figure()

      # plt.scatter(self.k[:,0,0,int(self.nx/2),1],gamma_est)
      # #plt.xscale('log')
      # plt.title(' Ni gamma estimate static')
      # plt.xlabel('t')
      # plt.savefig(pp, format='pdf')
 
  
  
      plt.figure()

      #GAMMA VS K
      
      # print self.k.shape,self.k.shape
      # print self.gamma.shape,self.gamma.shape

      for j in range(self.nrun):
         for i in range(self.nmodes[j]):
            plt.errorbar(self.k[j][0,i,int(self.nx/2),1].flatten(),
                         self.gamma[j][0,i,int(self.nx/2),0].flatten(),
                         yerr = self.gamma[j][0,i,self.nx/2,1].flatten(),
                         fmt ='ro',c=colors[j]) 
    
      #plt.yscale('log')
      plt.yscale('log') 
      plt.xscale('log')
      plt.title('gamma vs k, mode 0 all runs Ni ')
      plt.savefig(pp, format='pdf')

      
      # w vs k
      plt.figure()
      for j in range(self.nrun):
         for i in range(self.nmodes[j]):
            plt.scatter(self.k[j][0,i,self.nx/2,1].flatten(),
                        self.freq[j][0,i,self.nx/2,0].flatten(),c=colors[j])
      
      plt.title('omega vs k, all modes all runs Ni ')
      plt.yscale('log') 
      plt.xscale('log')
      plt.savefig(pp, format='pdf')
            
      plt.close() 
      pp.close()

      return gamma_w
