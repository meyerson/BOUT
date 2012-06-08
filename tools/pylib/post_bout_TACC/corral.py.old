#! /usr/bin/env python
# note - these commands are only run by default in interactive mode

import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')

import post_bout
from ListDict import ListDictKey, ListDictFilt
from read_inp import parse_inp, read_inp, read_log
from basic_info import weighted_avg_and_std
from read_cxx import read_cxx

import os
import numpy as np
import pickle
import json

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages


from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.lib import colors

import urllib2


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
         a = post_bout.save(path=path,IConly=1) #re post-process a run
         
      
   elif cached == False: #if all the ind. simulation pkl files are in place skip this part
  
      a = post_bout.save(path=current) #save to current dir
      cached = True
  
   if done:
      alldata = []
      alldb =[]
      print 'last_one: '
      for i,val in enumerate(runs):
         print val
         array,db = post_bout.read(path=val)
         alldata.append(array)
         alldb.append(db)
         
         #build the end database
         
         
      #remove the read in pickle
      alldb = sum(alldb,[])
      alldata = np.array(alldata)
     
      result = LinRes(alldata,alldb)
  
      #s = subset(alldata,alldb,'dz',[1.0])
      #result.show2()
      return result
   
#data[run#]['fields']['fieldNAME']['modes'][mode#][nx x nt] array

#i want to have a collection of modes
#that is populated by looping over all runs and their modes
#and is NOT indexed, but can

class LinRes(object):
   def __init__(self,data,alldb):
      self.raw = data
      self.db = alldb

      #self.modekeys = data[0]['fields']['Ni']['modes'][0].keys()
      print len(alldb)

      self.keys = alldb[0].keys()
      self.avekeys = data[0]['fields']['Ni']['ave'].keys()
      
      self.nrun = len(data) #number of runs
      
      self.nx = data[0]['meta']['Ni0']['v'].shape[0]
      self.ny = data[0]['meta']['Ni0']['v'].shape[1]
      self.nz = data[0]['meta']['MZ']['v']-1
      self.nt = int(data[0]['meta']['NOUT']['v']+1)
      self.Rxy = data[0]['meta']['Rxy']['v']

      

      self.meta = data[0]['meta']
     # a.raw[0]['meta']['Rxy']
# np.array(ListDictKey(alldb,'Rxy'))

      #self.dt = [data[i]['meta']['dt']['v'] for i in range(self.nrun)]
    

      #self.dt = alldb

      #self.nmodes = [len(data[i]['fields']['Ni']['modes']) for i in range(self.nrun)]
     
      self.dt = np.array(ListDictKey(alldb,'dt')) 
      self.nfields = len(data[0]['fields'])

      self.field = np.array(ListDictKey(alldb,'field'))
      self.k = np.array(ListDictKey(alldb,'k'))
      self.mn = np.array(ListDictKey(alldb,'mn'))
      self.phase = np.array(ListDictKey(alldb,'phase'))
      self.amp= np.array(ListDictKey(alldb,'amp'))
      self.amp_n=np.array(ListDictKey(alldb,'amp_n'))
      self.dc= []
      self.freq = np.array(ListDictKey(alldb,'k'))
      self.gamma = np.array(ListDictKey(alldb,'gamma'))
      self.freq  = np.array(ListDictKey(alldb,'freq'))
      
      self.IC = np.array(ListDictKey(alldb,'IC'))
      self.dz = np.array(ListDictKey(alldb,'dz'))
      self.nmodes = self.dz.size

   
      self.MN = np.float32(self.mn) 
      self.MN[:,1] = self.mn[:,1]/self.dz
      self.cxx = read_cxx(boutcxx='2fluid.cxx') #for now hardwired to read q3_simp.cxx in the current dir

   def show2(self):
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      pp = PdfPages('output.pdf')
      
      from reportlab.lib.pagesizes import letter
      from reportlab.pdfgen import canvas
      from reportlab.lib.units import inch

      #dir(canvas)

      c = canvas.Canvas('report.pdf', pagesize=letter)
      width, height = letter

      #c.drawString(100,750,"Welcome to Reportlab!")
      #c.drawString(3*inch, -3*inch, "Hello World")
      for i,elem in enumerate(self.meta):
         print i
         c.drawString(100,750-12*i,
                      str(elem)+str(self.meta[elem]))
      c.showPage()
      
      c.save()

      self.printmeta(pp)

      fig = Figure(figsize=(6,6))
      canvas = FigureCanvas(fig)
      
      ax =fig.add_subplot(221)
      ax.set_title(' Ni spectrum at t=0, all x run 0',fontsize=10)
      ax.set_xlabel(u'k \u03c1',fontsize=10)
      ax.grid(True,linestyle='-',color='.75')
      ax.scatter(self.k[:,0,int(self.nx/2)],self.amp[:,0,int(self.nx/2)])
 
      ax =fig.add_subplot(2,2,2)
      ax.scatter(self.k[:,0,int(self.nx/2)],self.amp[:,0,int(self.nx/2)])
      ax.set_xscale('log')
      ax.set_title(' Ni spectrum at t=0, all x run 0')
      ax.set_xlabel(u'k \u03c1')
      fig.savefig(pp,format='pdf')

      plt.figure()
      plt.scatter(self.mn[:,0],self.mn[:,1],s = 5*self.mn.max()*
                  self.amp[:,0,int(self.nx/2)]/self.amp[:,0,int(self.nx/2)].mean())
      plt.annotate(str(list(self.mn[0,:])),tuple(self.mn[0,:]+.2))
      plt.title('n-m spectrum at t=0')
      plt.xlabel('n')
      plt.ylabel('m')
      plt.grid(True,linestyle='-',color='.75')
      plt.savefig(pp, format='pdf')

      # NM 2D spectrum
      plt.figure()
      plt.scatter(self.mn[:,0],self.mn[:,1],s = 5*self.mn.max()*
                  self.amp[:,15,int(self.nx/2)]/self.amp[:,15,int(self.nx/2)].mean())
      plt.annotate(str(list(self.mn[0,:])),tuple(self.mn[0,:]+.2))
      plt.title('n-m spectrum at t=T/2')
      plt.xlabel('n')
      plt.ylabel('m')
      plt.grid(True,linestyle='-',color='.75')
      plt.savefig(pp, format='pdf')

      #2D true NM spectrum with color code and boxes around spectral res regions log scale
      plt.figure()
      i = 0
      for j in list(set(self.dz).union()):    #looping over runs, over unique 'dz' key values
         s = subset(self.raw,self.db,'dz',[j])  #subset where dz = j
         plt.scatter(s.MN[:,1],s.MN[:,0],c=colors[i])
         i+=1
      plt.title(' Ni spectrum at t=0, all x')
      plt.ylabel('M -parallel')
      plt.xlabel('N -  axisymmteric')
      plt.xscale('log')
      plt.grid(True,linestyle='-',color='.75')
      plt.savefig(pp, format='pdf')

      self.plotmodes(pp)
      
      #pick anothe field
      all_fields = list(set(self.field).union())
      self.plotmodes(pp,field=all_fields[1])
      #self.plotmodes(pp,math='gamma',ylim=0,yscale='linear',
      #               clip = 5)
      self.plotmodes(pp,comp='gamma',yscale='linear',
                     xaxis = 'k',xscale='log')
 
      self.plotmodes(pp,comp='gamma',yscale='linear',
                     xaxis = 'k',xscale='log',xrange=5)
      self.plotgamma(pp,xscale='log')
      
      self.plotfreq(pp,xscale='log')

      self.plotradeigen(pp,yscale='linear')
      self.plotradeigen(pp,yscale='linear')
      #self.plotradeigen(pp,field='rho',yscale='log')
      
      #self.plotmodes(pp,field='Ni',comp='gamma')
      

      plt.close() 
      pp.close()

   def plotgamma(self,pp,field='Ni',yscale='log',clip=0,
                 xaxis='t',xscale='linear',xrange=1):
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      plt.figure()
      
      s = subset(self.raw,self.db,'field',[field]) #pick  field
      #s = subset(s.raw,s.db,'k',
      
      #xrange = s.nx/2-2
      
      xrange = [s.nx/2,s.nx/2+xrange]
                    
      gamma = s.gamma #nmodes x 2 x nx ndarray
      k = s.k ##nmodes x 2 x nx ndarray
      print xrange
      # plt.errorbar(k[:,1,xrange[0]:xrange[1]],
      #              gamma[:,0,xrange[0]:xrange[1]],
      #              yerr=gamma[:,1,xrange[0]:xrange[1]],
      #              fmt='b.')
      plt.errorbar(k[:,1,s.nx/2],
                   gamma[:,0,s.nx/2],
                   yerr=gamma[:,1,s.nx/2],
                   fmt='b.')
      plt.yscale('log')
      plt.xscale(xscale)
      plt.savefig(pp, format='pdf')

   def plotfreq(self,pp,field='Ni',clip=0,
                xaxis='t',xscale='linear',yscale='linear'):
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      plt.figure()
      
      s = subset(self.raw,self.db,'field',[field]) #pick  field
      #s = subset(s.raw,s.db,'k',

      gamma = s.freq #nmodes x 2 x nx ndarray
      k = s.k ##nmodes x 2 x nx ndarray

      plt.errorbar(k[:,1,s.nx/2],gamma[:,0,s.nx/2],
                   yerr=gamma[:,1,s.nx/2],fmt='b.')
      plt.yscale(yscale)
      plt.xscale(xscale)
      plt.savefig(pp, format='pdf')
   
   def plotmodes(self,pp,field='Ni',comp='amp',math='1',ylim=1,
                 yscale='log',clip=0,xaxis='t',xscale='linear',
                 xrange=1):

      Nplots = self.nrun
      
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      fig = Figure()
      plt.figure()
      canvas = FigureCanvas(fig)
      
      Modes = subset(self.raw,self.db,'field',[field]) #pick  field
      
      k = 0
      
      

      for j in list(set(Modes.dz).union()):  #
         s = subset(Modes.raw,Modes.db,'dz',[j]) #pick run
         xr = range(s.nx/2-xrange/2,s.nx/2+xrange/2+1)
         data = np.array(ListDictKey(s.db,comp)) #pick component   
         ax =fig.add_subplot(round(Nplots/3.0 + 1.0),3,k+1)  
         
         for i in range(s.nmodes):
            if math=='gamma':
               out = np.gradient(data[i,:,xr])[1]/data[i,:,xr]
            else:
               out = data[i,:,xr]
               
            print j,i
            #out = out[clip:,]
            
            if xaxis=='t':
               x = range(out.size)
               plt.scatter(x,out,c=colors[k]) 
               ax.scatter(x,out,c=colors[k])#,alpha = (1 +i)/s.nmodes) 
            else:
               x = np.array(ListDictKey(s.db,xaxis))[i,:,xr] 
               x #an N? by nx array
               #print x.shape, out.shape,x,out
               plt.scatter(x[:,1],out[:,0])#,c=colors[k]) 
               ax.scatter(x[:,1],out[:,0])#,c=colors[k])#,alpha = (1 +i)/s.nmodes) 

            #detect error (bar data
               print 'error bars:',x,out
   
             
         
         ax.set_yscale(yscale)
         ax.set_xscale(xscale)
         ax.set_title(str(j),fontsize=10)
         ax.set_xlabel(xaxis)

         if ylim:
            ax.set_ylim(data.min(),5*data.max())
         k+=1
            
      plt.title(' Ni amp, all runs',fontsize=10)
      plt.xlabel(xaxis)
      if ylim:
         plt.ylim(data.min(),10*data.max())
         
      plt.yscale(yscale)
      plt.xscale(xscale)

      fig.savefig(pp, format='pdf')  
      plt.savefig(pp, format='pdf')
            
      
      #except:
      #   print "Sorry you fail"
   

   def plotradeigen(self,pp,field='Ni',comp='amp',
                    yscale='linear',xscale = 'linear'):
      colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
      plt.figure()
      
      s = subset(self.raw,self.db,'field',[field])
      
      amp = s.amp


      for mode in range(s.nmodes): 
         print 'mode: ', mode
         plt.plot(s.Rxy[:,s.ny/2],s.amp[mode,s.nt-1,:])
      
      plt.yscale(yscale)
      plt.xscale(xscale)   

      plt.savefig(pp, format='pdf')        
   
   def printmeta(self,pp):

      PAGE_HEIGHT=defaultPageSize[1]
      styles = getSampleStyleSheet()
      Title = "BOUT++ Results"
      Author = "Dmitry Meyerson"
      URL = ""
      email = "dmitry.meyerson@gmail.com"
      Abstract = """This document highlights some results from BOUT++ simulation"""
      Elements=[]
      HeaderStyle = styles["Heading1"]
      ParaStyle = styles["Normal"]
      PreStyle = styles["Code"]
 
      def header(txt, style=HeaderStyle, klass=Paragraph, sep=0.3):
         s = Spacer(0.2*inch, sep*inch)
         para = klass(txt, style)
         sect = [s, para]
         result = KeepTogether(sect)
         return result
 
      def p(txt):
         return header(txt, style=ParaStyle, sep=0.1)
 
      def pre(txt):
         s = Spacer(0.1*inch, 0.1*inch)
         p = Preformatted(txt, PreStyle)
         precomps = [s,p]
         result = KeepTogether(precomps)
         return result
      
      def graphout(name, datain):
         xaxis=range(datain.size)
         drawing = Drawing(400, 200)
         # data = [
         #    ((1,1), (2,2), (2.5,1), (3,3), (4,5)),
         #    ((1,2), (2,3), (2.5,2), (3.5,5), (4,6))
         #    ]
         dataview = [tuple([(xaxis[i],datain[i]) for i in range(datain.size)])]
         lp = LinePlot()
         lp.x = 50
         lp.y = 50
         lp.height = 125
         lp.width = 300
         lp.data = dataview
         # lp.joinedLines = 1
         # lp.lines[0].symbol = makeMarker('FilledCircle')
         # lp.lines[1].symbol = makeMarker('Circle')
         # lp.lineLabelFormat = '%2.0f'
         # lp.strokeColor = colors.black
         # lp.xValueAxis.valueMin = min(xaxis)
         # lp.xValueAxis.valueMax = max(xaxis)
         # lp.xValueAxis.valueSteps = xaxis
         # lp.xValueAxis.labelTextFormat = '%2.1f'
         # lp.yValueAxis.valueMin = min(datain)
         # lp.yValueAxis.valueMax = max(datain)
         # lp.yValueAxis.valueSteps = [1, 2, 3, 5, 6]
         drawing.add(lp)
         return drawing
 
      def go():
         doc = SimpleDocTemplate('gfe.pdf')
         doc.build(Elements)
 
      mytitle = header(Title)
      myname = header(Author, sep=0.1, style=ParaStyle)
      mysite = header(URL, sep=0.1, style=ParaStyle)
      mymail = header(email, sep=0.1, style=ParaStyle)
      abstract_title = header("ABSTRACT")
      myabstract = p(Abstract)
      head_info = [mytitle, myname, mysite, mymail, abstract_title, myabstract]
      Elements.extend(head_info)
 
      meta_title = header("metadata")
      metasection = []
      metasection.append(meta_title)
      for i,elem in enumerate(self.meta):
         if type(self.meta[elem])== type({}):
            print elem #np.array(self.meta[elem]['v']).shape() 
            if np.array(self.meta[elem]['v']).shape == (self.nx,self.ny):
               datastr = str(self.meta[elem]['v'][:,self.ny/2])
               metasection.append(graphout('stuff',
                                           self.meta[elem]['v'][:,self.ny/2]))
            else:
               datastr = str(self.meta[elem]['v'])
            metasection.append(header(str(elem)+': '+datastr
                                         + ' '+ str(self.meta[elem]['u']), 
                                         sep=0.1, style=ParaStyle))
          
      
     
      src = KeepTogether(metasection)
      Elements.append(src)

      cxxtitle = header("Equations in CXX")
      cxxsection = []
      print self.cxx
      cxxsection.append(header(self.cxx, sep=0.1, style=ParaStyle))
      cxxsrc = KeepTogether(cxxsection)

      Elements.append(cxxsrc)
      # for i,elem in enumerate(self.cxx):
      #    if type(self.meta[elem])== type({}):
      #       print elem #np.array(self.meta[elem]['v']).shape() 
      #       if np.array(self.meta[elem]['v']).shape == (self.nx,self.ny):
      #          datastr = str(self.meta[elem]['v'][:,self.ny/2])
      #          metasection.append(graphout('stuff',
      #                                      self.meta[elem]['v'][:,self.ny/2]))
      #       else:
      #          datastr = str(self.meta[elem]['v'])
      #       metasection.append(header(str(elem)+': '+datastr
      #                                    + ' '+ str(self.meta[elem]['u']), 
      #                                    sep=0.1, style=ParaStyle))

      
      go()

class subset(LinRes):
   def __init__(self,data,alldb,key,valuelist):
      LinRes.__init__(self,data,ListDictFilt(alldb,key,valuelist))
      self.skey = key

def get_image(filename=''):
    """ Get a python logo image for this example """
    if not os.path.exists(filename):
        response = urllib2.urlopen(
            'http://www.python.org/community/logos/python-logo.png')
        f = open(filename, 'w')
        f.write(response.read())
        f.close()
