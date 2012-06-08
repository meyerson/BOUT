from pb_draw import LinResDraw,subset
from pb_corral import LinRes
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
#from matplotlib.ticker import FuncFormatter
#from matplotlib.ticker import ScalarFormatter 

from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.lib import colors

from replab_x_vs_y import RL_Plot


#uses  LinResDraw to make a pdf

class LinResPresent(LinResDraw):
   def __init__(self,alldb):
      LinResDraw.__init__(self,alldb)

   def show(self,filter =True,quick=False,pdfname='output2.pdf'):
       colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
       pp = PdfPages('output.pdf')
       
      #start by removing modes above the maxN threshold
       modelist =[]
       [modelist.append(list(self.modeid[p])) for p in range(self.nmodes) if self.mn[p][1] <= self.maxN[p] ]
       
       print modelist
       s = subset(self.db,'modeid',modelist)
      #print s.mn
      
       fig = Figure(figsize=(6,6))
       
      #a_amp = [s.amp[:,1,int(s.nx/2)]
       
       a_amp = np.array([s.amp[i][1,int(s.nx/2)] for i in range(s.nmodes)])
       
       
       plt.scatter(s.mn[:,1],s.mn[:,0],s = 5*s.mn.max()*
                  s._amp(1,int(s.nx/2))/s._amp(1,int(s.nx/2)).mean())
       plt.annotate(str(list(s.mn[0,:])),tuple(s.mn[0,:]+.2))
       plt.title('n-m spectrum at t=0')
       plt.xlabel('n')
       plt.ylabel('m')
       plt.grid(True,linestyle='-',color='.75')
       plt.savefig(pp, format='pdf')
       plt.close() 

      #2D true NM spectrum with color code and boxes around spectral res regions log scale
       plt.figure()
       i = 0
       for j in list(set(s.dz).union()):    #looping over runs, over unique 'dz' key values
           
           ss = subset(s.db,'dz',[j])  #subset where dz = j
           plt.scatter(ss.MN[:,1],ss.MN[:,0],c=colors[i])
           plt.annotate(str(j),(ss.MN[0,1],ss.MN[0,0]))
           i+=1
         
       plt.title(' Ni spectrum at t=0, all x')
       plt.ylabel('M -parallel')
       plt.xlabel('N -  axisymmteric')
       plt.xscale('log')
       plt.grid(True,linestyle='-',color='.75')
     
       plt.savefig(pp, format='pdf')
       plt.close() 
      

       s.plotomega(pp,xscale='linear')
       s.plotmodes(pp,yscale='linear',comp='phase',linestyle='.')
       
       
       s.plotomega(pp,xscale='linear',field='rho')
       s.plotmodes(pp,yscale='linear',comp='phase',linestyle='.',field='rho')
 
       s.plotmodes(pp,yscale='linear')
       s.plotmodes(pp,yscale='linear',field='rho')


       s.plotmodes(pp,yscale='log',debug=True)
       
       modelist = []
       [modelist.append([1,p+1]) for p in range(7)]
        
       ss = subset(s.db,'mn',modelist)
       ss.plotmodes(pp,yscale='log',debug=True)
       ss.plotmodes(pp,yscale='linear',comp='phase')
      #ss.plotmodes(pp,yscale='log',comp='phase',clip=True)
       ss.plotmodes(pp,comp='phase',clip=True)


       if quick==True:
           pp.close()
           s.printmeta(pp)
         #plt.savefig(pp, format='pdf')
           return 0
      
       all_fields = list(set(s.field).union())


       s.plotomega(pp,xscale='log',overplot=True)
       s.plotgamma(pp,xscale='log',overplot=True)
  
       s.plotgamma(pp,yscale='log',xscale='log',overplot=True)
       s.plotgamma(pp,yscale='log',xscale='log',field='rho')
      
       s.plotfreq2(pp,xscale='log',yscale='linear',overplot=True)
       s.plotfreq2(pp,xscale='log',yscale='symlog',overplot=True)
       s.plotfreq2(pp,xscale='log',yscale='linear',field='rho')
       
       s.plotfreq2(pp,xscale='log',yscale='symlog')

       s.plotradeigen(pp,yscale='linear')
       s.plotradeigen(pp,field ='Vi',yscale='linear')
       s.plotradeigen(pp,field='rho',yscale='log')
  
       pp.close()
       s.printmeta(pp,filename=pdfname) #append a metadata header
