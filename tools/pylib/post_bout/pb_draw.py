#some standard analytic stuff to plot, if appending just overplot gam or omeg
from pb_corral import LinRes
from ListDict import ListDictKey, ListDictFilt
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
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

from replab_x_vs_y import RL_Plot
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator


class LinResDraw(LinRes):
    def __init__(self,alldb):
        LinRes.__init__(self,alldb)
       
        
    def plottheory(self,pp,m=1,canvas=None,comp='gamma',field='Ni'):
        if self.M is None:
            self.model()
      
        s = subset(self.db,'field',[field])
        
        modelist = []
        [modelist.append([m,n+1]) for n in range(7)]
   
        s = subset(s.db,'mn',modelist)
         
        allk = s.k[:,1,s.nx/2]
        ki = np.argsort(allk)
    
        ownpage = False
        if canvas is None:
            ownpage = True
        

        if ownpage: #if not an overplot
            fig1 = plt.figure()
            canvas = fig1.add_subplot(1,1,1) 
         
        label ='gamma analytic'

        if comp=='gamma':
            y = np.array(s.gammamax)[ki]
        else:
            y = np.array(s.omegamax)[ki]
            
        canvas.plot(allk[ki],y,'-',
                    label=label)
             
        canvas.annotate('theory',(allk[ki[0]],y[0]),fontsize = 8)

        if ownpage: #set scales if this is its own plot
            fig1.savefig(pp, format='pdf')
            plt.close(fig1)

    def plotomega(self,pp,canvas=None,field='Ni',yscale='linear',clip=0,
                  xaxis='t',xscale='linear',xrange=1,comp='gamma',
                  pltlegend='both',overplot=False,gridON=True,trans=True):
        colors = ['b.','r.','k.','c.','g.','y.','m.','b.','r.','k.','c.','g.','y.','m.']
        colordash = ['b','r','k','c','g','y','m','b','r','k','c','g','y','m']

        
        if canvas is None:
            ownpage = True
        else:
            ownpage = False
        
        if ownpage:    
            fig1 = plt.figure()
            fig1.subplots_adjust(bottom=0.12)
            fig1.subplots_adjust(top=0.80)
            fig1.subplots_adjust(right=0.83)
            fig1.subplots_adjust(left=0.17)
            canvas = fig1.add_subplot(1,1,1) 
            clonex = canvas.twinx()
            if np.any(self.trans) and trans:
                cloney = canvas.twiny()

       
        dzhandles = []
        parhandles =[]
        parlabels =[]
        dzlabels=[]
        
        m_shift = 1
        for q in np.array(range(1))+m_shift:
            s = subset(self.db,'field',[field]) #pick  field
            modelist = []
            [modelist.append([q,p+1]) for p in range(7)]
        print q,'in plotgamma'
        s = subset(s.db,'mn',modelist)
         
        xrange = s.nx/2-2
        
        xrange = [s.nx/2,s.nx/2+xrange]
      
 
        y = np.array(ListDictKey(s.db,comp))

         #y = s.gamma #nmodes x 2 x nx ndarray
        k = s.k ##nmodes x 2 x nx ndarray

         # plt.errorbar(k[:,1,s.nx/2],
         #              y[:,0,s.nx/2],
         #              yerr=y[:,2,s.nx/2],
         #              fmt=colors[q])
        parhandles.append(canvas.errorbar(k[:,1,s.nx/2],
                                       y[:,0,s.nx/2],
                                       yerr=y[:,1,s.nx/2],
                                       fmt=colors[q]))
        parlabels.append("m "+str(q))
         
         #loop over dz sets and connect with dotted line  . . .
        jj=0

        ymin_data =np.max(np.array(ListDictKey(s.db,comp)))
        ymax_data = 0 #for bookeeping

        for p in list(set(s.path).union()):
            print p, 'in plotomega'
            
            sub_s = subset(s.db,'path',[p])
            j = sub_s.dz[0]
            #print sub_s.amp.shape
            s_i = np.argsort(sub_s.mn[:,1]) #sort by 'local' m, global m is ok also
            #print s_i, sub_s.mn, sub_s.nx, jj
            y = np.array(ListDictKey(sub_s.db,comp))
            y_alt = 2.0*np.array(ListDictKey(sub_s.db,comp))

            k = sub_s.k ##
            if q == m_shift: #fix the parallel mode
               dzhandles.append(canvas.plot(k[s_i,1,sub_s.nx/2],
                                            y[s_i,0,sub_s.nx/2],color=colordash[jj],alpha=.5))
               #clonex.plot(k[s_i,1,sub_s.nx/2],
                           #y_alt[s_i,0,sub_s.nx/2],color=colordash[jj],alpha=.5)
               ymin_data = np.min([np.min(y[s_i,0,sub_s.nx/2]),ymin_data])
               ymax_data = np.max([np.max(y[s_i,0,sub_s.nx/2]),ymax_data])
               

               if np.any(sub_s.trans) and trans:
                   comp_r = comp+'_r'
                   y2 = np.array(ListDictKey(sub_s.db,comp_r))
                   #canvas.plot(k[s_i,1,sub_s.nx/2],
                    #           y2[s_i,0,sub_s.nx/2],'k.',ms = 3)
                   
                   cloney.plot(k[s_i,1,sub_s.nx/2],
                               y2[s_i,0,sub_s.nx/2],'k.',ms = 3)
                   
               print 'dzhandle color', jj
               #dzlabels.append("DZ: "+ str(2*j)+r'$\pi$')
               dzlabels.append(j) 
               
               if yscale=='log':
                  factor = 10
               else:
                  factor = 2
               print 'annotating'
               canvas.annotate(str(j),(k[s_i[0],1,sub_s.nx/2],y[s_i[0],0,sub_s.nx/2]),fontsize = 8)
               p = canvas.axvspan(k[s_i[0],1,sub_s.nx/2], k[s_i[-1],1,sub_s.nx/2], 
                               facecolor=colordash[jj], alpha=0.01)
               print 'done annotating'
            else:
               canvas.plot(k[s_i,1,sub_s.nx/2],
                        y[s_i,0,sub_s.nx/2],color=colordash[jj],alpha=.3)

            jj=jj+1
      
        dzhandles = np.array(dzhandles).flatten()
        dzlabels = np.array(dzlabels).flatten()
      
        dzlabels = list(set(dzlabels).union())
      
        dz_i = np.argsort(dzlabels)
        
        dzhandles = dzhandles[dz_i]
        dzlabels_cp = np.array(dzlabels)[dz_i]
      
        print type(dzlabels), np.size(dzlabels)
        for i in range(np.size(dzlabels)):
            dzlabels[i] = "DZ: "+ str(dzlabels_cp[i])#+r"$\pi$"
     

        parlabels = np.array(parlabels).flatten()
      
      #if pltlegend =='both':    #
      
        print 'legends'

      #l1 = legend(parhandles,parlabels,loc = 3,prop={'size':6})
      #l2 = legend(dzhandles,dzlabels,loc = 1,prop={'size':6})
      #plt.gca().add_artist(l1)
      
      # else:
      #    legend(dzhandles,dzlabels,loc=3,prop={'size':6})
        if overplot==True:
            self.plottheory(pp,canvas=canvas,comp=comp)
        
                
        #cloney.set_xlim(xmin,xmax)
        try:
            canvas.set_yscale(yscale)
            canvas.set_xscale(xscale)
            
            if yscale =='symlog':
                canvas.set_yscale(yscale,linthreshy=1e-9)
            if xscale =='symlog':
                canvas.set_xscale(xscale,linthreshy=1e-9)
            
            if gridON:
                canvas.grid()
        except:
            try:
                canvas.set_yscale('symlog')
            except:
                print 'scaling failed completely'
   
        [xmin, xmax, ymin, ymax] = canvas.axis()    

        try:
            #canvas.set_yscale(yscale)
            #canvas.grid(axis='x')
            #canvas.grid(axis='y')
        
            #clonex.plot([xmin,xmax],[20*ymin_data,20*ymax_data],alpha=0.001)
            #clonex.set_ylim(2*ymin,2*ymax)
            #print '[xmin, xmax, ymin, ymax] ',[xmin, xmax, ymin, ymax]
            #cloney.set_yscale(yscale)
           
            clonex.set_yscale(yscale) #must be called before limits are set 
              
            if np.any(s.trans) and trans:
                cloney.set_yscale(yscale)
                cloney.set_xscale(xscale)
            
            if yscale =='symlog':
                clonex.set_yscale(yscale,linthreshy=1e-9)
                  
                
            if xscale =='symlog' and np.any(s.trans) and trans:
                cloney.set_yscale(yscale,linthreshy=1e-9)
                cloney.set_xscale(xscale,linthreshy=1e-9)
            
            Ln_drive_scale = s.meta['w_Ln'][0]**-1
            #Ln_drive_scale = 2.1e3
            clonex.set_ylim(Ln_drive_scale*ymin, Ln_drive_scale*ymax)
            if np.any(s.trans) and trans:
                cloney.set_xlim(xmin, xmax)
                cloney.set_ylim(ymin,ymax) #because cloney shares the yaxis with canvas it may overide them, this fixes that
                cloney.set_xlabel(r'$k_{\perp} \rho_{ci}$',fontsize=18)
            #clonex.set_xscale(xscale)
        except:
            #canvas.set_xscale('symlog', linthreshx=0.1)  
            print 'extra axis FAIL'


        canvas.set_ylabel(r'$\frac{\gamma}{\omega_{ci}}$',fontsize=18,rotation='horizontal')
        #if yscale == 'linear':
        #canvas.yaxis.set_major_locator(ticker.LinearLocator(numticks=8))
    
        
        # minorLocator   = MultipleLocator(.005)
        # canvas.yaxis.set_minor_locator(minorLocator)
        #spawn another y label
        
       

        #clone = canvas.twinx()
        #s2 = np.sin(2*np.pi*t)
        #ax2.plot(x, s2, 'r.')
        
        
        #ion_acoust_str = r"$\frac{c_s}{L_{\partial_r n}}}$"
        
        clonex.set_ylabel(r'$\frac{\gamma}{\frac{c_s}{L_n}}$', color='k',fontsize=18,rotation='horizontal')
        canvas.set_xlabel(r'$k_{\zeta} \rho_{ci}$',fontsize=18)
        
        # for tl in clone.get_yticklabels():
        #     tl.set_color('r')

      #title = r'$\'+comp+'$ '+ 'computed from '+field'
      #title = 'r\$\'
      #title = r'$\frac{Ni}{Ni_0}$'
      #title = '\\$'+comp
      #title = title+'$ '+ ' computed from '+field
      #def totex(input):
         #return str(input)
        
        #cloney = canvas.twiny()
        

        
      
        if yscale == 'linear':
            formatter = ticker.ScalarFormatter()
            formatter.set_powerlimits((0, 0))  #force scientific notation
            canvas.yaxis.set_major_formatter(formatter)
            clonex.yaxis.set_major_formatter(formatter)
        
        # if self.trans and overplot:
        #     self.plotomega(pp,canvas=canvas,overplot=False,comp='gamma_r')
        
        title = comp+ ' computed from '+field
        #canvas.set_title(title,fontsize=14)
        fig1.suptitle(title,fontsize=14)
        
        if ownpage:    
            fig1.savefig(pp, format='pdf')
            plt.close(fig1)

    def plotfreq(self,pp,field='Ni',clip=0,
                 xaxis='t',xscale='linear',yscale='linear'):
      #colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        colors = ['b.','g.','r.','c.','m.','y.','k.','b.','g.','r.','c.','m.','y','k']
        plt.figure()
        
      #s = subset(self.db,'field',[field]) #pick  field

      
        
        for q in range(4):
            s = subset(self.db,'field',[field]) #pick  field across all dz sets
            modelist = []
            [modelist.append([q+1,p+1]) for p in range(5)]
            print q,'in plotgamma'
            s = subset(s.db,'mn',modelist)
      
      
            gamma = s.freq #nmodes x 2 x nx ndarray
            k = s.k ##nmodes x 2 x nx ndarray
         
            plt.errorbar(k[:,1,s.nx/2],gamma[:,0,s.nx/2],
                         yerr=gamma[:,1,s.nx/2],fmt=colors[q])
            plt.plot(k[:,1,s.nx/2],
                     gamma[:,0,s.nx/2],'k:',alpha=.3)
         
#loop over dz sets and connect with dotted line  . . .
            for j in list(set(s.dz).union()):
           # print j,len(s.mn)
                sub_s = subset(s.db,'dz',[j]) 
                gamma = sub_s.gamma 
                k = sub_s.k ##
                plt.plot(k[:,1,sub_s.nx/2],
                         gamma[:,0,sub_s.nx/2],'k:',alpha=.1)

        try:
            plt.yscale(yscale)
        except:
            print 'yscale fail'

       
        try:
            plt.xscale(yscale)
        except:
            plt.xscale('symlog')
        plt.xlabel(r'$k \rho_{ci}$',fontsize=14)
        plt.ylabel(r'$\frac{\omega}{\omega_{ci}}$',fontsize=14)
      #plt.title(r'$\frac{\omega}\{\omega_{ci}}$ '+ 'computed from'+field+ 'field',fontsize=10)
      
        plt.savefig(pp, format='pdf')
        plt.close()
      
    def plotgamma(self,pp,field='Ni',yscale='symlog',clip=0,
                  xaxis='t',xscale='linear',xrange=1,comp='gamma',
                  overplot=False,trans=True):
        self.plotomega(pp,field=field,yscale=yscale,clip=clip,
                     xaxis=xaxis,xscale=xscale,xrange=xrange,comp=comp,
                     overplot=overplot,trans=trans)
   
    def plotfreq2(self,pp,field='Ni',yscale='symlog',clip=0,
                  xaxis='t',xscale='linear',xrange=1,comp='freq',
                  overplot=False,trans=True):
        self.plotomega(pp,field=field,yscale=yscale,clip=clip,
                       xaxis=xaxis,xscale=xscale,xrange=xrange,comp=comp,
                       overplot=overplot,trans=trans)

    def plotmodes(self,pp,field='Ni',comp='amp',math='1',ylim=1,
                  yscale='symlog',clip=False,xaxis='t',xscale='linear',
                  xrange=1,debug=False,yaxis=r'$\frac{Ni}{Ni_0}$',
                  linestyle='-'):
      
        Nplots = self.nrun

        colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        fig1 = plt.figure()
        
        fig2 = plt.figure()
   

        Modes = subset(self.db,'field',[field]) #pick  field
     
  
        adj = fig2.subplots_adjust(hspace=0.4,wspace=0.4)
        fig2.suptitle('Dominant mode '+ comp+' for  '+ field)
        props = dict( alpha=0.8, edgecolors='none' )
        
        allcurves = fig1.add_subplot(1,1,1)
        fig1.suptitle('Dominant mode behavior for  '+ field)
        
        modenames = []
        k = 0
        
        for j in list(set(Modes.path).union()):  #
            s = subset(Modes.db,'path',[j]) #pick run
            dz = s.dz[0]
            xr = range(s.nx/2-xrange/2,s.nx/2+xrange/2+1)
            data = np.array(ListDictKey(s.db,comp)) #pick component should be ok for a fixed dz key
         
            data = data + 1e-32
            ax =fig2.add_subplot(round(Nplots/3.0 + 1.0),3,k+1)  
            ax.grid(True,linestyle='-',color='.75')
            handles=[]
         #modenames.append(str(j))

          #find the "biggest" mode for this dz
            d = data[:,s.nt[0]-1,:] #nmode X nx array
         #d = s.gamma[:,2,:]
            where = d == d.max()
            z = where.nonzero() #mode index and n index
            imax = z[0][0]
            #xi_max = z[1][0]
            xi_max = s.nx/2

            if debug and yscale=='log':
                gamma = np.array(ListDictKey(s.db,'gamma'))  #nmodes x 2 x nx
                
            for i in range(s.nmodes):
                if math=='gamma':
                    out = np.gradient(data[i,:,xr])[1]/data[i,:,xr]
                else:
                    out = data[i,2:,xi_max] #skip the first 2 points
                    
                if xaxis=='t':
               #print 'out.size', out.size, out.shape
                    x = range(out.size)
               #plt.plot(x,out.flatten(),c=colors[k])
                    handles.append(ax.plot(x,out.flatten(),
                                     c=cm.jet(1.*k))) 
       
                else:
                    x = np.array(ListDictKey(s.db,xaxis))[i,:,xr] 
               #x #an N? by nx array
                    print x[:,1], out[:,0]
                    plt.scatter(x[:,1],out[:,0])#,c=colors[k]) 
                    ax.scatter(x[:,1],out[:,0])#,c=colors[k])#,alpha = (1 +i)/s.nmodes) 
               
            #detect error (bar data
                    print 'error bars:',x,out
   
            formatter = ticker.ScalarFormatter()
            formatter.set_powerlimits((0, 0)) 
            ax.xaxis.set_major_formatter(formatter)        
         #ax.axis('tight')
            if yscale=='linear':
                ax.yaxis.set_major_formatter(formatter)
            else:
                try:
                    ax.set_yscale(yscale)  
                except:
                    print 'may get weird axis'
                    ax.set_yscale('symlog')

            artist.setp(ax.axes.get_xticklabels(), fontsize=6)
            artist.setp(ax.axes.get_yticklabels(), fontsize=8)
         #artist.setp(ax.axes.get_yscale(), fontsize=8)     
            ax.set_title(str(dz),fontsize=10)
            ax.set_xlabel(xaxis)
   
         #x = s.Rxy[imax,:,s.ny/2]
         
            t0=2
            if clip == True: 
                t0 = round(s.nt[0]/3)
            y = data[imax,t0:,xi_max]
            x = np.array(range(y.size))
         
        
            print imax,xi_max
         
         #label = str([round(elem,3) for elem in s.MN[imax]])+ str(s.mn[imax])+' at x= '+str(xi_max)+' ,'+str(round(s.gamma[imax][2,xi_max]/s.gamma[imax][0,xi_max],3))+'%  '+str(round(s.gamma[imax][0,xi_max],4))

            label = str([round(elem,3) for elem in s.MN[imax]])+ str(s.mn[imax])+' at x= '+str(xi_max)+' ,'+str(round(s.gamma[imax,2,xi_max]/s.gamma[imax,0,xi_max],3))+'%  '+str(round(s.gamma[imax,0,xi_max],4))

            short_label = str(dz)
            print short_label
            allcurves.plot(x,y,'.',c=cm.jet(1.*k/len(x)),
                           label=label) 
         #print len(x), k*len(x)/(Nplots+2),s.nrun
            allcurves.annotate(short_label,(x[k*len(x)/(Nplots+1)],
                                            y[k*len(x)/(Nplots+1)]),fontsize = 8)
            
            
         #modenames.append(str([round(elem,3) for elem in s.MN[imax]])
                  #        +str(s.mn[imax])+' at x= '+str(xi_max)+' ,'+str(s.gamma[imax,2,xi_max]))
         
            if debug and yscale=='log':
                gam = gamma[imax,0,xi_max]
                f0 = gamma[imax,1,xi_max]
                allcurves.plot(x,f0*np.exp(gam*s.dt[imax]*x),'k:')

          
            k+=1
      
      # if ylim:
      #    allcurves.set_ylim(data[,xi_max].min(),5*data[:,xi_max].max())
        
        fig2.savefig(pp, format='pdf')   
        
        handles, labels = allcurves.get_legend_handles_labels()
        allcurves.legend(handles,labels,loc='best',prop={'size':6}) 
#allcurves.legend(modenames,loc='best',prop={'size':6})     
        allcurves.set_title(field+ ' '+comp+', all runs, '+yscale+' yscale',fontsize=10)
        allcurves.set_ylabel(yaxis)
        allcurves.set_xlabel(xaxis)
     
         
        if yscale=='linear':
            allcurves.yaxis.set_major_formatter(formatter)
        else:
            try:
                allcurves.set_yscale(yscale) 
            except:
                print 'may get weird axis scaling'
            if yscale=='log':
                allcurves.axis('tight')
            #allcurves.set_ylim(data.min(),data.max())
      #allcurves.set_yscale(yscale,nonposy='mask')

      #plt.xscale(xscale)
      
      
      #plt.legend(modenames,loc='best')

        fig1.savefig(pp, format='pdf')  
        plt.close(fig1)
        plt.close(fig2)
    
      #except:
      #   print "Sorry you fail"
   

    def plotradeigen(self,pp,field='Ni',comp='amp',
                     yscale='linear',xscale = 'linear'):
        
        Nplots  = self.nrun
        colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        fig1 = plt.figure()
        
        fig2 = plt.figure()
        adj = fig2.subplots_adjust(hspace=0.4,wspace=0.4)
        
      #canvas = FigureCanvas(fig)
        
        Modes = subset(self.db,'field',[field])
        
        k = 0
        fig2.suptitle('Dominant mode behavior for  '+ field)
        props = dict( alpha=0.8, edgecolors='none' )    
        
        allcurves = fig1.add_subplot(1,1,1)
        fig1.suptitle('Dominant mode behavior for  '+ field)
        
        modeleg = []


        for p in list(set(Modes.path).union()):
            print p
            s = subset(Modes.db,'path',[p]) #pick run
         #data = np.array(ListDictKey(s.db,comp)) #pick component   
            j =  s.dz[0]
            ax =fig2.add_subplot(round(Nplots/3.0 + 1.0),3,k+1)  
            ax.grid(True,linestyle='-',color='.75')
            data = np.array(ListDictKey(s.db,comp)) #pick component 
            handles=[]
            
         #find the "biggest" mode for this dz
            d = data[:,s.nt[0]-1,:] #nmode X nx array
            where = d == d.max()
            z = where.nonzero() #mode index and n index
            imax = z[0][0]
            modeleg.append(str([round(elem,3) for elem in s.MN[imax]])+str(s.mn[imax]))
         
   #str(s.k[:,1,:][z])
         

            for i in range(s.nmodes): 
            #out = mode[mode.ny/2,:]
            #print i,s.Rxynorm.shape,s.ny
                x = s.Rxynorm[i,:,s.ny/2]
                y = data[i,s.nt[0]-1,:]
            
                handles.append(ax.plot(x,y,c=cm.jet(1.*k/len(x))))

            formatter = ticker.ScalarFormatter()
            formatter.set_powerlimits((0, 0)) 
            ax.xaxis.set_major_formatter(formatter)        

            if yscale=='linear':
                ax.yaxis.set_major_formatter(formatter)
            else:
                ax.set_yscale(yscale)
            
            artist.setp(ax.axes.get_xticklabels(), fontsize=6)
            artist.setp(ax.axes.get_yticklabels(), fontsize=8)
         #artist.setp(ax.axes.get_yscale(), fontsize=8)
            ax.set_title(str(j),fontsize=10)
         
         
            x = s.Rxynorm[imax,:,s.ny/2]
            y = data[imax,s.nt[0]-1,:]
         #allcurves.plot(x,y,c= colors[k])
            allcurves.plot(x,y,c=cm.jet(1.*k/len(x)))
            k= k+1
         
        fig2.savefig(pp, format='pdf')  
        if yscale=='linear':
            allcurves.yaxis.set_major_formatter(formatter)
        else:
            allcurves.set_yscale(yscale)
        try:
            allcurves.set_xscale(xscale)
        except:
            allcurves.set_xscale('symlog')

      #allcurves.xaxis.set_major_formatter(ticker.NullFormatter())
        allcurves.legend(modeleg,loc='best',prop={'size':6})
        allcurves.set_xlabel(r'$\frac{x}{\rho_{ci}}$')
        allcurves.set_ylabel(r'$\frac{Ni}{Ni_0}$')
        fig1.savefig(pp, format='pdf')        
        plt.close(fig1)
        plt.close(fig2)
 
               
    def plotmodes2(self,pp,field='Ni',comp='amp',math='1',ylim=1,
                  yscale='symlog',clip=0,xaxis='t',xscale='linear',
                  xrange=1,debug=False):
        
       Nplots = self.nrun
       Modes = subset(self.db,'field',[field]) #pick  field
       colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
       
       fig = Figure()
       plt.figure()
       
       canvas = FigureCanvas(fig)
       k = 0
       nrow = round(Nplots/3.0 + 1.0)
       ncol = 3
      #nrow = round(Nplots/3.0 + 1.0)
      #ncol = round(Nplots/3.0 + 1.0)
       f, axarr = plt.subplots(int(nrow),int(ncol))
       
       for p in list(set(Modes.path).union()):  # 
           s = subset(Modes.db,'path',[p]) #pick run
           j= s.dz[0]
           xr = range(s.nx/2-xrange/2,s.nx/2+xrange/2+1)
           data = np.array(ListDictKey(s.db,comp)) #pick component   
         #ax =fig.add_subplot(round(Nplots/3.0 + 1.0),3,k+1)  
         
           for i in range(s.nmodes):
               out = data[i,:,xr]
               
           print j,i
           if xaxis=='t':   
               x = range(out.size)
               #plt.scatter(x,out.flatten(),c=colors[k]) 
               plt.scatter(x,out.flatten(), 
                           c=cm.jet(1.*k/len(x))) 
               #axarr[j%(ncol),j/ncol].scatter(x,out.flatten(),c=colors[k])#,alpha = (1 +i)/s.nmodes) 
               axarr[j/ncol,j%(ncol)].scatter(x,out.flatten(),c=cm.jet(1.*k/len(x)))#
               
           else:
               x = np.array(ListDictKey(s.db,xaxis))[i,:,xr] 
               #x #an N? by nx array
               print x[:,1], out[:,0]
               plt.scatter(x[:,1],out[:,0])#,c=colors[k]) 
               axarr[j%(col),j/col].scatter(x[:,1],out[:,0])#,c=colors[k])#,alpha = (1 +i)/s.nmodes) 

            #detect error (bar data
               print 'error bars:',x,out
   
             
               
           axarr[j/ncol,j%(ncol)].set_yscale(yscale)
           axarr[j/ncol,j%(ncol)].set_xscale(xscale)
           axarr[j/ncol,j%(ncol)].set_title(str(j),fontsize=10)
           axarr[j/ncol,j%(ncol)].set_xlabel(xaxis)
         
           plt.setp([a.get_xticklabels() for a in axarr[0,:]], visible=False)
           plt.setp([a.get_yticklabels() for a in axarr[:,ncol-1]], visible=False)

           if ylim:
               axarr[j%(ncol),j/ncol].set_ylim(data.min(),5*data.max())
           k+=1
            
       plt.title(field+ ' '+comp+', all runs, '+yscale+' yscale',fontsize=10)
       plt.xlabel(xaxis)
       if ylim:
           plt.ylim(data.min(),10*data.max())
         
      
       plt.yscale(yscale,nonposy='mask')
       plt.xscale(xscale)

       fig.savefig(pp, format='pdf')  
       plt.savefig(pp, format='pdf')
            
       plt.close()

       return 0
          
    def plotMacroDep(self,pp,field='Ni',yscale='symlog',clip=0,
                     xaxis='t',xscale='linear',xrange=1):
        colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        plt.figure()
   
    def printmeta(self,pp,filename='output2.pdf'):

        import os
        from pyPdf import PdfFileWriter, PdfFileReader
        
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
 
        def header(txt, style=HeaderStyle, klass=Paragraph, sep=0.3): #return styled text with a space
            s = Spacer(0.2*inch, sep*inch)
            para = klass(txt, style)
            sect = [s, para]
            result = KeepTogether(sect)
            return result
 
        def p(txt): #wrapper for header
            return header(txt, style=ParaStyle, sep=0.1)
 
        def pre(txt): #return styled text with a space
            s = Spacer(0.1*inch, 0.1*inch)
            p = Preformatted(txt, PreStyle)
            precomps = [s,p]
            result = KeepTogether(precomps)
            return result
      
        def graphout(name, datain,xaxis=None):
            if xaxis is None:
                xaxis=range(datain.size)
            if xlabel is None:
                xlabel = ''
            if ylabel is None:
                ylabel = ''
      
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
            lp.xValueAxis.xLabelFormat       = '{mmm} {yy}'  
            lp.lineLabels.fontSize           = 6  
            lp.lineLabels.boxStrokeWidth     = 0.5  
            lp.lineLabels.visible            = 1  
            lp.lineLabels.boxAnchor          = 'c'  
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
            doc = SimpleDocTemplate('meta.pdf')
            doc.build(Elements)
 
        mytitle = header(Title)
        myname = header(Author, sep=0.1, style=ParaStyle)
        mysite = header(URL, sep=0.1, style=ParaStyle)
        mymail = header(email, sep=0.1, style=ParaStyle)
        abstract_title = header("ABSTRACT")
        myabstract = p(Abstract)
        head_info = [mytitle, myname, mysite, mymail, abstract_title, myabstract]
        Elements.extend(head_info)
 
        meta_title = header("metadata", sep=0)
        metasection = []
        metasection.append(meta_title)
      
        for i,elem in enumerate(self.meta):
            #if type(self.meta[elem]) != type(np.array([])):
                
            if type(self.meta[elem])== type({}):
                data = np.array(self.meta[elem]['v'])
                unit_label = str(self.meta[elem]['u'])
            else:
                data = np.array(self.meta[elem])
                unit_label = ''

            xaxis = self.meta['Rxy']['v'][:,self.ny/2]
            if data.shape == (self.nx,self.ny):
                datastr = data[:,self.ny/2]
            #metasection.append(graphout('stuff',datastr,xaxis=xaxis))  
            #metasection.append(RL_Plot(datastr,xaxis))
                metasection.append(RL_Plot(datastr,xaxis,linelabel=str(elem)))
            #metasection.append(RL_Plot(datastr,xaxis,xlabel='xlabel'))
            elif data.shape == self.nx:
                datastr = data
            #metasection.append(graphout('stuff',datastr,xaxis=xaxis))
                metasection.append(RL_Plot(datastr,xaxis,linelabel=str(elem)))
            else:
                print elem, data
                metasection.append(header(str(elem)+': '+str(data)+ ' '+ unit_label, 
                                          sep=0.1, style=ParaStyle))
                   
     
        src = KeepTogether(metasection)
        Elements.append(src)

        cxxtitle = header("Equations in CXX")
        cxxsection = []
      #print self.cxx
        cxxsection.append(header(self.cxx[0], sep=0.1, style=ParaStyle))
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
        
        output = PdfFileWriter()
        metapdf = PdfFileReader(file("meta.pdf", "rb"))
        mainpdf = PdfFileReader(file("output.pdf", "rb"))

        for i in range (0, metapdf.getNumPages() ):
            output.addPage(metapdf.getPage(i))
      
        for i in range (0, mainpdf.getNumPages() ):
            output.addPage(mainpdf.getPage(i)) 

        outputFile = filename
        outputStream = file(outputFile, "wb")
        output.write(outputStream)
        outputStream.close()
        print "Consolidation complete."


class subset(LinResDraw):
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


