#some standard analytic stuff to plot, if appending just overplot gam or omeg
from pb_corral import LinRes,subset
import numpy as np
import matplotlib.pyplot as plt


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
    
        if canvas is None:
            ownpage = True

        if ownpage: #if not an overplot
            fig1 = plt.figure()
            canvas = fig1.add_subplot(1,1,1) 
         
        print s.gammamax
        label ='gamma analytic'

        canvas.plot(allk[ki],np.array(s.gammamax)[ki],'-',
                    label=label)

        if comp=='gamma':
            canvas.plot(allk[ki],np.array(s.gammamax)[ki],'-',
                        label=label)
        else:
            canvas.plot(allk[ki],np.array(s.omegamax)[ki],'-',
                        label=label)
            
        print allk[ki]    
        if ownpage: #set scales if this is its own plot
            fig1.savefig(pp, format='pdf')
            plt.close(fig1)
