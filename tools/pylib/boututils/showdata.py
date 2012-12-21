# Display animations of data, similar to the IDL routine
# of the same name
#
# Ben Dudson, University of York, July 2009
#
#

#try:
try:
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    import subprocess 
    import os
    
except:
    print "ERROR: Showdata needs numpy, matplotlib and gobject modules"
    raise

try:
    from wx import *
    widget = "wx"
except:
    print "failed wx import"
    


def showdata(data, scale=True, loop=False,movie=False):
    try:
        print "Trying to import GTK..."
        import gobject
        widget = "gtk"
    except:
        print "failed\nTrying to import WX..."
   
    
   
    if widget == "gtk":
        matplotlib.use('GTKAgg')
    else:
        matplotlib.use('WXAgg') # do this before importing pylab

 


    """Animate a dataset, with time as first dimension
    
    2D - Show a line plot
    3D - Show a surface plot

    scale = True  Use the same scale for all times
    loop  = False Loop the dataset
    """

    if movie:
        return savemovie(data)

    size = data.shape
    ndims = len(size)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.set_autoscale(true)
                        
    # w,h = fig.figaspect(1);
   
    if ndims == 2:
        # Animate a line plot
        
        if widget == "gtk":
            # GTK method (default)
            def animate():
                line, = ax.plot(data[0,:])
                if scale:
                    ax.set_ylim([np.min(data), np.max(data)])
                while True:
                    for i in np.arange(size[0]):
                        line.set_ydata(data[i,:])
                        if not scale:
                            ax.set_ylim([np.min(data[i,:]), np.max(data[i,:])])
                        fig.canvas.draw()
                        yield True
                    if not loop: break
            
            gobject.idle_add(lambda iter=animate(): iter.next())
        else:
            # WX widgets method
            line, = ax.plot(data[0,:])
            def animate(idleevent):
                if scale:
                    # Get data range
                    ax.set_ylim([np.min(data), np.max(data)])

                if animate.i == size[0]:
                    wx.GetApp().Bind(wx.EVT_IDLE, None)
                    return False
              
                line.set_ydata(data[animate.i,:])
                if not scale:
                    ax.set_ylim([np.min(data[i,:]), np.max(data[i,:])])
                fig.canvas.draw_idle()
                animate.i += 1
            animate.i = 0
            wx.GetApp().Bind(wx.EVT_IDLE, animate)
        
        plt.show()
        
    elif ndims == 3:
        # Animate a contour plot

        if widget == "gtk":
            def animate():
                cmap = None
                m = plt.imshow(data[0,:,:], interpolation='bilinear', cmap=cmap, animated=True,aspect='auto')
                #c = plt.contour(data[0,:,:],8,colors='k')
                while True:
                    for i in np.arange(size[0]):
                        m.set_data(data[i,:,:])
                        #c = plt.contour(data[0,:,:],8,colors='k')
                        fig.canvas.draw()
                        yield True
                    if not loop: break

            gobject.idle_add(lambda iter=animate(): iter.next())
        else:
            # WX widgets method
            cmap = None
            m = plt.imshow(data[0,:,:], interpolation='bilinear', cmap=cmap, animated=True,aspect='auto')
            def animateContour(idleevent):
                if animateContour.i == size[0]:
                    wx.GetApp().Bind(wx.EVT_IDLE, None)
                    return False

                m.set_data(data[animateContour.i,:,:])
                fig.canvas.draw_idle()
                animateContour.i += 1
            
            animateContour.i = 0

            wx.GetApp().Bind(wx.EVT_IDLE, animateContour)
        
        plt.show()
    else:
      print "Sorry can't handle this number of dimensions"
        

def savemovie(data,data2=None,moviename='output.avi',norm=True,
              overcontour=True,aspect='auto',meta=None,mxg=2,
              cache='/tmp/'):
    size = data.shape
    ndims = len(size)
    print 'Saving pictures -  this make take a while'
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if meta != None:
        r= meta['Rxy']['v'][:,5]
        #DZ = meta['dz']
        DZ = 1
            


    #matplotlib.use('Agg')
    
    files = []
    if ndims == 2:
    
        line, = ax.plot(data[0,:])
        if scale:
            # Get data range
            ax.set_ylim([np.min(data), np.max(data)])
          # while True:
            
            for i in np.arange(size[0]):
                print i
                line.set_ydata(data[i,:])
                if not scale:
                    ax.set_ylim([np.min(data[i,:]), np.max(data[i,:])])
                #fig.canvas.draw()
                filename = str('%03d' % i) + '.png'
                plt.savefig(filename, dpi=100)
                files.append(filename)
          
    
    elif ndims == 3: #typical nt x (nx or ny) X nz

        nt,nx,nz = data.shape
        
        if norm: 
            # amp = abs(data).max(1).max(1)
            # fourDamp = np.repeat(amp,nx*ny)
            # fourDamp = fourDamp.reshape(nt,nx,ny)
            # data_n = data/fourDamp
            data_n = normalize(data)
        else:
            data_n = data

        if data2 != None:
            data_c = data2
            if norm:
                data_c = normalize(data_c)
        else:
            data_c = data_n

        cmap = None
        #m = plt.imshow(data_n[0,:,:], interpolation='bilinear', cmap=cmap, animated=True,aspect=aspect)
        if meta != None:
            dzeta = (2*DZ*np.pi)/(nz-1)
            x = np.outer(r,np.cos(dzeta*np.arange(0,nz)))
            y = np.outer(r,np.sin(dzeta*np.arange(0,nz)))
            mxg = meta['MXG']['v']
            x = x[mxg:-mxg,:]
            y = y[mxg:-mxg,:]
            data_n = data_n[:,mxg:-mxg,:]
            data_c = data_c[:,mxg:-mxg,:]
        else:
            x = np.arange(nx)
            y = np.arange(nz)
            data_n = np.transpose(data_n,[0,2,1])
            data_c = np.transpose(data_c,[0,2,1])
           
            print x.shape,y.shape,data_n.shape
        
        os.system("rm "+cache+"*png")
        m = plt.contourf(x,y,data_n[0,:,:],30,cmap=cmap)
        # if overcontour:
        c = plt.contour(x,y,data_c[0,:,:],8,colors='k')
        for i in np.arange(size[0]):
            print i
            #m.set_data(data_n[i,:,:])
            m = plt.contourf(x,y,data_n[i,:,:],30,cmap=cmap)
            
            #c = plt.contour(data[i,:,:],8,colors='k')
            #c.set_data(data[i,:,:])
            if overcontour:
               
                for coll in c.collections:
                    try:
                        plt.gca().collections.remove(coll)
                    except:
                        print 'not in this collection'
                c = plt.contour(x,y,data_c[i,:,:],8,colors='k')

           
                
        
            #fig.canvas.draw()         
            filename = cache+str('%03d' % i) + '.png'
            plt.savefig(filename, dpi=100)
            files.append(filename)
            #plt.clf()
        #plt.show()
    else:
      print "Sorry can't handle this number of dimensions"  
    
    print 'Making movie animation.mpg - this make take a while'
    command = ('mencoder',
               'mf://'+cache+'*.png',
               '-mf',
               'type=png:w=800:h=600:fps=10',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               moviename)
    
    subprocess.check_call(command)
    
    print files
    #cleanup = ('rm',files)
    os.system("rm "+cache+"*png")
    #matplotlib.use('pdf')
    #return 0

def test():
    x = np.arange(0, 2*np.pi, 0.01)
    t = np.arange(0, 2*np.pi, 0.1)
    nt = len(t)
    nx = len(x)

    # Animated line plots
    data = np.zeros([nt,nx])
    for i in np.arange(nt):
        data[i,:] = np.sin(x+i/10.0) * np.sin(10*x-i/5.0)
    
    showdata(data, loop=True)

    # Animated 2D plot
    y = x
    ny = len(y)

    data2d = np.zeros([nt,nx,ny])

    for i in np.arange(ny):
        data2d[:,:,i] = data * np.sin(y[i])

    showdata(data2d, loop=True)


def normalize(data):
    nt,nx,ny = data.shape
    amp = abs(data).max(1).max(1)
    fourDamp = np.repeat(amp,nx*ny)
    fourDamp = fourDamp.reshape(nt,nx,ny)
    data_n = data/fourDamp
    
    return data_n
