# Display animations of data, similar to the IDL routine
# of the same name
#
# Ben Dudson, University of York, July 2009
#
#



try:
    import matplotlib
    matplotlib.use('Agg') # do this before importing pylab
    #matplotlib.pyplot.switch_backend('Agg')
    import numpy as np
    import matplotlib.pyplot as plt
    import subprocess  
    
except ImportError:
    print "ERROR: Showdata needs numpy, matplotlib and gobject modules"
    raise

def showdata_tacc(data, scale=True, loop=False,movie=True):
    """Animate a dataset, with time as first dimension
    
    2D - Show a line plot
    3D - Show a surface plot

    scale = True  Use the same scale for all times
    loop  = False Loop the dataset
    """

    size = data.shape
    ndims = len(size)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ndims == 2:
        # Animate a line plot
     
        def animate():
            line, = ax.plot(data[0,:])
            if scale:
                # Get data range
                ax.set_ylim([np.min(data), np.max(data)])
            while True:
                for i in np.arange(size[0]):
                    line.set_ydata(data[i,:])
                    if not scale:
                        ax.set_ylim([np.min(data[i,:]), np.max(data[i,:])])
                    fig.canvas.draw()
                    filename = str('%03d' % i) + '.png'
                    plt.savefig(filename, dpi=100)
                    yield True
                if not loop: break
      
        plt.show()
        
    elif ndims == 3:
        # Animate a contour plot

        if widget == "gtk":
            def animate():
                cmap = None
                m = plt.imshow(data[0,:,:], interpolation='bilinear', cmap=cmap, animated=True)
                while True:
                    for i in np.arange(size[0]):
                        m.set_data(data[i,:,:])
                        fig.canvas.draw()
                        yield True
                    if not loop: break

            gobject.idle_add(lambda iter=animate(): iter.next())
        else:
            # WX widgets method
            cmap = None
            m = plt.imshow(data[0,:,:], interpolation='bilinear', cmap=cmap, animated=True)
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


