#basic_info return some statistical averages and harmonic info
import numpy as np
import math

def basic_info(data,meta,rescale=True,rotate=False,user_peak=0):
    
    print 'in basic_info'
    #from . import read_grid,parse_inp,read_inp,show
 
    dims = data.shape
    ndims = len(dims)
    
    
    
    if ndims ==4:
        nt,nx,ny,nz = data.shape
        print nt,nx,ny
    else:
        print "something with dimesions"

    dc = data.mean(1).mean(1).mean(1) # there MUST be a way to indicate all axis at once
    amp = abs(data).max(1).max(1).max(1)

    if rescale:
        amp_o = amp - dc
        fourDamp = np.repeat(amp_o,nx*ny*nz)
        fourDamp = fourDamp.reshape(nt,nx,ny,nz)
        dc_n = dc/amp_o
        data_n = data/fourDamp
        
        ave = {'data':data,'amp':amp,'dc':dc,'amp_o':amp_o}
        
    else:
        print "no rescaling"
        ave = {'amp':amp,'dc':dc}
    if rotate:
        print 'rotate stuff'
        # will need to provide some grid geometry to do this one
    else:        print 'or not'

    #let's identify the dominant modes, for now at every [t,x] slice
    #if the data set is too large we can average over x
    peaks = fft_info(data,user_peak,meta=meta) #Nt X Nx X (# of loc. max) list of dict
     
    #print peaks[0]['gamma']
    return peaks,ave

def fft_info(data,user_peak,dimension=[3,4],rescale=False,wavelet=False,show=False,meta=0):
    import numpy as np
    import math
 
    print 'in fft_inf0'

    dims = data.shape
    ndims = len(dims)
    
    if ndims==4:
        nt,nx,ny,nz = data.shape
        print data.shape
    else:
        print "something with dimesions"
        
    #dt, k labels for the revelant dimensions 
    

    dt = meta['dt']['v']
    ky_max = ny/2 
    kz_max = nz/2 
    amp = abs(data).max(2).max(2)
 
    print 'dt: ',dt

    fft_data = np.fft.fft2(data)[:,:,0:ky_max,0:kz_max] #by default the last 2 dimensions
 
    power = fft_data.conj() * fft_data
    cross_pow = (fft_data * (np.roll(fft_data,1,axis=0)).conj())/(ny*nz)
    
    if rescale:
        fft_data_n = np.fft.fft2(data_n)[:,:,0:ky_max,0:kz_max]
        pow_n = np.sqrt((fft_data_n.conj() * fft_data_n).real)
     
    peaks = [[[] for i in range(nx)] for j in range(nt)] #a list of dictionaries

    peak_hist =  [[0 for i in range(kz_max)] for j in range(ky_max)] #a 2d bin array

  #for now using a lame 2x loop method
    
    
    if user_peak != 0:
        for mem in user_peak:
            peak_hist[mem[0]][mem[1]] = abs(power.mean(0).mean(0)[mem[0],mem[1]])
            print mem
    else:
        for t in range(nt):
            for x in range(nx):
                peaks[t][x] = local_maxima(power[t,x,:,:],0)
                for p in peaks[t][x] : #looping over each returned peak at some fixed t,x pair
                    peak_hist[p['y_i']][p['z_i']] +=1 #average across t and x, at least exclude pad
                
    

    #this array is usefull for determining what the dominant modes are
    #but we want to retain the option of observing how the amplitude
    #of any give harmonic varies in space
    peak_hist = np.array(peak_hist)
    
    
     #let's find the top N overall powerfull harmnonics
    net_peak = local_maxima(peak_hist,user_peak)


    print 'net_peak: ',net_peak,user_peak != 0
    dom_mode = [{'amp':[],'amp_n':[],'phase':[],'freq':[],'gamma':[]} for x in net_peak]
    
    L_z = meta['L_z']
    L_y = meta['lpar']
    rho_s = meta['rho_s']['v']

    print 'L_z,Ly: ' , L_z,L_y

    #if user provides the harmo    nic info overide the found peaks


    #thi is where all the good stuff is picked up
    i=0
    
    
    #look at each mode and pull out some usefull linear measures
    for p in net_peak:
        print i,p['y_i'],p['z_i']
        dom_mode[i]['amp'] = (np.sqrt(power[:,:,p['y_i'],p['z_i']])/(kz_max*ky_max)).real
        dom_mode[i]['amp_n'] = (np.sqrt(power[:,:,p['y_i'],p['z_i']])/(kz_max*ky_max*amp)).real
        dom_mode[i]['phase'] = (np.angle(cross_pow[:,:,p['y_i'],p['z_i']],deg=False)).real
        #phase = dom_mode[i]['phase'][:,2:nx-3]/dt 
        #amp_n = dom_mode[i]['amp_n'][:,2:nx-3]
        phase = dom_mode[i]['phase']/dt 
        amp_n = dom_mode[i]['amp_n'] #nt x nx
       
        amp = dom_mode[i]['amp']#[:,2:nx-3] 
        
        #let come up with a good weight
        #amp is nt by nx

        #let just look over the nx range
        gamma_est2 = np.diff(amp,axis = 0)/(amp[0:-1,:]*dt)
        gamma_w =  np.diff(gamma_est2,axis=0)
        
        gamma_i = np.abs(gamma_w).argmin(0) #index of the minimum for any given run
        for j in range(nx):
            gamma_w[0:max([gamma_i[j],nt/3]),j] = np.average(gamma_w)*100000.0
           

        
        dom_mode[i]['freq'] = weighted_avg_and_std(phase[2:,:],weights=1/gamma_w)# ave(dom_mode[i]['phase'])/(time window for ave)
        dom_mode[i]['gamma'] =  weighted_avg_and_std(gamma_est2[1:,:],weights=1/gamma_w)# deriv(dom_mode[i]['amp'])/(dt amp)
        dom_mode[i]['k'] = [[p['y_i'],p['z_i']],
                            [2*math.pi*rho_s*float(p['y_i'])/L_y,2*math.pi*rho_s*p['z_i']/L_z]]
        
        i+=1

    
    return dom_mode

#return a 2d array fof boolean values, a very simple boolian filter
def local_maxima(array2d,user_peak,index=False,count=3):

    from operator import itemgetter, attrgetter
    
    if user_peak == 0:
        where = ((array2d >= np.roll(array2d,1,0)) &
                 (array2d >= np.roll(array2d,-1,0)) &
                 (array2d >= np.roll(array2d,0,1)) &
                 (array2d >= np.roll(array2d,0,-1)) &
                 (array2d >= array2d.max()/5.0) &
                 (array2d >= array2d.mean()))
    else:
        where = array2d

    #ignore the lesser local maxima
    
    peaks = zip(where.nonzero()[0],where.nonzero()[1],array2d[where.nonzero()])
    
    peaks = sorted(peaks,key=itemgetter(2),reverse=True)
   
    if len(peaks) > count and user_peak==0:
            peaks = peaks[0:count]
            
    keys = ['y_i','z_i','amp']
    
    peaks = [dict(zip(keys,peaks[x])) for x in range(len(peaks))]
    
    return peaks
    #return np.array(peak_dic)

def weighted_avg_and_std(values, weights):
    """
    Returns the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    print 'taking an average'

    # N dim in -> N-1 dim out
    print values.shape
    print weights.shape

    
    if len(values.shape) == 2:
         average = np.average(values,0)#, weights=weights)
         variance = (np.inner(weights.transpose(),((values-average)**2).transpose()).diagonal())/weights.sum(0)
    else:     
        average = np.average(values, weights=weights)
        variance = np.dot(weights, (values-average)**2)/weights.sum()  # Fast and numerically precise

    return [average,variance]
