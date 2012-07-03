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
    peaks_db = fft_info(data,user_peak,meta=meta) #Nt X Nx X (# of loc. max) list of dict
     
    #print peaks[0]['gamma']
    return peaks_db,ave

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
    dz = meta['dz']
    #IC = meta['IC']
    ky_max = ny/2 
    kz_max = nz/2 
    amp = abs(data).max(2).max(2) #nt x nx
 
    print 'dt: ',dt

    #print data[0,2,:,:]

    IC = amp[0,:].max() #intial condition, set
    print IC

    fft_data = np.fft.fft2(data)[:,:,0:ky_max,0:kz_max] #by default the last 2 dimensions
 
    power = fft_data.conj() * fft_data
    #print power[0].max(), (IC*(ky_max)*(kz_max))**2
    
    cross_pow = (fft_data * (np.roll(fft_data,1,axis=0)).conj())/(ny*nz)
    
    if rescale:
        fft_data_n = np.fft.fft2(data_n)[:,:,0:ky_max,0:kz_max]
        pow_n = np.sqrt((fft_data_n.conj() * fft_data_n).real)
     
    peaks = [[[] for i in range(nx)] for j in range(nt)] #a list of dictionaries
    peaks_db = []

    peak_hist =  [[0 for i in range(kz_max)] for j in range(ky_max)] #a 2d bin array

  #for now using a lame 2x loop method
    
    
    if user_peak != 0:
        for mem in user_peak:
            print mem
            peak_hist[int(mem[0])][int(mem[1])] = abs(power.mean(0).mean(0)[int(mem[0]),int(mem[1])])
        
        #floor =  ((IC*(kz_max*ky_max))**2)/10000
           
    else:
        for t in range(nt):
            for x in range(nx):
                peaks[t][x] = local_maxima(power[t,x,:,:],0,floor = (IC*(kz_max*ky_max))**2)
                for p in peaks[t][x] : #looping over each returned peakset at some fixed t,x pair
                    peak_hist[p['y_i']][p['z_i']] +=1 #average across t and x, at least exclude pad
        floor = 0        
                  

    #this array is usefull for determining what the dominant modes are
    #but we want to retain the option of observing how the amplitude
    #of any give harmonic varies in space
    
    peak_hist = np.array(peak_hist) 
    
     #let's find the top N overall powerfull harmnonics
    net_peak = local_maxima(peak_hist,user_peak,bug=False)


    print 'net_peak: ',net_peak,user_peak != 0
    #dom_mode = [{'amp':[],'amp_n':[],'phase':[],'freq':[],'gamma':[]} for x in net_peak]
    dom_mode_db = []
    
    
    
    Bp = meta['Bpxy']['v'][:,ny/2]
    B = meta['Bxy']['v'][:,ny/2]
    

    rho_s = meta['rho_s']['v']
    
    L_z = meta['L_z']/rho_s
    #L_z = 
    L_y = meta['lpar'] #already normalized earlier in read_inp.py
    L_norm = meta['lbNorm']/rho_s

    hthe0_n = 1e2*meta['hthe0']['v']/rho_s #no x dep
    hthe0_n_x = L_y/(2*np.pi) #no x dep

    print 'L_z,Ly: ' , L_z,L_y

    #if user provides the harmo    nic info overide the found peaks


    #thi is where all the good stuff is picked up

    #look at each mode annd pull out some usefull linear measures
    for i,p in enumerate(net_peak):
        print i,p['y_i'],p['z_i']
        amp =  (np.sqrt(power[:,:,p['y_i'],p['z_i']])/(kz_max*ky_max)).real
        
       
        #phase = (np.angle(cross_pow[:,:,p['y_i'],p['z_i']],deg=False)).real
        phase = -np.array(np.gradient((np.angle(fft_data[:,:,p['y_i'],p['z_i']],deg=False)).real)[0]) #nt x nx
        
        #loop over radaii
        phasenew = []
        from scipy.interpolate import interp2d,interp1d
        from scipy import interp
        for phase_r in np.transpose(phase):
            jumps = np.where(abs(phase) > np.pi/32)
            print jumps
            if len(jumps[0]) != 0:
                        
                all_pts = np.array(range(0,nt))
                good_pts = (np.where(abs(phase_r) < np.pi/3))[0] 
            #print good_pts[0],good_pts
                #f = interp1d(good_pts,phase_r[good_pts],fill_value=.001)
                #print max(all_pts), max(good_pts)
                #phasenew.append(f(all_pts))
                phase_r = (interp(all_pts,good_pts,phase_r[good_pts]))
            phasenew.append(phase_r)

        phase = np.transpose(phasenew)/dt    
        # phase_f =np.roll(phase,-1,axis=0) #forward
        # phase_b =np.roll(phase,1,axis=0) #back

        #jumps = np.where(abs(phase) > np.pi/8)
        
#         for loc in np.transpose(np.array(jumps)):
#             #loc is the (t,r) index of a jump
            
        #alli = np.where(abs(phase) == abs(phase))

# #nearest points
#         edge1 = (jumps[0]-1,jumps[1]) #t location
#         edge2 = (jumps[0]+1,jumps[1]) # will have as many or less index pairs thatn jumps

#         edge1i = np.where(edge1[0] > 0)
#         edge1 = (edge1[0][edge1i],edge1[1][edge1i])
 
#         edge2i = np.where(edge2[0] < nx-1)
#         edge2 = (edge2[0][edge2i],edge2[1][edge2i])
        
#         #construct pairs of vals that do not 

#         edge = []
#         #loop over the gaps
#         for loc in np.transpose(np.array(edge1)): #will iterate over the 1st var npairs x 2
#             if loc not in np.transpose(np.array(edge1)):
#                 edge.append[loc]


            
        #edge2i = np.setdiff1d(edge2[0][edge2i],jumps[0])

        
        
        # if len(jumps[0]) != 0:
        #     print 'jumps', len(jumps[0])
        #     phase[jumps] = .0001
            #for elem in jumps2:
             #   phase[elem[0]+1:elem[1]-1] = (phase[elem[0]]+phase[elem[1])/2.0


          
        #phase = phase/dt  
        #for x in range(nx):
            
        #print 'jumps',jumps.shape,phase.shape
       

        amp_n = (np.sqrt(power[:,:,p['y_i'],p['z_i']])/(kz_max*ky_max*amp)).real
        #amp_n = dom_mode[i]['amp_n'] #nt x nx
   
        #let just look over the nx range
        lnamp = np.log(amp)

        t = dt*np.array(range(nt)) #dt matters obviouslyww
        r = np.polyfit(t[nt/2:],lnamp[nt/2:,2:-2],1,full=True)

        gamma_est = r[0][0] #nx
        f0 = np.exp(r[0][1]) #nx
        res = r[1]
        pad =[0,0]       
        gamma_est = np.concatenate([pad,gamma_est,pad])
        f0 = np.concatenate([pad,f0,pad])
        res = np.concatenate([pad,res,pad])
        
                #sig = res/np.sqrt((x['nt']-2))
        sig = np.sqrt(res/(nt-2))
                #sig0 = sig*np.sqrt(1/(x['nt'])+ ) # who cares
        sig1 = sig*np.sqrt(1.0/(nt * t.var()))
        nt = np.array(nt)
        print 'shapes ', nt.shape, nt, lnamp.shape
        res = 1 - res/(nt*lnamp.var(0)) #nx 
        res[0:2] = 0
        res[-2:] = 0

        gamma=[gamma_est,sig1,f0,res]    
                
                
        
        # gamma_est2 = np.gradient(amp)[0]/(amp[:,:]*dt)
        # gamma_w =  np.gradient(gamma_est2)[0]
        
        # gamma_i = np.abs(gamma_w).argmin(0) #index of the minimum for any given run
        # for j in range(nx):
        #     gamma_w[0:max([gamma_i[j],nt/3]),j] = np.average(gamma_w)*100000.0
           

        freq =  weighted_avg_and_std(
            phase[-10:,:],weights=np.ones(phase[-10:,:].shape)) 
        
       
        # gamma = weighted_avg_and_std(
        #         gamma_est2[-5:,:],weights=np.ones(gamma_est2[-5:,:].shape))  
        
        
        k = [[p['y_i'],p['z_i']],
             [2*math.pi*rho_s*float(p['y_i'])/L_y,2*math.pi*p['z_i']/L_z]]
        #L_y is normalized
        
        #simple k def, works in drift-instability fine
        # k = [[p['y_i'],p['z_i']],
        #      [(B/Bp)**-1*2*math.pi*float(p['y_i'])/(L_y),(B/Bp)*2*math.pi*p['z_i']/L_z]]

        k_r = [[p['y_i'],p['z_i']],
             [(Bp/B)*2*math.pi*float(p['y_i'])/(L_y),
              (B/Bp)*2*math.pi*p['z_i']/L_z]]
      
        #what I think is the most general one, works in drift-instability again 
        # seems to work for Bz only helimak, now trying Bp = Bt
        # k = [[p['y_i'],p['z_i']],
        #      [((Bp/B)*float(p['y_i'])/(hthe0_n)) +
        #       2*np.pi*p['z_i']*np.sqrt(1-(Bp/B)**2)/L_z,
        #       2*math.pi*p['z_i']/L_norm - 
        #       float(p['y_i'])*np.sqrt(1-(Bp/B)**2)/(hthe0_n)]]
        # k = [[p['y_i'],p['z_i']],
        #      [((Bp/B)*float(p['y_i'])/(hthe0_n)),
        #       2*math.pi*p['z_i']/L_norm]]
        #BOTH SEEM TO PRODOCE SAME RESULTS

        # k = [[p['y_i'],p['z_i']],
        #      [(float(p['y_i'])/(hthe0_n_x)),
        #       2*math.pi*float(p['z_i'])/L_norm]]

        dom_mode_db.append({'modeid':i,'k':k[1],'gamma':gamma,'freq': freq,
                            'amp': amp,'amp_n':amp_n,'phase':phase,'mn':k[0],
                            'nt':nt,'k_r':k_r[1]})

    
    return dom_mode_db

#return a 2d array fof boolean values, a very simple boolian filter
def local_maxima(array2d,user_peak,index=False,count=4,floor=0,bug=False):

    from operator import itemgetter, attrgetter
    
    if user_peak == 0:
        where = ((array2d >= np.roll(array2d,1,0)) &
                 (array2d >= np.roll(array2d,-1,0)) &
                 (array2d >= np.roll(array2d,0,1)) &
                 (array2d >= np.roll(array2d,0,-1)) &
                 (array2d >= array2d.max()/5.0) &
                 (array2d > floor*np.ones(array2d.shape)) &
                 (array2d >= array2d.mean()))
    else: #some simpler filter if user indicated some modes
        where = array2d > floor

    #ignore the lesser local maxima, throw out anything with a ZERO
    if bug==True:    
        print array2d,array2d[where.nonzero()],where.nonzero()[0]
    
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
    
    if len(values.shape) == 2:
         average = np.average(values,0)#, weights=weights)
         variance = (np.inner(weights.transpose(),((values-average)**2).transpose()).diagonal())/weights.sum(0)
    else:     
        average = np.average(values, weights=weights)
        variance = np.dot(weights, (values-average)**2)/weights.sum()  # Fast and numerically precise

    return [average,variance]
