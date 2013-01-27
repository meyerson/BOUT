#basic_info return some statistical averages and harmonic info
import numpy as np
import math

def normed(data):
    
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
   
    amp_o = amp - dc
    fourDamp = np.repeat(amp_o,nx*ny*nz)
    fourDamp = fourDamp.reshape(nt,nx,ny,nz)
 
    data_n = data/fourDamp
    
    
    return data_n,amp

    
