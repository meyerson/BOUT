import sys
import pickle
import os
import json
 
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages   
    #print list(sys.modules.keys())

from multiprocessing import Queue,Pool
import multiprocessing
import math

try:
    boutpath = os.environ['BOUT_TOP']
    pylibpath = boutpath+'/tools/pylib'
    pbpath = pylibpath+'/post_bout'
    boutdatapath = pylibpath+'/boutdata'
    boututilpath = pylibpath+'/boututils'
    
    allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    # sys.path.append('/home/crxyosphere/BOUT/tools/pylib/post_bout')
    #sys.path.append(allpath)
    [sys.path.append(elem) for elem in allpath]
    print sys.path
        

    #import gobject
    import numpy as np
    
    print 'in post_bout/post_bout.py'
    from ordereddict import OrderedDict
    from scipy.interpolate import interp2d,interp1d,griddata
    #import scipy.interpolate.RectBivariateSpline as fastgrid
    from scipy.ndimage import interpolation as interp

    from boutdata import collect

    
except ImportError:
    print 'in post_bout/rotate.py'

    print "can't find the modules I need, you fail"
    sys.exit() #no point in going on
            
    #import some bout specific modules
try:
    from boutdata import collect
    
except:
    print 'in post_bout/post_bout.py'
    print "can't find bout related modules, you fail"
    sys.exit() #no point in going on, shoot yourself 
    


#i feel like I shouldn't have to include these line in virtue of __init__.py
from read_inp import parse_inp, read_inp, read_log,metadata
from read_grid import read_grid
from basic_info import basic_info, fft_info
from corral import corral


def rotate(data,meta,rt=False,spectrum=False,view=False,save=False):
    dims = data.shape
    ndims = len(dims)
    
    print 'rotate_mp'
    print 'dims: ', dims
    TWOPI = 2*np.pi
    #one can simply pass a slice or a full 4D simuation array 
    # later choice is made one can futher  chose to try to rotate it for all t
    #   and r or specify r and t to save memory
    try:
        if ndims ==4:
            nt,nx,ny,nz = dims
            print nt,nx,ny
        if ndims == 2:
            ny,nz =dims
            nx = meta['nx']
            #nt = meta['NOUT']+1
            nt = 0
            data = np.array([[data]])
    except:
        print "something with dimesions"
        
    n_zz = nz-1

    Bpxy = meta['Bpxy']['v']
    Bxy = meta['Bxy']['v']
    Btxy = meta['Btxy']['v']
    Zxy = meta['Zxy']['v']

    rho_s = meta['rho_s']['v'] #cm
    
    L_z = (meta['L_z']/rho_s)[nx/2]  #meta is setup in read_inp, L_z is just 2pi*R
    L_y = (meta['lpar'])[nx/2] #already normalized earlier in read_inp.py , B/Bp * 2p*hthe_n
    L_norm = (meta['lbNorm']/rho_s)[nx/2] #L_z * B/Bp normed

     
    hthe0_n = 1e2*meta['hthe0']['v']/rho_s #no x dep, minor rad
    hthe0_n_x = L_y/(2*np.pi) #with x dep

    # TWOPI = 2*np.pi
    
    dz = L_z/n_zz
    dy = L_y/ny #stupid way to do this better to look at hthe
    #dy =  1e2*Bxy[:,ny/2]/Bpxy[*]meta['dlthe']/rho_s

    # z = indgen(n_zz)*dz
    # y = indgen(ny)*dy
    z = np.array(range(0,nz))*dz
    y = np.array(range(0,ny))*dy

    
    cos_theta = np.max(Btxy,1)/np.max(Bxy,1) #nx sized array
    sin_theta = np.max(Bpxy,1)/np.max(Bxy,1) #
    
    print 'cos_theta,sin_theta:', cos_theta,sin_theta

    tan_phi_c = cos_theta 
    tan_phi_i = dz/dy * tan_phi_c #nx
    print tan_phi_i,sin_theta,dz,dy,L_z,L_y
    #new_z = z*sin_theta
    new_z = np.outer(z,sin_theta) #nz x nx 2d array 
    new_Lz = L_z*sin_theta #nx 
     #new_Lz = Lz,sin_theta)
    dz2 = dz*sin_theta #1d array
  
    z = np.outer(TWOPI/new_Lz ,np.array(range(0,n_zz))) #nz x nz
    
    ky = TWOPI/L_y * np.array(range(ny)) #ny
    
    N_turns = np.round(np.max((Zxy*Btxy)/(TWOPI*Bpxy)))


    #print 'zxy: ',Zxy
    print 'N_turns: ',N_turns

    
    sz = 1.0
    sy = 1.0
    dz= dz/sz
    dy = dy/sy
    
    #MZ = sz*N_turns*n_zz
    #want to keep the end points
    MZ = sz*nz 
    MY = sy*(ny)
    
#original index coordinates
    i2d = np.repeat(np.array(range(ny)),nz)
    i2d = i2d.reshape(ny,nz)

    j2d = np.repeat(np.array(range(nz)),ny)
    j2d = j2d.reshape(nz,ny)
    j2d = j2d.transpose()
#index coordinates from old to hi-res unrotated version,
#so all i and j values must be less than ny and nz 
#respectively
    ssy =(MY-1)/float(ny-1)
    ssz = (MZ-1)/float(nz-1)
    
    ii = np.repeat(np.array(range(MY))/float(ssy),MZ)
    ii = ii.reshape(MY,MZ)
    jj = np.repeat(np.array(range(MZ))/float(ssz),MY)
    jj = jj.reshape(MZ,MY)
    jj = jj.transpose()
     
#index coordinates from high-res unrotated to hd-rotated,
#so all i and j values must be LT MY,MZ 
    #ii2 = np.mod(((N_turns+2)*ny)+ii- (jj* tan_phi_i[2]),ny-1)
    #tan_phi_i[2] = .15
    offset_i = nz*tan_phi_i[2] 
    ii2 = offset_i + (ny - offset_i - 1)/(ny-1) * ii- (jj* tan_phi_i[2])
   
    
   
    output = np.zeros((nt,nx,MY,MZ))
    #can we speed this up?
    
    #let's use multiprocessing to speed things up
    
   

    
    def rotor(t,x,out_q,p_id):            
        hd_image = griddata((i2d.flatten(),j2d.flatten()), 
                            (data[t,x,:,:]).flatten(), 
                            (ii,jj), method='linear')
        rot_image = griddata((ii.flatten(),jj.flatten()), 
                             hd_image.flatten(), 
                             (ii2,jj), 
                             method='linear')

        print 't,x,p_id',t,x,p_id
        out = {'data':rot_image,'t':t,'x':x}
        
        out_q.put(out)

    out_q = Queue()
    nprocs = 8 # have to determine the maximum number later, careful here
    Nchunks = int(math.ceil(nt / float(nprocs))) #number of cells / proc
    procs = []
    
    #break up the work radially
    
    for x in range(nx): 
        #break up the work in time
        T = nt
        for chunk in range(Nchunks):
            if T < nprocs:
                nprocs = T
            else:
                nprocs = 8

            out_q = Queue()
            procs = []  
            for i in range(nprocs):
                t = i
                p = multiprocessing.Process(
                    target=rotor,
                    args=(T-t-1,x,out_q,i)) #into the past for some reason
                procs.append(p)
                p.start()
            
            T = T - nprocs
            for i in range(nprocs):
                temp = out_q.get()
                output[temp['t'],temp['x'],:,:] = temp['data']

            for p in procs:
                p.join()
   
    if not view:   
        return output      
    else:

        tile_image = data[0,2,:,:]
        print tile_image.shape
        
        new_image0 = griddata((i2d.flatten(),j2d.flatten()), 
                            (data[0,2,:,:]).flatten(), 
                            (ii,jj), method='linear')
        
        new_image1 = output[0,2,:,:]
        N_zeta2 = np.round(N_turns*n_zz)     #assume that N_zeta Ny are even
        Nyy = np.round(ny)
        
        pp = PdfPages('rotate.pdf')
        fig = plt.figure()
        img00 = fig.add_subplot(2,2,1)
        img01 = fig.add_subplot(2,2,2)
        img11 = fig.add_subplot(2,2,3)
        img10 = fig.add_subplot(2,2,4)


        cmap = None
        img00.contourf(i2d,j2d,tile_image,interpolation='bilinear', cmap=cmap)
        #img00.plot(i2d.flatten(), j2d.flatten(), 'k.', ms=1)
        #img00.plot(ii.flatten(), jj.flatten(), 'r.', ms=1)
        img00.axis('tight')

        img01.contourf(ii,jj,new_image0,interpolation='bilinear', cmap=cmap)

        img01.plot(ii2[0,:].flatten(), jj[0,:].flatten(), 'k.', ms=1)
        img01.plot(ii2[-1,:].flatten(), jj[-1,:].flatten(), 'k.', ms=1)

        img01.plot(ii2[:,110].flatten(), jj[:,110].flatten(), 'k*', ms=1)
        img01.plot(ii2[10:20,110].flatten(), jj[10:20,110].flatten(), 'k-', ms=1)


        img01.axis('tight')
        img11.contourf(ii,jj,new_image1,interpolation='bilinear', cmap=cmap)
        img11.plot(ii[10,:].flatten(), jj[10,:].flatten(), 'k.', ms=1)
        img11.plot(ii[:,110].flatten(), jj[:,110].flatten(), 'k*', ms=1)
        img11.axis('tight')

        img10.contourf(ii[-5:,0:10],jj[-5:,0:10],new_image0[-5:,0:10],interpolation='bilinear', cmap=cmap)
        img10.plot(ii2[-5:,0:10].flatten(), jj[-5:,0:10].flatten(), 'r-', ms=1)


        ky_max = ny/2 
        kz_max = nz/2 
        fft_data = np.fft.fft2(data)[:,:,0:ky_max,0:kz_max]
        power = fft_data.conj() * fft_data

        fft_data_new = np.fft.fft2(new_image1)[0:sy*ky_max-1,0:sz*kz_max-1]
        power_new = fft_data_new.conj() * fft_data_new

        fig.savefig(pp, format='pdf')
        plt.close(fig)
        fig = plt.figure()
        img00 = fig.add_subplot(2,2,1)
        img01 = fig.add_subplot(2,2,2)
        img11 = fig.add_subplot(2,2,3)
        img10 = fig.add_subplot(2,2,4)


        img00.contourf(np.abs(power[0,2,0:10,0:10]),interpolation='bilinear', cmap=cmap)
        img01.contourf(np.abs(power_new[0:10,0:10]),interpolation='bilinear', cmap=cmap)

        img00.contour(np.log(np.abs(power[0,2,0:10,0:10])),interpolation='bilinear', cmap=cmap)
        img01.contour(np.log(np.abs(power_new[0:10,0:10])),interpolation='bilinear', cmap=cmap)

        img11.contourf(np.log(np.abs(power[0,2,:,:])),interpolation='bilinear', cmap=cmap)
        img10.contourf(np.log(np.abs(power_new)),interpolation='bilinear', cmap=cmap)

        img00.axis('tight')
        img01.axis('tight')


        plt.close(fig)
        fig.savefig(pp, format='pdf') 
        pp.close()


        #plot old and new side by side

        return 0

