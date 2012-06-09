#eventually the default path shoudl really be the current directory
#import post_bout
#post_bout.show()
#from . import read_grid,parse_inp, read_inp,basic_info
#import sqlite3 as sql

import sys
import pickle
import os
import json
    
    #print list(sys.modules.keys())

try:
    boutpath = os.environ['BOUT_TOP']
    pylibpath = boutpath+'/pylib'
    pbpath = pylibpath+'/post_bout'
    boutdatapath = pylibpath+'/boutdata'
    boututilpath = pylibpath+'/boututils'
    
    allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
    #sys.path.append(allpath)
    [sys.path.append(elem) for elem in allpath]
    
        

    #import gobject
    import numpy as np
    from ordereddict import OrderedDict
    from scipy.interpolate import interp2d,interp1d
    from boutdata import collect
except ImportError:
    print 'in post_bout/post_bout.py'

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


def save(path='/home/cryosphere/BOUT/examples/Q3/data_short',
         savemovie=False,IConly=0): 
    #lets collect the data
    print 'path :', path
    print 'in post_bout/post_bout.save'
    boutpath = path

    print path
    print sys.path
 

    meta = metadata(path=path)
    
    #ICmodes = 
    print 'ok got meta'
    #return meta
    output = OrderedDict()
    data = OrderedDict()

   
    all_modes = []
    print meta['evolved']['v']

    for i,active in enumerate(meta['evolved']['v']): #loop over the fields
        print path, active
        data[active] = collect(active,path=path)
        if savemovie:
            movie(data,meta)
        

    minIC = np.array([np.max(data[x][0,:,:,:]) for x in meta['evolved']['v']])
    minIC = min(minIC[np.array(meta['IC']).nonzero()])  
    ICmodes =[]

    if(meta['ys_opt']['v']==2 and  meta['zs_opt']['v']==2): #if BOUT.inp IC indicates a single mode
        ICmodes.append([meta['ys_mode']['v'],meta['zs_mode']['v']])
    elif (meta['ys_opt']['v']==3 and meta['zs_opt']['v']==2): 
        [ICmodes.append([p+1,meta['zs_mode']['v']]) for p in range(5)]
    elif(meta['ys_opt']['v']==2 and meta['zs_opt']['v']==3):
        [ICmodes.append([meta['ys_mode']['v'],p+1]) for p in range(8)]
    elif(meta['ys_opt']['v']==3 and meta['zs_opt']['v']==3):
        [ICmodes.append([p+1,q+1]) for p in range(7) for q in range(8)]          
    elif(meta['ys_opt']['v']==0 and meta['zs_opt']['v']==3):
        [ICmodes.append([0,p+1]) for p in range(4)]   
  
    #some attemps to automate peak finding
    for i,active in enumerate(meta['evolved']['v']):
        modes_db,ave = basic_info(data[active],meta) #basic_info does not  correct for incomplete runs
       #print modes[0]['gamma']
        output[active] = {'ave':ave}
        #for x in output[active]['modes']: #keep track of relevant harmonics
        print 'debug:', len(modes_db),type(modes_db),len(modes_db)
        for x in modes_db: 
            print x['mn'],all_modes
            if x['mn'] not in  all_modes: 
                #don't append if the mode is weak . .
                ##minIC = np.array(meta['IC'])
                #minIC = min(minIC[minIC.nonzero()])
                print 'minIC: ',minIC, x['amp'][-10:,:].max()
                #if x['amp'][-10:,:].max()>=(minIC/1.0):
                all_modes.append(x['mn'])

    if IConly or len(all_modes)==0:
        all_modes = ICmodes
    else:
        z = []
        [z.append(x) for x in all_modes]
        [z.append(x) for x in ICmodes]
        all_modes = z
        print 'len(all_modes):',len(all_modes)
        
        
    #rebuild with all the relevant modes
    output = OrderedDict()

    output['meta'] = meta
    output['ave'] = {}
    allmodes_db =[]


    
    for i,active in enumerate(meta['evolved']['v']): 
        print 'once again', active
        print all_modes
        print meta['ys_opt']['v'], meta['zs_opt']['v']
        modes_db,ave = basic_info(data[active],meta,
                                  user_peak = all_modes)     
        
        for j,x in enumerate(modes_db):
            x['field'] = active
            x['dz']=meta['dz']
            x['IC']=meta['IC']
            x['L']=meta['L']
            x['Rxy']= meta['Rxy']['v']
            x['nfields'] = len(meta['evolved']['v'])
            x['meta'] = meta
            x['nx']= meta['nx']
            x['ny'] = meta['ny']
            x['nz'] = meta['MZ']['v']-1
            x['dt'] = meta['dt']['v']
            x['Rxynorm'] = 100*x['Rxy']/meta['rho_s']['v']
            x['rho_s'] = meta['rho_s']['v']
            
            ntt = x['nt'] #actual number of steps
            nt = meta['NOUT']['v'] #requested number of steps
            print j,' out of ',len(modes_db)
    #for incomplete runs    
            if nt > ntt: #if the run did not finish then rescale to keep data dims the same
                
                xx = np.array(range(0,ntt)) #dep var to interp FROM
                xxnew = np.linspace(0,ntt-1,nt+1) #new dep var grid, nt points up to ntt-1
                
                print ntt/nt,x['dt']*float(ntt/nt),nt,ntt
                x['dt'] = x['dt']*float(ntt/nt) #smaller dt if we are interp tp more points . . 
                
                print x['dt']
                x['nt'] = nt
                amp = x['amp'] #2d array ntt x nx
                phase = x['phase']
                ampnew=[] #will be a nx x nt 2d array
                phasenew =[]
                for mode_r in np.transpose(amp):
                    f = interp1d(xx,mode_r)
                    ampnew.append(f(xxnew))
                
                for mode_r in np.transpose(phase):
                    f = interp1d(xx,mode_r)
                    phasenew.append(f(xxnew))
                 
                x['amp'] = np.transpose(ampnew)
                x['phase'] = np.transpose(phasenew)
                print x['phase'].shape

        output['ave'][active] = {'ave':ave}
        allmodes_db.append(modes_db)
 

    allmodes_db = sum(allmodes_db,[])

    filename_db = path+'/post_bout.db'
    
    
    pickle_db = open(filename_db,'wb')

    pickle.dump(allmodes_db,pickle_db)

    pickle_db.close()

    
    # let pickle the results
    # filename = path+'/post_bout.jsn'
    # json_out = open(filename,'wb')
    # pickle.dump(output,json_out)
    # json_out.close()

    print filename_db
    
    
def read(path='.',filename='post_bout.pkl',trysave=True,
         gamma_updt=False):
    print 'in post_bout.read()'

    filepath = path+'/'+filename
    filepath_db = path+'/post_bout.db'
    datafile = path+'/BOUT.dmp.0.nc'


    print 'is file present: ', os.path.isfile(filepath_db)
    if trysave and not(os.path.isfile
                       (filepath_db)) and os.path.isfile(datafile):
        save(path=path)
  
        
        
    if os.path.isfile(filepath_db):
        #pkl_file = open(filepath, 'rb')
        pkl_file_db = open(filepath_db, 'rb')
        #output = pickle.load(pkl_file)
        output_db = pickle.load(pkl_file_db)
        
        meta = metadata(path=path)
        #update the meta data incase read_inp was altered . . .
        for i,x in enumerate(output_db):
            x['dz']=meta['dz'] 
            x['IC']=meta['IC']
            x['L']=meta['L']
            x['Rxy']= meta['Rxy']['v']
            x['nfields'] = len(meta['evolved']['v'])
            x['meta'] = meta
            x['nx'] = meta['nx']
            x['ny'] = meta['ny']
            x['nz'] = meta['MZ']['v']-1
            x['path'] = path
            #print x['mn']
            x['MN'] = list([1,2])
            x['MN'][0] = x['mn'][0]
            x['MN'][1] = x['mn'][1]/x['dz']
            x['modeid'] = [x['dz'],x['modeid']]
           # print x['mn']
            
            #x['MN'] = x['mn']
            #x['MN'][1] = x['MN'][1]/x['dz']
            
            x['Rxynorm'] = 100*x['Rxy']/meta['rho_s']['v']
            if 'dt' not in x:
                x['dt'] = meta['dt']['v'] #this one gets modified no read . 
            if gamma_updt: #experimental, recalculates gamma from amp, but for a while amp was 
                #wrong for incomplete runs s0 . . . dont' run for now
                lnamp = np.log(x['amp']) #nt x nx
                t = x['dt']*np.array(range(int(x['nt']))) 
                r = np.polyfit(t[x['nt']/2:],lnamp[x['nt']/2:,2:-2],1,full=True)
                gamma_est = r[0][0] #nx
                f0 = np.exp(r[0][1]) #nx
                res = r[1]
                pad =[0,0]       
                gamma_est = np.concatenate([pad,gamma_est,pad])
                f0 = np.concatenate([pad,f0,pad])
                res = np.concatenate([pad,res,pad])
               
                #sig = res/np.sqrt((x['nt']-2))
                sig = np.sqrt(res/(x['nt']-2))
                #sig0 = sig*np.sqrt(1/(x['nt'])+ ) # who cares
                sig1 = sig*np.sqrt(1.0/(x['nt'] * t.var()))
                res = 1 - res/(x['nt']*lnamp.var(0)) #nx 
                res[0:2] = 0
                res[-2:] = 0

                x['gamma'] = [gamma_est,f0,sig1,res]
    
        pkl_file_db.close()
        return output_db
    else:
        return 0

def movie(data,meta,name='output.avi'):
    #based on  Josh Lifton 2004 movie_demo.py file and showdata 

    
    #import modules for creating a movie  
    try:
        import matplotlib
        matplotlib.use('Agg')
        #import boututils
        import subprocess    
        import matplotlib.pyplot as plt
        
    except:
        print "can' find stuff to def movie function"

    not_found_msg = """
The mencoder command was not found;
mencoder is used by this script to make an avi file from a set of pngs.
It is typically not installed by default on linux distros because of
legal restrictions, but it is widely available.
"""

    try:
        subprocess.check_call(['mencoder'])
    except subprocess.CalledProcessError:
        print "mencoder command was found"
        pass # mencoder is found, but returns non-zero exit as expected
    # This is a quick and dirty check; it leaves some spurious output
    # for the user to puzzle over.
    except OSError:
        print not_found_msg
        print 'quitting movie'
        return 1

    
