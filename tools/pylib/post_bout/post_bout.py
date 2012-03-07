#eventually the default path shoudl really be the current directory
#import post_bout
#post_bout.show()
#from . import read_grid,parse_inp, read_inp,basic_info
import sys
import pickle
import os
import json
    
    #print list(sys.modules.keys())

try:
    sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
    import matplotlib
    import gobject
    import numpy as np
    from ordereddict import OrderedDict
except ImportError:
    print "can't find the modules I need, you fail"
    sys.exit() #no point in going on
            
    #import some bout specific modules
try:
    import boutdata
    import boututils
except:
    print "can't find bout related modules, you fail"
    sys.exit() #no point in going on, shoot yourself 
    
    #import modules for creating pdfs    
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except:
    print "can' find matplotlib"

#i feel like I shouldn't have to include these line in virtue of __init__.py
from read_inp import parse_inp, read_inp, read_log,metadata
from read_grid import read_grid
from basic_info import basic_info, fft_info
from corral import corral


def show(path='/home/cryosphere/BOUT/examples/Q3/data_short', var='Ni',runcase='runcase.sh'): 
    
    #lets collect the data
    #from . import read_grid,parse_inp, read_inp,basic_info


    #the purpose of this function is to pull out previously 
    #collect linear response data and present all information in  

    #a single plot    
    #collect

    #figure out what to collect from runcase
    



    boutpath = path
    
    print var 
    print path
    print sys.path
   
    ni = boutdata.collect(var,path=path)

    if len(ni.shape) ==4:
        nt,nx,ny,nz = ni.shape
    else:
        print "something with dimesions"

    data = ni[:,:,ny/2,:]
    pp = PdfPages('output.pdf')
    plt.figure()
    cs = plt.contour(data[nt/2,:,:])
    plt.clabel(cs, inline=1, fontsize=10)
    plt.title('Simplest default with labels')
    plt.savefig(pp, format='pdf')

    plt.figure()
    CS = plt.contour(data[0,:,:],6,
                    colors='k', # negative contours will be dashed by default
                    )
    plt.clabel(CS, fontsize=9, inline=1)
    plt.title('Single color - negative contours dashed')
    plt.savefig(pp, format='pdf')
   
    
    a = read_inp(path=path)
 
    b = parse_inp(a)
   
    #print b

   #now we also want to get some information from the grid file

    gridfile = b['[main]']['grid']
    #print b
    #print gridfile
    f = read_grid(gridfile)

    #we can try to bundle some key peices of data into a single meta data
    #container
    

    evolved = []
    #loop over the top level keys and look for any evolve subkeys, see what is evolved
    for section in b.keys():
        if 'evolve' in b[section]:
            if b[section]['evolve']=='true' or b[section]['evolve']=='True':
                evolved.append(section.strip('[]'))
                print evolved
            
    

    meta = OrderedDict()
    meta['dt'] = b['[main]']['TIMESTEP']
    meta['evolved'] = evolved
    meta['DZ'] = b['[main]']['ZMAX']#-b['[main]']['ZMIN']
    

    #now put everything you need into an ordered dictionary
    AA = b['[2fluid]']['AA']
    
    Ni0 = f.variables['Ni0'][:]*1.e14
    bmag = f.variables['bmag'][:]*1.e4 #to cgs
    Ni_x = f.variables['Ni_x'][:]*1.e14 # cm^-3
    rho_s = f.variables['Te_x'][:]/bmag # in cm

    # and so on just follow post_bout.pro, create a sep. function

    #find a nice way to print
    
   
    
    #loop over evoled fields ?
    #find out which are evolved to begin with
    
    output = OrderedDict()
    
    data = OrderedDict()

    all_modes = []
    for i,active in enumerate(meta['evolved']):
       data[active] = boutdata.collect(active,path=path)
       modes,ave = basic_info(data[active],meta)
       output[active] = {'modes':modes,'ave':ave}
       for x in output[active]['modes']:
          if x['k'][0] not in  all_modes:
              all_modes.append(x['k'][0])

    #rebuild with all the relevant modes
    output = OrderedDict()


    for i,active in enumerate(meta['evolved']):
        modes,ave = basic_info(data[active],meta,user_peak = all_modes)
        output[active] = {'modes':modes,'ave':ave}
    
    

    # let's make sure that if there any peaks that only appear
    # in one field can compared across all field
    
    
    
    
    # let pickle the results
    pickled_out = open('post_bout.pkl','wb')
    pickle.dump(output,pickled_out)
    pickled_out.close()

    
   
    plt.figure()
    #ax = plt.add_subplot(2,1,1)
    #ax.set_yscale('log')

    plt.semilogy(modes[0]['amp'].max(1))
    plt.title('simple plot - amp of the 1st dominant mode')
    plt.savefig(pp, format='pdf')

    plt.figure()
    plt.plot(ave['amp'])
    plt.plot(modes[0]['amp'].max(1))
    #plt.semilogy(modes[0]['amp'].max(1))
    plt.title('A_dom and A_ave')
    plt.savefig(pp,format='pdf')

    plt.figure()
    plt.semilogy(ave['amp'])
    plt.semilogy(modes[0]['amp'].max(1))
    #plt.semilogy(modes[0]['amp'].max(1))
    plt.title('A_dom and A_ave')
    plt.savefig(pp,format='pdf')


    plt.figure()
    plt.plot(modes[0]['phase'][3:,nx/2])
    #plt.semilogy(modes[1]['amp'].max(1))
    #plt.semilogy(modes[2]['amp'].max(1))
    #plt.semilogy(modes[0]['amp'].max(1))
    plt.title('phase')
    plt.savefig(pp,format='pdf')

    plt.figure()
    plt.semilogy(modes[0]['amp'][nt/2,:])
    plt.title('A[r] of the 1st mode at nt/2')
    plt.savefig(pp,format='pdf')
    
    #ax = plt.add_subplot(2,1,1)
    #ax.set_yscale('log')

    

    plt.close() 
    pp.close()

    return output

def save(path='/home/cryosphere/BOUT/examples/Q3/data_short'): 
    #lets collect the data
    print 'in post_bout.save'
    boutpath = path

    print path
    print sys.path
 

    meta = metadata(path=path)
   
    print 'ok got meta'
    #return meta
    output = OrderedDict()
    data = OrderedDict()

   
    all_modes = []
    print meta['evolved']['v']

    for i,active in enumerate(meta['evolved']['v']): #loop over the fields
        print path, active
        data[active] = boutdata.collect(active,path=path)
        modes,ave = basic_info(data[active],meta)
       #print modes[0]['gamma']
        output[active] = {'modes':modes,'ave':ave}
        for x in output[active]['modes']: #keep track of relevant harmonics
            if x['k'][0] not in  all_modes: 
                all_modes.append(x['k'][0])

    
    #rebuild with all the relevant modes
    output = OrderedDict()

    output['meta'] = meta
    output['fields'] = {}

    for i,active in enumerate(meta['evolved']['v']): 
        print 'once again', active
        modes,ave = basic_info(data[active],meta,user_peak = all_modes) 
        output['fields'][active] = {'modes':modes,'ave':ave}
        #print modes[0]['gamma']

    # let's make sure that if there any peaks that only appear
    # in one field can compared across all field
    
    
    
    
    # let pickle the results
    filename = path+'/post_bout.jsn'
    json_out = open(filename,'wb')
    pickle.dump(output,json_out)
    json_out.close()

    print filename
    
    
def read(path='.',filename='post_bout.jsn'):
    print 'in post_bout.read()'

    filepath = path+'/'+filename
    print 'is file present: ', os.path.isfile(filepath)
    pkl_file = open(filepath, 'rb')
    #pkl_file = open('/tmp/data_.05/post_bout.pkl', 'rb')
 
    output = pickle.load(pkl_file)
    print type(output)
    
    pkl_file.close()
    
    return output



# def convert_to_dict(object):  #takes a L000000ng time 
#     if isinstance(object, collections.Iterable):
#         final = type(obj)()
#         for i,elem in enumerate(object):
#             print i,elem, type(object).__name__ 
#             if type(object).__name__ == 'ndarray':
#                 print object.shape
#             if type(object).__name__ == 'dict' or  type(object).__name__ == 'OrderedDict':
#                 final[elem] = convert_to_dict(object[elem])
#             else:
#                 final[i] = convert_to_dict(object[i])    
#         return final
#     else:
#         if object.__class__.__name__ == 'ValUnit':
#             d = {}
#             d = { '__class__':obj.__class__.__name__, 
#                   '__module__':obj.__module__,
#                   }
#             d.update(object.__dict__)
            
#             return d
#         else:
#             return objec
        
    

# def convert_to_dict(object):  #takes a L000000ng time 
    
#     for elem in object['meta']:
        
        
  
