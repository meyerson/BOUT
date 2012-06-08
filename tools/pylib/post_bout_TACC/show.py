#eventually the default path shoudl really be the current directory
#import post_bout
#post_bout.show()
def show(path='/home/cryosphere/BOUT/examples/Q3/data_short', var='Ni'): 

    from . import read_grid,parse_inp, read_inp,basic_info
    import sys
    import pickle
    
    #print list(sys.modules.keys())

    try:
        sys.path.append('/home/cryosphere/BOUT/tools/pylib')
        sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
        sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
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

    
    
    #lets collect the data
    
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
