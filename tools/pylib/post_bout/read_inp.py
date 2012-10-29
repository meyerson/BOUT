from read_grid import read_grid
from ordereddict import OrderedDict
import numpy as np

def read_inp(path='',boutinp='BOUT.inp'):
   
   
   boutfile = path+'/'+boutinp
   boutinp = open(boutfile,'r').readlines()
   

   # start by stripping out all comments 
   # look at the 1st character of all list elements
   # for now use a gross loop, vectorize later
   boutlist = []
  
   for i,val in enumerate(boutinp):
      if val[0] != '#' and val.isspace() == False:
         boutlist.append(val.split("#")[0])

   return boutlist

def parse_inp(boutlist):

   import re
   from ordereddict import OrderedDict
   
   if not boutlist:
      return 0

 
   
   #boutdic={} unordered standard dict
   boutdic = OrderedDict()
   
   #regex is messy see http://docs.python.org/howto/regex.html#regex-howto
   pattern = "\[\S*\]" #pattern for a section start [All],[Ni], etc
   
   pattern = re.compile(pattern)
         
   boutdic['[main]'] = {}
   current='[main]'

   for i,val in enumerate(boutlist):
      print i,val
      result =pattern.match(val)
      #while the current value is not a new section name add everything to the current section
      
      if result is None:
         #print val
         key,value = val.split("=")
         value = value.replace('\"','')
         print current, key,value

         
         boutdic[current][key.strip()] = value.strip()
      else:
         boutdic[result.group()] = {}
         current = result.group()
      
   return boutdic

def read_log(path='.',logname='status.log'):
   
   print 'in read_log'
   import re
   from ordereddict import OrderedDict
   
   #logfile = path+'/'+logname
   logfile = logname
   print logfile
   logcase = open(logfile,'r').readlines()
   
   # start by stripping out all comments 
   # look at the 1st character of all list elements
   # for now use a gross loop, vectorize later
   loglist = []
  
   for i,val in enumerate(logcase):
      if val[0] != '#' and val.isspace() == False:
         loglist.append(val.split("#")[0])
         
   if not loglist:
       return 0

   logdict = OrderedDict()
   logdict['runs'] = []
   #print len(loglist)
   print loglist
   #print loglist[len(loglist)-1] == 'last one\n'
   
   # last = loglist.pop().rstrip()
   
   # logdict['done'] = last == 'done'

  
   logdict['current'] = loglist.pop().rstrip()
   for i,val in enumerate(loglist):
      print val
      logdict['runs'].append(val.rstrip())
      
   logdict['runs'].append(logdict['current'])
  

   #print logdict
   return logdict
   
def metadata(inpfile='BOUT.inp',path ='.',v=False):    
    filepath = path+'/'+inpfile
    print filepath
    inp = read_inp(path=path,boutinp=inpfile)
    inp = parse_inp(inp)
    
    try:
       gridname = inp['[main]']['grid']
       IC = read_grid(gridname) #have to be an ansoulte file path for now
       print 'IC: ',type(IC)
       # print IC.variables
       # print gridname
    except:
       print gridname
       print 'Fail to load the grid file'
       

    #print gridname
    #print len(IC)
    #print IC
       
    evolved = []
    ICscale = []
   
    fieldkeys = ['[Ni]','[jpar]','[Te]','[Ti]','[Vi]','[rho]']
    
    defaultIC = float(inp['[All]'].get('scale',0.0))

    for section in inp.keys(): #loop over section keys 
       #print section
       if section in fieldkeys: #pick the relevant sections
          print section
          #print type(inp[section].get('evolve','True'))
          print (inp[section].get('evolve','True')).lower().strip()
          if inp[section].get('evolve','True').lower().strip() == 'true':
              print 'ok reading'
              evolved.append(section.strip('[]'))
              ICscale.append(float(inp[section].get('scale',defaultIC)))
             
           
   

                
    meta = OrderedDict()
    
    class ValUnit(object):
       def __init__(self,value=0,units=''):
          self.u = units
          self.v = value
       def todict(self):
          return {'u':self.u,'v':self.v}

   

    #def decode_valunit(d):
       
    
    def ToFloat(metaString):
       try:
          return float(metaString)
       except ValueError:
          return metaString
    
    meta['evolved'] = ValUnit(evolved,'')
    meta['IC']= np.array(ICscale)
    d = {}

    # read meta data from .inp file
    for section in inp.keys():
        if (('evolve' not in inp[section]) and ('first' not in inp[section])): #hacky way to exclude some less relevant metadata
           for elem in inp[section].keys():
              meta[elem] = ValUnit(ToFloat(inp[section][elem]))
              d[elem] = np.array(ToFloat(inp[section][elem]))
              
    #read in some values from the grid(IC) and scale them as needed
    norms = {'Ni0':ValUnit(1.e14,'cm^-3'),'bmag':ValUnit(1.0e4,'gauss'),
             'Ni_x':ValUnit(1.e14,'cm^-3'),
             'Te_x':ValUnit(1.0,'eV'),'Ti_x':ValUnit(1.0,'eV'),'Rxy':ValUnit(1,'m'),
             'Bxy':ValUnit(1.0e4,'gauss'),'Bpxy':ValUnit(1.0e4,'gauss'),
             'Btxy':ValUnit(1.0e4,'gauss'),'Zxy':ValUnit(1,'m'),
             'dlthe':ValUnit(1,'m'),'dx':ValUnit(1,'m'),'hthe0':ValUnit(1,'m')}

    for elem in norms.keys():
       #print 'elem: ',elem
       meta[elem] = ValUnit(IC.variables[elem][:]*norms[elem].v,norms[elem].u)
       d[elem] = np.array(IC.variables[elem][:]*norms[elem].v)
    
  

    #if case some values are missing   
    default = {'bmag':1,'Ni_x':1,'NOUT':100,'TIMESTEP':1,
               'MZ':32,'AA':1,'Zeff':ValUnit(1,''),'ZZ':1,
               'zlowpass':0.0}
    diff = set(default.keys()).difference(set(d.keys()))
       
    for elem in diff:
       meta[elem] = default[elem]
       d[elem] = np.array(default[elem])

    nx,ny  = d['Rxy'].shape
    #compute some quantities that are usefull
        
    print meta['AA'].v
    

    if meta['zlowpass'].v != 0:
       print meta['MZ'].v, meta['zlowpass'].v
       meta['maxZ'] = np.floor(meta['MZ'].v*meta['zlowpass'].v)
    else:
       meta['maxZ'] = 5

    meta['nx'] = nx
    meta['ny']= ny
    meta['dt'] = meta['TIMESTEP'] 
    
    meta['rho_s'] = ValUnit(1.02e2*np.sqrt(d['AA']*d['Te_x'])/(d['ZZ']* d['bmag']),'cm')   # ion gyrorad at T_e, in cm 
    meta['rho_i'] = ValUnit(1.02e2*np.sqrt(d['AA']*d['Ti_x'])/(d['ZZ']* d['bmag']),'cm') 
    meta['rho_e'] = ValUnit(2.38*np.sqrt(d['Te_x'])/(d['bmag']),'cm') 
    
    meta['fmei']  = ValUnit(1./1836.2/d['AA'])   
    
    meta['lambda_ei'] = 24.-np.log(np.sqrt(d['Ni_x'])/d['Te_x']) ;
    meta['lambda_ii'] = 23.-np.log(d['ZZ']**3 * np.sqrt(2.*d['Ni_x'])/(d['Ti_x']**1.5)) #

    meta['wci']       = 1.0*9.58e3*d['ZZ']*d['bmag']/d['AA'] # ion gyrofrteq
    meta['wpi']       = 1.32e3*d['ZZ']*np.sqrt(d['Ni_x']/d['AA']) # ion plasma freq 

    meta['wce']       = 1.78e7*d['bmag'] #electron gyrofreq
    meta['wpe']       = 5.64e4*np.sqrt(d['Ni_x'])#electron plasma freq
    
    meta['v_the']    = 4.19e7*np.sqrt(d['Te_x'])#cm/s
    meta['v_thi']    = 9.79e5*np.sqrt(d['Ti_x']/d['AA']) #cm/s
    meta['c_s']      = 9.79e5*np.sqrt(5.0/3.0 * d['ZZ'] * d['Te_x']/d['AA'])#
    meta['v_A']     = 2.18e11*np.sqrt(1.0/(d['AA'] * d['Ni_x']))
    
    meta['nueix']     = 2.91e-6*d['Ni_x']*meta['lambda_ei']/d['Te_x']**1.5 #
    meta['nuiix']     = 4.78e-8*d['ZZ']**4.*d['Ni_x']*meta['lambda_ii']/d['Ti_x']**1.5/np.sqrt(d['AA']) #
    meta['nu_hat']    = meta['Zeff'].v*meta['nueix']/meta['wci'] 
    
    meta['L_d']      = 7.43e2*np.sqrt(d['Te_x']/d['Ni_x'])
    meta['L_i_inrt']  = 2.28e7*np.sqrt(d['AA']/d['Ni_x'])/ d['ZZ'] #ion inertial length in cm
    meta['L_e_inrt']  = 5.31e5*np.sqrt(d['Ni_x']) #elec inertial length in cm
    
    meta['Ve_x'] = 4.19e7*d['Te_x']

    
    meta['R0'] =  (d['Rxy'].max()+d['Rxy'].min())/2.0 
    
 
    print d['Rxy'].mean(1) 
    print d['ZMAX']
    print  d['ZMIN'] 
    meta['L_z'] = 1e2 * 2*np.pi * d['Rxy'].mean(1) *(d['ZMAX'] - d['ZMIN']) # in cm toroidal range
    meta['dz'] = (d['ZMAX'] - d['ZMIN'])
 
    #meta['lbNorm']=meta['L_z']*(d['Bpxy']/d['Bxy']).mean(1)     #-binormal coord range [cm]
    meta['lbNorm']=meta['L_z']*(d['Bxy']/d['Bpxy']).mean(1)
    
    #meta['zPerp']=np.array(meta['lbNorm']).mean*np.array(range(d['MZ']))/(d['MZ']-1) 
  #let's calculate some profile properties
    dx = np.gradient(d['Rxy'])[0]
    meta['L'] = 1.0*1e2*dx*(meta['Ni0'].v)/np.gradient(meta['Ni0'].v)[0]/meta['rho_s'].v
    
    meta['w_Ln']     =  meta['c_s']/(np.min(abs(meta['L']))*meta['wci'] *meta['rho_s'].v) #normed to wci

    AA = meta['AA'].v
    ZZ = d['ZZ']
    Te_x = d['Te_x']
    Ti_x = d['Ti_x']
    fmei = meta['fmei'].v
    
    meta['lpar'] =1e2*((d['Bxy']/d['Bpxy'])*d['dlthe']).sum(1)/meta['rho_s'].v #-[normed], average over flux surfaces, parallel length
    #meta['lpar']=1e2*(d['Bxy']/d['Bpxy']).mean(1)*d['dlthe'].mean(1) #function of x
    meta['sig_par'] = 1.0/(fmei*0.51*meta['nu_hat'])
    #meta['heat_nue'] = ((2*np.pi/meta['lpar'])**2)/(fmei*meta['nu_hat'])
    #kz_e = kz_i*(rho_e/rho_i)
    # kz_s = kz_i*(rho_s/rho_i)
    # kz_i = (TWOPI/L_z)*(indgen((*current_str).fft.nz+1))*rho_i
    
    # knorm = (TWOPI/lbNorm)*(indgen((*current_str).fft.nz+1))*rho_s

    # for now just translate
    for elem in meta:
       if type(meta[elem]).__name__ =='ValUnit':
          meta[elem] = {'u':meta[elem].u,'v':meta[elem].v}
    
       
    print 'meta: ', type(meta)
    return meta

    # meta['DZ'] =inp['[main]']['ZMAX']#-b['[main]']['ZMIN']
    # AA = inp['[2fluid]']['AA']
    # Ni0 = IC.variables['Ni0'][:]*1.e14
    # bmag = IC.variables['bmag'][:]*1.e4 #to cgs
    # Ni_x = IC.variables['Ni_x'][:]*1.e14 # cm^-3
    # Te_x
    
    # rho_s = 1.02e2*sqrt(AA.v*Te_x.v)/ZZ.v/bmag.v
    # rho_i
    # rho_e 
    
    

    

   #for i,val in enumerate(boutlist):
