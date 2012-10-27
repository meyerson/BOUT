from read_grid import read_grid
from ordereddict import OrderedDict
import numpy as np
import string
import re


def findlowpass(cxxstring):
    ##p1="lowPass\(.*\)"
    p1="lowPass\(.*\,(.*)\)"
    maxN = np.array(re.findall(p1, cxxstring)) 
    #print substrings

    #p2=("[0-9]")

    #maxN = np.array([re.findall(p2,elem) for elem in substrings]).flatten()
    
    if maxN.size == 0:
        return 100
    else:
        output = int(min(maxN))     
        if output ==0:
            return 20
        else:
            return output
def read_cxx(path='.',boutcxx='q3_simp.cxx',evolved=''):
    
    
    #print path, boutcxx
    boutcxx = path+'/'+boutcxx
   #boutcxx = open(boutcxx,'r').readlines()
    f = open(boutcxx,'r')
    boutcxx = f.read()
    f.close()

   
   
   # start by stripping out all comments 
   # look at the 1st character of all list elements
   # for now use a gross loop, vectorize later

    start =string.find(boutcxx,'/*')
    end =string.find(boutcxx,'*/')+2

    s = boutcxx[0:start]
    for i in range(string.count(boutcxx,'/*')):
        start =string.find(boutcxx,'/*',end)
        s=s+boutcxx[end+1:start-1]
        
        end =string.find(boutcxx,'*/',end)+2
        
        s=s+boutcxx[end+1:] 

    section_0 = string.find(s,'int physics_run(BoutReal t)')
    section_1 = string.find(s,'return',section_0)
    s =  s[section_0:section_1]
   
    
    tmp = open('./read_cxx.tmp','w')
    tmp.write(s)
    tmp.close()
    tmp = open('./read_cxx.tmp','r')
    
    cxxlist = ''

    for line in tmp:
        if line[0] != '//' and line.isspace()==False:
            cxxlist=cxxlist+line.split("//")[0]
            
    return cxxlist
