import sys
sys.path.append('/home/cryosphere/BOUT/tools/pylib')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
import numpy as np


def ListDictKey(input,key):
   #given a key and a list of dictionaries this method returns an ordered
   #list of requested key values
   output = []
   for x in input:
      try :
         #print x[key]
         output.append(x[key])
      except:
         print 'Key not found'
         return 1

   return output


def ListDictFilt(input,key,valuelist):
   #given a key,value pair and a list of dictionaries this 
   #method returns an ordered list of dictionaries where (dict(key)==value) = True
   #http://stackoverflow.com/questions/5762643/how-to-filter-list-of-dictionaries-with-matching-values-for-a-given-key
   try:
      x = copyf(input,key,valuelist)
      return x
   except:
      return []

def copyf(dictlist, key, valuelist):
      return [dictio for dictio in dictlist if dictio[key] in valuelist]
