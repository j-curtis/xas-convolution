#!/usr/bin/env python
#This is a simple python script that will flip the sign and center the xps spectral function
#Jonathan Curtis
#08/10/15

import numpy as np
import sys 

def main():
  #We count the arguments passed
  argc = len(sys.argv)
  if not argc == 3:
    print "Implementation error. Usage is 'program_name' 'inputxps.dat' 'outputxps.dat' "
    return 1
  
  #If this is not the case, we can continue with processing the files
  xpsin = sys.argv[1]
  xpsout = sys.argv[2]
  
  #Now we load the data using numpy loadtxt
  data = np.loadtxt(str(xpsin))
  
  print "data array has shape ",data.shape
  
  #Now we locate the global maximum
  location = data[0,0]  #We start at the lowest frequency
  value = 0.0
  for i in xrange(data.shape[0]):
    #We record whenever a value is larger
    #We record both the value and the location
    if data[i,2] > value:
      value = data[i,2]
      location = data[i,0]  
  
  print "Maximum of ",value," occurs at ",location
  
  #Now we subtract this and flip the frequencies 
  data[:,0] = -(data[:,0] - location)
  
  #And we print the data
  np.savetxt(xpsout,data,delimiter = " ")
  
  #And we close
  return 0
  
if __name__ == "__main__":
  main()
  
