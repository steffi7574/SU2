#!/usr/bin/env python
#
# Script for comparing a flow solution with the reference solution
# USAGE: python3 compareflow.py <Reference Direction> <Number of Time Steps Of the Reference>
#

import sys
from numpy import *


def getdata(filename):
# returns the data of the second line of the file

  # open file for read only
  ufile=open(filename,'r')

  # parse the second line of the file
  ufile.readline()        # skip the first line (header)
  line=ufile.readline()   # read the sond line
  linelist = line.split() # split the string by separator " "


  datalist=[]
  for datastring in linelist:
    datalist.append(float(datastring))  

    ufile.close()
  return datalist

# parse command line arguments
args = sys.argv
if (len(args) != 3):
  print("\nERROR while parsing argument list")
  print("Usage: python3 compareflow.py <ReferenceDirection> <NumberOfTimeStepOfReference>\n")
  stop
else:
  refdir = str(args[1])  # Set the reference directory
  Iter   = int(args[2])  # Set the number of time step that you want to compare

# define some constants
EPS = 1E-14

# Get the file names
filenameref = refdir + '/restart_flow_' + '%05i'%(Iter) + '.dat'
filenamecur = 'restart_flow_' + '%05i'%(Iter+1)  + '.dat'
print("Comparing reference ", filenameref, " with ", filenamecur)

# read data from reference and current file
dataref     = getdata(filenameref)
datacur     = getdata(filenamecur)

# compute the relative error and its norm
relerr = []
relerrnorm = 0.0
print(" Reference           Current             ERROR ")
for i in range(len(dataref)):
  #print( "%1.12E" %(dataref[i]))
  #print( "%1.12E" %(datacur[i]))

  # compute the relative error
  if ( abs(dataref[i]) > EPS ):
    relerr.append( ( datacur[i] - dataref[i] ) / dataref[i] )
  else:
    relerr.append( datacur[i] - dataref[i] )
  
  # print the relative error and the values
  print( "% 1.12E % 1.12E % 1.12E " %(dataref[i], datacur[i], relerr[i]) )

  # add to the error norm
  relerrnorm += relerr[i]**2
relerrnorm = sqrt(relerrnorm)

# print the norm of the relative error
print( "\n ERROR NORM: %1.12E \n" %(relerrnorm) )


