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

  # parse the line of the file
  for i, line in enumerate(ufile):
    if i == LINENUMBER-1:
      ufileline = ufile.readline()   # read the line
      linelist  = ufileline.split()  # split the string by separator " "
    elif i > LINENUMBER -1:
      break

  datalist=[]
  for datastring in linelist:
    datalist.append(float(datastring))  

    ufile.close()
  return datalist

# parse command line arguments
args = sys.argv
if (len(args) == 3):
  refdir     = str(args[1])  # Set the reference directory
  IterRef    = int(args[2])  # Set the number of time step of the reference that you want to compare
  IterComp   = IterRef 
elif (len(args) == 4):
  refdir     = str(args[1])  # Set the reference directory
  IterRef    = int(args[2])  # Set the number of time step of the reference
  IterComp   = int(args[3])  # Set the number of time step that you want to compare to
else:
  print("\nERROR while parsing argument list !")
  print("Usage: python3 compareflow.py <ReferenceDirection> <NumberOfTimeStepOfReference> <NumberOfTimeStopToCompare>\n")
  stop

# define some constants
LINENUMBER= 201
EPS = 1E-14

# Get the file names
filenameref = refdir + '/restart_flow_' + '%05i'%(IterRef) + '.dat'
filenamecur = 'restart_flow_' + '%05i'%(IterComp)  + '.dat'
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


