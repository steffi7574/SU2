#!/usr/bin/env python
#
# Script for comparing a flow solution with the reference solution
# USAGE: python3 compareflow.py <Reference Direction> <Number of Time Steps Of the Reference>
#

import sys
from numpy import *


def getdata(ufile, linenumber):
# returns the data of line <linenumber> of the file

  # parse the line of the file
  for i, line in enumerate(ufile):
    if i == linenumber:
      #ufileline = ufile.readline()   # read the line
      linelist  = line.split()  # split the string by separator " "
    elif i > linenumber-1:
      break

  datalist=[]
  for datastring in linelist:
    datalist.append(float(datastring))  

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

# Get the file names
filenameref = refdir + '/restart_flow_' + '%05i'%(IterRef) + '.dat'
filenamecur = 'restart_flow_' + '%05i'%(IterComp)  + '.dat'
 
# open the files for read only
ufileref=open(filenameref,'r')
ufilecur=open(filenamecur,'r')

# define some constants
EPS = 1E-14

# initialize
avgerr  = 0.0
maxerr  = -1.0

# compute the rel. error for each line
#for linenumber in range(201,203):
for linenumber, line in enumerate(ufileref):

    # skip the first lines since they contain identifiers
    if (linenumber==0):
        ufilecurline = ufilecur.readline()   # skip the line
        continue

    # parse the line of the reference file
    linelistref  = line.split()  # split the string by separator " "
    if (linelistref[0] == 'EXT_ITER='):
        break
    # get data of that line
    dataref=[]
    for datastring in linelistref:
        dataref.append(float(datastring))  

    # parse line of current file
    ufilecurline = ufilecur.readline()   # read the line
    linelistcur  = ufilecurline.split()  # split the string by separator " "
    # get data of that line
    datacur=[]
    for datastring in linelistcur:
        datacur.append(float(datastring))  


    # read data from reference and current file
    #dataref     = getdata(ufileref, linenumber)
    #datacur     = getdata(ufilecur, linenumber)
    
    # compute the relative error and its norm
    relerr = []
    relerrnorm = 0.0
    #print(" Reference           Current             rel. ERROR ")
    for i in range(len(dataref)):

          # compute the relative error
          if ( abs(dataref[i]) > EPS ):
            relerr.append( ( datacur[i] - dataref[i] ) / dataref[i] )
          else:
            relerr.append( datacur[i] - dataref[i] )
          
          # print the relative error and the values
          #print( "% 1.12E % 1.12E % 1.12E " %(dataref[i], datacur[i], relerr[i]) )
        
          # add to the error norm
          relerrnorm += relerr[i]**2

          if (relerr[i] > maxerr):
              maxerr  = relerr[i]
              maxts   = dataref[0]
              print(" LINE %06d % 1.12E % 1.12E % 1.12E " %(linenumber, dataref[i], datacur[i], relerr[i])  )

    relerrnorm = sqrt(relerrnorm)
    avgerr += relerrnorm
    #if (relerrnorm > maxerr):
    #    maxerr  = relerrnorm
    #    maxts   = dataref[0]
    #    print(" LINE ", linenumber, ", Rel. error norm: %1.12E" %(relerrnorm) )

    # print the norm of the relative error
    #print(" LINE ", linenumber, ", Rel. error norm: %1.12E" %(relerrnorm) )

#OUTPUT
print("\nComparing reference ", filenameref, " with ", filenamecur, "\n")
print("Scanned lines: ", linenumber)
print("TOTAL SUM REL. ERROR:  %1.12E " %(avgerr))
print("AVERAGE REL. ERROR:    %1.12E \n" %(avgerr/linenumber))
print("MAXIMUM REL. ERROR:    %1.12E " %(maxerr), " at line ", maxts, "\n")

ufileref.close()
ufilecur.close()

