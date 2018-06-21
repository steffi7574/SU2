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
filenameref = refdir + '/flow_' + '%05i'%(IterRef) + '.dat'
filenamecur = 'flow_' + '%05i'%(IterComp)  + '.dat'
 
# open the files for read only
ufileref=open(filenameref,'r')
ufilecur=open(filenamecur,'r')

# define some constants
EPS = 1E-14

err = 0.0
nenner = 0.0
maxdiff = 0.0

# compute the rel. error for each line
for linenumber, line in enumerate(ufileref):

    # skip the first lines since they contain identifiers
    if (linenumber==0 or linenumber == 1 or linenumber == 2):
        ufilecurline = ufilecur.readline()   # skip the line
        continue

    # parse the line of the reference file
    linelistref  = line.split()  # split the string by separator " "
    if (len(linelistref) < 5):
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


    # add to the error //Pressure!
    for i in range(2,3):
        #print(datacur, dataref, i)
        diff    = abs(datacur[i] - dataref[i])
        test   = abs(diff/dataref[i])
        err    += diff**2
        nenner += dataref[i]**2
        if (test > maxdiff):
            maxdiff = test
            maxline = linenumber
            print(datacur[i], dataref[i], maxdiff, linenumber)
          
# relative errornorm
relerrnorm = sqrt(err / nenner)

#OUTPUT
print("\nComparing reference ", filenameref, " with ", filenamecur, "\n")
print("Scanned lines: ", linenumber)

#print("TOTAL SUM REL. ERROR:  %1.12E " %(avgerr))
print("MAX. REL. ERROR:    %1.12E at line " %(maxdiff), (maxline))
print("\nREL. ERROR:         %1.12E \n" %(relerrnorm))

ufileref.close()
ufilecur.close()

