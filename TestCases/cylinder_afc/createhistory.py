#!/usr/bin/env python
#
# Script for putting together the history.dat file, collecting data from all temporal processors. 
# USAGE: python3 createhistory.py <FileNameWithoutLast"_000*.dat"> <NumberOfTemporalProcessors>
#

import sys
from numpy import *

# parse command line arguments
args = sys.argv
if (len(args) == 3):
  filename   = str(args[1])  # Get the filenames without the last "_000*.dat"
  npt        = int(args[2])  # Get the number of temporal processors
else:
  print("\nERROR while parsing argument list !")
  print("Usage: python3 createhistory.py <FileNameWithoutLast'_000*.dat'> <NumberOfTemporalProcessors\n")
  stop

# Open the output filename
outfile = open(filename + '.dat', 'w')

# Iterate over all local history files
for inpt in range(npt):
    
    # Open the local history file and read it
    currfile = open(filename + '_%05i'%(inpt) + '.dat','r') 
    currfilelines = currfile.readlines()
    currfilesize  = size(currfilelines)

    # Sort the data with respect to first column
    currfilelines.sort()

    # If first local file: copy the header
    if (inpt == 0):
        for lineid in range(currfilesize - 3, currfilesize):
            header = currfilelines[lineid]
            outfile.write(header)
    
    # Copy the rest to the global history file
    for lineid in range(currfilesize - 3):
        outfile.write(currfilelines[lineid])

    # close local history file
    currfile.close()


# Close the global history file
outfile.close()

