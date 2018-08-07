#!/usr/bin/env python 

## \file parallel_computation.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
#  \author T. Economon, T. Lukaczyk, F. Palacios
#  \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os, sys
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-x", "--n_process_space", dest="npx", 
                      help="number of processes in space", metavar="PROCESSES_SPACE")
    parser.add_option("-t", "--n_process_time", dest="npt",
                     help="number of processes in time", metavar="PROCESSES_TIME")
    parser.add_option("-o", "--folder", dest="folder", help="Folder to run the computation", metavar="FOLDER")
                       
    (options, args)=parser.parse_args()
    options.npx         = int (options.npx)
    options.npt         = int (options.npt)

    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")
    
    parallel_time_computation( options.filename    ,
                          options.npx  ,
                          options.npt  ,
                          options.folder)
        
#: def main()


# -------------------------------------------------------------------
#  CFD Solution
# -------------------------------------------------------------------

def parallel_time_computation( filename           , 
                               npx, 
                               npt, 
                               folder):
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART      = npt*npx
    config.BRAID_NPROC_TIME = npt
    config.NZONES = 1
    
    # State
    state = SU2.io.State()
    if config.RESTART_SOL == 'YES':
        state.FILES.DIRECT= config.SOLUTION_FLOW_FILENAME

    state.FILES.MESH = config.MESH_FILENAME    
    
    funcs = SU2.eval.aerodynamics(config, state, folder )
    
#: parallel_computation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

