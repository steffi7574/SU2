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

import os, sys, copy, subprocess
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
                       
    (options, args)=parser.parse_args()
  
    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")
    
    batch_job( options.filename)
        
#: def main()


batch_base_command = "sbatch"

batch_args_mapping_slurm = {"NODES" : "--nodes",
                            "NAME"  : "--job-name",
                            "ERROR" : "--error",
                            "OUTPUT" : "--output",
                            "TIME"  : "--time", 
                            "CONSTRAINT" : "--constraint",
                            "TASKS_PER_NODE" : "--ntasks-per-node",
                            "MEMORY": "--mem"}

# use slurm argument names
batch_args_mapping = batch_args_mapping_slurm 

default_batch_args = {batch_args_mapping["NODES"]         : 1,
                      batch_args_mapping["NAME"]          : "default",
                      batch_args_mapping["OUTPUT"]        : "default-%j.out",
                      batch_args_mapping["ERROR"]         : "default-%j.err",
                      batch_args_mapping["TASKS_PER_NODE"]: 16,
                      batch_args_mapping["MEMORY"]        : 30000,
                      batch_args_mapping["TIME"]          : "01:00:00",
                      batch_args_mapping["CONSTRAINT"]    : "XEON_E5_2670"}

run_base_command = "python ../SU2_xbraid/SU2_PY/parallel_time_computation.py"

default_run_args    = {"-f"  : "default.cfg",
                        "-o"  : "TEST", 
                        "-x": 8,
                        "-t": 10}

# -------------------------------------------------------------------
#  CFD Solution
# -------------------------------------------------------------------

def batch_job( filename ):

    # Read the base config file
    config = SU2.io.Config(filename)

    # Total number of processes
    config.NUMBER_PART      = 16

    # number of time processes
    config.BRAID_NPROC_TIME = 2

    # modify more config options here

    # assemble a jobname
    jobname = "npx" + str(config.NUMBER_PART/config.BRAID_NPROC_TIME) + "npt" + str(config.BRAID_NPROC_TIME)

    # submit the job
    submit_job(config, jobname, "02:00:00")


def submit_job(config, jobname, time_limit):

    config.dump(jobname+".cfg")

    run_args   = copy.deepcopy(default_run_args)
    batch_args = copy.deepcopy(default_batch_args)

    # set number of space processes
    run_args["-x"] = config.NUMBER_PART/config.BRAID_NPROC_TIME

    # set number of time processes
    run_args["-t"] = config.BRAID_NPROC_TIME

    # set the folder to execute
    run_args["-o"] = jobname

    # set the config file
    run_args["-f"] = jobname+".cfg"
    
    batch_args[batch_args_mapping["NODES"]]           = config.NUMBER_PART/batch_args[batch_args_mapping["TASKS_PER_NODE"]]
    batch_args[batch_args_mapping["NAME"]]            = jobname
    batch_args[batch_args_mapping["ERROR"]]           = jobname+"%j.err"
    batch_args[batch_args_mapping["OUTPUT"]]          = jobname+"%j.out"
    batch_args[batch_args_mapping["TIME"]]            = time_limit
    

    command =  assemble_command(run_base_command, run_args, " ")
    assemble_batch_script(jobname+".batch", command, batch_args)

    subprocess.call("sbatch " + jobname + ".batch", shell=True)

def assemble_command(command, args, sep):

    ex_command = command

    for arg, value in args.iteritems():
        ex_command = ex_command + str(" ") + arg + sep + str(value)

    return ex_command

def assemble_batch_script(name, run_command, args):
    
    outfile = open(name, 'w')

    outfile.write("#!/usr/bin/bash\n")

    for arg,value in args.iteritems():
        outfile.write("#SBATCH " + arg + "=" + str(value) + "\n")

    outfile.write("export SU2_MPI_COMMAND=\"mpirun -n %i %s\"\n")
    outfile.write(run_command)
    outfile.close()



#: parallel_computation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

