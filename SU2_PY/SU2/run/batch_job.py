import os, sys, shutil, copy
import subprocess
from ..io import Config

def submit_job(config, jobname, time_limit, run_base_command, run_args, batch_base_command, batch_args, batch_args_mapping):

    config.dump(jobname+".cfg")
    
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
