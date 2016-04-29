# run SU2 with submitting a job
#mpirun -np 4 ../../bin/SU2_CFD lam_airfoil_shedding.cfg

# submit a job on Elwetritsch with 100 cores of the same model and reserve 4000mb memory per node (or per core) for it
bsub -R "rusage[mem=4000]" -a openmpi -n 100 -R "same[model]" mpirun ../../bin/SU2_CFD lam_airfoil_shedding.cfg
