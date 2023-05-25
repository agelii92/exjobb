simulation_setup.txt: modify this document to design the simulation parameters and  values

prepare.py: This script will setup a sub-directory 'runs_a-b' where a,b are taken from simulation_setup.txt. In the sub-directory, for i in [a,b], run_i.gen files will be generated, with modifications that can be custom made by modifying setup_simulation.py. Additionally, for each .gen file a run_i_parameters.csv file will be generated including the parameters that define the properties of the run_i.gen file. For example, if setup_simulation.py is modified to create .gen files that randomize a sample density, the .csv will include the randomized density values. Finally, a davinci_runs_a-b.sh file is generated which will execute the cretin simulation when sent to the que system slurm via bash command 'sbatch davinci_runs_a-b.sh'

get_database.sh: This script will create another sub-directory 'database' containing a directory 'files' with run_i.hdf5 files containing the simulaton results. Additionally, it is possible to create a virtual database of these files by running the script 'create_virtual_database.py' which will be copied to the 'database' directory.

create_database.py: does the same as 'get_database.sh' but on the login node. For large simulation datasets, use the slurm que system instead and ignore this script.

functions.py: contain helper functions used by 'create_database.py' and 'create_virtual_database.py'

davinci.sh: template for sending simulations to slurm system, modyfied by 'setup_simulation.py'

template.gen: template for creating .gen files, used by 'setup_simulation.py'
