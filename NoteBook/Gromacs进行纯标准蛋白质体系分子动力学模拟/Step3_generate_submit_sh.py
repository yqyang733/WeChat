import os
import shutil

class config:

    def __init__(self):

        self.queue = "quick"
        self.name = "MD"

def submit_sh(queue, name):
    
    submitsh = open(os.path.join(".", "job.sh"), "w")
    submitsh.write(
'''#!/bin/bash
#SBATCH -J {0}
#SBATCH -p {1}
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Decide the software version
source /public/software/profile.d/apps_gromacs_2023.2.sh

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#############################
mkdir em
cd em
if [ ! -f em.gro ]; then
    gmx grompp -f ../mdp/em.mdp -c ../solv_ions.gro -p ../topol.top -r ../solv_ions.gro -o em.tpr -maxwarn 2
    gmx mdrun -s em.tpr -deffnm em -ntmpi 1 -nb gpu -gpu_id 0
fi
###########################
mkdir ../nvt
cd ../nvt
if [ ! -f nvt.gro ]; then
    gmx grompp -f ../mdp/nvt.mdp -c ../em/em.gro -p ../topol.top -o nvt.tpr -r ../em/em.gro -maxwarn 4 -n ../index.ndx
    gmx mdrun -s nvt.tpr -deffnm nvt -ntmpi 1 -nb gpu -bonded gpu -pme gpu -gpu_id 0
fi
##########################
mkdir ../npt
cd ../npt
if [ ! -f npt.gro ]; then
    gmx grompp -f ../mdp/npt.mdp -c ../nvt/nvt.gro -t ../nvt/nvt.cpt -p ../topol.top -o npt.tpr -r ../nvt/nvt.gro -maxwarn 4 -n ../index.ndx
    gmx mdrun -s npt.tpr -deffnm npt -ntmpi 1 -nb gpu -bonded gpu -gpu_id 0 -pme gpu
fi
###################################
mkdir ../prod
cd ../prod
if [ ! -f prod.gro ]; then
    gmx grompp -f ../mdp/prod.mdp -c ../npt/npt.gro -t ../npt/npt.cpt -p ../topol.top -o prod.tpr -r ../npt/npt.gro -maxwarn 4 -n ../index.ndx
    gmx mdrun -s prod.tpr -deffnm prod -dhdl dhdl -ntmpi 1 -nb gpu -bonded gpu -gpu_id 0 -pme gpu  -nsteps 50000000
fi
###################################
cd ..                                                                                                                                        
'''.format(name, queue)
    )

def run():

    settings = config()
    submit_sh(settings.queue, settings.name)                                                                                                                                                         

def main():

    run()
    
if __name__=="__main__":
    main() 