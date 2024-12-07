import os
import sys

def ions_mdp():
    
    ions_mdp = open("ions.mdp", "w")
    ions_mdp.write(
'''; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
'''
    )

def run(input_pdb):

    ions_mdp()

    run_sh = open("build.sh", "w")
    run_sh.write(
'''gmx pdb2gmx -f {0} -o build.gro -water tip3p -ignh << -EOF
    2
-EOF

gmx editconf -f build.gro -o newbox.gro -bt cubic -d 0.8
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc 0.15 << -EOF
    SOL
-EOF

gmx make_ndx -f solv_ions.gro -o index.ndx << -EOF
    q
-EOF
'''.format(input_pdb)
    )

    command = "sh build.sh"
    os.system(command)

def main():

    input_pdb = sys.argv[1]
    run(input_pdb)
    
if __name__=="__main__":
    main() 