import os
import time
import shutil

class config:

    def __init__(self):

        self.title = ""
        self.queue = "single"
        self.system_pdb = os.path.join("complex.pdb")
        self.system_psf = os.path.join("complex.psf")
        self.vmd_path = os.path.join("~/../../software/apps/vmd/1.9.3/vmd")
        self.ff_path = os.path.join("/public/home/yqyang/file/vegf-toppar")
        self.npt_step = "100000"                       # timestep: 2fs
        self.prod_step = "500000"                      # timestep: 2fs  
        self.constraints = "protein and noh"

def submit(title, queue):
    
    mk_submit = open("do_md.sh", "w")
    mk_submit.write(
'''#!/bin/bash
#SBATCH -J {0}
#SBATCH -p {1}
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Decide the software version
source /public/software/profile.d/apps_namd_3.0alpha9.sh

NAMD="namd3 +p1 +devices 0"

## min
cd prod
base_em=em.config
$NAMD $base_em > $base_em.log

## nvt
base_nvt=nvt.config
$NAMD $base_nvt > $base_nvt.log

## npt
base_nptstep=npt.config
$NAMD $base_nptstep > $base_nptstep.log

## prod
base_prodstep=prod.config
$NAMD $base_prodstep > $base_prodstep.log
'''.format(title, queue)
    )
    mk_submit.close()

def md_pbc_box(vmd_path, pdb, psf):
    mk_pbcbox = open("mk_pbcbox.tcl", "w")
    mk_pbcbox.write(
"""
#!/bin/bash
# vmd -dispdev text -e mk_pbcbox.tcl

package require psfgen
psfcontext reset
mol load psf {1} pdb {0}
set everyone [atomselect top all]
set minmax [measure minmax $everyone]
foreach {{min max}} $minmax {{ break }}
foreach {{xmin ymin zmin}} $min {{ break }}
foreach {{xmax ymax zmax}} $max {{ break }}

set file [open "PBCBOX.dat" w]
puts $file "cellBasisVector1 [ expr $xmax - $xmin ] 0 0 "
puts $file "cellBasisVector2 0 [ expr $ymax - $ymin ] 0 "
puts $file "cellBasisVector3 0 0 [ expr $zmax - $zmin ] "
puts $file "cellOrigin [ expr ($xmax + $xmin)/2 ] [ expr ($ymax + $ymin)/2 ] [ expr ($zmax + $zmin)/2 ] "

exit
""".format(pdb, psf)
    )
    mk_pbcbox.close()
    cmd = vmd_path + " -dispdev text -e mk_pbcbox.tcl"
    os.system(cmd)
    time.sleep(1)

def position_constraints(pdb, vmd_path, constraints):

    mk_constraints = open("constraints.tcl", "w")
    mk_constraints.write(
'''
mol new {0} type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "{1}"]
$sel set beta 1
$all writepdb constraints.pdb

quit
'''.format(pdb, constraints)
    )
    mk_constraints.close()
    cmd = vmd_path + " -dispdev text -e constraints.tcl"
    os.system(cmd)
    time.sleep(1)

def generate_em_config(pdb, psf, ff_path, pbcbox):
    
    mk_em_config = open("em.config", "w")
    mk_em_config.write(
'''
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# Minimization
# namd3 +p1 em.config > em.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          {1}
coordinates        {0}
set outputbase     com

firsttimestep      0

# open names all, later will control
set ITEMP 310
set FTEMP 310
# if you do not want to open this option, assign 0
set INPUTNAME   0                      ;# use the former outputName, for restarting a simulation
set PSWITCH     1                      ;# whether to use langevinPiston pressure control
set FIXPDB      0
set CONSPDB     0
set CONSSCALE   1                      ;# default; initial value if you want to change
set parpath     {2}

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${{parpath}}/addition.prm
parameters          ${{parpath}}/par_all36m_prot.prm
parameters          ${{parpath}}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if {{ $INPUTNAME != 0 }} {{
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
}} else {{
    # Periodic Boundary Conditionsc
    temperature $ITEMP
    source {3}
}}

## Force-Field Parameters
exclude             scaled1-4;         # non-bonded exclusion policy to use "none,1-2,1-3,1-4,or scaled1-4"
                                    # 1-2: all atoms pairs that are bonded are going to be ignored
                                    # 1-3: 3 consecutively bonded are excluded
                                    # scaled1-4: include all the 1-3, and modified 1-4 interactions
                                    # electrostatic scaled by 1-4scaling factor 1.0
                                    # vdW special 1-4 parameters in charmm parameter file.
1-4scaling          1.0

# CUT-OFFS
switching                on
switchdist              10.0
cutoff                  12.0
pairlistdist            13.5

PME                     yes
PMEGridSpacing          1.0
PMETolerance            10e-6
PMEInterpOrder          4

wrapWater               on;                # wrap water to central cell
wrapAll                 on;                # wrap other molecules too
wrapNearest             off;               # use for non-rectangular cells (wrap to the nearest image)

# SPACE PARTITIONING
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  2
longSplitting           C2

# RESPA PROPAGATOR
# timestep                1.0
timestep                2.0
useSettle               on
fullElectFrequency      2
nonbondedFreq           1

# SHAKE
rigidbonds              all
rigidtolerance          0.000001
rigiditerations         400

# COM
ComMotion               no

# vdw
vdwForceSwitching       on

# Constant Temperature Control
if {{ $ITEMP == $FTEMP }} {{
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
}} else {{
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if {{ $FTEMP > $ITEMP }} {{
        reassignIncr 10
    }} else {{
        reassignIncr -10
    }}
    reassignTemp $ITEMP
    reassignHold $FTEMP
}}

# Constant Pressure Control (variable volume)
if {{ $PSWITCH != 0 }} {{
    # if running G-actin remove/comment out these 3 lines
    # by Junfan
    # CONSTANT-P, not in tutorial
    useGroupPressure        yes;           # use a hydrogen-group based pseudo-molecular viral to calcualte pressure and
                                        # has less fluctuation, is needed for rigid bonds (rigidBonds/SHAKE)
    useFlexibleCell         no;            # yes for anisotropic system like membrane
    useConstantRatio        no;            # keeps the ratio of the unit cell in the x-y plane constant A=B
    #    useConstatntArea     yes;
    langevinPiston          on
    langevinPistonTarget    1.01325
    langevinPistonPeriod    100;         # 100? 2000?
    langevinPistonDecay     50;         # 50?
    langevinPistonTemp      $FTEMP
    StrainRate              0.0 0.0 0.0
}}

# Output
outputname $outputbase-em;

#@ equilibration work flow. have to put in the end!
# run one step to get into scripting mode
minimize                0

# turn off until later
langevinPiston          off

# min all atoms
minimize                10000
'''.format(pdb, psf, ff_path, pbcbox)
    )
    mk_em_config.close()

def generate_nvt_config(pdb, psf, ff_path, time, dcdfrequency, constraints):

    mk_nvt_config = open("nvt.config", "w")
    mk_nvt_config.write(
'''
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# NVT
# namd3 +p1 +devices 0 nvt.config 2 > nvt.config.log  

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          {1}
coordinates        {0}
set outputbase     com

firsttimestep      0

set                 ITEMP 310
set                 FTEMP 310
set                 INPUTNAME   0                      ;# use the former outputName, for restarting a simulation
set                 PSWITCH     0                      ;# whether to use langevinPiston pressure control
set                 CONSSCALE   1                      ;# default; initial value if you want to change
set                 CONSPDB     {5}
set parpath         {2}

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${{parpath}}/addition.prm
parameters          ${{parpath}}/par_all36m_prot.prm
parameters          ${{parpath}}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if {{ $INPUTNAME != 0 }} {{
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
}} else {{
    bincoordinates      ${{outputbase}}-em.coor
    binvelocities       ${{outputbase}}-em.vel
    extendedSystem      ${{outputbase}}-em.xsc
}}

## Force-Field Parameters
exclude             scaled1-4;         # non-bonded exclusion policy to use "none,1-2,1-3,1-4,or scaled1-4"
                                    # 1-2: all atoms pairs that are bonded are going to be ignored
                                    # 1-3: 3 consecutively bonded are excluded
                                    # scaled1-4: include all the 1-3, and modified 1-4 interactions
                                    # electrostatic scaled by 1-4scaling factor 1.0
                                    # vdW special 1-4 parameters in charmm parameter file.
1-4scaling              1.0

# CUT-OFFS
switching                on
switchdist              10.0
cutoff                  12.0
pairlistdist            13.5

PME                     yes
PMEGridSpacing          1.0
PMETolerance            10e-6
PMEInterpOrder          4

wrapWater               on;                # wrap water to central cell
wrapAll                 on;                # wrap other molecules too
wrapNearest             off;               # use for non-rectangular cells (wrap to the nearest image)

# SPACE PARTITIONING
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  2
longSplitting           C2

# RESPA PROPAGATOR
# timestep                1.0
timestep                2.0
useSettle               on
fullElectFrequency      2
nonbondedFreq           1

# SHAKE
rigidbonds              all
rigidtolerance          0.000001
rigiditerations         400

# vdw
vdwForceSwitching       on

# Constant Temperature Control
if {{ $ITEMP == $FTEMP }} {{
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
}} else {{
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if {{ $FTEMP > $ITEMP }} {{
        reassignIncr 10
    }} else {{
        reassignIncr -10
    }}
    reassignTemp $ITEMP
    reassignHold $FTEMP
}}

# according to P. Blood use "no" for first NPT run
# then use "yes" for all NPT runs afterward
COMmotion no

#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Output
outputname $outputbase-nvt;

# 500steps = every 1ps
restartfreq         50000
dcdfreq             {4}
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

CUDASOAintegrate        on

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if {{ $CONSPDB != 0 }} {{
    Constraints          yes
    ConsRef              $CONSPDB.pdb
    ConsKFile            $CONSPDB.pdb
    ConskCol             B
    constraintScaling    $CONSSCALE
}}

# NVT
langevinPiston          off
run                     {3}
'''.format(pdb, psf, ff_path, time, dcdfrequency, constraints)
    )
    mk_nvt_config.close()

def generate_nptstep_config(pdb, psf, ff_path, time, dcdfrequency, constraints):
    
    mk_npt_config = open("npt.config", "w")
    mk_npt_config.write(
'''
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# namd3 +p1 +devices 0 npt.config > npt.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          {1}
coordinates        {0}
set outputbase     com

firsttimestep      0

set                 ITEMP 310
set                 FTEMP 310
set                 INPUTNAME   0                      ;# use the former outputName, for restarting a simulation
set                 PSWITCH     0                      ;# whether to use langevinPiston pressure control
set                 CONSSCALE   1                      ;# default; initial value if you want to change
set                 CONSPDB     {5}
set parpath         {2}

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${{parpath}}/addition.prm
parameters          ${{parpath}}/par_all36m_prot.prm
parameters          ${{parpath}}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if {{ $INPUTNAME != 0 }} {{
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
}} else {{
    bincoordinates      ${{outputbase}}-nvt.coor
    binvelocities       ${{outputbase}}-nvt.vel
    extendedSystem      ${{outputbase}}-nvt.xsc
}}

## Force-Field Parameters
exclude             scaled1-4;         # non-bonded exclusion policy to use "none,1-2,1-3,1-4,or scaled1-4"
                                    # 1-2: all atoms pairs that are bonded are going to be ignored
                                    # 1-3: 3 consecutively bonded are excluded
                                    # scaled1-4: include all the 1-3, and modified 1-4 interactions
                                    # electrostatic scaled by 1-4scaling factor 1.0
                                    # vdW special 1-4 parameters in charmm parameter file.
1-4scaling              1.0

# CUT-OFFS
switching                on
switchdist              10.0
cutoff                  12.0
pairlistdist            13.5

PME                     yes
PMEGridSpacing          1.0
PMETolerance            10e-6
PMEInterpOrder          4

wrapWater               on;                # wrap water to central cell
wrapAll                 on;                # wrap other molecules too
wrapNearest             off;               # use for non-rectangular cells (wrap to the nearest image)

# SPACE PARTITIONING
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  2
longSplitting           C2

# RESPA PROPAGATOR
# timestep                1.0
timestep                2.0
useSettle               on
fullElectFrequency      2
nonbondedFreq           1

# SHAKE
rigidbonds              all
rigidtolerance          0.000001
rigiditerations         400

# vdw
vdwForceSwitching       on

# Constant Temperature Control
if {{ $ITEMP == $FTEMP }} {{
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
}} else {{
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if {{ $FTEMP > $ITEMP }} {{
        reassignIncr 10
    }} else {{
        reassignIncr -10
    }}
    reassignTemp $ITEMP
    reassignHold $FTEMP
}}

# according to P. Blood use "no" for first NPT run
# then use "yes" for all NPT runs afterward
COMmotion no

#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Output
outputname $outputbase-nptstep;

# 500steps = every 1ps
restartfreq         50000
dcdfreq             {4}
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

CUDASOAintegrate        on

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if {{ $CONSPDB != 0 }} {{
        Constraints          yes
        ConsRef              $CONSPDB.pdb
        ConsKFile            $CONSPDB.pdb
        ConskCol             B
        constraintScaling    $CONSSCALE
}}

set PSWITCH 1
# Constant Pressure Control (variable volume)
if {{ $PSWITCH != 0 }} {{
    # if running G-actin remove/comment out these 3 lines
    # by Junfan
    # CONSTANT-P, not in tutorial
    useGroupPressure        yes;           # use a hydrogen-group based pseudo-molecular viral to calcualte pressure and
                                        # has less fluctuation, is needed for rigid bonds (rigidBonds/SHAKE)
    useFlexibleCell         no;            # yes for anisotropic system like membrane
    useConstantRatio        no;            # keeps the ratio of the unit cell in the x-y plane constant A=B
    #    useConstatntArea     yes;
    langevinPiston          on
    langevinPistonTarget    1.01325
    langevinPistonPeriod    100;         # 100? 2000?
    langevinPistonDecay     100;         # 50?
    langevinPistonTemp      $FTEMP
    #StrainRate              0.0 0.0 0.0
}}

run                     {3}
'''.format(pdb, psf, ff_path, time, dcdfrequency, constraints)
    )
    mk_npt_config.close()

def generate_prodstep_config(pdb, psf, ff_path, time, dcdfrequency):

    mk_prod_config = open("prod.config", "w")
    mk_prod_config.write(
'''
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# namd3 +p1 +devices 0 prod.config > prod.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          {1}
coordinates        {0}     ;# or reports error
set outputbase     com               ;# consistent with equil
firsttimestep      0

#############################################################

set ITEMP 310
set FTEMP 310
set INPUTNAME       0                   ;# restart
set PSWITCH         1                   ;# whether to use langevinPiston pressure control
set FIXPDB          0
set CONSPDB         0
set CONSSCALE       0                   ;# default:1
set parpath         {2}

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${{parpath}}/addition.prm
parameters          ${{parpath}}/par_all36m_prot.prm
parameters          ${{parpath}}/toppar_water_ions_namd.str
mergeCrossterms yes

if {{ $INPUTNAME != 0 }} {{
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
}} else {{
    # from equil. use the former outputName
    bincoordinates          $outputbase-nptstep.coor
    binvelocities           $outputbase-nptstep.vel
    extendedSystem      $outputbase-nptstep.xsc
}}

## Force-Field Parameters
exclude             scaled1-4;         # non-bonded exclusion policy to use "none,1-2,1-3,1-4,or scaled1-4"
                                    # 1-2: all atoms pairs that are bonded are going to be ignored
                                    # 1-3: 3 consecutively bonded are excluded
                                    # scaled1-4: include all the 1-3, and modified 1-4 interactions
                                    # electrostatic scaled by 1-4scaling factor 1.0
                                    # vdW special 1-4 parameters in charmm parameter file.
1-4scaling          1.0

# CUT-OFFS
switching                on
switchdist              10.0
cutoff                  12.0
pairlistdist            13.5

PME                     yes
PMEGridSpacing          1.0
PMETolerance            10e-6
PMEInterpOrder          4

wrapWater               on;                # wrap water to central cell
wrapAll                 on;                # wrap other molecules too
wrapNearest             off;               # use for non-rectangular cells (wrap to the nearest image)

# SPACE PARTITIONING
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  2
longSplitting           C2

# RESPA PROPAGATOR
# timestep                1.0
timestep                2.0
useSettle               on
fullElectFrequency      2
nonbondedFreq           1

# SHAKE
rigidbonds              all
rigidtolerance          0.000001
rigiditerations         400

# vdw
vdwForceSwitching       on

# Constant Temperature Control
if {{ $ITEMP == $FTEMP }} {{
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
}} else {{
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if {{ $FTEMP > $ITEMP }} {{
        reassignIncr 10
    }} else {{
        reassignIncr -10
    }}
    reassignTemp $ITEMP
    reassignHold $FTEMP
}}

# Constant Pressure Control (variable volume)
if {{ $PSWITCH != 0 }} {{
    # if running G-actin remove/comment out these 3 lines
    # by Junfan
    # CONSTANT-P, not in tutorial
    useGroupPressure        yes;           # use a hydrogen-group based pseudo-molecular viral to calcualte pressure and
                                        # has less fluctuation, is needed for rigid bonds (rigidBonds/SHAKE)
    useFlexibleCell         no;            # yes for anisotropic system like membrane
    useConstantRatio        no;            # keeps the ratio of the unit cell in the x-y plane constant A=B
    #    useConstatntArea     yes;
    langevinPiston          on
    langevinPistonTarget    1.01325
    langevinPistonPeriod    100;         # 100? 2000?
    langevinPistonDecay     50;         # 50?
    langevinPistonTemp      $FTEMP
    StrainRate              0.0 0.0 0.0
}}

# according to P. Blood use "no" for first NPT run
# then use "yes" for all NPT runs afterward
COMmotion yes

# Fixed atoms
# port first, h2o 2nd, 1 means not move
if {{ $FIXPDB != 0 }} {{
    fixedAtoms      yes
    fixedAtomsForces yes
    fixedAtomsFile  $FIXPDB.pdb
    fixedAtomsCol   B                   ;# beta
}}

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if {{ $CONSPDB != 0 }} {{
    Constraints          yes
    ConsRef              $CONSPDB.pdb
    ConsKFile            $CONSPDB.pdb
    ConskCol             B
    constraintScaling    $CONSSCALE
}}

CUDASOAintegrate         on

# Output
outputName          $outputbase-prodstep

restartfreq         50000     ;# 500steps = every 1ps. name=default
dcdfreq             {4}
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

run                 {3}
'''.format(pdb, psf, ff_path, time, dcdfrequency)
    )
    mk_prod_config.close()

def run():
    import sys
    settings = config()
    settings.title = sys.argv[1]
    submit(settings.title, settings.queue)
    md_pbc_box(settings.vmd_path, "complex.pdb", "complex.psf")
    position_constraints("complex.pdb", settings.vmd_path, settings.constraints)
    os.remove("mk_pbcbox.tcl")
    os.remove("constraints.tcl")

    if os.path.exists(os.path.join(".", "common")):
        shutil.rmtree(os.path.join(".", "common"))
        os.makedirs(os.path.join(".", "common"))
    else:
        os.makedirs(os.path.join(".", "common"))
    shutil.move("PBCBOX.dat", os.path.join(".", "common", "PBCBOX.dat"))
    shutil.move("constraints.pdb", os.path.join(".", "common", "constraints.pdb"))
    shutil.move("merged.pdb", os.path.join(".", "common", "merged.pdb"))
    shutil.move("merged.psf", os.path.join(".", "common", "merged.psf"))
    shutil.move("solvated.log", os.path.join(".", "common", "solvated.log"))
    shutil.move("solvated.pdb", os.path.join(".", "common", "solvated.pdb"))
    shutil.move("solvated.psf", os.path.join(".", "common", "solvated.psf"))
    shutil.move("complex.pdb", os.path.join(".", "common", "complex.pdb"))
    shutil.move("complex.psf", os.path.join(".", "common", "complex.psf"))

    if os.path.exists(os.path.join(".", "prod")):
        shutil.rmtree(os.path.join(".", "prod"))
        os.makedirs(os.path.join(".", "prod"))
    else:
        os.makedirs(os.path.join(".", "prod"))
    os.chdir(os.path.join(".", "prod"))
    generate_em_config("../common/complex.pdb", "../common/complex.psf", settings.ff_path, os.path.join("..", "common", "PBCBOX.dat"))
    generate_nvt_config("../common/complex.pdb", "../common/complex.psf", settings.ff_path, "20000", "50000", "../common/constraints")
    generate_nptstep_config("../common/complex.pdb", "../common/complex.psf", settings.ff_path, settings.npt_step, "50000", "../common/constraints")
    generate_prodstep_config("../common/complex.pdb", "../common/complex.psf", settings.ff_path, settings.prod_step, "50000")
    os.chdir(os.path.join(".."))

def main():
    run()

if __name__=="__main__":
    main() 
