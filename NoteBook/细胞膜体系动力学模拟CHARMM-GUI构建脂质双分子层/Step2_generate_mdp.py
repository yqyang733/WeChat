import os
import shutil

class config:

    def __init__(self):

        self.time = 1000           # unit: ns
        self.frame_steps = 500000

def mdp_file():

    if os.path.exists(os.path.join(".", "mdp")):
        shutil.rmtree(os.path.join(".", "mdp"))
        os.makedirs(os.path.join(".", "mdp"))
    else:
        os.makedirs(os.path.join(".", "mdp"))

def step6_0_minimizationmdp():

    step6_0_mdp = open(os.path.join(".", "mdp", "step6.0_minimization.mdp"), "w")
    step6_0_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=4000.0 -DPOSRES_FC_SC=2000.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 
integrator              = steep
emtol                   = 1000.0
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS
nstxout                = 0          ; save coordinates to .trr every 250 steps
nstvout                = 0          ; don't save velocities to .trr
nstfout                = 0          ; don't save forces to .trr

nstxout-compressed     = 0      ; xtc compressed trajectory output every 500 steps
compressed-x-precision = 1000
nstlog                 = 500        ; update log file every 500 steps
nstenergy              = 500        ; save energies every 500 steps
nstcalcenergy          = 500                                           
'''
    )

def step6_1_equilibrationmdp():

    step6_1_mdp = open(os.path.join(".", "mdp", "step6.1_equilibration.mdp"), "w")
    step6_1_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=4000.0 -DPOSRES_FC_SC=2000.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 
integrator              = md
dt                      = 0.001
nsteps                  = 125000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 10000
nstlog                  = 10000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
;
gen-vel                 = yes
gen-temp                = 310.15
gen-seed                = -1
'''
    )

def step6_2_equilibrationmdp():

    step6_2_mdp = open(os.path.join(".", "mdp", "step6.2_equilibration.mdp"), "w")
    step6_2_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=2000.0 -DPOSRES_FC_SC=1000.0 -DPOSRES_FC_LIPID=400.0 -DDIHRES -DDIHRES_FC=400.0 
integrator              = md
dt                      = 0.001
nsteps                  = 125000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 10000
nstlog                  = 10000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
'''
    )

def step6_3_equilibrationmdp():

    step6_3_mdp = open(os.path.join(".", "mdp", "step6.3_equilibration.mdp"), "w")
    step6_3_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=500.0 -DPOSRES_FC_LIPID=400.0 -DDIHRES -DDIHRES_FC=200.0 
integrator              = md
dt                      = 0.001
nsteps                  = 125000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 10000
nstlog                  = 10000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
'''
    )

def step6_4_equilibrationmdp():

    step6_4_mdp = open(os.path.join(".", "mdp", "step6.4_equilibration.mdp"), "w")
    step6_4_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=500.0 -DPOSRES_FC_SC=200.0 -DPOSRES_FC_LIPID=200.0 -DDIHRES -DDIHRES_FC=200.0 
integrator              = md
dt                      = 0.002
nsteps                  = 250000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 10000
nstlog                  = 10000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU_MEMB SOLV
'''
    )

def step6_5_equilibrationmdp():

    step6_5_mdp = open(os.path.join(".", "mdp", "step6.5_equilibration.mdp"), "w")
    step6_5_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=200.0 -DPOSRES_FC_SC=50.0 -DPOSRES_FC_LIPID=40.0 -DDIHRES -DDIHRES_FC=100.0 
integrator              = md
dt                      = 0.002
nsteps                  = 250000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 10000
nstlog                  = 10000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
'''
    )

def step6_6_equilibrationmdp():

    step6_6_mdp = open(os.path.join(".", "mdp", "step6.6_equilibration.mdp"), "w")
    step6_6_mdp.write(
'''define                  = -DPOSRES -DPOSRES_FC_BB=50.0 -DPOSRES_FC_SC=0.0 -DPOSRES_FC_LIPID=0.0 -DDIHRES -DDIHRES_FC=0.0 
integrator              = md
dt                      = 0.002
nsteps                  = 250000
;nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed     = 0      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = berendsen
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
'''
    )

def step7_productionmdp(time, frame_steps):

    nsteps = int(time*1000000/2)

    prod_mdp = open(os.path.join(".", "mdp", "step7_production.mdp"), "w")
    prod_mdp.write(
'''integrator              = md
dt                      = 0.002
nsteps                  = {0}
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 1000
nstenergy               = {1}
nstlog                  = {1}
nstxout-compressed     = {1}      ; xtc compressed trajectory output every 5000 steps
compressed-x-precision = 1000       ; precision with which to write to the compressed trajectory file
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = V-rescale
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
'''.format(nsteps, frame_steps)
    )

def run():

    settings = config()
    mdp_file()
    step6_0_minimizationmdp()
    step6_1_equilibrationmdp()
    step6_2_equilibrationmdp()
    step6_3_equilibrationmdp()
    step6_4_equilibrationmdp()
    step6_5_equilibrationmdp()
    step6_6_equilibrationmdp()
    step7_productionmdp(settings.time, settings.frame_steps)                                                                           

def main():

    run()

if __name__=="__main__":
    main()