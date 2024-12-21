# 纯标准蛋白质体系分子动力学模拟：vmd建模+NAMD模拟
在生物学和化学研究中，分子动力学（Molecular Dynamics, MD）模拟已成为一种不可或缺的工具，用于研究蛋白质的结构与动态行为。通过 MD 模拟，研究者可以以原子分辨率探索蛋白质的构象变化、配体结合机制以及蛋白质功能的内在动力学基础。尤其是针对纯标准蛋白质体系的模拟，这种方法为解析蛋白质在天然状态下的动态特性提供了关键帮助。  

VMD（Visual Molecular Dynamics） 和 NAMD（Nanoscale Molecular Dynamics） 是 MD 模拟中的两款主流工具，分别承担建模和模拟计算的任务。VMD 提供强大的分子可视化功能，可以轻松完成体系构建与初步检查；而 NAMD 以其高性能和良好的可扩展性，在大规模分子体系的模拟中表现出色。通过这两种工具的结合，研究者能够高效完成从模型构建到动力学计算的整个流程。本文主要简介vmd建模+NAMD模拟进行纯标准蛋白质体系分子动力学模拟的流程。  
## 建模前准备
获取体系的pdb文件，并将体系中的所有蛋白链分别保存为单独的pdb文件。如此示例体系由一个单链蛋白和一条单链多肽组成。因此将此体系pdb文件分别保存为蛋白REC.pdb和多肽PEP.pdb。  
## VMD对体系建模生成pdb和psf文件
（1）首先需要安装VMD软件。  
（2）使用VMD的psfgen模块建模蛋白。命令如下所示：  
```shell
mkdir common
cd common
cp ../REC.pdb .
cp ../PEP.pdb .
cat > pipline.tcl << EOF
package require psfgen   # 加载psfgen模块
psfcontext reset

topology top_all36_prot.rtf   # 加载力场
topology toppar_water_ions.str

pdbalias residue HIS HSE   # 无法识别HIS，可以使用该方式，也可以先确定质子化状态将组氨酸名称改为HSD，HSE和HSP再进行建模
alias atom ILE CD1 CD   # 对一些原子名字进行重命名，否则无法识别
alias atom SER HG HG1
alias atom CYS HG HG1

segment REC {pdb REC.pdb}   # 创建并加载蛋白
coordpdb REC.pdb REC

segment PEP {pdb PEP.pdb}   # 创建并加载多肽
coordpdb PEP.pdb PEP

guesscoord   # 对一些缺失原子，比如氢原子坐标补全

writepdb merged.pdb   # 保存为pdb文件
writepsf merged.psf   # 保存为psf文件

exit
EOF
vmd/1.9.3/vmd -dispdev text -e pipline.tcl
rm pipline.tcl
```
（3）确定PBC盒子尺寸。    
在后续加盒子加水时候可以使用如下-t选项添加长方体盒子也可以使用-minmax指定特定的盒子尺寸，如下所示。  
```shell
solvate merged.psf merged.pdb -t 10 -o solvated   # 指定盒子边界距离蛋白所有原子的最近距离是10埃，一般是创建的长方体盒子。
solvate merged.psf merged.pdb -minmax {{0 -5 2} {115 110 118}} -o solvated  # 通过指定xyz对应的最小值和最大值确定盒子
```
使用第一种方式创建长方体盒子的优点在于创建的盒子较小，使得体系的原子数目不多，模拟效率高。缺陷在于如果长时间模拟让溶质方向转动比较大可能导致其首尾撞在一起计算出错，并且也可能导致PBC处理的问题。所以我个人一般使用第二种方式创建正方体盒子。可使用下述脚本计算体系正方体盒子对应的-minmax参数，读取pdb的所有原子坐标计算其几何中心，确定pdb文件xyz方向的最大边长作为盒子长宽高，通过盒子几何中心和盒子长宽高确定盒子xyz方向对应的最小值和最大值。    
```shell
cat > do.py << EOF
def set_watbox(file_in):
    with open(file_in) as f:
        f1 = f.readlines()
    x = []
    y = []
    z = []
    for i in f1:
        if i.startswith("ATOM"):
            x.append(float(i[30:38]))
            y.append(float(i[38:46]))
            z.append(float(i[46:54]))
    center = ((max(x)+min(x))/2, (max(y)+min(y))/2, (max(z)+min(z))/2)
    x_com = max(x) - min(x)
    y_com = max(y) - min(y)
    z_com = max(z) - min(z)
    radius = max(x_com, y_com, z_com)/2 + 12
    x_min = int(center[0] - radius)
    x_max = int(center[0] + radius)
    y_min = int(center[1] - radius)
    y_max = int(center[1] + radius)
    z_min = int(center[2] - radius)
    z_max = int(center[2] + radius)
    print("{{{{{0} {1} {2}}} {{{3} {4} {5}}}}}".format(x_min, y_min, z_min, x_max, y_max, z_max))
    return x_min, x_max, y_min, y_max, z_min, z_max

def main():
    import sys

    set_watbox(sys.argv[1])

if __name__=="__main__":
    main()
EOF
box_size=`python do.py merged.pdb`
rm do.py
```
（4）加水加离子。  
使用VMD的pbctools模块，solvate模块和autoionize模块对体系加水加离子生成建模完成的complex.pdb文件和complex.psf文件。   
```shell
box_size=`python do.py merged.pdb`
cat > pipline.tcl << EOF
package require pbctools
psfcontext reset
mol load psf merged.psf pdb merged.pdb
package require solvate
solvate merged.psf merged.pdb -minmax ${box_size} -o solvated
mol delete all
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o complex
pbc box -center centerofmass
exit
EOF
vmd/1.9.3/vmd -dispdev text -e pipline.tcl
rm pipline.tcl
```
（5）生成PBC盒子信息文件。  
读取加水加离子之后的pdb文件输出PBC盒子信息，命令如下所示：  
```shell
cat > mk_pbcbox.tcl << EOF
#!/bin/bash
# vmd -dispdev text -e mk_pbcbox.tcl
package require psfgen
psfcontext reset
mol load psf complex.psf pdb complex.pdb
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
EOF
vmd/1.9.3/vmd -dispdev text -e mk_pbcbox.tcl
rm mk_pbcbox.tcl
```
（6）生成位置限制文件。如下所示：  
```shell
cat > constraints.tcl << EOF
mol new complex.pdb type pdb waitfor all
set all [atomselect top "all"]
$all set beta 0
set sel [atomselect top "protein and noh"]
$sel set beta 1
$all writepdb constraints.pdb
quit
EOF
vmd/1.9.3/vmd -dispdev text -e constraints.tcl
rm constraints.tcl
```
## 准备NAMD模拟的配置文件
（1）准备em，nvt，npt，prod所需的配置文件em.config，nvt.config，npt.config和prod.config。将它们存在prod文件夹中。  
```shell
mkdir prod
cd prod
```
（2）em.config文件如下：  
```text
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# Minimization
# namd3 +p1 em.config > em.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          ../common/complex.psf
coordinates        ../common/complex.pdb
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
set parpath     ./toppar

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${parpath}/par_all36m_prot.prm
parameters          ${parpath}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if { $INPUTNAME != 0 } {
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
} else {
    # Periodic Boundary Conditionsc
    temperature $ITEMP
    source ../common/PBCBOX.dat
}

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
if { $ITEMP == $FTEMP } {
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
} else {
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if { $FTEMP > $ITEMP } {
        reassignIncr 10
    } else {
        reassignIncr -10
    }
    reassignTemp $ITEMP
    reassignHold $FTEMP
}

# Constant Pressure Control (variable volume)
if { $PSWITCH != 0 } {
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
}

# Output
outputname $outputbase-em;

#@ equilibration work flow. have to put in the end!
# run one step to get into scripting mode
minimize                0

# turn off until later
langevinPiston          off

# min all atoms
minimize                10000
```
（3）nvt.config文件如下：  
```text
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# NVT
# namd3 +p1 +devices 0 nvt.config 2 > nvt.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          ../common/complex.psf
coordinates        ../common/complex.pdb
set outputbase     com

firsttimestep      0

set                 ITEMP 310
set                 FTEMP 310
set                 INPUTNAME   0                      ;# use the former outputName, for restarting a simulation
set                 PSWITCH     0                      ;# whether to use langevinPiston pressure control
set                 CONSSCALE   1                      ;# default; initial value if you want to change
set                 CONSPDB     ../common/constraints
set parpath         ./toppar

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${parpath}/par_all36m_prot.prm
parameters          ${parpath}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if { $INPUTNAME != 0 } {
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
} else {
    bincoordinates      ${outputbase}-em.coor
    binvelocities       ${outputbase}-em.vel
    extendedSystem      ${outputbase}-em.xsc
}

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
if { $ITEMP == $FTEMP } {
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
} else {
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if { $FTEMP > $ITEMP } {
        reassignIncr 10
    } else {
        reassignIncr -10
    }
    reassignTemp $ITEMP
    reassignHold $FTEMP
}

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
dcdfreq             50000
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

CUDASOAintegrate        on

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if { $CONSPDB != 0 } {
    Constraints          yes
    ConsRef              $CONSPDB.pdb
    ConsKFile            $CONSPDB.pdb
    ConskCol             B
    constraintScaling    $CONSSCALE
}

# NVT
langevinPiston          off
run                     20000
```
（4）npt.config文件如下：  
```text
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# namd3 +p1 +devices 0 npt.config > npt.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          ../common/complex.psf
coordinates        ../common/complex.pdb
set outputbase     com

firsttimestep      0

set                 ITEMP 310
set                 FTEMP 310
set                 INPUTNAME   0                      ;# use the former outputName, for restarting a simulation
set                 PSWITCH     0                      ;# whether to use langevinPiston pressure control
set                 CONSSCALE   1                      ;# default; initial value if you want to change
set                 CONSPDB     ../common/constraints
set parpath         ./toppar

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${parpath}/par_all36m_prot.prm
parameters          ${parpath}/toppar_water_ions_namd.str
mergeCrossterms yes

# restart or PBC
if { $INPUTNAME != 0 } {
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
} else {
    bincoordinates      ${outputbase}-nvt.coor
    binvelocities       ${outputbase}-nvt.vel
    extendedSystem      ${outputbase}-nvt.xsc
}

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
if { $ITEMP == $FTEMP } {
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
} else {
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if { $FTEMP > $ITEMP } {
        reassignIncr 10
    } else {
        reassignIncr -10
    }
    reassignTemp $ITEMP
    reassignHold $FTEMP
}

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
dcdfreq             50000
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

CUDASOAintegrate        on

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if { $CONSPDB != 0 } {
        Constraints          yes
        ConsRef              $CONSPDB.pdb
        ConsKFile            $CONSPDB.pdb
        ConskCol             B
        constraintScaling    $CONSSCALE
}

set PSWITCH 1
# Constant Pressure Control (variable volume)
if { $PSWITCH != 0 } {
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
}

run                     100000
```
（5）prod.config文件如下：  
```text
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# namd3 +p1 +devices 0 prod.config > prod.config.log

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          ../common/complex.psf
coordinates        ../common/complex.pdb     ;# or reports error
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
set parpath         ./toppar

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ${parpath}/par_all36m_prot.prm
parameters          ${parpath}/toppar_water_ions_namd.str
mergeCrossterms yes

if { $INPUTNAME != 0 } {
    # restart
    BinVelocities $INPUTNAME.restart.vel.old
    BinCoordinates $INPUTNAME.restart.coor.old
    ExtendedSystem $INPUTNAME.restart.xsc.old
} else {
    # from equil. use the former outputName
    bincoordinates          $outputbase-nptstep.coor
    binvelocities           $outputbase-nptstep.vel
    extendedSystem      $outputbase-nptstep.xsc
}

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
if { $ITEMP == $FTEMP } {
    langevin                   on;         # do langevin dynamics
    langevinDamping             1;         # damping coefficient (gamma) of 1/ps
                                        # 5/ps by Junfan
    langevinTemp           $FTEMP;
    langevinHydrogen          off;         # don't couple langevin bath to hydrogens
} else {
    reassignFreq 1000;                     # use this to reassign velocity every 1000 steps
    if { $FTEMP > $ITEMP } {
        reassignIncr 10
    } else {
        reassignIncr -10
    }
    reassignTemp $ITEMP
    reassignHold $FTEMP
}

# Constant Pressure Control (variable volume)
if { $PSWITCH != 0 } {
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
}

# according to P. Blood use "no" for first NPT run
# then use "yes" for all NPT runs afterward
COMmotion yes

# Fixed atoms
# port first, h2o 2nd, 1 means not move
if { $FIXPDB != 0 } {
    fixedAtoms      yes
    fixedAtomsForces yes
    fixedAtomsFile  $FIXPDB.pdb
    fixedAtomsCol   B                   ;# beta
}

# Positional restraints
# Write out a separate pdb file in which the B values for
# the backbone, the non-hydrogen nucleotide atoms, the ion,
# and the water oxygens within 2.5 A of magnesium are set to 2
if { $CONSPDB != 0 } {
    Constraints          yes
    ConsRef              $CONSPDB.pdb
    ConsKFile            $CONSPDB.pdb
    ConskCol             B
    constraintScaling    $CONSSCALE
}

CUDASOAintegrate         on

# Output
outputName          $outputbase-prodstep

restartfreq         50000     ;# 500steps = every 1ps. name=default
dcdfreq             50000
xstFreq             50000
outputEnergies      50000
outputPressure      50000
outputTiming        50000

run                 500000   # only 1ns
```
## em，nvt，npt，prod
使用NAMD进行em，nvt，npt，prod，命令如下所示：  
```shell
source namd_3.0alpha9.sh

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
```
## 一键式流程化脚本
同样，该流程化脚本是我在做项目时候针对特定的项目写定的成功运行的脚本。但并不一定适用于所有体系的一键化建模。该脚本只是提供一个样例作为记录，针对具体的项目，需要具体问题具体分析。可在该样例脚本中进行修改使用即可。   
```shell
build_md(){
    cat > pipline.tcl << EOF
package require psfgen
psfcontext reset
topology toppar/top_all36_prot.rtf
topology toppar/toppar_water_ions.str
pdbalias residue HIS HSE
alias atom ILE CD1 CD
alias atom SER HG HG1
alias atom CYS HG HG1
segment REC {pdb REC.pdb}
coordpdb REC.pdb REC
segment PEP {pdb PEP.pdb}
coordpdb PEP.pdb PEP
guesscoord
writepdb merged.pdb
writepsf merged.psf
exit
EOF
    vmd/1.9.3/vmd -dispdev text -e pipline.tcl
    rm pipline.tcl
    rm do.py
    cat > do.py << EOF
def set_watbox(file_in):
    with open(file_in) as f:
        f1 = f.readlines()
    x = []
    y = []
    z = []
    for i in f1:
        if i.startswith("ATOM"):
            x.append(float(i[30:38]))
            y.append(float(i[38:46]))
            z.append(float(i[46:54]))
    center = ((max(x)+min(x))/2, (max(y)+min(y))/2, (max(z)+min(z))/2)
    x_com = max(x) - min(x)
    y_com = max(y) - min(y)
    z_com = max(z) - min(z)
    radius = max(x_com, y_com, z_com)/2 + 12
    x_min = int(center[0] - radius)
    x_max = int(center[0] + radius)
    y_min = int(center[1] - radius)
    y_max = int(center[1] + radius)
    z_min = int(center[2] - radius)
    z_max = int(center[2] + radius)
    print("{{{{{0} {1} {2}}} {{{3} {4} {5}}}}}".format(x_min, y_min, z_min, x_max, y_max, z_max))
    return x_min, x_max, y_min, y_max, z_min, z_max

def main():
    import sys

    set_watbox(sys.argv[1])

if __name__=="__main__":
    main()
EOF
    box_size=`python do.py merged.pdb`
    rm do.py
    cat > pipline.tcl << EOF
package require pbctools
psfcontext reset
mol load psf merged.psf pdb merged.pdb
package require solvate
solvate merged.psf merged.pdb -minmax ${box_size} -o solvated
mol delete all
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o complex
pbc box -center centerofmass
exit
EOF
    vmd/1.9.3/vmd -dispdev text -e pipline.tcl
    rm pipline.tcl
    python mk_md_namd.py
}
build_md
```
mk_md_namd.py的内容如下所示：   
```python
import os
import time
import shutil

class config:

    def __init__(self):

        self.title = ""
        self.queue = ""
        self.system_pdb = os.path.join("complex.pdb")
        self.system_psf = os.path.join("complex.psf")
        self.vmd_path = os.path.join("/vmd/1.9.3/vmd")
        self.ff_path = os.path.join("./toppar")
        self.npt_step = "100000"                       # timestep: 2fs
        self.prod_step = "500000"                      # timestep: 2fs  
        self.constraints = "protein and noh"

def submit():
    
    mk_submit = open("do_md.sh", "w")
    mk_submit.write(
'''#!/bin/bash
source namd_3.0alpha9.sh

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
''')
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
    submit()
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
```
## 参考
1. [workflow.sh](./纯标准蛋白质体系分子动力学模拟vmd建模NAMD模拟/workflow.sh)  
2. [mk_md_namd.py](./纯标准蛋白质体系分子动力学模拟vmd建模NAMD模拟/mk_md_namd.py)  