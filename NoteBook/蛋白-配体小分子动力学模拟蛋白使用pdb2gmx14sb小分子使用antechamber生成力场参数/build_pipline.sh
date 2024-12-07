antechamber -i LIG.mol2 -fi mol2 -o LIG-1.mol2 -fo mol2 -c bcc -nc 0 -pf y
parmchk2 -i LIG-1.mol2 -f mol2 -o LIG.frcmod

cat << EOL > tleap.in
source oldff/leaprc.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2
loadamberparams frcmod.ionsjc_tip3p
loadamberparams LIG.frcmod
LIG = loadmol2 LIG-1.mol2
saveamberparm LIG LIG.prmtop LIG.inpcrd
quit
EOL

tleap -f tleap.in

cat <<EOL > do1.py
import parmed as pmd

parm = pmd.load_file('LIG.prmtop', 'LIG.inpcrd')
parm.save('gromacs.top', format='gromacs')
parm.save('gromacs.gro')
EOL

python do1.py

sed -i 's/HSD/HIS/g' complex.pdb
sed -i 's/HSE/HIS/g' complex.pdb

echo 1|gmx pdb2gmx -f complex.pdb -o build.gro -water tip3p -ignh

cat <<EOL > do2.py
def merge_gro_files(gro_file1, gro_file2, output_file):
    with open(gro_file1, 'r') as file1, open(gro_file2, 'r') as file2:
        # 读取第一个文件的内容
        lines1 = file1.readlines()
        # 读取第二个文件的内容
        lines2 = file2.readlines()

        # 合并注释行和总原子数
        header = lines1[0]
        atom_count1 = int(lines1[1].strip())
        atom_count2 = int(lines2[1].strip())
        total_atoms = atom_count1 + atom_count2
        
        # 获取坐标数据并去掉最后一行的box尺寸信息
        coords1 = lines1[2:-1]
        coords2 = lines2[2:-1]
        
        # 获取box尺寸信息
        box_size = lines1[-1].strip()  # 假设两个文件的box尺寸相同

        # 写入合并后的文件
        with open(output_file, 'w') as output:
            output.write(header)  # 写入注释行
            output.write(f'{total_atoms}\n')  # 写入总原子数
            output.writelines(coords1)  # 写入第一个文件的坐标
            output.writelines(coords2)  # 写入第二个文件的坐标
            output.write(f'{box_size}\n')  # 写入box尺寸信息

# 使用示例
merge_gro_files('build.gro', 'gromacs.gro', 'build1.gro')
EOL

python do2.py

head -n -20 topol.top > Protein_chain_A.itp
sed -i '1,27d' Protein_chain_A.itp

head -n -7 gromacs.top > LIG.itp
sed -i '1,19d' LIG.itp

echo 0 |gmx genrestr -f gromacs.gro -o posre_LIG.itp -fc 1000 1000 1000

echo -e "\n; Include Position restraint file\n#ifdef POSRES\n#include \"posre_LIG.itp\"\n#endif" >> LIG.itp

rm topol.top
cat <<EOL >> topol.top
;
;       File 'topol.top' was generated
;       By user: yqyang (1032)
;       On host: tc6000
;       At date: Thu May 16 15:05:39 2024
;
;       This is a standalone topology file
;
;       Created by:
;                           :-) GROMACS - gmx pdb2gmx, 2023.2 (-:
;
;       Executable:   /public/software/apps/gromacs/2023.2/bin/gmx
;       Data prefix:  /public/software/apps/gromacs/2023.2
;       Working dir:  /public/home/yqyang/SerpinB9/SerpinB9-SMS/SMS-peps-MMGBSA/pep_1/dup1
;       Command line:
;         gmx pdb2gmx -f complex.pdb -o build.gro -water tip3p -ignh
;       Force field data was read from:
;       /public/home/yqyang/software/forcefield
;
;       Note:
;       This might be a non-standard force field location. When you use this topology, the
;       force field must either be present in the current directory, or the location
;       specified in the GMXLIB path variable or with the 'include' mdp file option.
;

; Include forcefield parameters
#include "/public/home/yqyang/software/forcefield/amber14sb_OL15.ff/forcefield.itp"

; Include chain topologies
#include "LIG.itp"
#include "Protein_chain_A.itp"

; Include water topology
#include "/public/home/yqyang/software/forcefield/amber14sb_OL15.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "/public/home/yqyang/software/forcefield/amber14sb_OL15.ff/ions.itp"

[ system ]
; Name
Protein in water

[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG     1
EOL

gmx editconf -f build1.gro -o newbox.gro -bt cubic -d 0.8
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f ~/file/gmx_file/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo 15|gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
echo -e "1|13\nname 24 SOLU\nname 23 SOLV\nq\n"|gmx make_ndx -f solv_ions.gro -o index.ndx

cp ~/file/gmx_file/Step* .
python Step2_generate_mdp.py
python Step3_generate_submit_sh.py
sbatch job.sh
