# 蛋白-配体小分子动力学模拟：蛋白使用pdb2gmx charmm36，小分子使用cgenff生成力场参数
在现代生物分子模拟中，蛋白-配体相互作用的研究为药物发现与设计提供了重要的理论依据。分子动力学（MD）模拟作为一种强有力的工具，通过精确地模拟分子的运动和相互作用，能够揭示蛋白质和小分子之间的结合机制，帮助我们深入理解其生物学功能。然而，进行蛋白-配体小分子动力学模拟时，准确的力场参数设置至关重要，特别是在蛋白质和小分子的模拟中。本文我们将介绍如何使用GROMACS工具进行蛋白-配体小分子动力学模拟。对于蛋白质部分，采用pdb2gmx与CHARMM36力场进行模拟，以确保模拟的高准确性和可靠性；对于小分子部分，则使用cgenff生成力场参数，以适应小分子的多样性和复杂性。
![](蛋白-配体小分子动力学模拟蛋白使用pdb2gmxcharmm36小分子使用cgenff生成力场参数/蛋白-配体小分子动力学模拟蛋白使用pdb2gmxcharmm36小分子使用cgenff生成力场参数_2024-12-10-23-51-27.png)  
## 使用charmmgui的Ligand Reader & Modeler模块准备配体小分子参数文件
（1）准备配体小分子mol2文件，这里仅说一下我的个人习惯，我一般将小分子放在pymol中进行编辑（包括加氢，保证单双键的正确性，保证结构的正确性），然后将其保存成mol2格式文件。  
（2）将（1）中的配体mol2文件上传charmmgui的Ligand Reader & Modeler模块生成参数（下载整个压缩包）。
  
## 使用charmmgui的Force Field Converter模块将力场转为gmx使用的力场文件
（3）将（2）中的ligandrm.psf，ligandrm.crd，lig.rtf和lig.prm文件上传
## 生成蛋白的gro文件和top文件
## 合并蛋白配体的gro文件和top文件并加水加离子
## em，nvt，npt，md
## 一键式流程化脚本
同样，该流程化脚本是我在做项目时候针对特定的项目写定的成功运行的脚本。但并不一定适用于所有体系的一键化建模。该脚本只是提供一个样例作为记录，针对具体的项目，需要具体问题具体分析。可在该样例脚本中进行修改使用即可。   
```shell
echo 2 | gmx pdb2gmx -f complex.pdb -o build.pdb -water tip3p -ignh   # 选择charmm36力场
head -n -2 build.pdb > build1.pdb
cat step3_input.pdb >> build1.pdb
cat << EOL > do.py
import os

class topol_file():

    def __init__(self, top):

        self.title = list()
        self.forcefield = list()
        self.watertop = list()
        self.ionstop = list()
        self.system = list()
        self.molecules = list()
        
        identifier = {"; Include forcefield parameters":1, "; Include water topology":2, "; Include topology for ions":3, "[ system ]":4, "[ molecules ]":5,}
        index = 0
        flag = 0
        with open(top) as f:
            lines = f.readlines()

        while index < len(lines):
            line = lines[index]
            
            if line.strip() in identifier.keys():
                flag = identifier[line.strip()]
            if flag == 0:
                self.title.append(line)
            elif flag == 1:
                self.forcefield.append(line)
            elif flag == 2:
                self.watertop.append(line)                  
            elif flag == 3:
                self.ionstop.append(line) 
            elif flag == 4:
                self.system.append(line)
            elif flag == 5:
                self.molecules.append(line)

            index += 1

    def lig_add(self):

        rt = open("topol.top1", "w")
        rt.write("".join(self.title))
        rt.write("".join(self.forcefield) + "#include \"LIG.itp\"\n\n")
        rt.write("".join(self.watertop))
        rt.write("".join(self.ionstop))
        rt.write("".join(self.system))
        rt.write("".join(self.molecules)+"LIG       1\n")
        os.remove("topol.top")
        os.renames("topol.top1", "topol.top")

topol_file("topol.top").lig_add()
EOL
python do.py
gmx editconf -f build1.pdb -o newbox.gro -bt cubic -d 0.8
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f ~/file/gmx_file/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo 15 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc 0.15
echo -e "1|13\nname 19 SOLU\n14|15|16\nname 20 SOLV\nq\n"|gmx make_ndx -f solv_ions.gro -o index.ndx   # 设置SOLU和SOLV

python Step2_generate_mdp.py
python Step3_generate_submit_sh.py
sh job.sh
```
## 参考
1. [build_pipline.sh](./蛋白-配体小分子动力学模拟蛋白使用pdb2gmxcharmm36小分子使用cgenff生成力场参数/build_pipline.sh)  
2. [Step2_generate_mdp.py](./Gromacs进行纯标准蛋白质体系分子动力学模拟/Step2_generate_mdp.py)  
3. [Step3_generate_submit_sh.py](./Gromacs进行纯标准蛋白质体系分子动力学模拟/Step3_generate_submit_sh.py)  