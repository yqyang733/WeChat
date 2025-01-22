# NAMD进行共价体系蛋白配体动力学模拟，共价小分子配体参数由SwissParam生成
## Temp
```shell
$ ~/../../software/apps/vmd/1.9.3/vmd -dispdev text   # 启动vmd交互界面执行命令

vmd > package require psfgen   # 加载psfgen包建模
vmd > psfcontext reset

vmd > topology /public/home/yqyang/file/vegf-toppar/top_all36_prot.rtf   # 加载力场文件
vmd > topology /public/home/yqyang/file/vegf-toppar/toppar_water_ions.str

vmd > pdbalias residue HIS HSE
vmd > alias atom ILE CD1 CD
vmd > alias atom SER HG HG1
vmd > alias atom CYS HG HG1

vmd > segment PRO {pdb Protein.pdb}   # 加载蛋白文件
vmd > coordpdb Protein.pdb PRO

vmd > topology post.rtf   # 加载共价配体文件
vmd > segment LIG {pdb post.pdb}
vmd > coordpdb post.pdb LIG

vmd > topology reaction.str   # 构建共价键
vmd > patch REAC PRO:72 LIG:1

vmd > regenerate angles dihedrals   # 重新生成键角，二面角
vmd > guesscoord

vmd > writepsf merged.pdb
vmd > writepdb merged.psf

vmd > quit
```

确定PBC盒子尺寸。  
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

echo ${box_size}
```
加水加离子。  
```shell
$ ~/../../software/apps/vmd/1.9.3/vmd -dispdev text

vmd > package require pbctools   # 加载pbctools软件包

vmd > mol load psf merged.psf pdb merged.pdb   # 加载整个分子

vmd > package require solvate   # 加载solvate软件包进行溶剂化
vmd > solvate merged.psf merged.pdb -minmax {{-42 -35 -23} {26 32 44}} -o solvated
vmd > mol delete all

vmd > package require autoionize   # 加载autoionize软件包进行离子化
vmd > autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o complex
vmd > pbc box -center centerofmass

vmd > exit
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
set minmax [measure minmax \$everyone]
foreach {min max} \$minmax { break }
foreach {xmin ymin zmin} \$min { break }
foreach {xmax ymax zmax} \$max { break }
set file [open "PBCBOX.dat" w]
puts \$file "cellBasisVector1 [ expr \$xmax - \$xmin ] 0 0 "
puts \$file "cellBasisVector2 0 [ expr \$ymax - \$ymin ] 0 "
puts \$file "cellBasisVector3 0 0 [ expr \$zmax - \$zmin ] "
puts \$file "cellOrigin [ expr (\$xmax + \$xmin)/2 ] [ expr (\$ymax + \$ymin)/2 ] [ expr (\$zmax + \$zmin)/2 ] "
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
\$all set beta 0
set sel [atomselect top "protein and noh"]
\$sel set beta 1
\$all writepdb constraints.pdb
quit
EOF
vmd/1.9.3/vmd -dispdev text -e constraints.tcl
rm constraints.tcl
```

（7）准备键，键角，二面角等参数文件：  
与标准残基或单独配体小分子模拟不同，配体与标准残基形成共价键之后，以共价键（对应的原子类型）为中心对应原子类型的键，键角，二面角参数是没有的，需要用户以 parameters .par 的形式自行添加进跑模拟的config文件中。所以需要找出哪些键类型，键角类型和二面角类型参数没有。  

如下图所示，首先找出涉及共价键的键类型，键角类型和二面角类型都有哪些。如下图红圈框出的所示，涉及共价键的所有原子包括共价键向两侧延申两个原子后所包含的所有原子。先确定共价键，然后共价键向外两侧延申一个原子即可得到键角，然后键角再向外两侧延申一个原子即可得到二面角。如下图键，键角，二面角的类型在原始的.par参数文件中是没有的，需要将这些参数补进去。    
![](NAMD进行共价体系蛋白配体动力学模拟共价小分子配体参数由SwissParam生成/NAMD进行共价体系蛋白配体动力学模拟共价小分子配体参数由SwissParam生成_2025-01-22-17-01-34.png)  

但是共价反应产生的力场中仅延伸到CYS残基的CB原子，并没有延申到CA原子。所以根据上述已有的文件是没办法获得含有CA原子的二面角参数的。所以这里将分子结构mol2文件多保存一些然后使用SwissParam基于MMFF在产生一次力场参数。输入结构如下图所示：  
