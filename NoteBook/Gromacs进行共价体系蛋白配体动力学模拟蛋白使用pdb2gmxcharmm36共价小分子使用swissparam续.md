# Gromacs进行共价体系蛋白配体动力学模拟：蛋白使用pdb2gmx charmm36，共价小分子使用swissparam（续）
**写在前面：** 这篇是文章[]()的后续。文章[]()中设置共价的方式是在[ intermolecular_interactions ]字段使用分子间限制（距离限制，角度限制，二面角限制）模拟共价作用。[ intermolecular_interactions ]字段中不能使用成化学键的bond type (type 1)，所以选用了距离限制（type 6）代替。所以真实模拟时候确实能模拟近似共价的效果。本篇文章则不使用[ intermolecular_interactions ]字段，而是将蛋白配体itp合并，将共价相关信息直接加入[ bonds ]，[ angles ]，[ dihedrals ]中进行更加准确的模拟。  

本文在文章[]()的基础上仍以PDBid：5VBM为例将重点介绍如何在GROMACS平台上进行共价体系的蛋白-配体动力学模拟，利用CHARMM36为蛋白质提供强大的力场支持，同时使用SwissParam为共价小分子配体生成定制的力场参数。另外，本文重点介绍共价相关的文件准备，其他步骤可参考文章[]()。  

## 使用SwissParam产生共价配体92V的MMFF力场参数itp文件
参考文章[]()中的 使用SwissParam产生共价配体92V的MMFF力场参数itp文件。  
## pdb2gmx产生共价反应后蛋白的拓扑与力场
参考文章[]()中的 pdb2gmx产生共价反应后蛋白的拓扑与力场。  
## 修改拓扑文件构建共价结构和力场参数（键，键角，二面角等参数）
本文不同之处主要在此部分。  