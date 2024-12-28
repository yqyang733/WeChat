# 蛋白-配体小分子动力学模拟：蛋白使用pdb2gmx charmm36，小分子使用swissparm生成的MATCH力场参数
## 使用NAMD进行模拟
## 使用Gromacs进行模拟
## Temp
直接网站上下载的只有prm，pdb，rtf，par等文件，没有psf文件，不足以直接使用charmm-gui中的Force Field Converter模块将其准换成gmx兼容的itp文件格式。解决方式是使用vmd的psfgen模块根据这些文件生成psf文件。然后使用Force Field Converter进行转换。  

## 