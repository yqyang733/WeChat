# gmx后处理|单链体系周期性边界条件处理生成平滑轨迹
写在前面，在gmx跑轨迹时候个人觉得比较麻烦一步就是原生轨迹的PBC处理了。由于模拟时候的PBC条件设置，所以原生轨迹往往会出现周期跳跃的情况，原生轨迹断断续续，跳来跳去，非常不利于观察。特别是涉及多链情况时候周期性处理不当就会得到破碎的轨迹。所以这里根据之前做过的一些项目为各种情况下的PBC处理尽量找一些通用的流程并剖析处理过程。本文仅介绍单链的PBC处理，后续也会更新介绍其他多种情况下的PBC处理流程。  

在分子动力学模拟中，周期性边界条件（PBC）是用来模拟无限大系统的常用技术，它通过在模拟盒中不断复制粒子，从而避免了边界效应对系统行为的干扰。然而，在处理单链体系时，尤其是链状分子（如蛋白质、RNA、聚合物等）时，周期性边界往往会导致链末端的接续问题，进而影响轨迹的平滑性和后续分析的准确性。尤其是当链条在模拟过程中穿越了边界，产生了不连续的轨迹数据，这些“断裂”会影响进一步的分析，比如计算二级结构、构象变化或相互作用等。  

本文将重点介绍如何在GROMACS中处理单链体系的周期性边界条件，并通过后处理方法生成平滑的轨迹。这一过程不仅能够消除由周期性边界引起的伪影，还能有效改善轨迹的连续性，从而提高分析结果的可靠性。掌握这种后处理技巧，对于分子动力学模拟中的链状分子体系研究具有重要价值，能够确保研究人员获得更加真实和精确的模拟数据，进而为分子设计和药物研发等领域提供有力支持。  

![](gmx后处理单链体系周期性边界条件处理生成平滑轨迹/gmx后处理单链体系周期性边界条件处理生成平滑轨迹_2025-01-12-17-02-47.png)  
## PBC后处理的必要性
原生轨迹如下所示：  
![](gmx后处理单链体系周期性边界条件处理生成平滑轨迹/gmx后处理单链体系周期性边界条件处理生成平滑轨迹_2025-01-12-17-10-47.gif)    

PBC处理之后生成平滑轨迹效果如下：  
![](gmx后处理单链体系周期性边界条件处理生成平滑轨迹/gmx后处理单链体系周期性边界条件处理生成平滑轨迹_2025-01-12-17-24-47.gif)    
## gmx trjconv 一些参数说明

## 单链体系PBC梳理通用流程
```shell
mkdir analysis
cd analysis
mkdir pbc
cd pbc

cp ../../npt/npt.gro ../../prod/npt.gro
cp ../../index.ndx .

echo "[ atom ]" >> index.ndx
echo "1" >> index.ndx

gmx trjconv -f ../../prod/npt.gro -s ../../prod/prod.tpr -o new.pdb -n index.ndx   # 选择需要输出的group，拿npt.gro生成的pdb文件作为参考结构便于可视化轨迹。   
gmx trjconv -f ../../prod/prod.xtc -s ../../prod/prod.tpr -o md_pbcmol_new.xtc -pbc atom -ur compact -center -n index.ndx   # 首先选择group atom将1号原子放在盒子中心。然后选择输出整个system。
gmx trjconv -f md_pbcmol_new.xtc -s ../../prod/prod.tpr -o md_pbcwhole_new.xtc -pbc whole -n index.ndx   # 选择输出整个system。  
gmx trjconv -f md_pbcwhole_new.xtc -s ../../prod/prod.tpr -o md_pbcfit_all_new.xtc -fit rot+trans -n index.ndx   # 首先选择蛋白进行对齐。然后选择蛋白进行输出。  
rm md_pbcmol_new.xtc md_pbcwhole_new.xtc
```
## 最终效果