# 分子碎片化工具MacFrag

在药物发现和化学研究中，分子碎片化是一个至关重要的步骤。它不仅能帮助我们理解复杂分子的结构，还能为虚拟筛选、骨架跃迁和药物设计提供丰富的片段库。然而，传统的分子分割方法往往存在效率低、片段多样性不足等问题，难以满足大规模数据处理的需求。  

今天，我们为大家介绍一款高效、灵活的分子碎片化工具——MacFrag。这款工具不仅显著提升了分割速度，还通过改进的BRICS规则支持环键断裂，生成更多新颖且高质量的片段。无论你是药物研发人员，还是化学信息学爱好者，MacFrag都将成为你探索化学空间、构建高质量片段库的得力助手。  

## MacFrag使用
这里使用方式是在Linux平台上运行。  
**1. 从代码仓库中下载MacFrag工具。** 从https://github.com/yydiao1025/MacFrag/tree/main中下载MacFrag_Linux.tar.xz文件并将文件上传Linux服务器中。  
**2. 软件工具解压缩。**  
```shell
tar xvf MacFrag_Linux.tar.xz
```
**3. 使用示例。**  
```shell
cd MacFrag.dist
./MacFrag -i /data/MacFrag/examp.smi -o /data/MacFrag/ -maxBlocks 6 -maxSR 8 -asMols False -minFragAtoms 1
```
**4. 参数说明。**  


**（5）实例探究maxBlocks参数。**  


maxBlocks用于控制生成的片段中包含的最大构建块。构建块：分子被分割后的最小结构单元（如苯环、酰胺键等）。​片段：由一个或多个构建块组成的子结构。  

示例：假设有一个分子 ​A，其构建块如下：​构建块1：苯环；​构建块2：酰胺键；​构建块3：羟基；​构建块4：甲基。如果maxBlocks设置为2，则片段可能包括：​苯环；​酰胺键；​羟基；​甲基；苯环+酰胺键；酰胺键+羟基；羟基+甲基。  

较小的maxBlocks：片段简单，计算速度快，适用于骨架跃迁；较大的maxBlocks：片段数量多复杂，计算速度慢，适用于虚拟筛选。调整 maxBlocks可以覆盖不同大小和复杂度的化学空间。    

## MacFrag文献

MacFrag 分割大规模分子，生成高质量且多样化的分子片段。  

MacFrag 是一种高效的分子分割方法，具有以下特点：

​基于改进的BRICS规则：优化了化学键的断裂规则。
​引入高效子图提取算法：显著提升了片段枚举的速度。
​性能优越：在 ChEMBL 数据集上的评估表明，MacFrag 比 RDKit 中的 BRICS 和改良版 molBLOCKS 更快。
​高质量片段：生成的片段更符合 ‘Rule of Three’，适合药物发现。


分子分割的重要性：通过片段替换可以优化先导化合物的结构。  
现有的方法：RECAP和BRICS是目前基于逆合成化学的两种经典分割算法。  
现有方法的局限性：（1）现有方法通常只断裂非环键，限制了片段的多样性。（2）现有方法倾向于将分子切割成最小的构建块，可能破坏药效团。（3）在处理大规模数据时存在性能瓶颈。  

MacFrag的技术细节：
**键断裂识别**  
1. **识别断裂键：** 基于BRICS规则，但是去除了仅断裂非环键的限制。**定义原子环境：** 使用SMARTS字符串定义了19种原子环境。**组合断裂键：** 将这些环境组合成49种可断裂的键。  
2. 用户可通过maxSR参数控制环键断裂。如果环的原子数目 ≤ maxSR 则不断裂该环。如果maxSR很大则于原始BRICS一致，不断裂任何环键。  
**子图和片段枚举**   
1. 将分子简化为图，使用 Simple 算法高效枚举子图。  
2. 提取片段并标记断裂位点，过滤冗余片段。  
3. 用户可通过 maxBlocks 参数控制片段大小。  

性能评估：  
1. **评估对象：** MacFrag，RDKit中的BRICS和改良版的molBLOCKs。
**计算性能评估**  
  
2. **性能对比：** 小分子（分子量＜500），MacFrag比BRICS快2.5倍，比molBLOCKs快11.8倍。大分子（分子量500-1000），MacFrag比BRICS快7.9倍，比molBLOCKs快104倍。  

**片段质量评估**  
1. 数据集：ChEMBL中1921745个分子（分子量＜1000）。  
2. 生成片段数量对比：MacFrag（10,336,743）；BRICS（20,733,058）；molBLOCKS（27,468,335）。  
3. 分子量（中位分子量）对比：MacFrag（288.3）；molBLOCKS（369.2），BRICS（389.4）。  
4. 符合 RO3 的片段比例：MacFrag（13%）；molBLOCKS（5%）；BRICS（6%）。  
MacFrag 通过环键断裂生成了更多独特的片段，提高了片段多样性。MacFrag 生成了 82,758 个独特的 RO3 合规片段。  

## 参考
文献：Diao Y, Hu F, Shen Z, et al. MacFrag: segmenting large-scale molecules to obtain diverse fragments with high qualities[J]. Bioinformatics, 2023, 39(1): btad012.   
代码：https://github.com/yydiao1025/MacFrag/tree/main  