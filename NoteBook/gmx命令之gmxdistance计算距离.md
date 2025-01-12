# gmx命令|gmx distance计算距离
## Temp
计算index中溶质组SOLU和膜组MEMB质心之间的距离   
```shell
gmx distance -s ../pbc/equil.pdb -f ../pbc/prod_pbc.xtc -select "com of group SOLU plus com of group MEMB" -oav -oall -n ../../index.ndx   # 计算index中溶质组SOLU和膜组MEMB质心之间的距离
```