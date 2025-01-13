# gmx命令|gmx cluster聚类

## temp
```shell
echo 22 22|gmx cluster -s ../../prod/prod.tpr -f ../pbc/md_pbcfit_all_new.xtc -g -dist -sz -clid -cl -method linkage -cutoff 0.2 -n ../../index.ndx
```