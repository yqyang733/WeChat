# gmx命令|gmx rms计算与参考结构之间的RMSD

## temp
```shell
# 一般情况下先对轨迹文件进行处理，变成流畅的动画格式的xtc文件再计算RMSD。
cp ../../index.ndx .
gmx rms -f ../pbc/md_pbcfit_all_new.xtc -s ../../prod/prod.tpr -o rms.xvg -n index.ndx
```