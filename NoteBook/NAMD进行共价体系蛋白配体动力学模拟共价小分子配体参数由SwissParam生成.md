# NAMD进行共价体系蛋白配体动力学模拟，共价小分子配体参数由SwissParam生成
## Temp
```shell
$ ~/../../software/apps/vmd/1.9.3/vmd -dispdev text
vmd > package require psfgen
vmd > psfcontext reset
vmd > topology /public/home/yqyang/file/vegf-toppar/top_all36_prot.rtf
vmd > topology /public/home/yqyang/file/vegf-toppar/toppar_water_ions.str
vmd > pdbalias residue HIS HSE
vmd > alias atom ILE CD1 CD
vmd > alias atom SER HG HG1
vmd > alias atom CYS HG HG1
vmd > segment PRO {pdb Protein.pdb}
vmd > coordpdb Protein.pdb PRO
vmd > topology post.rtf
vmd > segment LIG {pdb post.pdb}
vmd > coordpdb post.pdb LIG
vmd > topology reaction.str
vmd > patch REAC PRO:72 LIG:1
vmd > regenerate angles dihedrals
vmd > guesscoord
vmd > writepsf test.psf
vmd > writepdb test.pdb
vmd > quit
```
