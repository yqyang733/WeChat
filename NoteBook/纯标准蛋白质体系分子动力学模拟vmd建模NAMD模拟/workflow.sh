build_md(){
    cat > pipline.tcl << EOF
package require psfgen
psfcontext reset
topology /public/home/yqyang/file/vegf-toppar/top_all36_prot.rtf
topology /public/home/yqyang/file/vegf-toppar/toppar_water_ions.str
pdbalias residue HIS HSE
alias atom ILE CD1 CD
alias atom SER HG HG1
alias atom CYS HG HG1
segment REC {pdb REC.pdb}
coordpdb REC.pdb REC
segment PEP {pdb PEP.pdb}
coordpdb PEP.pdb PEP
guesscoord
writepdb merged.pdb
writepsf merged.psf
exit
EOF
    ~/../../software/apps/vmd/1.9.3/vmd -dispdev text -e pipline.tcl
    rm pipline.tcl
    rm do.py
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
    cat > pipline.tcl << EOF
package require pbctools
psfcontext reset
mol load psf merged.psf pdb merged.pdb
package require solvate
solvate merged.psf merged.pdb -minmax ${box_size} -o solvated
mol delete all
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o complex
pbc box -center centerofmass
exit
EOF
    ~/../../software/apps/vmd/1.9.3/vmd -dispdev text -e pipline.tcl
    rm pipline.tcl
    cp ../mk_md_namd.py .
    python mk_md_namd.py md_b
}
build_md