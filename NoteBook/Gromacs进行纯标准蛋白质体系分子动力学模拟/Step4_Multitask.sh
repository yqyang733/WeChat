generate_files(){
    ls *.pdb|awk -F "." '{print $1}' > lst
    for i in `cat lst`;do mkdir ${i};mv ${i}.pdb ${i};done
}

build_system(){
    for i in `cat lst`
    do
        cd ${i}
        python ../Step1_generate_gro_top_command.py ${i}.pdb   
        python ../Step2_generate_mdp.py    
        python ../Step3_generate_submit_sh.py 
        cd ..
    done
}

input=$*
generate_files
build_system