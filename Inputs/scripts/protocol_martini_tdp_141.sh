#!/bin/bash
MOL=tdp_ctd50_cg

#martinize2 -f tdp_ctd50_DS.pdb -ff martini30b32 -x tdp_ctd50_cg.pdb -dssp /usr/bin/dssp -elastic -el 0.3 -eu 0.8 -ef 500. -scfix -cys auto -p backbone -o tdp_ctd50_cgZ.top &> martinize2.out
#gmx_mpi editconf -f $MOL.pdb -o $MOL.gro -box 44.37 44.37 44.37 -noc &> $MOL.gengro.out -d 0.
# change box size at the end of gro file
#gmx_mpi grompp -f min-vac.mdp -c ${MOL}.gro -p ${MOL}.top -o min-vac.tpr &> min-vac.grompp.out
#gmx_mpi mdrun -deffnm min-vac -v &> min-vac.mdrun.out

#gmx_mpi solvate -cp min-vac5.gro -cs ../../CWD/WTD_24_233_WT_50_nacl100/water_CG_310K_1bar.gro -radius 0.21 -o ${MOL}_Wat.gro &> ${MOL}.solvate.out
#cp ${MOL}.top ${MOL}_Wat.top
#NWATERS=$(grep -c WN ${MOL}_Wat.gro)
#echo -e "\nWN  $NWATERS" >> ${MOL}_Wat.top
#gmx_mpi grompp -f min.mdp -c ${MOL}_Wat.gro -r ${MOL}_Wat.gro -p ${MOL}_Wat.top -o min.tpr &> min.grompp.out -maxwarn 1
#cp ${MOL}_Wat.top ${MOL}_WatNa.top
#echo WN | gmx_mpi genion -s min.tpr -neutral -nn 5357 -np 5263 -p ${MOL}_WatNa.top -o ${MOL}_WatNa.gro &> ${MOL}_WatNa.genion.out
# 945

#gmx_mpi grompp -f min.mdp -c ${MOL}_WatNa.gro -r ${MOL}_WatNa.gro -p ${MOL}_WatNa.top -o min.tpr -maxwarn 1 &> min.grompp.out
#gmx_mpi mdrun -deffnm min -v &> min.mdrun.out

gmx_mpi grompp -f eq-10fs.mdp -c min.gro -r min.gro -p ${MOL}_WatNa.top -o eq-10fs.tpr &> eq-10fs.grompp.out
gmx_mpi mdrun -deffnm eq-10fs -v &> eq-10fs.mdrun.out

gmx_mpi grompp -f eq-20fs.mdp -c eq-10fs.gro -r eq-10fs.gro -p ${MOL}_WatNa.top -o eq-20fs.tpr -maxwarn 1 &> eq-20fs.grompp.out
gmx_mpi mdrun -deffnm eq-20fs -v &> eq-20fs.mdrun.out

gmx_mpi grompp -f dynamic_works.mdp -c eq-20fs.gro -r eq-20fs.gro -p ${MOL}_WatNa.top -o pr00.tpr -maxwarn 1 &> pr00.grompp.out
gmx_mpi mdrun -deffnm pr00 -v &> pr00.mdrun.out



