#!/bin/bash

# 0. Build 50/30-system by using Gromacs insert-molecule:
# Monomer=wtd_H95_15ns
 mer30=md_D20F225_30
# gmx_mpi insert-molecules -ci ${Monomer}.pdb -nmol 30 -box 46. 46. 46. -seed 432567 -o $mer30.pdb
# Inserting is not ideal, many molecules will be positioned on the edge of box sticking out. One can them in VMD: pbc box to see the box, then Mouse->Move->Fragment
# Then the chains need to be separated by TER, MODEL, and ENDMDL removed in TextEditor. Then chains can be renamed in DS from all X to A,B,C,...

InpPDB=${mer30}_ds
# the molecule name:
MOL=md_D20F225_30_cg
# Edit ${MOL}.top by changing lines
#include "/home/luda/Documents/Prion/Modeling/martini3/martini_3.0.b.3.2/martini_v3.0.b.3.2_resc4.itp"
#include "/home/luda/Documents/Prion/Modeling/martini3/martini_3.0.b.3.2/martini_v3.0_solvents.itp"
#include "/home/luda/Documents/Prion/Modeling/martini3/martini_3.0.b.3.2/martini_v3.0_ions_.itp"
# Then edit all 30 itp-files molecule_29.itp by deleting long SS-string

# 1. Converting existing all-atom-pdb into coarse-grained=cg:
#martinize2 -f ${InpPDB}.pdb -ff martini30b32 -x ${MOL}.pdb -dssp /usr/bin/dssp -ef 500. -el 0.3 -eu 0.8 -elastic -scfix -cys auto -p backbone -o ${MOL}.top -maxwarn 30 &> ${InpPDB}.martinize.out
#gmx_mpi editconf -f $MOL.pdb -o $MOL.gro &> $MOL.gengro.out

# 2. Generating restraints file for minimizations
echo 2 q | gmx_mpi make_ndx -f $MOL.gro -o $MOL.ndx  # choice 2;q
echo 2 | gmx_mpi genrestr -f $MOL.gro -n $MOL.ndx -o posre_0.itp -fc 0 0 0 
echo 2 | gmx_mpi genrestr -f $MOL.gro -n $MOL.ndx -o posre_10.itp -fc 10 10 10 
echo 2 | gmx_mpi genrestr -f $MOL.gro -n $MOL.ndx -o posre_100.itp -fc 100 100 100 
echo 2 | gmx_mpi genrestr -f $MOL.gro -n $MOL.ndx -o posre_1000.itp -fc 1000 1000 1000 
echo 2 | gmx_mpi genrestr -f $MOL.gro -n $MOL.ndx -o posre_10000.itp -fc 10000 10000 10000 

# 3. Creating box around protein (NO -bt dodecahedron - moving molecules away from each other). Minimizing protein in vacuum 
gmx_mpi editconf -f $MOL.gro -o $MOL.gro -d 2.0  &> $MOL.genbox.out 
cp posre_10.itp posre.itp

cat > Min_noSol.mdp << _EOF_
;include             =  -I/ParmFiles
define               =  -DPOSRES -DFLEXIBLE
cutoff-scheme        =  group
nstxout              =  10000
nstcomm              =  1
nstcalcenergy        =  1
comm_mode            =  Linear ; remove com linear motion
integrator           =  steep
nstlist              =  1
ns_type              =  grid
rlist                =  1.4
rcoulomb             =  1.4
coulombtype          =  Reaction-Field
fourierspacing       =  0.135
pme_order            =  4
ewald_rtol           =  1e-5
rvdw                 =  1.4
epsilon_rf           = 78
pbc                  =  xyz
nsteps              =  500000
emtol               =  1.0
emstep              =  0.001
_EOF_

gmx_mpi grompp -f Min_noSol.mdp -c $MOL.gro -r $MOL.gro -p $MOL.top -o $MOL.tpr -maxwarn 1 &> ${MOL}_noSolMin.grompp.out
gmx_mpi mdrun -s $MOL.tpr -o ${MOL}_noSolMin.trr -e ${MOL}_noSolMin.edr -c ${MOL}_noSolMin.gro &> ${MOL}_noSolMin.mdrun.out

# 4. Adding solvent CG-water to box
cp $MOL.top ${MOL}_Wat.top
gmx_mpi editconf -f ${MOL}_noSolMin.gro -o ${MOL}_noSolMinBox.gro -d 2.0 &> ${MOL}_noSolMin.editconf.out  # -bt dodecahedron - again, not a good idea
gmx_mpi solvate -cp ${MOL}_noSolMinBox.gro -cs /home/luda/Documents/Prion/Modeling/martini3/water.gro -o ${MOL}_Wat.gro -p ${MOL}_Wat.top -radius 0.21 &> ${MOL}_Wat.solvate.out

# 5. Adding counter-ions into box
cp ${MOL}_Wat.top ${MOL}_WatNa.top
cat > Min_noSol2.mdp << _EOF_
;include             =  -I$CDIR/ParmFiles
define               =  -DPOSRES -DFLEXIBLE
nstxout              =  10000
nstcomm              =  1
nstcalcenergy        =  1
comm_mode            =  Linear ; remove com linear motion
cutoff-scheme        =  verlet ;group
integrator           =  steep
nstlist              =  10
ns_type              =  simple
rlist                =  1.5
rcoulomb             =  1.5
coulombtype          =  Reaction-Field
ewald_rtol           =  1e-5
rvdw                 =  1.5
pbc                  =  xyz
nsteps              =  500000
emtol               =  1.0
emstep              =  0.001
_EOF_
gmx_mpi grompp -f Min_noSol2.mdp -c ${MOL}_Wat.gro -r ${MOL}_Wat.gro -p ${MOL}_WatNa.top -o ${MOL}_WatNa.tpr -maxwarn 1 &> ${MOL}_WatNa.grompp.out
echo WN | gmx_mpi genion -s ${MOL}_WatNa.tpr -neutral -conc 0.05 -p ${MOL}_WatNa.top -o ${MOL}_WatNa.gro &> ${MOL}_WatNa.genion.out

# might need to edit ${MOL}_WatNa.top changing names of NA->TNA, CL->TCL; or /martini3/solvents.itp changing name of atom WN to W

echo "Ready for minimization"

#################################### Martinize install ###############################################################
# need to have python3 installed, made primary: 
# sudo update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1
# sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.6 2
# sudo update-alternatives --config python # choose 2
# sudo apt update
# sudo apt install python3-pip
# pip install vermouth
# martinize2 -h # test
##################################### New sm molecule FF  ############################################################
### copied cp -r /home/luda/apps/gromacs-2018.8/share/top/oplsaa.ff/ . here, renamed oplsaa.ff --> oplsaa_cbu.ff
### 
### cp PrP_C_Sc.top wtd_csc_cbu1.top
### modified wtd_csc_cbu1.top: added 
###   ; Include inositol topology
###   #include "oplsaa_cbu.ff/cbu_oplsaa.itp"
###   ...
###   oplsaa.ff/forcefield.itp --> oplsaa_cbu.ff/forcefield.itp
###   oplsaa.ff/spce.itp --> oplsaa_cbu.ff/spce.itp
###   oplsaa.ff/ions.itp --> oplsaa_cbu.ff/ions.itp
###   in: [ molecules ]
###   CBU                 6
###
### Converted all coordinates with editconf:
### gmx_mpi editconf -f new_all_inositol2.pdb -o new_all_inositol.gro
### Tested if works: 
### changed box size at the bottom of gro:    7.80000   9.60000   7.80000
### gmx_mpi grompp  -f Min_noSol.mdp -c wtd_csc_cbu1.gro -p wtd_csc_cbu1.top -o tmp.tpr 
### there is a warning about charge, but we do not worry about this right now.
######################################################################################################################
#gmx_mpi genion -s ${MOL}_WatNa.tpr -np 2 -nn 6 -pname NA -nname CL -rmin 0.7 -p ${MOL}_WatNa.top -o ${MOL}_WatNa.gro &> ${MOL}_WatNa.genion.out # choice 15-SOL



gmx_mpi insert-molecules -f bck_wtd_H95_2.gro -ci Gdn.pdb -nmol 10265 -o  bck_wtd_H95_Gdn.gro
