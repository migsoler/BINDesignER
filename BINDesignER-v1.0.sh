#!/bin/sh
##PBS -N TESTMIG 
##PBS -l nodes=1:ppn=20
##PBS -q qgpu
##PBS -l walltime=5000:00:00


### EXTERNAL INPUT FILES ######

# Start.pdb
# Start.gro and topology files: topol-0.top, topol_Protein_chain_A-0.itp, topol_Protein_chain_B-0.itp(, solvent-0.itp if organic solvent) 
# Important: Target = chain A ; Binder = chain B

# mdp FILES:
# minimization input files: minim_scmut.mdp, minim_overlap.mdp and minim.mdp, for the minimization of mutants. minim_scmut: minimization of mutated side-chain and the two attached residues using only the binding complex; minim_overlap: minimization of the mutated side chain and the surrounding solvent molecules; minim: global minimization.

# Input files for NVT and NPT MD simulations: md-NVT.mdp, md-NPT.mdp

# sf folder: folder with all scoring programs.

# scwrl4: program Scwrl4 for homology modelling.

# inputs for haddock (optional):
# scoring0.inp: initial haddock file
# inpA, inpC: input haddock file prepared by pieces to include new mutations 

# topology file of water molecules where the mass of H is modified to speed up the thermodynamics (optional, not active)
# tip3p-m.itp

##############################



##############################
# CHANGE PROGRAM PATHS HERE!!!
# environment variable definitions
#cd $PBS_O_WORKDIR

##############################

#--Function-Set-Up
setup () {
##### LOADING MODULES, IF NECECSSARY #######
#module load profile/chem intel/pe-xe-2017--binary intelmpi/2017--binary matheval/1.1.11--intelmpi--2017--binary gromacs/5.1.2
#module load mkl/2017--binary amber/16.0
source /etc/profile.d/modules.sh
#module purge
#source /u/sbp/igladich/opt/MGLTools-1.5.6/bin/mglenv.sh
#source ~/.bash_profile
#module load python/2.7.5
#module load numpy/1.9.2/python/2.7
#module load Gromacs/4.6.7s2
#module load gcc/4.8.2 cudatoolkit/6.0 openmpi/1.6.4/gcc
#module load gromacs/5.0
module load Gromacs/4.6.1d/intel/13.1
###########################################

#restart the simulation: yes (1) or no (0)
rest=1

#names of gromacs (depends on the version used)
trjconv="trjconv_mpi_d"
make_ndx="make_ndx_mpi_d"
grompp="grompp_mpi_d"
editconf="editconf_mpi_d"
mdrun="mdrun_mpi_d -v -ntomp 8"
#mdrund=mdrun_mpi_d
g_select="g_select_mpi_d"
pdb2gmx="pdb2gmx_mpi_d"


# Method of mutation acceptance: Monte Carlo (MC) or consensus criterion (consensus)
ACCMET=consensus

################## MC OPTIONS #####################################
# scoring function you want to use in MC:
# Irad=irad, Pie*Pisa=pisa, Bluues=bluues, Bach6=bach, Haddock2.1=hadd, FireDock=fired
# VINA=vina
SF_MC=irad

# Montecarlo Temperature for mutation process
Trep=(3)

# MC Replica exchange active (1) or inactive (0)
REXCH="0"
###################################################################

################## consensus criterion options ####################
# List of scoring functions to be used:
# Irad=irad, Pie*Pisa=pisa, Bluues=bluues, Bach6=bach, Haddock2.1=hadd, FireDock=fired
# VINA=vina
SF_CONS=(irad pisa bluues hadd bach fired)

# consensus threshold (how many SF must be evaluated positively to accept the mutation)
Tcons=3
###################################################################

# Run MD for the initial complex (yes, no) 
MD0=no

# Export trajectories from $bprint ps
bprint="100"

# Minimum distance to delete solvent molecules after mutation (nm)
disSOL=0.2

#list of residues that can be mutated
residue_list=(13 14 15 16 17 18 19 20 21 22 23 24 25) 
# Total number of Montecarlo steps or of mutations attemps
NMT=5

# specify the type of solvent: Water_and_ions, Water, SOL (for organic solvents).
typesol=Water_and_ions

# number of atoms of your solvent molecule (e.g. Water: 3)
nsolv=3

# List of amino acids used during mutation process
aalistL=(ALA ASP GLU GLY PHE HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR)

# initiating parameters
MM=0
acc=0
# Total number of residues for mutation
NA=${#residue_list[@]}
# Total number of different amino acids used
NAA=${#aalistL[@]}
# Total number of replicas in MonteCarlo mutation
Nrep=${#Trep[@]}
# Total number of scoring functions in consensus score option
Nsf=${#SF_CONS[@]}

# current folder
heref=`pwd`
cd ..
here0=`pwd`
cd ${heref}

# folder scrwl
scrwl=${here0}/scwrl4/Scwrl4

# save initial number of solvent molecules
SOL_orig=`grep SOL topol-0.top | awk '{print $2}'`

# loading environment variables and paths of the scoring functions
PIE="${here0}/sf/pie/bin/pie_score"
PIEpar="${here0}/sf/pie/bin/pie.params"
PISA="${here0}/sf/pisa/pisaEnergy_linux"
export PIEDOCK_HOME="${here0}/sf/pie"
PISApar="${here0}/sf/pisa/pisa.params"
IRAD="${here0}/sf/irad/irad"
source ${here0}/sf/VINA/MGLTools-1.5.6/bin/mglenv.sh
UTILITIES="${here0}/sf/VINA/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24"
VINA_FOLD="${here0}/sf/VINA/autodock_vina_1_1_2_linux_x86/bin"

# Setting up the list. First score in the list is always the one used in MC
if [ ${ACCMET} = consensus ]; then
    echo "Consensus method"
    for ((isf=0;isf<${Nsf};isf++)); do
	sf_list[${isf}]=${SF_CONS[${isf}]}
	echo "scoring function: ${sf_list[${isf}]}, ${isf} of ${Nsf}"
    done
else
    sf_list[0]=${SF_MC}
    Nsf=1
    echo "MonteCarlo method, scoring function: ${sf_list[0]}"
fi
}
#--END-Function-Set-Up-----------------------------------------

#--Function-initial_box
gromacs () {
####### here you can change it to adapt it to your hpc environment ###########

    # example 1. MD run in bash parallel if you have more than one replica
    if [ ${Nrep} -eq 1 ]; then
	it=0
        Teff=${Trep[${it}]}
        rundir=Teff_${Trep[${it}]}
        cd ${rundir}
	pwd
	${mdrun} -deffnm md-npt &> mdrun-npt_${MM}
	cd ..
    fi

    # syncro all parallel simulations
    if [ ${Nrep} -gt 1 ]; then
	for ((it=0;it<${Nrep};it++)); do
	    (
		Teff=${Trep[${it}]}
		rundir=Teff_${Trep[${it}]}
		cd ${rundir}
		pwd
		${mdrun} -deffnm md-npt &> mdrun-npt_${MM}
		cd ..
		touch DONE_${Teff}
	    ) &
	done

	((c=0))
	while [ $c -lt ${Nrep} ];do
	    sleep 1
	    ((c=0))
	    for ((iit=0;iit<${Nrep};iit++)); do
		Teff=${Trep[${iit}]}
		if [ -e DONE_$Teff ]; then
		    ((c++))
		fi
	    done
	    if [ -f EXIT ]; then
		exit
	    fi
	done
    fi
    rm DONE_*


# example 2. Multidir option, by using mpirun and more than one node.    
#export OMP_NUM_THREADS=32
#mpirun -x OMP_NUM_THREADS=10 -npernode 2 ${BIN_DIR}/mdrun_4_7_gpu_s -v -multidir Teff_[43] -pin on -deffnm md-npt &> mdrun-npt_${MM}
#mpirun -x OMP_NUM_THREADS=4 -np 2 ${mdrun} -v -multidir Teff_[42] -gpu_id "00" -deffnm md-npt &> mdrun-npt_${MM}

}

#--END-Function-Gromacs-----------------------------------------

#--Function-initial_box
initial_box () {

####### PREPARATION INITIAL STRUCTURE ################
####### here you can change it to adapt it to your initial files ###########

# initial cleaning
rm -f REMC.dat
rm -f evol*
rm -f iradsc* bluuesonly* hadd* ppisasc* FireDock* Bach-sint* VINA*
rm -fr Teff*
rm -f system*
rm -f output*
rm -f DONE_*
rm -f mdrun-npt_*
rm -f solvent*.pdb binder*.pdb target*.pdb box*.txt
rm -f \#*
if [[ $REXCH -eq 1 ]];then
touch "REMC.dat"
fi

# creation input file in_templ    
echo 'FFFF&!LLLL' > in_templ
echo 'q' >> in_templ

# CREATION OF INPUT FILES NECESSARIES FOR SOME SCORING FUNCTIONS
# loop of scoring functions
for ((isf=0;isf<${Nsf};isf++)); do
    SF=${sf_list[${isf}]}
# creation input file for VINA
if [ ${SF} = vina ]; then
    echo "receptor = receptor.pdbqt" > configure.txt
    echo "ligand = ligand.pdbqt" >> configure.txt
    echo "num_modes = 1" >> configure.txt
    echo "center_x = 0.26" >> configure.txt
    echo "center_y = -1.39" >> configure.txt
    echo "center_z = -5.28" >> configure.txt
    echo "size_x = 25" >> configure.txt
    echo "size_y = 25" >> configure.txt
    echo "size_z = 25" >> configure.txt
    echo "exhaustiveness = 1" >> configure.txt
    echo "energy_range = 2" >> configure.txt
    echo "cpu=4" >> configure.txt
fi
done

# Creation of gromacs input files 
#$pdb2gmx_mpi_d -f complex_0.pdb -p comp.top -o comp.gro < inp-pdb2g.in > gromacs.log

# Creation of index file    
echo "chain A" >select.txt # Target
echo "chain B" >>select.txt # Binder
echo "q" >>select.txt
$make_ndx -f Start.pdb -o fixbs2.ndx<select.txt>gromacs.log

# copying the initial topology files (if you have them)
sed "s/-0.itp/.itp/g" topol-0.top > topol.top
cp topol_Protein_chain_A-0.itp topol_Protein_chain_A.itp
cp topol_Protein_chain_B-0.itp topol_Protein_chain_B.itp
if [[ ${typesol} = SOL ]]; then
cp solvent-0.itp solvent.itp    
fi
#$editconf -f Start.pdb -o Start.gro

# creating the topology file of binder-target complex (if you dont have it)
sed "/Ion/d" topol-0.top > complex.top
sed -i "/SOL/d" complex.top
sed -i "/NA/d" complex.top
sed -i "/CL/d" complex.top
sed -i "/solvent/d" complex.top

#create initial pdb files from the original start.pdb
grep "ATOM" Start.pdb > temp
grep " A " temp > targ0.pdb 
grep " B " temp > bind0.pdb
grep -v " A " temp > temp2
grep -v " B " temp2 > SOLBOX0.pdb
rm temp*
#echo -e "non-Protein" | $editconf -f ${CONF} -n fixbs2.ndx -o SOLBOX0.pdb
#echo -e "target" | $editconf -f ${CONF} -n fixbs2.ndx -o prot0.pdb   
#echo -e "binder" | $editconf -f ${CONF} -n fixbs2.ndx -o bind0.pdb    

# saving the last residue number
LAA=`grep ATOM targ0.pdb|tail -1 | awk '{print $4}'`
LAB=`grep ATOM bind0.pdb|tail -1 | awk '{print $4}'`

}
#--END--INITIAL---BOX
##################################################################

    ####################################### END OF INPUT SET UP. FROM HERE, DO NOT TOUCH ##################################################    

#--Function-block analisis
block () {
    n_data=`wc -l $1 | awk '{print $1}'`

awk -v N=${n_data} -v nblock=4 '
	BEGIN { 
		nconf=int( N/nblock);
		j=1; conta=1; mediavera=0;

		for(i=1;i<=nblock;i++){
			media[i]=0;
		}
	}
	
	{	
		if( j <= nconf*conta  ){
			media[conta]= media[conta]+$1;
		}
		
		j++;
		if(j%nconf==0) {conta=conta+1;}	
		mediavera=mediavera+$1;
		
	}


	END { 	
		mediavera = mediavera/N;
		
		for(i=1;i<=nblock;i++){
			media[i]=media[i]/nconf;
			sigma2 = sigma2 + (media[i] - mediavera)**2    ; 
    		}
	
		sigma2=sigma2/nblock;

		errore = (sigma2 /nblock)**0.5;
	
		print errore, media[1]-media[nblock];

	}' $1 > temp
}
# END block ##############################################################

#--Function compress sequence
compress () {
###################################################################
#### To convert a sequence from 3 letter code to one letter code ##
###################################################################
oneaa () {
#local variables
local ii=0
SEQ=''
#for (( ii=0; ii<=$NA; ii++ )); do
        if [[ ${aaseq[${ii}]} = ALA ]]
        then
                tmp=a
        fi
        if [[ ${aaseq[${ii}]} = CYS ]]
        then
                tmp=c
        fi
        if [[ ${aaseq[${ii}]} = CYX ]]
        then
            tmp=c
        fi
        if [[ ${aaseq[${ii}]} = ASP ]]
        then
                tmp=d
        fi
        if [[ ${aaseq[${ii}]} = GLU ]]
        then
                tmp=e
        fi
        if [[ ${aaseq[${ii}]} = PHE ]]
        then
                tmp=f
        fi
        if [[ ${aaseq[${ii}]} = GLY ]]
        then
                tmp=g
        fi
        if [[ ${aaseq[${ii}]} = HIS ]]
        then
                tmp=h
        fi
        if [[ ${aaseq[${ii}]} = HIE ]]
        then
                tmp=h
        fi
        if [[ ${aaseq[${ii}]} = HID ]]
        then
                tmp=h
        fi
        if [[ ${aaseq[${ii}]} = HIP ]]
        then
                tmp=h
        fi
        if [[ ${aaseq[${ii}]} = ILE ]]
        then
                tmp=i
        fi
        if [[ ${aaseq[${ii}]} = LYS ]]
        then
                tmp=k
        fi
        if [[ ${aaseq[${ii}]} = LEU ]]
        then
                tmp=l
        fi
        if [[ ${aaseq[${ii}]} = MET ]]
        then
                tmp=m
        fi
        if [[ ${aaseq[${ii}]} = ASN ]]
        then
                tmp=n
        fi
        if [[ ${aaseq[${ii}]} = PRO ]]
        then
                tmp=p
        fi
        if [[ ${aaseq[${ii}]} = GLN ]]
        then
                tmp=q
        fi
        if [[ ${aaseq[${ii}]} = ARG ]]
        then
                tmp=r
        fi
        if [[ ${aaseq[${ii}]} = SER ]]
        then
                tmp=s
        fi
        if [[ ${aaseq[${ii}]} = THR ]]
        then
                tmp=t
        fi
        if [[ ${aaseq[${ii}]} = VAL ]]
        then
                tmp=v
        fi
        if [[ ${aaseq[${ii}]} = TRP ]]
        then
                tmp=w
        fi
        if [[ ${aaseq[${ii}]} = TYR ]]
        then
                tmp=y
        fi
        SEQ=$SEQ$tmp
#done
echo $SEQ
}


 while read line
 do 
  aaseq=($line)
  SEQ=$(oneaa) # get the one-letter code sequence
  echo $SEQ
 done<$1

}

####################################################################################

#--Function-Mutated-Box
mutated_box () {

# initial clean
rm index.ndx *overlapping.pdb indexsol.ndx tot.ndx newSOL.pdb
# copy the solvent file
cp  box_solvent_to_use.pdb  SOL.pdb

# prepare the pdb file of the system. Adding the solvent. 
head -n4 SOL.pdb  > primo
nlines=$(wc -l SOL.pdb | awk '{print $1-4}')
tail -n ${nlines} SOL.pdb > last
grep "ATOM" complex.pdb > secondo
cat primo secondo last > all.pdb
${editconf} -f all.pdb -o solvated_and_overlapping.pdb
rm -f primo secondo last tmp.ivan all.pdb

# remove the overlapping solvent molecules. At 2 A from heavy atoms    
echo -e "chain B \n q" | ${make_ndx} -f solvated_and_overlapping.pdb
${g_select} -f solvated_and_overlapping.pdb -n index.ndx -s solvated_and_overlapping.pdb -select "group ${typesol} and same residue as within ${disSOL} of (group chB and resnr ${resid})" -on indexsol.ndx

##################### new
echo "\"SideChain\" & \"chB\" & r ${Nmut}" >notover.txt #group sidechain of mutated residue
echo "name 0 overlap" >>notover.txt #changing name
echo "\"overlap\" | \"SideChain_&_chB_&_r_${resid}\"" >>notover.txt #merging sidechain of mutated resid and olverlapping solvent molecules
echo "\"System\" &! \"overlap_SideChain_&_chB_&_r_${resid}\"" >>notover.txt
echo "q" >>notover.txt

${make_ndx} -f solvated_and_overlapping.pdb -n indexsol.ndx index.ndx -o min_sol.ndx<notover.txt
sed -i "s/System_&_\!overlap_SideChain_&_chB_&_r_${resid}/to_block/g" min_sol.ndx

${editconf} -f solvated_and_overlapping.pdb -o solvated_and_overlapping.gro
# Partial minimization of selected solvent molecules and the side chain of mutated residue
${grompp} -f minim_overlap.mdp -n min_sol.ndx -c solvated_and_overlapping.gro -p topol.top -o em_overlap.tpr
${mdrun} -deffnm em_overlap
cp em_overlap.gro solvated_and_NOToverlapping.gro
##############################################################################################################

#### old one #####
# Obtain the numbers of the type of molecules
#first=`grep "${typesol}" pru | awk '{print $1}'`
#last=`grep "chB" pru | awk '{print $1}'`
#(( last=last+1 ))
#rm pru

# Creating the new index and deleting the solvent molecules which are placed less than disSOL nm of mutated residue
#sed "s/LLLL/${last}/g" in_templ | sed "s/FFFF/${first}/g" | ${make_ndx} -n index.ndx indexsol.ndx -o tot.ndx
#(( last=last+1 ))
#echo -e "${last}" | ${editconf} -f solvated_and_overlapping.pdb -n  tot.ndx -o newSOL.pdb

# Calculating the difference between the number of solvent molecules
#SOL_before=$(grep "SOL" SOL.pdb | wc -l |awk -v a=${nsolv} '{print $1/a}')
#echo ${SOL_before}
#SOL_after=$(grep "SOL" newSOL.pdb | wc -l |awk -v a=${nsolv} '{print $1/a}')
#echo ${SOL_after}
#DSOL=$((SOL_before-SOL_after))
#echo "DSOL= ${DSOL}"

# Changing the number in the topology file
#sed "s/${SOL_orig}/${SOL_after}/g" topol-0.top > topol.top

# creating the final pdb of the system
#head -n4 newSOL.pdb  > primo
#nlines=$(wc -l newSOL.pdb | awk '{print $1-4}')
#tail -n ${nlines} newSOL.pdb > last
#grep -v "TER" complex.pdb > secondo
#grep -v "REMARK" secondo > tmp.ivan
#grep -v "END" tmp.ivan > secondo
#cat primo secondo last > all.pdb
#${editconf} -f all.pdb -o solvated_and_NOToverlapping.pdb  
#rm -f primo secondo last tmp.ivan all.pdb

#Transforming to gro. The pdb does not have information on the box size
#${editconf} -f solvated_and_NOToverlapping.pdb -o solvated_and_NOToverlapping.gro
#sed -i '$ d' solvated_and_NOToverlapping.gro
#cat  solvated_and_NOToverlapping.gro box_dimension_to_use.txt > tmp.gro
#mv tmp.gro solvated_and_NOToverlapping.gro
# only if you want to accelerate the MD
#sed -i "s#amber99sb-ildn.ff/tip3p.itp#${heref}/tip3p-m.itp#g" system.top

####################################

#${make_ndx} -f solvated_and_NOToverlapping.gro -o fixbs2.ndx<select.txt>>gromacs.log 
cat  solvated_and_NOToverlapping.gro box_dimension_to_use.txt > tmp.gro
cp solvated_and_NOToverlapping.gro system.gro
cp topol.top system.top

# only if you want to fix your molecules
#atfix=`grep " N     35    TRP" system_Protein_chain_B.itp | awk '{print $1}'`
#sed "s/XXX/${atfix}/g" addrestB >> system_Protein_chain_B.itp
#cp -f topol_Protein_chain_A-0.itp system_Protein_chain_A.itp

rm -f \#*

}

#--END-Function-Mutated box

#####################################################################################################

# Function Scoring
scoring () {

HERE=`pwd`
rm -f scoring/*.pdb
#rm -f FireDock.ene 

# creation of pdb poses #######################################################################
# extract
echo -e "chA" | $trjconv -f npt-pbc.xtc -s md-npt.tpr -n fixbs2.ndx -o scoring/HER2.pdb -sep
echo -e "chB" | $trjconv -f npt-pbc.xtc -s md-npt.tpr -n fixbs2.ndx -o scoring/VHH.pdb -sep

cd scoring

# total number of poses
files=$(( `ls -l VHH*.pdb|wc -l`-1 ))

# creation file declist
rm -f declist
> declist

# formatting the pdb poses
for i in `seq 0 ${files}`
do
#sed -i "s/ HZ1 LYS/1HZ  LYS/g" HER2${i}.pdb
#sed -i "s/ HZ2 LYS/2HZ  LYS/g" HER2${i}.pdb
#sed -i "s/ HZ3 LYS/3HZ  LYS/g" HER2${i}.pdb
sed -i "s/CD  ILE/CD1 ILE/g" HER2${i}.pdb
#awk '$1=="ATOM"&&$6=="1" {print $4}' HER2${i}.pdb > pru
#first=`head -n 1 pru`
#sed -i "s/H1  ${first}/H   ${first}/g" HER2${i}.pdb
sed -i "s/OC1/O  /g" HER2${i}.pdb
sed -i "s/OC2/OXT/g" HER2${i}.pdb
#awk '$1=="ATOM"' HER2${i}.pdb | awk 'BEGIN {FS=""} {if($14!="H")print}' > tt.pdb
#cp tt.pdb HER2${i}.pdb
#awk '$1=="ATOM"' HER2${i}.pdb | awk 'BEGIN {FS=""} {if($13!="H")print}' > tt.pdb
#cp tt.pdb HER2${i}.pdb

#sed -i "s/ HZ1 LYS/1HZ  LYS/g" VHH${i}.pdb
#sed -i "s/ HZ2 LYS/2HZ  LYS/g" VHH${i}.pdb
#sed -i "s/ HZ3 LYS/3HZ  LYS/g" VHH${i}.pdb
sed -i "s/CD  ILE/CD1 ILE/g" VHH${i}.pdb
#awk '$1=="ATOM"&&$6=="1" {print $4}' VHH${i}.pdb > pru
#first=`head -n 1 pru`
#sed -i "s/H1  ${first}/H   ${first}/g" VHH${i}.pdb
sed -i "s/OC1/O  /g" VHH${i}.pdb
sed -i "s/OC2/OXT/g" VHH${i}.pdb
#awk '$1=="ATOM"' VHH${i}.pdb | awk 'BEGIN {FS=""} {if($14!="H")print}' > tt.pdb
#cp tt.pdb VHH${i}.pdb
#awk '$1=="ATOM"' VHH${i}.pdb | awk 'BEGIN {FS=""} {if($13!="H")print}' > tt.pdb
#cp tt.pdb VHH${i}.pdb

sed -i "s/HIE/HIS/g" HER2${i}.pdb
sed -i "s/HIE/HIS/g" VHH${i}.pdb
sed -i "s/HID/HIS/g" HER2${i}.pdb
sed -i "s/HID/HIS/g" VHH${i}.pdb
sed -i "s/HIP/HIS/g" HER2${i}.pdb
sed -i "s/HIP/HIS/g" VHH${i}.pdb

# creation of file list
echo "scr${i}.pdb" >> declist
(grep ATOM HER2${i}.pdb
 echo "TER" 
 grep ATOM VHH${i}.pdb
 echo "END") > scr${i}.pdb
done

# restarting output temporal file
echo $MM > tempsf

# loop of scoring functions
for ((isf=0;isf<${Nsf};isf++)); do
    SF=${sf_list[${isf}]}
#    echo "${SF}, SF ${isf}/${Nsf}" 
#IRAD ###
if [ ${SF} = irad ]; then
${IRAD} declist > /dev/null
awk '{print $2}' declist.irad.out > iradsc
block iradsc
scirad=0
scirad=` awk '{sum += $1} END {print sum/NR}' iradsc`
eirad=`awk '{print $1}' temp`
dirad=`awk '{print $2}' temp`
#sdirad=`awk '{x+=$2;y+=$2^2}END{print sqrt(y/NR-(x/NR)^2)}' iradsc`
dif=`echo "${scirad} ${scirad0}" | awk '{print $1 - $2}'`
contirad=`echo "${dif} <= 0.0" | bc -l`
cp iradsc ../iradsc${MM}
Enew=${scirad}
echo -e "${scirad} \n${eirad} \n${dirad}" >> tempsf
fi

#PIE*PISA ###
#echo "pisa"
if [ ${SF} = pisa ]; then
> orig.pie
> orig.pisa
for x in `cat declist`
	do ${PIE} ${x} A B ${PIEpar} >> orig.pie
	${PISA} ${x} A B ${PISApar} >> orig.pisa
	awk 'BEGIN{printf "\n"}' >> orig.pisa
done
paste orig.pie orig.pisa > pru
awk '{print $1*$2}' pru > ppisasc
block ppisasc

cp ppisasc ../ppisasc${MM}
scpisa=0
scpisa=` awk '{sum += $1} END {print sum/NR}' ppisasc`
#sdpisa=`awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' ppisasc`
episa=`awk '{print $1}' temp`
dpisa=`awk '{print $2}' temp`
dif=`echo "${scpisa} ${scpisa0}" | awk '{print $1 - $2}'` 
contpisa=`echo "${dif} <= 0.0" | bc -l`
Enew=${scpisa}
echo -e "${scpisa} \n${episa} \n${dpisa}" >> tempsf
fi

# BLUUES ######
#echo "bluues"
if [ ${SF} = bluues ]; then
> bluuesonly
for i in `seq 0 ${files}`
do
cp scr${i}.pdb 1complex.pdb
${here0}/sf/pdb2pqr-2.1.0/pdb2pqr --ff=AMBER 1complex.pdb tmp.pqr --chain >& /dev/null
grep " A " tmp.pqr > tmp_A.pqr
grep " B " tmp.pqr > tmp_B.pqr
#echo tmpc
${here0}/sf/bluues_new tmp.pqr tmpc  >& /dev/null
#echo tmpa
${here0}/sf/bluues_new tmp_A.pqr tmpa  >& /dev/null
#echo tmpb
${here0}/sf/bluues_new tmp_B.pqr tmpb  >& /dev/null
c=`(grep ^"Total    " tmpc.solv_nrg | awk '{print $3}')`
a=`(grep ^"Total    " tmpa.solv_nrg | awk '{print $3}')`
b=`(grep ^"Total    " tmpb.solv_nrg | awk '{print $3}')`
echo $c $a $b | awk '{print $1 - $2 - $3}' >> bluuesonly
done

block bluuesonly
cp bluuesonly ../bluuesonly${MM}
scbluues=0
scbluues=` awk '{sum += $1} END {print sum/NR}' bluuesonly`
#sdbluues=`awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' bluuesonly`
ebluues=`awk '{print $1}' temp`
dbluues=`awk '{print $2}' temp`
dif=`echo "${scbluues} ${scbluues0}" | awk '{print $1 - $2}'` 
contbluues=`echo "$dif <= 0.0" | bc -l`
Enew=${scbluues}
echo -e "${scbluues} \n${ebluues} \n${dbluues}" >> tempsf
fi

# HADDOCK #####
if [ ${SF} = hadd ]; then
#echo "hadd"
> Haddock.ene
rm -rf had
mkdir had
cp -rf ${here0}/sf/rescoring-scripts/* had/.
cp declist had/.
cd had
echo "\"had2.pdb\"" > filelist.list

for x in `cat declist`
do rm -rf had.pdb
rm -rf had2_conv.pdb
rm -rf had2_conv.psf
cp ../${x} had.pdb
cp -f ../../scoring.inp scoring.inp # You have to create .inp file before
./pdb_chain-segid had.pdb > had2.pdb
./run_scoring.csh
awk '$1=="energies:" {print $2}' had2_conv.psf | sed 's/,//' >> ../Haddock.ene
done
cd ..
block Haddock.ene
cp Haddock.ene ../Haddock${MM}.ene
schadd=0
schadd=` awk '{sum += $1} END {print sum/NR}' Haddock.ene`
#sdhadd=`awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' Haddock.ene`
ehadd=`awk '{print $1}' temp`
dhadd=`awk '{print $2}' temp`
dif=`echo "${schadd} ${schadd0}" | awk '{print $1 - $2}'` 
conthadd=`echo "$dif <= 0.0" | bc -l`
echo -e "${schadd} \n${ehadd} \n${dhadd}" >> tempsf
Enew=${schadd}
fi

# BACH-6SENSE #####
if [ ${SF} = bach ]; then
#echo "bach"
cp ${here0}/sf/BACH-SixthSense/BSS.* .
cp ${here0}/sf/BACH-SixthSense/ATOMIC_PARAMETERS_BSS .

./BSS.x -COMPUTE_ENE -STRICT_INTERFACE -PDBLIST declist -o complex.bss >& /dev/null
grep "E " complex.bss|awk '{print $4}' > Bach6-sint.dat
block Bach6-sint.dat
cp Bach6-sint.dat ../Bach6-sint${MM}.dat
scbach=0
scbach=` awk '{sum += $1} END {print sum/NR}' Bach6-sint.dat`
ebach=`awk '{print $1}' temp`
dbach=`awk '{print $2}' temp`
#sdbach=`awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' Bach6-sint.dat`
dif=`echo "${scbach} ${scbach0}" | awk '{print $1 - $2}'` 
contbach=`echo "$dif <= 0.0" | bc -l`
echo -e "${scbach} \n${ebach} \n${dbach}" >> tempsf
Enew=${scbach}
fi

# FIREDOCK #####
if [ ${SF} = fired ]; then
#echo "fired"
cp ${here0}/sf/FireDock/examples/1ACB/myex/1ACB.zdock.trans fake.trans
> FireDock.ene
for i in `seq 0 ${files}`
do
cp HER2${i}.pdb chainA.pdb
sed -i "s/ HZ1 LYS/1HZ  LYS/g" chainA.pdb
sed -i "s/ HZ2 LYS/2HZ  LYS/g" chainA.pdb
sed -i "s/ HZ3 LYS/3HZ  LYS/g" chainA.pdb
awk '$1=="ATOM"&&$6=="1" {print $4}' chainA.pdb > pru
first=`head -n 1 pru`
sed -i "s/H1  ${first}/H   ${first}/g" chainA.pdb

cp VHH${i}.pdb chainB.pdb
sed -i "s/ HZ1 LYS/1HZ  LYS/g" chainB.pdb
sed -i "s/ HZ2 LYS/2HZ  LYS/g" chainB.pdb
sed -i "s/ HZ3 LYS/3HZ  LYS/g" chainB.pdb
awk '$1=="ATOM"&&$6=="1" {print $4}' chainB.pdb > pru
first=`head -n 1 pru`
sed -i "s/H1  ${first}/H   ${first}/g" chainB.pdb
${here0}/sf/FireDock/buildFireDockParams.pl chainA.pdb chainB.pdb B B Default fake.trans out 0 0 0.8 0 FireDock.par
${here0}/sf/FireDock/runFireDock.pl FireDock.par > FireDock.log
tail -n1 out.unref | awk 'BEGIN{FS="|"} {print $6}' >> FireDock.ene
done
block FireDock.ene
cp FireDock.ene ../FireDock${MM}.ene
scfired=0
scfired=` awk '{sum += $1} END {print sum/NR}' FireDock.ene`
efired=`awk '{print $1}' temp`
dfired=`awk '{print $2}' temp`
#sdfired=`awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' FireDock.ene`
dif=`echo "${scfired} ${scfired0}" | awk '{print $1 - $2}'`
contfired=`echo "$dif <= 0.0" | bc -l`
echo -e "${scfired} \n${efired} \n${dfired}" >> tempsf
Enew=${scfired}
fi

if [ ${SF} = vina ]; then
##### CALCULATING BINDING ENERGIES WITH VINA #####
# VINA scoring ######################
rm VINAscore
touch VINAscore
cp ../configure.txt .
for i in `seq 0 ${files}`
do
  rm vina.log
  rm -f *.pdbqt
  cp VHH${i}.pdb ligand.pdb
  cp HER2${i}.pdb receptor.pdb
  python ${UTILITIES}/prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
  python ${UTILITIES}/prepare_ligand4.py -l ligand.pdb -Z -o ligand.pdbqt
  nice -19 ${VINA_FOLD}/vina --score_only --config configure.txt &> vina.log
  af=`grep "Affinity" vina.log |awk '{print $2}'`
  echo ${af} >>VINAscore
done
block VINAscore
scvina=0
scvina=`awk '{sum += $1} END {print sum/NR}' VINAscore`
evina=`awk '{print $1}' temp`
dvina=`awk '{print $2}' temp`
#sdvina=`awk '{x+=$2;y+=$2^2}END{print sqrt(y/NR-(x/NR)^2)}' iradsc`
dif=`echo "${scirad} ${scirad0}" | awk '{print $1 - $2}'`
contvina=`echo "${dif} <= 0.0" | bc -l`
cp VINAscore ../VINAsc${MM}
Enew=${scvina}
echo -e "${scvina} \n${evina} \n${dvina}" >> tempsf
fi

done

# column to row
awk 'BEGIN { ORS = " " } { print }' tempsf >> scorings.out
echo >> scorings.out

cd ..

}

######### END of FUNCTIONS ########################

############ MAIN PROGRAM #########################
######### DO NOT TOUCH FROM HERE DOWN  ############

#--Function-Set-Up
# environment variable definitions
setup

#inital coordinates
if [ ${rest} -eq 0 ]; then
initial_box
fi

### INITIAL STRUCTURE READY FOR RUN!!! RUN!!! ######
######################################################################


########################################################################
####### SECTION 1. INITIAL CALCULATIONS BEFORE RE-MONTECARLO SIMULATION
if [ ${rest} -eq 0 ]; then

#NPT OF INITIAL COMPLEX (OPTIONAL)
if [ ${MD0} = yes ]; then
$grompp -f md-NPT.mdp -c Start.gro -p topol.top -n fixbs2.ndx -o npt-start.tpr -maxwarn 5
${mdrun} -deffnm npt-start
#${mdrun} -v -deffnm npt-start -ntomp 8
#mpirun -x OMP_NUM_THREADS=10 -npernode 2 ${BIN_DIR}/mdrun_4_7_gpu_s -deffnm npt-start -v -pin on
fi

# Export stage ####################
echo -e  "chA \n 0" | $trjconv -f npt-start.gro -s npt-start.tpr -n fixbs2.ndx -pbc mol -ur compact -center -o npt-pbc.gro  # Makes molecules whole
echo -e  "chA \n non-Water" | $trjconv -f npt-start.xtc -s npt-start.tpr -n fixbs2.ndx -pbc mol -ur compact -center -b ${bprint} -o npt-pbc.xtc      # make molecules all  
rm \#*
echo -e "chA"| $trjconv -f npt-pbc.gro -s npt-start.tpr -n fixbs2.ndx -o target$MM.pdb # export binder
#for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/"$i  "/"$i A"/g target$MM.pdb; done
echo -e "chB"     | $trjconv -f npt-pbc.gro -s npt-start.tpr -n fixbs2.ndx -o binder$MM.pdb        # export target
#for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/"$i  "/"$i B"/g binder$MM.pdb; done
echo -e "${typesol}" | $trjconv -f npt-pbc.gro -s npt-start.tpr -n fixbs2.ndx -o solvent$MM.pdb    # export solvent

# creating system pdb file
(grep ATOM target$MM.pdb
          echo "TER" 
          grep ATOM binder$MM.pdb) > complex.pdb
 head -n4 solvent$MM.pdb  > primo
 nlines=$(wc -l solvent$MM.pdb | awk '{print $1-4}')
 tail -n ${nlines} solvent$MM.pdb > last
 cat primo complex.pdb last > system$MM.pdb

# CREATING THE FOLDER OF SCORING ########################################################
rm -rf scoring
mkdir scoring
cp block-mig.sh scoring/.
echo "Step " > tempsf
for ((isf=0;isf<${Nsf};isf++)); do
    SF=${sf_list[${isf}]}
    echo -e "sc${SF} \nErr${SF} \nDiff${SF} " >> tempsf
done
# column to row
awk 'BEGIN { ORS = " " } { print }' tempsf > scoring/scorings.out
echo >> scoring/scorings.out

# SCORING THE BINDING
for ((isf=0;isf<${Nsf};isf++)); do
    SF=${sf_list[${isf}]}
if [ ${SF} = hadd ]; then
    cp scoring0.inp scoring.inp
fi
done
cp npt-start.tpr md-npt.tpr
scoring
Eprev=${Enew}


##### 1F. PREPARING FILES FOR LOOP #####
## 1F.1 init parallel step by copying inside its folders
for ((it=0;it<${Nrep};it++)); do
  rundir=Teff_${Trep[${it}]}
  rm -rf ${rundir}
  mkdir ${rundir}
  cp ./* ${rundir}
  cp -rf scoring ${rundir}/.
  cd ${rundir}
  
  # Creating the BindingScore file
  echo "Step " > tempsf
  echo "${MM} " > tempsf2
  for ((isf=0;isf<${Nsf};isf++)); do
      SF=${sf_list[${isf}]}
      echo "Score${SF} " >> tempsf
      tt=sc${SF}
      echo "${!tt} " >> tempsf2
  done
  # column to row
  awk 'BEGIN { ORS = " " } { print }' tempsf > Edt_${Trep[${it}]}.dat
  echo >> Edt_${Trep[${it}]}.dat
  awk 'BEGIN { ORS = " " } { print }' tempsf2 >> Edt_${Trep[${it}]}.dat
  echo >> Edt_${Trep[${it}]}.dat
  #echo "${MM}  ${Enew}" > Edt_${Trep[${it}]}.dat
  #echo "${MM}  ${scirad} ${scpisa} ${scbluues} ${schadd} ${scbach} ${scfired}" > Edt_${Trep[${it}]}.dat
  cd ..
done

for ((it=0;it<${Nrep};it++)); do
  Teff=${Trep[${it}]}
  rundir=Teff_${Trep[${it}]}
  cd ${rundir}

  # saving geometries of step 0 (before entering in the loop)
  cp binder${MM}.pdb evol_binder${MM}.pdb 
  cp target${MM}.pdb evol_target${MM}.pdb
  cp solvent${MM}.pdb evol_solvent${MM}.pdb
  tail -n1 npt-start.gro >  box${MM}.txt
  
  tail -n1 npt-start.gro >  box_dimension_to_use.txt
  cp -f solvent${MM}.pdb box_solvent_to_use.pdb
  cp -f binder${MM}.pdb kk.pdb
  cp -f target${MM}.pdb substrate.pdb

  >output_${Teff}
  cd ..
done
# starting step
M1=1

fi  
####### SECTION 1 FINISHED ########


######## SECTION 1.B RESTARTING SIMULATION ###############
if [ ${rest} -eq 1 ]; then
    echo "RESTART"
    for ((iit=0;iit<${Nrep};iit++)); do
    rundir=Teff_${Trep[${iit}]}
    cd ${rundir}
    
    M0=`ls -lrt evol_binder*.pdb | tail -n1 | awk '{print $9}' | sed "s/evol_binder//;s/.pdb//"`
    (( M1=M0+1 ))
    cp -f evol_target${M0}.pdb substrate.pdb
    cp -f evol_binder${M0}.pdb kk.pdb
    cp -f evol_solvent${M0}.pdb box_solvent_to_use.pdb
    cp -f  box${M0}.txt  box_dimension_to_use.txt
    Eprev=`tail -1 Edt_${Trep[${iit}]}.dat | awk '{print $2}'`
    echo "${rundir}: Last step $M0 with score ${Eprev}"

    # restarting clean    
    rm -f \#*
    rm -f nvt.*
    rm -f primo last
    rm -f complex.pdb
    rm -f min.*
    rm -f system.*
    rm -f mdout.*
    rm -f ./#mdout.*
    rm -f mod
    rm -f md-npt.*

    cd ..
done
fi
##########################################################

##########################################################
####### SECTION 2. MUTATION LOOP  #######

# beginning iterative process 
for MM in `seq ${M1} ${NMT}`; do
# replica exchange process loop
  for ((iit=0;iit<${Nrep};iit++)); do
    Teff=${Trep[${iit}]}
    rundir=Teff_${Trep[${iit}]}
    rm -f DONE_${Teff}
    cd ${rundir}
    echo "step $MM folder ${rundir}"

  ##### 2A. MUTATION ALGORITHM #####

## 2A.1 GET SEQUENCE FROM PDB FILE 
    for (( i=0; i<${NA};i++ )); do
	j=`awk -v a=${residue_list[$i]} '{if ($6 ==a) print $4}' kk.pdb|tail -n1`
	sequence[${i}]=${j}
    done
    cp kk.pdb kk0.pdb
    # writing the sequence
    echo "SEQ OLD: "${sequence[*]} >>output_${Teff}

## 2A.2 CHOOSE THE RESIDUE NUMBER
    rn=`echo ${RANDOM} | awk '{print $1/32767.}'`
    ia=`echo ${rn} ${NA} | awk '{printf "%i\n",$1*$2}'`
    echo ${ia}"/" ${NA} >>output_${Teff}
    ndum=${residue_list[${ia}]}
## 2A.3 CHOOSE THE AMINO ACID 
    rn=`echo ${RANDOM} | awk '{print $1/32767.}'`
    na=`echo ${rn} ${NAA} | awk '{printf "%i\n",$1*$2}'`
    echo ${na}"/" ${NAA} >>output_${Teff}
    # writing the mutation chosen
    echo "${ia} ${na}" > mutation
    echo "MUTATION: "${sequence[${ia}]} ${ndum} ${aalistL[${na}]} >>output_${Teff}


##### 2B. SUBSTITUTION OF AMINO ACID #####

## 2B.1 PERFORM MUTATION IN PDB BY SUBSTITUTION OF THE RESIDUE NAME
    resid=${ndum}
    oldrestype=${sequence[${ia}]}
    newrestype=${aalistL[${na}]}
    rm -f kk1.pdb
    cp kk.pdb kk1.pdb
    # changing the name of the residue
    sed -i "s/${oldrestype} B ${resid}/${newrestype} B ${resid}/g" kk1.pdb
    sed -i "s/${oldrestype} B  ${resid}/${newrestype} B  ${resid}/g" kk1.pdb
    sed -i "s/${oldrestype} B   ${resid}/${newrestype} B   ${resid}/g" kk1.pdb
    # delete not backbone atoms of the mutated residue
    cat kk1.pdb|awk -v x=${resid} '{if ($6 != x) {print} else if (($3 == "O")||($3 == "N")||($3 == "CA")||($3 == "C")) {print}}' > pdb_mutated.pdb

## 2B.2 RECONSTRUCTION OF MUTATED RESIDUE USING SCWRL4

# set up the binder file format
    awk '$1=="ATOM"' pdb_mutated.pdb | awk 'BEGIN {FS=""} {if($14!="H")print}' > tt.pdb
    cp -f tt.pdb pdb_mutated.pdb
    awk '$1=="ATOM"' pdb_mutated.pdb | awk 'BEGIN {FS=""} {if($13!="H")print}' > tt.pdb
    cp -f tt.pdb pdb_mutated.pdb

    sed -i "s/OC1 ${LAB}/O   ${LAB}/g" pdb_mutated.pdb
    sed -i "s/OC2 ${LAB}/OXT ${LAB}/g" pdb_mutated.pdb
    sed -i "s/CD  ILE/CD1 ILE/g" pdb_mutated.pdb
    mv pdb_mutated.pdb kk.pdb


# set up the target file format
    #awk '$1=="ATOM"' substrate.pdb | awk 'BEGIN {FS=""} {if($14!="H")print}' > tt.pdb
    #cp -f tt.pdb substrate.pdb
    #awk '$1=="ATOM"' substrate.pdb | awk 'BEGIN {FS=""} {if($13!="H")print}' > tt.pdb
    #mv tt.pdb substrate.pdb

    #sed -i "s/OC1 ${LAA}/O   ${LAA}/g" substrate.pdb
    #sed -i "s/OC2 ${LAA}/OXT ${LAA}/g" substrate.pdb
    #sed -i "s/CD  ILE/CD1 ILE/g" substrate.pdb


    # SCWRL reconstruction of side chain ###########################
    nini=`awk -v a=${resid} '{if ($6 ==a && $3 =="N") print NR}' kk.pdb`
    (( nini=nini-1 ))
    head -n${nini} kk.pdb > start
    nfin=`awk -v a=${resid} '{if ($6 ==a && $3 =="O") print NR}' kk.pdb`
    tot=`wc -l kk.pdb | awk '{print $1}'`
    (( fin=tot-nfin ))
    tail -n${fin} kk.pdb > fin

    # creation of sequence file
    rm seq seq2
    grep " N " kk.pdb|awk '{print $4}' > seq
    eo=`grep " N " kk.pdb|head -1|awk '{print $6}'`
    compress seq > seq2
    fila=$((resid-eo+1))
    sed "${fila}s/^./\U&/" seq2 > seq
    # column to row
    awk 'BEGIN { ORS = " " } { print }' seq > seq2

    # executing scwrl
    $scrwl -i kk.pdb -f substrate.pdb -o output.pdb -s seq2 -h > logscrwl

    # extracting the geometry of mutant residue
    awk -v a=${resid} '{if ($6 ==a) print $0}' output.pdb > mut
    cat start mut fin > pdb_repaired.pdb
    rm start mut fin

########################################

## 2B.3 GROMACS FILE SET UP

# set up the pdb file format for gromacs. 1. Delete hydrogens. 2. Change atom names.
    awk '$1=="ATOM"' pdb_repaired.pdb | awk 'BEGIN {FS=""} {if($14!="H")print}' > tt.pdb
    cp -f tt.pdb pdb_repaired.pdb
    awk '$1=="ATOM"' pdb_repaired.pdb | awk 'BEGIN {FS=""} {if($13!="H")print}' > tt.pdb
    mv tt.pdb pdb_repaired.pdb
    echo -e "6 \n 6" | ${pdb2gmx} -f pdb_repaired.pdb -o pdb_repaired.gro -p binder.top

    # creation of topology file of the binder
    sed -i "/forcefield/d" binder.top
    sed -i '/\[ system \]/,$d' binder.top
    mv  binder.top topol_Protein_chain_B.itp

#   final transformation of pdb of binder and creation of complex pdb
#    sed -i "s/O   ${LAB}/OC1 ${LAB}/g" pdb_repaired.pdb
#    sed -i "s/OXT ${LAB}/OC2 ${LAB}/g" pdb_repaired.pdb
    ${editconf} -f pdb_repaired.gro -o pdb_repaired.pdb
    for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/"$i  "/"$i B"/g pdb_repaired.pdb; done

    (grep ATOM substrate.pdb
	echo "TER" 
        grep ATOM pdb_repaired.pdb
        echo "END") > complex.pdb

# creation of gro file of complex
    ${editconf} -f complex.pdb -o comp.gro


##### 2C. MINIMIZATION OF THE GEOMETRY AFTER MUTATION. ########
# BASED ON A PARTIAL MINIMIZATION ALGORITHM

## 2C.1 Making Index file to define the groups. We define the group scmut as complex target-binder except sidechain of mutated residue
#$grompp -f minim.mdp -c comp.gro -p complex.top -o em.tpr -maxwarn 5 #I need of this just to make the molecule whole             
#echo -e "1" | $trjconv -f comp.gro -s em.tpr -o index.pdb

    Nmut=${resid} 
    rm -f ssmut.txt
    rm -f ssmut.ndx
    echo "chain B" >ssmut.txt 
    echo "\"SideChain\" & \"chB\" & r ${Nmut}" >>ssmut.txt #group sidechain of mutated residue
    echo "\"System\" &! \"SideChain_&_chB_&_r_${Nmut}\"" >>ssmut.txt #group rest of the system
    echo "q" >>ssmut.txt
    #cp ssmut.txt ssmut${MM}.txt
    ${make_ndx} -f complex.pdb -o ssmut.ndx < ssmut.txt &> indexlog
    ngroup=`grep "System_&_\!SideChain_&_chB_&_r_${Nmut}" indexlog | awk '{print $1}'`
    echo -e "name ${ngroup} scmut \n q" | ${make_ndx} -n ssmut.ndx -o ssmut.ndx 

    ## 2C.2 First minimization. Only side chain of mutated residue
#    export OMP_NUM_THREADS=12
    sed -i '$ d' comp.gro
    echo "   20.0   20.0   20.0" > dim-min.txt
    cat comp.gro dim-min.txt > tt
    mv tt comp.gro
    ${grompp} -v -f minim_scmut.mdp -c comp.gro -p complex.top -o min.tpr -n ssmut.ndx -maxwarn 5
    ${mdrun} -deffnm min
    cp -f min.gro comp.gro

    echo "0 \n" | ${trjconv} -f comp.gro -s min.tpr -o complex.pdb     #complex
    rm -f min.*

    #function mutated box # Creating the water box.
    mutated_box

########################

    # saving complex.pdb file
    cp solvated_and_NOToverlapping.pdb system${MM}.pdb

    cp system.gro system${MM}.gro
    cp system.top system${MM}.top

#    $editconf -f system.gro -o system.box.gro -bt cubic -d $boxsize

# Our GROMACS file is ready!!!! ################################

## 2C.4 Last global minimization
    ${grompp} -f minim.mdp -c system.gro -p system.top -o min-all.tpr -maxwarn 5
    ${mdrun} -deffnm min-all
    mv min-all.gro system.gro
    rm index.ndx
    echo "0 \n" | ${trjconv} -f system.gro -s min-all.tpr -o min.pdb     #index
    ${make_ndx} -f min.pdb -o fixbs2.ndx<select.txt
    rm -f min-all.*

# Our file is already super minimized :)

# Molecular Dynamics simulation
#1. NVT equilibration
rm nvt.gro
${grompp} -f md-NVT.mdp -c system.gro -p system.top -n fixbs2.ndx -o nvt.tpr -maxwarn 5
${mdrun} -deffnm nvt

#2. NPT
${grompp} -f md-NPT.mdp -t nvt.cpt -c nvt.gro -p system.top -n fixbs2.ndx -o md-npt.tpr -maxwarn 5

cd ..
done

# call function gromacs: NPT MD run.
gromacs

for ((it=0;it<${Nrep};it++)); do
    #(
	Teff=${Trep[${it}]}
	rundir=Teff_${Trep[${it}]}
	rm -f DONE_${Teff}
	cd ${rundir}
	echo "step $MM folder ${rundir}"
	
	#tail -n1 md-npt.gro >  box_dimension_to_use.txt

  # clean
  rm -f npt-pbc.gro npt-pbc.xtc
  # loading binding scores of previous step	
  Eprev=`tail -1 Edt_${Trep[${it}]}.dat | awk '{print $2}'`
  Enew=0.0
  echo "Eprev = ${Eprev} ; Enew = ${Enew} ; Nsf = ${Nsf}"
  for ((isf=0;isf<${Nsf};isf++)); do
      SF=${sf_list[${isf}]}
      echo "isf=${isf} SF=${SF}"
      declare "sc${SF}=0"
      (( icol=isf+2 ))
      declare "sc${SF}0=`tail -1 Edt_${Teff}.dat | awk -v a=${icol} '{print $a}'`"
  done

  # IF MD simulation finished correctly
  if [ -e md-npt.gro ]; then
      echo -e  "chA \n 0" | $trjconv -f md-npt.gro -s md-npt.tpr -n fixbs2.ndx -pbc mol -ur compact -center -o npt-pbc.gro  # Makes molecules whole
      echo -e  "chA \n non-Water" | $trjconv -f md-npt.xtc -s md-npt.tpr -n fixbs2.ndx -pbc mol -ur compact -center -b ${bprint} -o npt-pbc_${MM}.xtc      # make molecules all  
      cp npt-pbc_${MM}.xtc npt-pbc.xtc
      echo -e "chA"| $trjconv -f npt-pbc.gro -s md-npt.tpr -n fixbs2.ndx -o target$MM.pdb #protein
#for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/"$i  "/"$i A"/g target$MM.pdb; done
      echo -e "chB"     | $trjconv -f npt-pbc.gro -s md-npt.tpr -n fixbs2.ndx -o binder$MM.pdb        #VHH
#for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/"$i  "/"$i B"/g binder$MM.pdb; done
      echo -e "${typesol}" | $trjconv -f npt-pbc.gro -s md-npt.tpr -n fixbs2.ndx -o solvent$MM.pdb    #solvent
      
# creating system pdb file
(grep ATOM target$MM.pdb
          echo "TER" 
          grep ATOM binder$MM.pdb) > complex.pdb
 head -n4 solvent$MM.pdb  > primo
 nlines=$(wc -l solvent$MM.pdb | awk '{print $1-4}')
 tail -n ${nlines} solvent$MM.pdb > last
 cat primo complex.pdb last > system$MM.pdb


# SCORING THE BINDING ######################
# preparing .inp file for haddock scoring. Evaluating histidines in VHHs
if [ ${SF} = haddock ]; then
>hieB
>hidB
>inpB
grep "N   HIS B" binder$MM.pdb | awk '{print $6}' >> hieB
grep "N   HIE B" binder$MM.pdb | awk '{print $6}' >> hieB
grep "N   HID B" binder$MM.pdb | awk '{print $6}' >> hidB
((c=0))
for x in `cat hieB`
do
    ((c++))
    echo "evaluate (\$toppar.hise_resid_2_${c}=${x})">>inpB
done
((c=0))
for x in `cat hidB`
do
    ((c++))
    echo "evaluate (\$toppar.hisd_resid_2_${c}=${x})">>inpB
done
cat inpA inpB inpC > scoring.inp
fi
#########################
scoring
#########################
# if MD simulation crushed
  else
      contall=0
      echo -e  "non-Water"      | $trjconv -f md-npt.xtc -s md-npt.tpr -n fixbs2.ndx -pbc whole       -o npt-pbc_${MM}.xtc      # make molecules all  
  #mv md-npt.xtc npt-pbc.xtc
  #tail -n1 system.gro >  box_dimension_to_use.txt
  fi

# save  
mv md-npt.tpr md-npt_${MM}.tpr
mv md-npt.gro md-npt_${MM}.gro
mv md-npt.log md-npt_${MM}.log


##### 2G. METROPOLIS ALGORITHM FOR THE MUTATION ########   

# Applying the metropolis criterion.
if [ ${ACCMET} = MC ]; then
    rn=`echo ${RANDOM} | awk '{print $1/32767.}'`
    acc=`echo ${Eprev} ${Enew} ${Teff} ${rn} | awk '{acc=0;if(exp(($1-$2)/$3)>$4)acc=1;print acc}'`
fi

# Applying the consensus criterion.
if [ ${ACCMET} = consensus ]; then	
    contall=0	
    for ((isf=0;isf<${Nsf};isf++)); do
	SF=${sf_list[${isf}]}
	tt=cont${SF}
	contall=`echo ${contall}+${!tt} | bc -l`
    done
    if [ ${contall} -ge ${Trep[${it}]} ]; then
       acc=1
    fi
fi
       
       # Recover the mutation data    
ia=`cat mutation | awk '{print $1}'`
na=`cat mutation | awk '{print $2}'`
ndum=${residue_list[${ia}]}
(( NO=MM-1 ))
for (( i=0; i<${NA};i++ )); do
    j=`awk -v a=${residue_list[$i]} '{if ($6 ==a) print $4}' kk0.pdb|tail -n1`
    sequence[${i}]=${j}
done
rm kk0.pdb

if [ ${acc} -eq 1 ]; then
	# Condition accepted. The new sequence is stored.
	cp -f solvent${MM}.pdb box_solvent_to_use.pdb
	tail -n1 md-npt_${MM}.gro >  box_dimension_to_use.txt
	cp -f binder${MM}.pdb kk.pdb
	cp -f target${MM}.pdb substrate.pdb
	tail -n1 md-npt_${MM}.gro >  box${MM}.txt
	cp binder${MM}.pdb evol_binder${MM}.pdb 
	cp target${MM}.pdb evol_target${MM}.pdb
	cp solvent${MM}.pdb evol_solvent${MM}.pdb
	echo "mut:${MM} ${sequence[${ia}]} ${ndum} => ${aalistL[${na}]} accepted ">>output_${Teff}
	if [ ${ACCMET} = MC ]; then
	echo "prevE:${Eprev} newE:${Enew} accepted ">>output_${Teff}
	Eprev=${Enew}
	fi
	if [ ${ACCMET} = consensus ]; then	
	    echo "consensus criterion ACCEPTED with ${contall}/${Nsf} positive BScores, passed the threshold ${Trep}">>output_${Teff}
	fi
	
	for ((isf=0;isf<${Nsf};isf++)); do
	    SF=${sf_list[${isf}]}
	    tt=sc${SF}
	    declare "sc${SF}0=${!tt}"
	done
else 
	# Condition not accepted. We come back to the old sequence.
	(( NO=MM-1 ))
	echo "mut:${MM} ${sequence[${ia}]} ${ndum} => ${aalistL[${na}]} rejected ">>output_${Teff}
	if [ ${ACCMET} = MC ]; then
	echo "prevE:${Eprev} newE:${Enew} REJECTED ">>output_${Teff}
	fi
	if [ ${ACCMET} = consensus ]; then	
	echo "consensus criterion REJECTED with ${contall}/${Nsf} positive BScores, NOT passed the threshold ${Trep}">>output_${Teff}
	fi
	cp -f evol_target${NO}.pdb substrate.pdb
	cp -f evol_binder${NO}.pdb kk.pdb
	cp -f evol_target${NO}.pdb evol_target${MM}.pdb
	cp -f evol_binder${NO}.pdb evol_binder${MM}.pdb
	cp -f evol_solvent${NO}.pdb evol_solvent${MM}.pdb
	cp -f box${NO}.txt box${MM}.txt
	
    fi

# WRITING THE FINAL BINDING SCORE
  echo "${MM} " > tempsf2
  for ((isf=0;isf<${Nsf};isf++)); do
      SF=${sf_list[${isf}]}
      tt=sc${SF}0
      echo "${!tt} " >> tempsf2
  done
  # column to row
  awk 'BEGIN { ORS = " " } { print }' tempsf2 >> Edt_${Trep[${it}]}.dat
  echo >> Edt_${Trep[${it}]}.dat

# clean process ###############################
    rm -f \#* 
    rm -f nvt.*
    rm -f primo last
    rm -f complex.pdb
    rm -f min.*
    rm -f system.* 
    rm -f mdout.*
    rm -f ./#mdout.*
    rm -f mod
    rm -f md-npt.*
    
    cd ..
# FINISHING THE LOOP
    #touch DONE_${Teff}
   # ) &
  done
#----------------checked--------------------
# syncro all parallel simulations
#if [ ${Nrep} -gt 1 ]; then
#  ((c=0))
#  while [ $c -lt ${Nrep} ];do
#    sleep 1
#    ((c=0))
#    for ((iit=0;iit<${Nrep};iit++)); do
#      Teff=${Trep[${iit}]}
#      if [ -e DONE_$Teff ]; then
#        ((c++))
#      fi
#    done
#    if [ -f EXIT ]; then
#      exit
#    fi
#  done
#fi

# restarting variable acceptance to 0
  acc=0
# replica exchange process
 ####################################################### 
  # now attempt swapping the state in two temperatures ##
  ######################################################
  #### Selecting 2 temperatures ####
  ##################################
 
  if [[ $REXCH -eq 1 ]];then

  index1=0
  index2=0
  FLOOR=0
  while [ "$index1" -le $FLOOR ] && [ "$index2" -le $FLOOR ]; do
    index1=$RANDOM
    index2=$RANDOM
    let "index1 %= ${Nrep}"  # Scales $index down within ${Nrep} (the number of temperatures).
    let "index2 %= ${Nrep}"
    if [ "$index1" -eq $index2 ]; then
       index1=0
       index2=0
    fi
  done
  Teff1=${Trep[${index1}]}
  rundir1=Teff_${Teff1}
  Teff2=${Trep[${index2}]}
  rundir2=Teff_${Teff2}



  ####################################################
  #### Taking final energies for replica exchange #### 
  ####################################################

  iex1=`tail -n1 ${rundir1}/Edt_${Teff1}.dat | awk '{print $2}'`
  iex2=`tail -n1 ${rundir2}/Edt_${Teff2}.dat | awk '{print $2}'`

  rn=`echo $RANDOM | awk '{print $1/32767.}'`
  compare_result=` echo $rn $iex1 $iex2 $Teff1 $Teff2 | awk '{r=0;factor=exp(($2-$3)/(1/$4-1/$5));if($1<factor)r=1;print r}'`
  
  if [[ $compare_result -eq 1 ]];then
# exchange binder (kk.pdb), target (substrate.pdb) and energies (EdT.dat)
    cp -f ${rundir1}/kk.pdb ${rundir1}/dumm.pdb
    cp -f ${rundir2}/kk.pdb ${rundir1}/kk.pdb
    mv -f ${rundir1}/dumm.pdb ${rundir2}/kk.pdb
    cp -f ${rundir1}/substrate.pdb ${rundir1}/dumm.pdb
    cp -f ${rundir2}/substrate.pdb ${rundir1}/substrate.pdb
    mv -f ${rundir1}/dumm.pdb ${rundir2}/substrate.pdb
    cp -f ${rundir1}/box_solvent_to_use.pdb ${rundir1}/dumm.pdb
    cp -f ${rundir2}/box_solvent_to_use.pdb ${rundir1}/box_solvent_to_use.pdb
    mv -f ${rundir1}/dumm.pdb ${rundir2}/box_solvent_to_use.pdb
    cp -f ${rundir1}/box_dimension_to_use.txt ${rundir1}/dumm.txt
    cp -f ${rundir2}/box_dimension_to_use.txt ${rundir1}/box_dimension_to_use.txt
    mv -f ${rundir1}/dumm.txt ${rundir2}/box_dimension_to_use.txt

    tail -n1 ${rundir1}/Edt_${Teff1}.dat > ${rundir1}/dumm.dat
    tail -n1 ${rundir2}/Edt_${Teff2}.dat > ${rundir2}/dumm.dat
    sed -i "$ d" ${rundir1}/Edt_${Teff1}.dat
    sed -i "$ d" ${rundir2}/Edt_${Teff2}.dat
    cat ${rundir1}/dumm.dat >> ${rundir2}/Edt_${Teff2}.dat
    cat ${rundir2}/dumm.dat >> ${rundir1}/Edt_${Teff1}.dat
    
# exchanging the permanent files of the system 
    cp -f ${rundir1}/kk.pdb ${rundir1}/evol_binder${MM}.pdb
    cp -f ${rundir1}/substrate.pdb ${rundir1}/evol_target${MM}.pdb
    cp -f ${rundir2}/kk.pdb ${rundir2}/evol_binder${MM}.pdb
    cp -f ${rundir2}/substrate.pdb ${rundir2}/evol_target${MM}.pdb
    cp ${rundir1}/evol_solvent${MM}.pdb ${rundir2}/evol_solvent${MM}f.pdb
    cp ${rundir2}/evol_solvent${MM}.pdb ${rundir1}/evol_solvent${MM}f.pdb
    mv ${rundir1}/evol_solvent${MM}f.pdb ${rundir1}/evol_solvent${MM}.pdb
    mv ${rundir2}/evol_solvent${MM}f.pdb ${rundir2}/evol_solvent${MM}.pdb
    cp ${rundir1}/box${MM}.pdb ${rundir2}/box${MM}f.pdb
    cp ${rundir2}/box${MM}.pdb ${rundir1}/box${MM}f.pdb
    mv ${rundir1}/box${MM}f.pdb ${rundir1}/box${MM}.pdb
    mv ${rundir2}/box${MM}f.pdb ${rundir2}/box${MM}.pdb

  fi
  echo "MUT=${MM}, i1=${index1}, T1=${Teff1}, E1=${iex1}; i2=${index2}, T2=${Teff2}, E2=${iex2}, Random=${rn}, Accept=${compare_result}" >>"REMC.dat"
  if [ -f EXIT ]; then
    exit
  fi
  fi
done
