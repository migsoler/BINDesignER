# BINDesignER (now PARCE)
# (note: the new version of the code, called now PARCE, is located in https://github.com/PARCE-project) 

ALGORITHM FOR THE DESIGN OF PROTEIN BINDERS (PEPTIDES, NANOBODIES, ...) FOR THE MOLECULAR RECOGNITION OF PROTEIN TARGETS OR DRUGS.

YOU CAN FIND TWO EXAMPLES IN THIS FOLDER: PEPTIDE BINDER OF DRUG TARGET IN METHANOL / PEPTIDE BINDER OF PROTEIN TARGET IN WATER

Before starting the simulation, you must install and configure the external programs and prepare the external input files.

### EXTERNAL PROGRAMS ######
GROMACS. A priori, you can use any version of Gromacs since you must define the names of the gromacs executables at the beginning of the design algorithm. Moreover, the commands for executing the MD simulations should be adapted to your work station in the function "gromacs".

SCWRL4. You can download the program for free at http://dunbrack.fccc.edu/scwrl4/ and install it in the folder scwrl4.

Scoring functions. They are located inside of the sf folder. If you add other scoring functions, they should be located in this folder. Available scoring functions: Irad, Pie*Pisa, Bluues, Haddock, Bach6, FireDock, AutoDock Vina. Paths of scoring functions that must be configured: Haddock and AutoDock Vina.
In Haddock: edit the paths in rescoring_scripts/run_scoring.csh and follow the README file in the folder haddock (haddock/haddock2.1/haddock_configure.csh).
In Vina: edit the path in VINA/MGLTools-1.5.6/bin/mglenv.sh

### EXTERNAL INPUT FILES ######
- Start.pdb 
- Start.gro and topology files: topol-0.top, topol_Protein_chain_A-0.itp (target topology), topol_Protein_chain_B-0.itp (binder topology) (, solvent-0.itp if organic solvent) 
Important: Target = chain A ; Binder = chain B
Follow the topology templates in the provided examples if you do not know how to configure the gromacs topology files.  

- Minimization input files: minim_scmut.mdp, minim_overlap.mdp and minim.mdp, for the minimization of mutants. minim_scmut: minimization of mutated side-chain and the two attached residues using only the binding complex; minim_overlap: minimization of the mutated side chain and the surrounding solvent molecules; minim: global minimization.

- Input files for NVT and NPT MD simulations: md-NVT.mdp, md-NPT.mdp

- Inputs for haddock (optional):
scoring0.inp: initial haddock file. Be aware that you must define the type of histidine residues that you have in your system in scoring0.inp. HIE=hise, HID=hisd
inpA, inpC: input haddock file prepared by pieces to include new mutations 

-topology file of water molecules where the mass of H is modified to speed up the thermodynamics (optional, not active):
tip3p-m.itp

##############################


