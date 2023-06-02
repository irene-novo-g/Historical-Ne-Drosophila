#!/bin/bash
#$ -cwd

#script_estimNe.sh <d> <EXP> <MS> <mig>
#To be used only in a machine with /state/partition1 directory

#######################################################################################
#####                    OUTLINE -  variation in Nº bottles                       #####
#######################################################################################
#                                                                                     #
#                  >>START EXP                                                        #
#                  .                 ________ (1.1)                                   #
#                  .                /        *                                        #
#                  .               /                                                  #
# ________________ .              /__________ (1.2)                                   #
#  BASE.P :  (0.1)*|             /           *                                        #
#         :        |      (1.0) /                                                     #
#         :        |__________*/                                                      #
#         :        .                                                                  #
#         :        .   _________  (2)                                                 #
#         :        .  /         *                                                     #
#         :        . /                                                                #
#         :_______ ./                          0. None Experiment (mant. BP and lines)#
#  LINES  :  (0.2)*.                           1. Bottleneck-Expansion Experiment (A) #
#         :        .                           2. Line-Expansion Experiment (B)       #
#         :_______ .                           3. Crossing-Lines Experiment (C)       #
#         :  L1   |.____________  (3)                                                 #
#         :_______|.  MIXED L   *              * outputs: data.ped and data.map       #
#            L2    .                                                                  #
#                                                                                     #
#######################################################################################

rm script_estimNe.sh.*

#Check number of arguments
if [ $# -ne 4 ]  
then
	echo "Usage: $0 <d> <EXP> <MS> <mig>" 
	exit 1
fi

#Set arguments
d=$1
EXP=$2
MS=$3
mig=$4

#Variables
LAMB=0.5
NCRO=1000
CROX=0
AVEs=0.2
AVEh=0.283
#mig=0.5
LINES=100
SEED=1234
NP=10000

#################################### DIRECTORIES #####################################

#Working directory
WDIR=$PWD 

#Output directory

if [[ "$MS" == 0 ]] ; then
	mkdir -p $WDIR/estimNe_results/BPcircular.n100.exp$EXP.m$mig.L$LAMB.K$NCRO.s$AVEs.h$AVEh/
	DIR="estimNe_results/BPcircular.n100.exp$EXP.m$mig.L$LAMB.K$NCRO.s$AVEs.h$AVEh"
fi
if [[ "$MS" == 1 ]] ; then
	mkdir -p $WDIR/estimNe_results/BPoffspool.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh/
	DIR="estimNe_results/BPoffspool.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh"
fi
if [[ "$MS" == 2 ]] ; then
	mkdir -p $WDIR/estimNe_results/BPpanmixia.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh/
	DIR="estimNe_results/BPpanmixia.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh"
fi

#Scratch directory
mkdir -p /state/partition1/noeliaNe$d/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaNe$d"

############################ TRANSFER OF FILES TO SCRATCH ############################

#Copy all files in scratch directory
cp seedfile $SCDIR/$SLURM_JOBID/
cp estimNe $SCDIR/$SLURM_JOBID/
cp POPFILE_N$NP.L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/popfile
cp DATAFILE_N$NP.L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/datafile

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

############################# SIMULATION - estimNecm #################################

for ((i=1; i<=$LINES; i++))
do

seed=$(( $SEED + $i - 1 ))

START=$(date +%s)
time ./estimNe>>out<<@
0
$seed
$EXP	Experiment: mant(0.1; 0.2), A(1.0; 1.1; 1.2), B(2), C(3), ALL(99)
$MS	Mating system (0 circular, 1 offspring pool, 2 complete panmixia)
$mig	Migration between bottles
50	NIND per Bottle (max 80)
30	Number Bottles - Base Population (max 60)
209	Generations before experiments (+1)
83	Formation of lines (needed for expB, expC and mant.0.2)
1	Number Bottles - Lines (needed for expB, expC and mant.0.2) (max 60)
1	if experiment=0.1, maintain Pb bottles isolated after formation of lines (yes 0, no 1)
1	(expA)     Nº Bottles bottleneck (max 60) 
216	(expA)     Iniciate recovery at generation (+1) 
218	(expA)     Stop recovery at generation
221	(expA 1.2) Mantain pop.size until generation 
60	(expA 1.1) Nº Bottles expansion (max 60) 
220	(expA 1.1) Expansion until generation
224	(expA 1.1) Mantain pop.size until generation -
30	(expB)     Nº Bottles expansion (max 60) 
212	(expB)     Expansion until generation 
216	(expB)     Mantain pop.size until generation
218	(expC)     Mantain mixed line until generation
0.5	Lenght per chromosome in Morgans (99=free)
2	Number of Chromosomes (min 1, max 50)
$CROX	Chromosome X (no 0, yes 1)
$NCRO	NCRO (min 1, max 2000)
30	NLOCI (2-30)
$LAMB	Lambda (neutral model)
100	Sample size for data_files
1	Replicates
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "line$i  estimNe took 		$DIFF seconds" >> timefile

### TRANSFER OF FILES TO MAIN DIRECTORY 

if [[ "$EXP" != 0.2 && "$EXP" != 2 && "$EXP" != 3 ]] ; then
	mv BP.map BP.r$i.map; mv BP.ped BP.r$i.ped
	mv BPbot.map BPbot.r$i.map; mv BPbot.ped BPbot.r$i.ped
	mv BPbot_exp.map BPbot_exp.r$i.map; mv BPbot_exp.ped BPbot_exp.r$i.ped
	mv BPred.map BPred.r$i.map; mv BPred.ped BPred.r$i.ped
	mv genfile_BP.dat genfile_BP.r$i.dat

	mkdir -p $WDIR/$DIR/expA_BP
	cp -r $SCDIR/$SLURM_JOBID/genfile_BP* $WDIR/$DIR/expA_BP/
	cp -r $SCDIR/$SLURM_JOBID/BP* $WDIR/$DIR/expA_BP/
fi
if [[ "$EXP" == 99 || "$EXP" == 0.2 || "$EXP" == 2 ]] ; then
	mv LINE.map LINE.r$i.map; mv LINE.ped LINE.r$i.ped
	mv LINEexp.map LINEexp.r$i.map; mv LINEexp.ped LINEexp.r$i.ped
	mv genfile_Lines.dat genfile_Lines.r$i.dat

	mkdir -p $WDIR/$DIR/expB_Lines
	cp -r $SCDIR/$SLURM_JOBID/genfile_Lines* $WDIR/$DIR/expB_Lines/
	cp -r $SCDIR/$SLURM_JOBID/LINE* $WDIR/$DIR/expB_Lines/
fi
if [[ "$EXP" == 99 || "$EXP" == 3 ]] ; then
	mv CrossL.map CrossL.r$i.map
	mv CrossL.ped CrossL.r$i.ped
	mv genfile_Cross.Line.dat genfile_Cross.Line.r$i.dat

	mkdir -p $WDIR/$DIR/expC_CrossL
	cp -r $SCDIR/$SLURM_JOBID/genfile_Cross.L* $WDIR/$DIR/expC_CrossL/
	cp -r $SCDIR/$SLURM_JOBID/CrossL* $WDIR/$DIR/expC_CrossL/
fi
if [[ "$EXP" == 99 ]] ; then
	mv genfile_Fst_Lines.dat genfile_Fst_Lines.r$i.dat
	cp -r $SCDIR/$SLURM_JOBID/genfile_Fst_Lines* $WDIR/$DIR/expB_Lines/
fi

cp -r $SCDIR/$SLURM_JOBID/timefile $WDIR/$DIR/
done

cp -r $SCDIR/$SLURM_JOBID/summary_outline.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/dfilename.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/out $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/seedfile $WDIR

############################# CLEANING OF SCRATCH ###################################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*

