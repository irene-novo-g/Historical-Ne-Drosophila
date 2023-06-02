#!/bin/bash
#$ -cwd

rm script_natBOTchrom_neutral.sh.*

#script_natBOTchrom_neutral.sh <d> <LAMB> <NCRO> <BETA> <AVEs> <AVEh>
#To be used only in a machine with /state/partition1 directory

#Check number of arguments
if [ $# -ne 6 ]  
then
	echo "Usage: $0 <d> <LAMB> <NCRO> <BETA> <AVEs> <AVEh>" 
	exit 1
fi

#Set arguments
d=$1
LAMB=$2
NCRO=$3
BETA=$4
AVEs=$5
AVEh=$6
 

#Working directory
WDIR=$PWD 

#Scratch directory
mkdir -p /state/partition1/noeliaNAT$d/$SLURM_JOBID/

#Copy all files in scratch directory
cp seedfile /state/partition1/noeliaNAT$d/$SLURM_JOBID/
cp naturalvBOTchrom_neutral /state/partition1/noeliaNAT$d/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd /state/partition1/noeliaNAT$d/$SLURM_JOBID

START=$(date +%s)
time ./naturalvBOTchrom_neutral>>out<<@
0
-99
10000	N
0	PS(99=random)
0.5	Lenght of chromosomes in Morgans (99=free)
2	Number of Chromosomes (min 1, max 50)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_a
0.0	Lambda_L
0.0	absolute effect of lethal (QT): normal (aL,aL)
0	random proportion(0) or large mutants (1)
1.0	Psi
$BETA	beta_s
$BETA	beta_a
$AVEs	ave |s|
$AVEs	ave |a|
0.0	PP_s
0.0	PP_a
2	dom model (0=cnt; 1:Deng, 2:CK94 gamma)
$AVEh	h_s (mod 0), k_s (mod 1)
$AVEh	ave h_s (mod 2)
$AVEh	h_a (mod 0), k_a (mod 1)
$AVEh	ave h_a (mod 2)
99	rho (99:a=s)
1	neutral
0	Vs
1	multi(1), add(2)
10000	generations
5000	gen/block
0	GENBOT
@

cat popfile >> POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh
cat datafile >> DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh
cat genfile >> GENFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "naturalv took 		$DIFF seconds" >> timefile

#Copy output files to main directory
cp -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/
cp -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/
cp -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/GENFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/
cp -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/timefile $WDIR/
cp -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/dfilename $WDIR/

#Cleaning of scratch
rm -r /state/partition1/noeliaNAT$d/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
