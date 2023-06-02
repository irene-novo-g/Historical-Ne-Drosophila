#!/bin/bash
#$ -cwd

#script_GONE_scratch_Xchro.sh <d> <FOLDER> <FILE>
#To be used only in a machine with /state/partition1 directory

rm script_GONE_scratch_Xchro.sh.*

######################################################################################

#Check number of arguments
if [ $# -ne 3 ]  
then
	echo "Usage: $0 <d> <FOLDER> <FILE>" 
	exit 1
fi

### Set arguments
d=$1
FOLDER=$2  #BPcircular; BPoffspool; BPpanmixia
FILE=$3  #BP; BPred; BPbot; BPbot_exp; (LINE; LINEexp; CrossL) (map and ped files)

################################ SIMULATION VARIABLES ################################

LAMB=0.125
NCRO=250
AVEs=0.2
AVEh=0.283
mig=0.5
EXP=99
LINES=100
CROX=1

#################################### DIRECTORIES #####################################

#Working directory
WDIR=$PWD 

#Output directory

if [[ "$FOLDER" == BPcircular ]] ; then
mkdir -p $WDIR/estimNe_results/$FOLDER.n100.exp$EXP.m$mig.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX/GONE/$FOLDER.$FILE/
DIR="estimNe_results/$FOLDER.n100.exp$EXP.m$mig.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX"
fi
if [[ "$FOLDER" == BPoffspool ]] ; then
mkdir -p $WDIR/estimNe_results/$FOLDER.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX/GONE/$FOLDER.$FILE/
DIR="estimNe_results/$FOLDER.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX"
fi
if [[ "$FOLDER" == BPpanmixia ]] ; then
mkdir -p $WDIR/estimNe_results/$FOLDER.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX/GONE/$FOLDER.$FILE/
DIR="estimNe_results/$FOLDER.n100.exp$EXP.L$LAMB.K$NCRO.s$AVEs.h$AVEh.X$CROX"
fi


#Scratch directory
mkdir -p /state/partition1/noeliaGONE$d/PROGRAMMES/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaGONE$d"

############################ TRANSFER OF FILES TO SCRATCH ############################

#Copy all (GONE) files in scratch directory
cp $WDIR/GONE/PROGRAMMES/GONE $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/PROGRAMMES/GONEaverage $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/PROGRAMMES/GONEparallel.sh $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/PROGRAMMES/LD_SNP_REAL3 $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/PROGRAMMES/MANAGE_CHROMOSOMES2  $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/PROGRAMMES/SUMM_REP_CHROM3 $SCDIR/$SLURM_JOBID/PROGRAMMES/
cp $WDIR/GONE/INPUT_PARAMETERS_FILE_Xchro $SCDIR/$SLURM_JOBID/INPUT_PARAMETERS_FILE

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

#map and ped files are transferred inside the loop

############################### LOOP OF REPLICATES ###################################

for ((rep=1; rep<=$LINES; rep++))
do

if [[ "$FILE" == BP* ]] ; then
cp $WDIR/$DIR/expA_BP/$FILE.r$rep.map $SCDIR/$SLURM_JOBID/$FILE.map
cp $WDIR/$DIR/expA_BP/$FILE.r$rep.ped $SCDIR/$SLURM_JOBID/$FILE.ped
fi

if [[ "$FILE" == LINE* ]] ; then
cp $WDIR/$DIR/expB_Lines/$FILE.r$rep.map $SCDIR/$SLURM_JOBID/$FILE.map
cp $WDIR/$DIR/expB_Lines/$FILE.r$rep.ped $SCDIR/$SLURM_JOBID/$FILE.ped
fi

if [[ "$FILE" == CrossL ]] ; then
cp $WDIR/$DIR/expC_CrossL/$FILE.r$rep.map $SCDIR/$SLURM_JOBID/$FILE.map
cp $WDIR/$DIR/expC_CrossL/$FILE.r$rep.ped $SCDIR/$SLURM_JOBID/$FILE.ped
fi

# **************** #

### Take input parameters from file INPUT_PARAMETERS_FILE

source INPUT_PARAMETERS_FILE

###################### FILES NEEDED ########################

### data.ped
### data.map

### EXECUTABLES FILES NEEDED IN DIRECTORY PROGRAMMES:

### MANAGE_CHROMOSOMES2
### LD_SNP_REAL3
### SUMM_REP_CHROM3
### GONE (needs gcc/7.2.0)
### GONEaverages
### GONEparallel.sh

################### Remove previous output files ##################

if [ -f "OUTPUT_$FILE" ]
then
rm OUTPUT_$FILE
fi

if [ -f "Ne_$FILE" ]
then
rm Ne_$FILE
fi

################### Create temporary directory ##################

if [ -d "TEMPORARY_FILES" ]
then
rm -r TEMPORARY_FILES
fi

mkdir TEMPORARY_FILES

################### Obtain sample size, number of chromosomes, number of SNPs ##################

cp $FILE.map data.map
cp $FILE.ped data.ped

tr '\t' ' ' < data.map > KK1
cut -d ' ' -f1 < KK1 > KK2

grep -w "" -c data.ped > NIND

tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1 > NCHR

SAM=$(grep -w "" -c $FILE.ped)

NCHR=$(tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1)

for((i=1;i<=$NCHR;i++))
do
grep -wc "$i" < KK2 > NCHR$i
done

if [ -f "SNP_CHROM" ]
then
rm SNP_CHROM
fi

for((i=1;i<=$NCHR;i++))
do
cat NCHR$i >> SNP_CHROM
done

rm KK*
rm $FILE.map
rm $FILE.ped

 ################### Divide ped and map files into chromosomes ##################

echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES"
echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES" > timefile

num=$RANDOM
echo "$num" > seedfile

./PROGRAMMES/MANAGE_CHROMOSOMES2>>out<<@
-99
$maxNSNP
@

rm NCHR*
rm NIND
rm SNP_CHROM
###mv checkfile TEMPORARY_FILES/

################### LOOP CHROMOSOMES ##################
### Analysis of linkage disequilibrium in windows of genetic
### distances between pairs of SNPs for each chromosome

if [ $maxNCHROM != -99 ]
then
NCHR=$maxNCHROM
fi

echo "RUNNING ANALYSIS OF CHROMOSOMES ..."
echo "RUNNING ANALYSIS OF CHROMOSOMES" >> timefile

options_for_LD="$SAM $MAF $PHASE $NGEN $NBIN $ZERO $DIST $cMMb"

if [ $threads -eq -99 ]
then
threads=$(getconf _NPROCESSORS_ONLN)
fi

START=$(date +%s)

cp chromosome* TEMPORARY_FILES/

###### LD_SNP_REAL3 #######

### Obtains values of c, d2, etc. for pairs of SNPs in bins for each chromosome
for ((n=1; n<=$NCHR; n++)); do echo $n; done | xargs -I % -P $threads bash -c "./PROGRAMMES/LD_SNP_REAL3 % $options_for_LD"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "CHROMOSOME ANALYSES took $DIFF seconds"
echo "CHROMOSOME ANALYSES took $DIFF seconds" >> timefile

######################## SUMM_REP_CHROM3 #########################
### Combination of all data gathered from chromosomes into a single output file

### Adds results from all chromosomes

for ((n=1; n<=$NCHR; n++))
do
cat outfileLD$n >> CHROM
echo "CHROMOSOME $n" >> OUTPUT
sed '2,3d' outfileLD$n > temp
mv temp outfileLD$n
cat parameters$n >> OUTPUT
done

mv outfileLD* TEMPORARY_FILES/
rm parameters*

./PROGRAMMES/SUMM_REP_CHROM3>>out<<@
$NGEN	NGEN
$NBIN	NBIN
$NCHR	NCHR
@

mv chrom* TEMPORARY_FILES/

echo "TOTAL NUMBER OF SNPs" >> OUTPUT_$FILE
cat nsnp >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
echo "HARDY-WEINBERG DEVIATION" >> OUTPUT_$FILE
cat outfileHWD >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
cat OUTPUT >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
echo "INPUT FOR GONE" >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
cat outfileLD >> OUTPUT_$FILE

rm nsnp
rm OUTPUT
rm CHROM

############################# GONE.cpp ##########################
### Obtain estimates of temporal Ne from GONE

echo "Running GONE"
echo "Running GONE" >> timefile
START=$(date +%s)

./PROGRAMMES/GONEparallel.sh -hc $hc outfileLD $REPS

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "GONE run took $DIFF seconds"
echo "GONE run took $DIFF seconds" >> timefile

echo "END OF ANALYSES"
echo "END OF ANALYSES" >> timefile

mv outfileLD_Ne_estimates Output_Ne_$FILE
mv outfileLD_d2_sample Output_d2_$FILE
rm outfileLD
rm data.ped
rm data.map
rm out
mv outfileLD_TEMP TEMPORARY_FILES/

###################################################################

# **************** #

### TRANSFER OF FILES TO MAIN DIRECTORY 
cp -r $SCDIR/$SLURM_JOBID/OUTPUT_$FILE $WDIR/$DIR/GONE/$FOLDER.$FILE/OUTPUT_$FILE.r$rep
cp -r $SCDIR/$SLURM_JOBID/Output_d2_$FILE $WDIR/$DIR/GONE/$FOLDER.$FILE/Output_d2_$FILE.r$rep
cp -r $SCDIR/$SLURM_JOBID/Output_Ne_$FILE $WDIR/$DIR/GONE/$FOLDER.$FILE/Output_Ne_$FILE.r$rep

cp -r $SCDIR/$SLURM_JOBID/timefile $WDIR/$DIR/GONE/$FOLDER.$FILE/

### REMOVE FILES 
rm OUTPUT_$FILE.r$rep
rm Output_d2_$FILE.r$rep
rm Output_Ne_$FILE.r$rep
rm -r TEMPORARY_FILES

done

############################# CLEANING OF SCRATCH ###################################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*

