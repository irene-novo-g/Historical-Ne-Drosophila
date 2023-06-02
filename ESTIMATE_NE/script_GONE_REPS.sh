#script_GONE_REPS.sh
#!/bin/bash
#$ -cwd

rm script_GONE_REPS.sh.*

########################################################

#Check number of arguments
if [ $# -ne 3 ]  
then
	echo "Usage: $0 <FILE> <REPLICATES> <d>" 
	exit 1
fi

### Set arguments

FILE=$1  ### Data file name
REPLICATES=$2  ### Number of replicates
d=$3     ### Directory number

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

module load gcc/7.2.0

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

#Working directory
WDIR=$PWD 

#Output directory
if [ -d "$WDIR/OUTPUT_GONE_$FILE" ]
then
rm -r $WDIR/OUTPUT_GONE_$FILE
fi

mkdir -p $WDIR/OUTPUT_GONE_$FILE

#Scratch directory
mkdir -p /state/partition1/GONE_REAL_REPS$d/$SLURM_JOBID/
mkdir -p /state/partition1/GONE_REAL_REPS$d/$SLURM_JOBID/PROGRAMMES
SCDIR="/state/partition1/GONE_REAL_REPS$d"

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

###################### TRANSFER TO SCRATCH ########################

if [ $maxNCHROM != -99 ]
then
NCHR=$maxNCHROM
fi

cp data/$FILE.map SCDIR/$SLURM_JOBID/
cp data/$FILE.ped SCDIR/$SLURM_JOBID/

cp INPUT_PARAMETERS_FILE SCDIR/$SLURM_JOBID/
cp PROGRAMMES/MANAGE_CHROMOSOMES2 SCDIR/$SLURM_JOBID/PROGRAMMES/
cp PROGRAMMES/LD_SNP_REAL3 SCDIR/$SLURM_JOBID/PROGRAMMES/
cp PROGRAMMES/SUMM_REP_CHROM3 SCDIR/$SLURM_JOBID/PROGRAMMES/
cp PROGRAMMES/GONE SCDIR/$SLURM_JOBID/PROGRAMMES/
cp PROGRAMMES/GONEaverage SCDIR/$SLURM_JOBID/PROGRAMMES/
cp PROGRAMMES/GONEparallel.sh SCDIR/$SLURM_JOBID/PROGRAMMES/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd SCDIR/$SLURM_JOBID

cp -r SCDIR/$SLURM_JOBID/INPUT_PARAMETERS_FILE $WDIR/OUTPUT_GONE_$FILE/

################################ REPLICATES #############################

for ((r=1; r<=$REPLICATES; r++))
do

echo "REPLICATE $r"
echo "REPLICATE $r" > timefile

################### Obtain sample size, number of chromosomes, number of SNPs ##################

cp $FILE.map data.map
cp $FILE.ped data.ped

tr '\t' ' ' < data.map > KK1
cut -d ' ' -f1 < KK1 > KK2

grep -w "" -c data.ped > NIND

tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1 > NCHR

SAM=$(grep -w "" -c data.ped)

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
###mv checkfile $WDIR/OUTPUT_GONE_$FILE/

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

#cp chromosome* $WDIR/OUTPUT_GONE_$FILE/

###### LD_SNP_REAL3 #######

### Obtains values of c, d2, etc. for pairs of SNPs in bins for each chromosome
for ((n=1; n<=$NCHR; n++)); do echo $n; done | xargs -I % -P $threads bash -c "./PROGRAMMES/LD_SNP_REAL3 % $options_for_LD"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "CHROMOSOME ANALYSES took $DIFF seconds"
echo "CHROMOSOME ANALYSES took $DIFF seconds" >> timefile

cp -r SCDIR/$SLURM_JOBID/LD* $WDIR/OUTPUT_GONE_$FILE/

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

#mv outfileLD* $WDIR/OUTPUT_GONE_$FILE/
rm parameters*

./PROGRAMMES/SUMM_REP_CHROM3>>out<<@
$NGEN	NGEN
$NBIN	NBIN
$NCHR	NCHR
@

#mv chrom* $WDIR/OUTPUT_GONE_$FILE/

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
#rm data.ped
#rm data.map
rm out
#mv outfileLD_TEMP $WDIR/OUTPUT_GONE_$FILE/

cp -r SCDIR/$SLURM_JOBID/Output_Ne_$FILE $WDIR/OUTPUT_GONE_$FILE/Ne_$FILE.$r
#cp -r SCDIR/$SLURM_JOBID/Output_d2_$FILE $WDIR/OUTPUT_GONE_$FILE/d2_$FILE.$r
#cp -r SCDIR/$SLURM_JOBID/OUTPUT_$FILE $WDIR/OUTPUT_GONE_$FILE/OUTPUT_$FILE.$r

rm chromosome*

echo "REPLICATE $r" >> GONE_d2
cat Output_d2_$FILE >> GONE_d2
cp SCDIR/$SLURM_JOBID/GONE_d2 $WDIR/OUTPUT_GONE_$FILE/GONE_d2

echo "REPLICATE $r" >> GONE_STATS 
cat OUTPUT_$FILE >> GONE_STATS 
cp SCDIR/$SLURM_JOBID/GONE_STATS $WDIR/OUTPUT_GONE_$FILE/GONE_STATS 

tail -n +3 Output_Ne_$FILE > temp
mv temp Output_Ne_$FILE

if [ $r -eq 1 ]
then

cut -f1 Output_Ne_$FILE > generations
cut -f2 Output_Ne_$FILE > Ne.all

else

cut -f2 Output_Ne_$FILE > Ne.rep
paste Ne.all Ne.rep > temp2
mv temp2 Ne.all
rm Ne.rep

awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' Ne.all > Ne.mean
awk '{ A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ; for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }' Ne.all > Ne.sd

paste generations Ne.mean Ne.sd > Ne.summary
sed  -i '1i Generations\tMean_Ne\tStandard_Deviation_Ne' Ne.summary
cp SCDIR/$SLURM_JOBID/Ne.summary $WDIR/OUTPUT_GONE_$FILE/Ne.summary 

fi

paste generations Ne.all > GONE_Ne_all
cp SCDIR/$SLURM_JOBID/GONE_Ne_all $WDIR/OUTPUT_GONE_$FILE/GONE_Ne_all 




done

######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf SCDIR/$SLURM_JOBID/ 2> /dev/null

###################################################################


