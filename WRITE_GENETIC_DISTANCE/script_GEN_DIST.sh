#script_GEN_DIST
#!/bin/bash
#$ -cwd

#rm script_GEN_DIST.sh.*

################## POP AS ARGUMENT ###########################

#Check number of arguments
if [ $# -ne 2 ]  
then
	echo "Usage: $0 <POP> <n>" 
	exit 1
fi

#Set arguments
POP=$1
n=$2

########################################################

module load gsl/2.1

#Working directory
WDIR=$PWD 

#Scratch directory
mkdir -p /state/partition1/IreneGEN_DIST$n/$SLURM_JOBID/
SCDIR="SCDIR"

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

###################### TRANSFER TO SCRATCH ########################

## MAKE SURE CHROMOSOMES IN reference.map AND INPUT MAPS ARE IN SAME CHROMOSOME CODE
## CHROMOSOMES NEED TO BE INTEGER NUMBERS

cp WRITE_GEN_DIST SCDIR/$SLURM_JOBID/
cp ${POP}_no_gen_dist_mod.map SCDIR/$SLURM_JOBID/
cp ${POP}_no_gen_dist_chrX_mod.map SCDIR/$SLURM_JOBID/
cp reference_autosomes.map SCDIR/$SLURM_JOBID/
cp reference_chrX.map SCDIR/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd SCDIR/$SLURM_JOBID

###################### AUTOSOMES #########################

awk '{print $1,$3,$4}' ${POP}_no_gen_dist_mod.map > data_no_gen_dist.map #remove name of SNPs ($2), hard to read
rm ${POP}_no_gen_dist_mod.map

awk '{print $1,$3,$4}' reference_autosomes.map > dataref.map #remove name of SNPs ($2), hard to read
rm reference_autosomes.map

START=$(date +%s)

./WRITE_GEN_DIST

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "WRITE_GEN_DIST took 		$DIFF seconds" >> timefile

sed -i '$ d' lost_list		# remove last line (writes an extra non-existing SNP)
sed -i '$ d' out.map		# remove last line (writes an extra non-existing SNP)

cp -r SCDIR/$SLURM_JOBID/timefile $WDIR/$POP.timefile
cp -r SCDIR/$SLURM_JOBID/lost_list $WDIR/$POP.lost_list
cp -r SCDIR/$SLURM_JOBID/out.map $WDIR/$POP.map

rm lost_list
rm out.map
rm dataref.map

###################### CHR X #########################

awk '{print $1,$3,$4}' ${POP}_no_gen_dist_chrX_mod.map > data_no_gen_dist.map #remove name of SNPs ($2), hard to read
rm ${POP}_no_gen_dist_chrX_mod.map

awk '{print $1,$3,$4}' reference_chrX.map > dataref.map #remove name of SNPs ($2), hard to read
rm reference_chrX.map

START=$(date +%s)

./WRITE_GEN_DIST

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "WRITE_GEN_DIST CHR X took 		$DIFF seconds" >> timefile

sed -i '$ d' lost_list		# remove last line (writes an extra non-existing SNP)
sed -i '$ d' out.map		# remove last line (writes an extra non-existing SNP)

cp -r SCDIR/$SLURM_JOBID/timefile $WDIR/$POP.timefile
cp -r SCDIR/$SLURM_JOBID/lost_list $WDIR/$POP.chrX.lost_list
cp -r SCDIR/$SLURM_JOBID/out.map $WDIR/$POP.chrX.map

######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf SCDIR/$SLURM_JOBID/ 2> /dev/null
