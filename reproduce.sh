#!/bin/bash
set -e

######### Run strainline #########

srcpath=/root/capsule/code/src
threads=16

# small example dataset
echo "Running the example dataset..."
input=/root/capsule/data/example/reads.fa
$srcpath/strainline.sh -i $input -o out -p pb -k 20 -t $threads


## Note: Running the following codes on 'Code Ocean' may be time-consuming because of limited computational resources ##
## For all datasets, we provide the final output haplotyps in 'result/' ##

if false
then
# simulated data, pacbio clr, 5-HIV
input=/root/capsule/data/simulated/clr/5-HIV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p pb -k 50 --maxGD 0.02 --maxLD 0.01 --maxOH 50 -t $threads

## For other different sequencing coverage data (500x,1000x,2000x,5000x,10000x)
# input=/root/capsule/data/simulated/clr/5-HIV/diff_coverage/*/reads.fa.gz
# $srcpath/strainline.sh -i $input -o out -p pb -k 50 --maxGD 0.01 --maxLD 0.01 --maxOH 50 -t $threads

# simulated data, pacbio clr, 6-polio
input=/root/capsule/data/simulated/clr/6-polio/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p pb -k 100 --maxOH 50 -t $threads

# simulated data, pacbio clr, 10-HCV
input=/root/capsule/data/simulated/clr/10-HCV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p pb -k 100 --maxGD 0.02 --maxOH 50 -t $threads

# simulated data, pacbio clr, 15-ZIKV
input=/root/capsule/data/simulated/clr/15-ZIKV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p pb -k 500 --minIdentity 0.995 --maxGD 0.005 --minAbun 0.01 -t $threads

# simulated data, pacbio clr, 5-SARSCov2
input=/root/capsule/data/simulated/clr/5-SARSCov2/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p pb -k 100 --maxLD 0.01 --minAbun 0.085 --iter 3 -t $threads




# simulated data, ont, 5-HIV
input=/root/capsule/data/simulated/ont/5-HIV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 20 --maxGD 0.02 --maxLD 0.01 --maxOH 50 -t $threads

# simulated data, ont , 6-polio
input=/root/capsule/data/simulated/ont/6-polio/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 100 --maxOH 50 --maxLD 0.01 -t $threads

# simulated data, ont , 10-HCV
input=/root/capsule/data/simulated/ont/10-HCV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 100 --maxGD 0.02 --maxLD 0.01 --maxOH 50 -t $threads

# simulated data, ont , 15-ZIKV
input=/root/capsule/data/simulated/ont/15-ZIKV/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 100 --minIdentity 0.995 --minAbun 0.01 -t $threads

# simulated data, ont , 5-SARSCov2
input=/root/capsule/data/simulated/ont/5-SARSCov2/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 100 --maxLD 0.01 --minAbun 0.085 -t $threads



# real data, ont, 5-PVY
input=/root/capsule/data/real/5-PVY/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 100 --maxGD 0.02 --maxLD 0.01 --minOvlpLen 400 --minSeedLen 2000 --minAbun 0.03 --maxOH 20 -t $threads

# real data, ont, SARSCov2 [SRA:SRP250446]
input=/root/capsule/data/real/SARSCov2/SRP250446/reads.fa.gz
$srcpath/strainline.sh -i $input -o out -p ont -k 50 --rmMisassembly True --minAbun 0.05 --minIdentity 0.98 --maxLD 0.01 -t $threads

fi



######### Run other tools #########

## Run canu --version, Canu snapshot (), installed on Nov4,2019 via conda ##
## use parameters recommended for metagenome assembly ##
## For PacBio CLR
# canu -p out -d out  genomeSize=GENOMESIZE minReadLength=500 minOverlapLength=200 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200   -pacbio-raw  reads.fa
## For ONT
# canu -p out -d out  genomeSize=GENOMESIZE minReadLength=500 minOverlapLength=200 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200   -nanopore-raw  reads.fa

## Run wtdbg2 Version: 0.0 (19830203) ##
# wtdbg2 -x rs -g GENOMESIZE -t 8 -i $fq -fo out
# wtpoa-cns -t 8 -i out.ctg.lay.gz -fo out.ctg.fa

echo The output assemblies of simulated datasets are located here:
ls /root/capsule/result/simulated/*/*/haplotypes.*.fa
echo
echo the output assemblies of real datasets are located here:
ls /root/capsule/result/real/*/haplotypes.*.fa
echo
echo The output assemblies of simulated datasets[5-HIV, various sequencing coverages] are located here:
ls /root/capsule/result/simulated/*/5-HIV/*/haplotypes.*.fa
echo

######### Evaluation #########
## For assembly evaluation, we used the following commands ##
# quast-5.1.0rc1/metaquast.py -r ref.fa --min-contig 500 -o out --unique-mapping -t 16 contigs.fa
# cat out/runs_per_reference/ref/report.tsv
