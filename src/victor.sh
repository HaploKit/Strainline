#!/bin/bash



#input parameters
input_fa=$1
threads=$3
outdir=$4

min_trimmed_len=1000 #3000 for SARS-CoV-2 datasets
basepath=`dirname $0`

##############################################
######## Step1: read error correction ########
##############################################
if [ ! -d $outdir ];then
  mkdir $outdir
else
  echo $outdir 'already exists, please use a new one, exiting...'
  exit
fi

cd $outdir || exit
ln -fs $input_fa reads.fasta
fasta2DAM reads.dam reads
DBsplit -s256 -x$min_trimmed_len reads.dam # -x: Trimmed DB has reads >= this threshold.
mkdir -p tmp
HPC.daligner reads.dam -P./tmp -T$threads| bash
daccord -t$threads reads.las reads.dam >reads.corrected.fa

##############################################
########### Step2: read clustering ###########
##############################################


#reformat fasta
seqkit sort --by-length ../../reads.corrected.fa  -w 0 -r|perl -ne 'if (/^>/){print ">r$.\n";}else{print;}' >corrected.fa
python $basepath/reformat_fa.py reads.corrected.fa corrected.fa

in_fa=corrected.fa
outdir='.'
topk=500
platform=pb
threads=46
min_ovlp_len=1000
min_identity=0.995
o=30
r=0.8
max_ovlps=10000
min_sread_len=3000

python $basepath/clustering.py $in_fa $outdir $topk $platform $threads $min_ovlp_len $min_identity $o $r $max_ovlps $min_sread_len

echo 'done'




##############################################
#### Step3: redundant haplotypes removal #####
##############################################


##############################################
### Step4: low frequent haplotypes removal ###
##############################################


