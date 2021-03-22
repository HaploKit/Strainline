#!/bin/bash

set -e

##############################################
### Step4: low frequent haplotypes removal ###
##############################################
old_hap_fa=$1
corrected_reads_fa=$2 #corrected read fasta
min_abun=$3 #0.02

percent_identity=97
threads=36
export min_abundance=$min_abun

mkdir -p filter_by_abun2
cd filter_by_abun2 || exit
#fa_read=../corrected.0.fa #corrected reads

for i in {1..2};
do
  if [[ $i == 1 ]];then
    cat $old_hap_fa | perl -ne 'BEGIN{$i=1;}if(/^>/){print ">contig_$i\n";$i+=1;}else{print;}' >haps.fa
  else
    cat haplotypes.final.fa | perl -ne 'BEGIN{$i=1;}if(/^>/){print ">contig_$i\n";$i+=1;}else{print;}' >haps.fa
  fi
  minimap2 -a --secondary=no -t $threads haps.fa $corrected_reads_fa | samtools view -F 3584 -b -t $threads | samtools sort - >haps.bam

  jgi_summarize_bam_contig_depths haps.bam --percentIdentity $percent_identity --outputDepth haps.depth

  perl -e 'open A,"haps.depth";<A>;my$alldepth=0;while(<A>){my@a=split;$alldepth+=$a[2];}close A; \
  open A,"haps.depth";<A>;while(<A>){my@a=split;my$d=$a[2]/$alldepth;print "$a[0]\t$a[1]\t$a[2]\t$d\n";}close A; ' |
    sort -k4nr >haps.depth.sort

  perl -e ' my%id2seq;$/=">";open A,"haps.fa";<A>;while(<A>){chomp;my@a=split;$id2seq{$a[0]}=$a[1];}close A;
  $/="\n"; open A,"haps.depth.sort";my$k=0;while(<A>){chomp;my@a=split;next if $a[-1]<$ENV{"min_abundance"};
  my$abun=sprintf "%.3f",$a[-1];my$cov=sprintf "%.0f",$a[-2]; $k+=1;print ">hap$k $cov"."x freq=$abun\n$id2seq{$a[0]}\n";}close A; ' >haplotypes.final.fa
done

cp haplotypes.final.fa ../victor.fa

echo 'All steps finished successfully.'
