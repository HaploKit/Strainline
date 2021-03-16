
#Prints a help message
function print_help() {
  echo "Usage: $0 [options] -i reads.fasta -o ./ -p sequencingPlatform"
  echo ""
  echo "Full-length De Novo Viral Haplotype Reconstruction from Noisy Long Reads"
  echo ""
  echo "Author: Xiao Luo"
  echo "Date:   Mar 2021"
  echo ""
  echo "	Input:"
  echo "	reads.fasta:                      fasta file of input long reads."
  echo "	./:                               directory where to output the results."
  echo "	sequencingPlatform:               long read sequencing platform: PacBio (-p pb) or Oxford Nanopore (-p ont)"
  echo ""
  echo "	Options:"
  echo "	--minTrimmedLen INT:              Minimum trimmed read length. (default: 1000)"
  echo "	--topk INT:                       Choose top k seed reads. (default: 200)"
  echo "	--minOvlpLen INT:                 Minimum read overlap length. (default: 1000)"
  echo "	--minIdentity FLOAT:              Minimum identity of overlaps. (default: 0.99)"
  echo "	--minSeedLen INT:                 Minimum seed read length. (default: 3000)"
  echo "	--iter INT:                       Number of iterations for contig extension. (default: 2)"
  echo "	--minDiv FLOAT:                   Minimum global divergence for merging haplotypes. (default: 0.01)"
#  echo "	--perIdentity INT:                Percent identity for haplcomputation. (default: 2)"
  echo "	--minAbun FLOAT:                  Minimum abundance for filtering haplotypes (default: 0.02)"
  echo "	--threads INT, -t INT:            Number of processes to run in parallel (default: number of cores)."
  echo "	--help, -h:                       Print this help message."
  exit 1
}

#Set options to default values
threads=$(nproc)
outdir="."

min_trimmed_len=1000 #3000 for SARS-CoV-2 datasets

topk=500
platform="pb"
min_ovlp_len=1000
min_identity=0.995
o=30
r=0.8
max_ovlps=10000
min_sread_len=3000

iter=2
min_divergence=0.01

percent_identity=97
min_abun=0.02 #TODO