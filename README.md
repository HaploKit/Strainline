# Strainline
## Description

Haplotype-resolved assembly of highly diverse virus genomes is critical in prevention, control and 
treatment of viral diseases. Current methods either only handle accurate short reads, or collapse 
haplotype-specific variations. Here, we present **Strainline**, a novel approach for full-length 
*de novo* viral haplotype reconstruction from noisy long reads. 
Benchmarking experiments on both simulated and real data demonstrate 
that our approach drastically outperforms others in various aspects such as haplotype coverage, 
assembly contiguity and accuracy. To our knowledge, our approach is the first one which can accurately 
and completely reconstruct viral haplotypes at full length and provide abundance estimates from 
error-prone long reads alone.

## Installation and dependencies

Strainline relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- [daccord](https://github.com/gt1/daccord)
- [samtools](http://www.htslib.org/)
- [spoa](https://github.com/rvaser/spoa)
- `jgi_summarize_bam_contig_depths` program from [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- Python3


To run Strainline, firstly, it is recommended to install the dependencies through [Conda](https://docs.conda.io/en/latest/).
Also, [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) and [DALIGNER](https://github.com/thegenemyers/DALIGNER) 
are required before running `daccord`.
```
conda create -n strainline
conda activate strainline
conda install -c bioconda minimap2 samtools dazz_db daligner metabat2
```
Then, install `daccord` and `spoa` following the instructions in their corresponding github pages.
Make sure that `daccord -h` and `spoa -h` can work successfully.

## Running and options
The input read file is required and the format should be FASTA. Other parameters are optional.
Please run `strainline.sh -h` to get details of optional parameters setting.
Before running Strainline, please read through the following basic parameter settings,
which may be helpful to achieve better assemblies. 
```
Usage: strainline.sh [options] -i reads.fasta -o out/ -p sequencingPlatform

Input:
	reads.fasta:                      fasta file of input long reads.
	out/:                             directory where to output the results.
	sequencingPlatform:               long read sequencing platform: PacBio (-p pb) or Oxford Nanopore (-p ont)

Options:
	--minTrimmedLen INT:              Minimum trimmed read length. (default: 1000)
	--topk INT, -k INT:               Choose top k seed reads. (default: 100)
	--minOvlpLen INT:                 Minimum read overlap length. (default: 1000)
	--minIdentity FLOAT:              Minimum identity of overlaps. (default: 0.99)
	--minSeedLen INT:                 Minimum seed read length. (default: 3000)
	--iter INT:                       Number of iterations for contig extension. (default: 2)
	--minDiv FLOAT:                   Minimum global divergence for merging haplotypes. (default: 0.01)
	--minAbun FLOAT:                  Minimum abundance for filtering haplotypes (default: 0.02)
	--rmMisassembly BOOL:             Break contigs at potential misassembled positions (default: False)
	--threads INT, -t INT:            Number of processes to run in parallel (default: 8).
	--help, -h:                       Print this help message.
```


## Examples

One can test the `strainline.sh` program using the small PacBio CLR reads file `example/reads.fa`.
- PacBio CLR reads
```
cd example
../src/strainline.sh -i reads.fa -o out -p pb -k 20 -t 32
```

- ONT reads
```
../src/strainline.sh -i reads.fa -o out -p ont -t 32
```


## Citation