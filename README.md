# Victor
## Description
Full-length de novo viral haplotype reconstruction from noisy long reads



## Installation and dependencies

Victor relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- [daccord](https://github.com/gt1/daccord)
- [samtools](http://www.htslib.org/)
- [spoa](https://github.com/rvaser/spoa)
- `jgi_summarize_bam_contig_depths` program from [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- Python3


To run Victor, firstly, it is recommended to install the dependencies through [Conda](https://docs.conda.io/en/latest/).
Also, [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) and [DALIGNER](https://github.com/thegenemyers/DALIGNER) 
are required before running `daccord`.
```
conda create -n victor
conda activate victor
conda install -c bioconda minimap2 samtools dazz_db daligner metabat2
```
Then, install `daccord` and `spoa` following the instructions in their corresponding github pages.
Make sure that `daccord -h` and `spoa -h` can work successfully.

## Running and options


## Examples

- PacBio CLR reads
```
    python victor.py -i reads.fa -t 8 -p pb 
```

- ONT reads
```
    python victor.py -i reads.fa -t 8 -p ont 
```


## Citation