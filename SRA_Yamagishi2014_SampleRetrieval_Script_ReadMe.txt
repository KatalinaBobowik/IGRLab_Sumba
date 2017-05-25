02.04.17
This script was created to download all of the SRA files from the paper “Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum” (SRA study: DRP000987, downloaded from the SRA Run selector: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=DRP000987&go=go). In total, there are 116 files to be downloaded. 

The following steps to make the script were as follows:

1. Downloaded SRA toolkit from Github (sratoolkit.2.8.2-1-centos_linux64): https://github.com/ncbi/sra-tools/wiki/Downloads

2. Used the following command from https://edwards.sdsu.edu/research/fastq-dump/ to download sequences from SRA and convert to fastq format:
    parallel fastq-dump --outdir fastq_parallel --gzip --read-filter pass --dumpbase --split-files ::: $i
    
Parallel (parallel-latest.tar.bz2) was downloaded from the GNU website and then uploaded to the server from the following website:
http://ftp.gnu.org/gnu/parallel/.
