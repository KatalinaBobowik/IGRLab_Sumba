#02.04.17
#This script was created to download all of the SRA files from the paper “Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum” (SRA study: DRP000987, downloaded from the SRA Run selector: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=DRP000987&go=go). In total, there are 116 files to be downloaded.

# Used the following command from https://edwards.sdsu.edu/research/fastq-dump/ to download sequences from SRA and convert to fastq format:
# If, on the off chance your script is interrupted or you run out of time on the server, this loop allows you to skip over the files you've already downloaded
for i in `cat /home/users/ntu/kbobowik/SRR_Acc_List-1.txt`; do
    if [ ! -f /home/users/ntu/kbobowik/scratch/fastq_parallel/${i}_pass_1.fastq.gz ]; then
	echo $i
	parallel fastq-dump --outdir ~/scratch/fastq_parallel --gzip --read-filter pass --dumpbase --split-files ::: $i
	rm /home/users/ntu/kbobowik/ncbi/public/sra/*
    fi
done

# The following script was also executed to get files downloaded without any flags (stored in a separate, "noFlags" folder). Run on 24.04.2017
for i in `cat /home/users/ntu/kbobowik/SRR_Acc_List-1.txt`; do
    if [ ! -f /home/users/ntu/kbobowik/scratch/fastq_noFlags/${i}.fastq.gz ]; then
	echo $i
	parallel fastq-dump --gzip --outdir ~/scratch/fastq_noFlags ::: $i
	rm /home/users/ntu/kbobowik/ncbi/public/sra/*
    fi
done

# In parallel to speed things up...
parallel fastq-dump --gzip --outdir ~/scratch/fastq_noFlags2 ::: `cat /home/users/ntu/kbobowik/SRR_Acc_List-1.txt`
