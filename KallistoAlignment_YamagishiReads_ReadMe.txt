File created on 21.04.2017 by Katalina Bobowik.

The objective of this script is to align all 122 reads from the paper “Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum” (SRA study: DRP000987, downloaded from the SRA Run selector: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=DRP000987&go=go).
The pseudoaligner Kallisto was used to align all of the 122 files (separately), stored in the folder /home/users/ntu/kbobowik/scratch/fastq_parallel. 
The manual for the Kallisto pseudoaligner can be found here: https://pachterlab.github.io/kallisto/manual
Human cDNA (downloaded from ENSEMBL, release 88) and Pfalciparum cDNA (downloaded from Sanger UK) were concatenated to give a combined genome.


A typical Kallisto pipeline requires the following:

1) Build and index:
	Kallisto quantifies read files directly without the need for read alignment, so no reference genome is required. However, Kallisto does perform a "pseudoaligment" procedure 
	which requires processing a transcriptome file to create a “transcriptome index”. The command for producing a transcriptome file is:
	kallisto index -i transcripts.idx transcripts.fasta.gz
	
	Where -i is a required field (--index=STRING) and transcripts.idx is the new index you will be creating.
	transcripts.fasta.gz is the FASTA formatted transcript file/transcriptome of the target sequences.
	
	Note: The human GRCh38.p10 cDNA file (Homo_sapiens.GRCh38.cdna.all.fa.gz) was downloaded from ensembl.org (http://www.ensembl.org/info/data/ftp/index.html)
	
2) Next, "kallisto quant" runs the quantification algorithm:
	For paired end data:
	
		kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
		
	Where -o output is a mandatory flag, which indicates the directory to save the output to. 
	This is followed by the two paired end/mate pair reads. Fragment length is automatically inferred.
		
	And for single-end data:
		
		kallisto quant -i index -o output --single -l 200 -s 20 file1.fastq.gz file2.fastq.gz file3.fastq.gz
		
	
	"--single" must be specified for single-end reads, as well as "-l 200", which is the fragment length (NOT the read length),
	and "-s 20", which is the standard deviation of the fragment length. Multiple fastq files can be input into one line.
	
	Note: The paper by Yamagishi et al 2014 uses single-end read data, so only this second command is used. 
