10.04.17

This script is a workflow for RNASeq analysis using samples obtained from the paper “Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum” (SRA study: DRP000987, downloaded from the SRA Run selector: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=DRP000987&go=go).
This script has 6 main components:

1) Read processing from BWA alignment: utilise built-in annotation file and summarize mapped reads to NCBI RefSeq gene
2) Count loading and annotation
3) Filtering and normalisation
4) Quality Control and data exploration
5) Differential expression
6) GO analysis


Reads were downloaded from the SRA Run Selector with the following command: --gzip --read-filter pass --dumpbase --split-files. In order to make sure that there is a positive linear relationship between filtered and unfiltered reads, TPM (mapped with Kallisto) for samples with the highest, middle, and lowest library size were plotted against each other.
Three samples with the highest, middle, and lowest library sizes (DRR006389, DRR006380, DRR006487, respectively) were mapped with
Kallisto pseudoaligner for filtered and unfiltered files (i.e., filtered being an SRA fastq dump download with the flags specified above and unfiltered being no flags at all).
Plotting the tpm of filtered vs unfiltered TPM from both libraries showed a positive linear relationship between each, indicating that the filtered library is sufficient to use.

Quality of reads from the SRA Run selector (above) were visualised by plotting the total number of mapped and unmapped reads.
Many reads were of poor quality so the following command was executed in Samtools:
samtools view -h -q 1 -F 4 -F 256 $file > ~/scratch/aligned_reads/filtered/${file}
Where -f 4 filters out reads which did not map to the genome, -f 256 filters not primarily aligned reads, and -q 1 filters out 
reads that have a mapping quality <1.

In order to obtain information on sex and age, supplementary table 11 from the Yamagishi et al study was downloaded, as well as the SRA Run table summary
from the SRA Run Selector. This information was used to create a matrix with patient names (used in Yamagishi et al study), accession numbers in SRA Run,
sex, and age.
Note: one patient name in the Yamagishi et al supplementary 11 table was not the same as in the SRA Run table summary (i.e., one says "malaria7#09" and one says "malaria7#009").

Subread for featureCounts was downloaded from SourceForge, version 1.5.2 and human HG38 annotation was downloaded from the UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables.
Counts for all samples were performed using featureCounts via command line rather than inside R.


Kallisto pseudoaligner was used to quantify transcript abundance. Tximport was then used to import and summarize the transcript-level estimates computed by Kallisto.
Files were manipulated in command line to remove decimals from ensembl IDs. The following BASH script was used:

for file in `ls -d /home/users/ntu/kbobowik/scratch/Kallisto/Human_PF_Combined/*/`; do
	echo $(basename $file)
	awk -F "\t" '{gsub(/\..*$/,"",$1)}1' OFS="\t" ${file}/abundance.tsv > ${file}/noDecimal_abundance.tsv
done
