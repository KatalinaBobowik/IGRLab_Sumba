# File created on 24.04.2017 by KSB


## Create genome index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ~/scratch/STAR/genomeIndex/hg38_PF_Combined --genomeFastaFiles ~/genomes/Combined_Hg38_Pf3D7.fasta --limitGenomeGenerateRAM 55996781270

## First Pass
for j in `cat ~/SRR_Acc_List-1.txt`; do
    echo ${j}_firstPass
    # Map the reads
    STAR --genomeDir ~/scratch/STAR/genomeIndex/hg38_PF_Combined --readFilesIn ~/scratch/fastq_parallel/${j}_pass_1.fastq.gz --outFileNamePrefix ~/scratch/STAR/alignedReads_SAM/STAR_aligned_HumanPFCombined_${j} --readFilesCommand zcat --runThreadN 12 --sjdbOverhang 99

## Second Pass- Create genome index and add in SJ.out.tab info from first pass
    echo ${j}_secondPass
    STAR --runMode genomeGenerate --genomeDir ~/scratch/STAR/genomeIndex/hg38_PF_secondpass --genomeFastaFiles ~/genomes/Combined_Hg38_Pf3D7.fasta --sjdbFileChrStartEnd ~/scratch/STAR/alignedReads_SAM/STAR_aligned_HumanPFCombined_${j}SJ.out.tab --limitGenomeGenerateRAM 55996781270 --sjdbOverhang 99 --runThreadN 12

## Map the reads
    echo ${j}_MappingReads
    STAR --genomeDir ~/scratch/STAR/genomeIndex/hg38_PF_secondpass --readFilesIn ~/scratch/fastq_parallel/DRR006375_pass_1.fastq.gz --outFileNamePrefix ~/scratch/STAR/alignedReads_SAM/STAR_aligned_HumanPFCombined_Second_DRR006375 --runThreadN 12 --readFilesCommand zcat

done




