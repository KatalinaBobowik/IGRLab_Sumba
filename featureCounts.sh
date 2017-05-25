### File created on 26.04.17
### Objective of file creation is to get the feature counts of all of the Yamagishi et al. reads.

for file in /home/users/ntu/kbobowik/scratch/aligned_reads/filtered/*.sam; do
    featureCounts -T 5 -a ~/genomes/genes.gtf -t exon -g gene_id -o ~/scratch/counts/NCBI/illumina_iGenomes/`basename $file .sam`_counts.txt $file
    cat `basename $file .sam`_counts.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > trimmed_GeneLengthCount_`basename $file .sam`_counts.txt
done

