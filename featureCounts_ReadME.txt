File created on 26.04.2017.
The objective of this script is to download featureCounts via SourceForge (command line) in order to 
avoid recounting the counts from the Yamagishi sample reads every time R is loaded. Additionally, getting featureCounts via the command
line is faster and it is clear which annotation is being used.

Subread for featureCounts was downloaded from SourceForge, version 1.5.2.

Human HG38 annotation was downloaded from 

Five different annotations were downloaded:
1. NCBI annotation from the UCSC table browser (file = hg_ucsc_NCBI.gtf (created via the table browser); http://genome.ucsc.edu/cgi-bin/hgTables)
2. NCBI annotation from illumina iGenomes (file = gene.gtf; entire package downloaded from https://support.illumina.com/sequencing/sequencing_software/igenome.html)
3. Gencode v26 annotation (file= gencode.v26.annotation.gtf; downloaded from https://www.gencodegenes.org/releases/current.html)
4. Ensembl annotation (file=Homo_sapiens.GRCh38.88.gtf.gz; downloaded from ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/)
5. inbuilt annotation = hg38_RefSeq_exon.txt stored in subread package under annotation

All of the above files can be found on the NSCC server under the path: /home/users/ntu/kbobowik/genomes/annotation, with the exception
of the inbuilt annotation which is in: /home/users/ntu/kbobowik/BASH/subread-1.5.2-source/annotation.


The following options were selected to generate the annotation file from the UCSC table:

Clade : Mammal
Genome : Human
Assembly : Dec 2013 GRCh38/hg38
Group : Genes and Gene Predictions
Track : NCBI RefSeq
Table : RefSeq All (ncbiRefSeq)
region : genome
Output fromat : GTF - gene transfer format
Output file : hg_ucsc_NCBI.gtf

The following command was used to execute featureCounts:
featureCounts -T 5 -a ~/genomes/hg_ucsc_NCBI.gtf -t exon -g gene_id -o ~/scratch/counts/`basename $file .sam`_counts.txt $file

Where -T is the number of threads, -a is the UCSC annotation file, -t is the specified feature type (exon by default), -g is the 
attribute type used to group features (eg.  exons) into meta-features (gene_id is the default when using gtf annotation), and -o
is the output file. 

Files were then cleaned up by taking the GeneID, length, and count.

More information on subread featureCounts can be found here: http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf