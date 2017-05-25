### File created on 21.04.2017
### output can be analyzed with the sleuth tool

### Download human cDNA from ENSEMBL and Pfalciparum cDNA from Sanger using wget, then combined:
wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/CURRENT/Pfalciparum.cdna.fasta.gz
cat Homo_sapiens.GRCh38.cdna.all.fa.gz Pfalciparum.cdna.fasta.gz > Homo_sapiens.GRCh38.PF.cdna.combined.fa.gz

### Build an index from HG38 cDNA transcripts
kallisto index -i ~/scratch/Kallisto/Human_PF_Combined/HG38_PF_combined_transcripts.idx ~/genomes/Homo_sapiens.GRCh38.PF.cdna.combined.fa.gz

### Quantification for single-end reads. Note: the fragment lenth is unknown, so as per teh google groups site (https://groups.google.com/forum/#!topic/kallisto-sleuth-users/h5LeAlWS33w), 200bp is used. For more accuracy, and still with the sake of time considered, 100 bootstraps are used.
for j in `cat ~/SRR_Acc_List-1.txt`; do
    if [ ! -d /home/users/ntu/kbobowik/scratch/Kallisto/Kallisto_filteredSRA/sample_${j} ]; then
    echo $j
    kallisto quant -i ~/scratch/Kallisto/Human_PF_Combined/HG38_PF_combined_transcripts.idx -t 4 -o ~/scratch/Kallisto/Human_PF_Combined/sample_${j} -b 100 --single -l 200 -s 20 ~/scratch/fastq_parallel/${j}_pass_1.fastq.gz
    fi
done

### Remove decimals from ensembl IDs
for file in `ls -d /home/users/ntu/kbobowik/scratch/Kallisto/Human_PF_Combined/*/`; do
echo $(basename $file)
awk -F "\t" '{gsub(/\..*$/,"",$1)}1' OFS="\t" ${file}/abundance.tsv > ${file}/noDecimal_abundance.tsv
done

