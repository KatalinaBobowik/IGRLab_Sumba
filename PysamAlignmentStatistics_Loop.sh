### File created on 17.04.17
### Purpose of this file is to get summary statistics of the 122 aligned samples from the Yamagishi et al study (aligned using BWA). Statistics are gathered below for Plasmodium Falciparum alignment and human alignemnt files.

### Loop through all files in human reads folder
for i in `cat ~/SRR_Acc_List-1.txt`; do
    if [ ! -f /home/users/ntu/kbobowik/Sumba/Human_BWA_AlignmentStats/humanAlignment_${i}_YamagishiReads_PysamSummaryStatistics.txt ]; then
	#keep track of which file you're on
	echo $i
	python /home/users/ntu/kbobowik/Scripts/Pysam_AlignmentStatistics_Loop.py /home/users/ntu/kbobowik/scratch/aligned_reads/aligned_human_${i}.sam > /home/users/ntu/kbobowik/Sumba/Human_BWA_AlignmentStats/humanAlignment_${i}_YamagishiReads_PysamSummaryStatistics.txt
	fi
done


### Loop through all files in P. Falciparum reads folder
for i in `cat ~/SRR_Acc_List-1.txt`; do
    if [ ! -f /home/users/ntu/kbobowik/Sumba/PF_BWA_AlignmentStats/PF_Alignment_${i}_YamagishiReads_PysamSummaryStatistics.txt ]; then
	#keep track of which file you're on
	echo $i
	python /home/users/ntu/kbobowik/Scripts/Pysam_AlignmentStatistics_Loop.py /home/users/ntu/kbobowik/scratch/aligned_reads_pfalciparum/aligned_PF_${i}.sam > /home/users/ntu/kbobowik/Sumba/PF_BWA_AlignmentStats/PF_Alignment_${i}_YamagishiReads_PysamSummaryStatistics.txt
	fi
done


### For simplified human output. concatenated into one file
for i in `cat ~/SRR_Acc_List-1.txt`; do
    if [ -f /home/users/ntu/kbobowik/Sumba/Human_BWA_AlignmentStats/humanAlignment_${i}_YamagishiReads_PysamSummaryStatistics.txt ]; then
	#keep track of which file you're on
	echo $i
	a=`python /home/users/ntu/kbobowik/Scripts/Pysam_AlignmentStatistics_YamagishiReads_Altered.py /home/users/ntu/kbobowik/scratch/aligned_reads/aligned_human_${i}.sam`
	b="$i $a"
	echo $b >> /home/users/ntu/kbobowik/Sumba/PF_BWA_AlignmentStats/humanAlignment_all_YamagishiReads_PysamSummaryStatistics.txt
	fi
done
