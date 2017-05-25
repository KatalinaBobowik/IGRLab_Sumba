17.04.17

This script was created to get stats for each aligned human read in the folder /home/users/ntu/kbobowik/scratch/aligned_reads/aligned_human_${i}.sam
on the NSCC server. 

Requirements: Pysam installtion

File should produce summary statistics of all alignemts in a text file (output: /home/users/ntu/kbobowik/Sumba/Human_BWA_AlignmentStats/humanAlignment_${i}_YamagishiReads_PysamSummaryStatistics.txt) 

The following information was gathered for human and malaria files:
1) total reads
2) total mapped reads
3) total unmapped reads
4) number of qcfail fails
5) mapping quality <= 30",
6) number of mismatches

One more file was run which concatenated all files together for human alignment statistics (saved as /home/users/ntu/kbobowik/Sumba/PF_BWA_AlignmentStats/humanAlignment_all_YamagishiReads_PysamSummaryStatistics.txt)
This gathered the same information with the exception of te number of mismatches. Additionally, percentages were left out so that only one file with "clean" numbers was produced. This is so that 
downstream manipulation can be conducted with R.