File created on 17.04.2017 by Katalina Bobowik

Objective: Get statistics on .sam files aligned to the P. Falciparum genome.
This script is made to collect data from the BWA alignment files aligned_PF_${file}.sam. The code for this script was taken from the book “Bioinformatics Data Skills”, by Vincent Buffalo (page 391). This script modifies the original script by omitting unnecessary commands, as well as implementing mapping quality and obtaining the number of unmapped reads and unmapped genes.
*Note, $file is the .sam file aligned to the P. Falciparum genome (120 in total)

To run the script, type the following three commands into the terminal:

#The first command calls upon Python, the second command is the Python script, and the third command is the bam/sam file

1) python 2) ~/Documents/Singapore_StemCells/Scripts/Sumba/Python/Pysam_AlignmentStatistics_PF.py 3) ~/Documents/Singapore_StemCells/Script_Output/BWA_OutputFiles/Sumba/PFalciparum_Transcriptome/aligned_PF_${file}.sam

This is then used to loop over every file (120 in total).

Four libraries were also used for this script: 
sys (which enables python to interact with the command line), 
Pysam (an interface for reading and writing SAM files), 
the Counter function from the collections library, 
and the csv function to write .csv files

