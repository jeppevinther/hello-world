#######
#Final reanalysis of the pilot m7G experiment using optimised Bowtie2 settings and getFreq script
#
#
#single read sequencing, 75 bp (NextSeq platform)
#data received on 04.02.16

#samples 
#m7G analysis
# index	NaBH4	ARP	Selc	Genotype
#     1				WT
#     2	+			WT
#     3		+		WT
#     4	+	+		WT
#     5			+	WT
#     6		+	+	WT
#     7	+	+	+	WT
#     9				RsmG KO
#     10	+			RsmG KO
#     11		+		RsmG KO
#     12	+	+		RsmG KO
#     13			+	RsmG KO
#     14		+	+	RsmG KO
#     15	+	+	+	RsmG KO


#Make figures for the first figure.using index 1,2 and 9,10.
#Estimate background mutation frenquency in same data.

#Make supplementary figure 1: demonstrate rates of mutation and RT-stop and compare with background
#

#05.02.16
#download data
cd /seqdata/krogh/linedp/ftp.dna.ku.dk
mkdir 050216
cd 050216
nice wget --http-user=NmRkYjk4 --http-password=ZjVhYjYw --timeout 10 --tries 10 -m http://185.45.23.199/pickup/NmRkYjk4/HMG5MBGXX_casava_1_mismatches.tar.gz
cd /seqdata/krogh/linedp/ftp.dna.ku.dk/050216/185.45.23.199/pickup/NmRkYjk4
#uncompress
tar zxvf HMG5MBGXX_casava_1_mismatches.tar.gz

#Put the data in isilon /binf-isilon/vintherlab/jvinther/050216/fastq
cp -r /seqdata/krogh/linedp/ftp.dna.ku.dk/050216/185.45.23.199/pickup/NmRkYjk4 /binf-isilon/vintherlab/jvinther/050216/fastq

mkdir /binf-isilon/vintherlab/jvinther/050216/scratch

#Combine files and trim adapters with Cutadapt and do quality trimming for nextseq data
#data is in /binf-isilon/vintherlab/jvinther/050216/fastq/NmRkYjk4/N127

cd /binf-isilon/vintherlab/jvinther/050216/scratch
mkdir /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/


#for the WT data

for lab_number in 1 2 3 4 5 6 7
do
FILES=(/binf-isilon/vintherlab/jvinther/050216/fastq/NmRkYjk4/N127/m7G-Wt-"$lab_number"*R1*.fastq.gz)
if (( ${#FILES[@]}>0 ));
then
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data
mkdir $lab_number
cd $lab_number
nice zcat /binf-isilon/vintherlab/jvinther/050216/fastq/NmRkYjk4/N127/m7G-Wt-"$lab_number"*R1*.fastq.gz | nice cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2>$lab_number.cutadapt.error | gzip > "$lab_number".fastq.gz & 
fi
done
wait


# For the RsmG data
for lab_number in 9 10 11 12 13 14 15
do
FILES=(/binf-isilon/vintherlab/jvinther/050216/fastq/NmRkYjk4/N127/m7G-RsmG-"$lab_number"*R1*.fastq.gz)
if (( ${#FILES[@]}>0 ));
then
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data
mkdir $lab_number
cd $lab_number
nice zcat /binf-isilon/vintherlab/jvinther/050216/fastq/NmRkYjk4/N127/m7G-RsmG-"$lab_number"*R1*.fastq.gz | nice cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2>$lab_number.cutadapt.error | gzip > "$lab_number".fastq.gz & 
fi
done
wait

##########################RNAprobR

#use guide from http://people.binf.ku.dk/lukasz/Guide.html
#Download scripts and decompress them in scripts directory
mkdir /data/jvinther/scripts/RNAprobBash_scripts
cd /data/jvinther/scripts/RNAprobBash_scripts
wget http://people.binf.ku.dk/~jvinther/data/rna_probing/RNAprobBash.tar.gz
cd RNAprobBash_scripts
tar -zxf RNAprobBash.tar.gz

#Add the scripts location to the path
PATH=$PATH:/data/jvinther/scripts/RNAprobBash_scripts

#########################remove 7n with awk script (RNAprobr Bash shell)

for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
preprocessing.sh -b NNNNNNN -t 9 -1 "$lab_number".fastq.gz -o /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
gzip read1.fastq
done


########################Map with high sensitivity to E coli

#Coli mapping
for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
nice bowtie2 --local -N 1 -D 20 -R 3 -L 15 --norc -p32 --quiet -x /binf-isilon/vintherlab/jvinther/Sequences/Transcriptome_scratch/Coli_mRNA_rRNA/Coli_mRNA_rRNA -U read1.fastq.gz  2>bowtie2.error | gzip > tx_mapped.sam.gz &
done


######################
#Filter for barcodes with N in them
# strategy remove lines from barcode. txt file
# in join function for barcode collapse the lines with "missing barcodes will be removed

for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
grep -P '^.*\t[^N]{7}' barcodes.txt > barcodes_filtered.txt   #bash grep does not support \t
done



########################Collapse on barcodes using collapse.sh script

PATH=$PATH: /home/jvinther/scripts

for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
collapse.sh tx_mapped.sam.gz barcodes_filtered.txt > tx_mapped_debarcoded.sam.gz 
done

#strict collapse
for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
collapse_strict.sh tx_mapped.sam.gz barcodes_filtered.txt > tx_mapped_sdebarcoded.sam.gz 
done


##########################Making sorted bamfiles
mkdir /binf-isilon/vintherlab/jvinther/050216/scratch/output

for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do 
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
nice samtools view -u -S tx_mapped_debarcoded.sam.gz|samtools sort > /binf-isilon/vintherlab/jvinther/050216/scratch/output/$lab_number.bam
done
wait

#Strict
mkdir /binf-isilon/vintherlab/jvinther/050216/scratch/soutput

for lab_number in 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do 
cd /binf-isilon/vintherlab/jvinther/050216/scratch/trimmed_data/$lab_number
nice samtools view -u -S tx_mapped_debarcoded.sam.gz|samtools sort > /binf-isilon/vintherlab/jvinther/050216/scratch/soutput/$lab_number.bam
done
wait



#######################
#Coli: 
#WT comparison: 2 (treated) vs 1 (control)
cd /binf-isilon/vintherlab/jvinther/050216/scratch/output
mkdir wt
cd wt
echo "" > bam_file.txt

for lab_number in 1 2  #control then treated
do
echo /binf-isilon/vintherlab/jvinther/050216/scratch/output/$lab_number.bam >> bam_file.txt
done
wait

cd /binf-isilon/vintherlab/jvinther/050216/scratch/output/wt
nice samtools mpileup -A -d 300000 -f /binf-isilon/vintherlab/jvinther/Sequences/Transcriptome_scratch/Coli_mRNA_rRNA/Coli_mRNA_rRNA.fa -b bam_file.txt > wt_coli.mpileup

#STRICT
#Coli: 
#WT comparison: 2 (treated) vs 1 (control)
cd /binf-isilon/vintherlab/jvinther/050216/scratch/soutput
mkdir wt
cd wt
echo "" > bam_file.txt

for lab_number in 1 2  #control then treated
do
echo /binf-isilon/vintherlab/jvinther/050216/scratch/soutput/$lab_number.bam >> bam_file.txt
done
wait

cd /binf-isilon/vintherlab/jvinther/050216/scratch/soutput/wt
nice samtools mpileup -A -d 300000 -f /binf-isilon/vintherlab/jvinther/Sequences/Transcriptome_scratch/Coli_mRNA_rRNA/Coli_mRNA_rRNA.fa -b bam_file.txt > swt_coli.mpileup

#######################
#Coli: 
#RSMG comparison: 10 (treated) vs 9 (control)
cd /binf-isilon/vintherlab/jvinther/050216/scratch/output
mkdir rsmg
cd rsmg
echo "" > bam_file.txt
for lab_number in 9 10  #control then treated
do
echo /binf-isilon/vintherlab/jvinther/050216/scratch/output/$lab_number.bam >> bam_file.txt
done
wait
cd /binf-isilon/vintherlab/jvinther/050216/scratch/output/rsmg
nice samtools mpileup -A -d 300000 -f /binf-isilon/vintherlab/jvinther/Sequences/Transcriptome_scratch/Coli_mRNA_rRNA/Coli_mRNA_rRNA.fa -b bam_file.txt > rsmg_coli.mpileup

#STRICT
#Coli: 
#RSMG comparison: 10 (treated) vs 9 (control)
cd /binf-isilon/vintherlab/jvinther/050216/scratch/soutput
mkdir rsmg
cd rsmg
echo "" > bam_file.txt
for lab_number in 9 10  #control then treated
do
echo /binf-isilon/vintherlab/jvinther/050216/scratch/soutput/$lab_number.bam >> bam_file.txt
done
wait
cd /binf-isilon/vintherlab/jvinther/050216/scratch/soutput/rsmg
nice samtools mpileup -A -d 300000 -f /binf-isilon/vintherlab/jvinther/Sequences/Transcriptome_scratch/Coli_mRNA_rRNA/Coli_mRNA_rRNA.fa -b bam_file.txt > srsmg_coli.mpileup


#################################
#Parse mpileup with getFreq2000 
R
source("/home/jvinther/scripts/mutationmaps/getFreq2000.R")
setwd("/binf-isilon/vintherlab/jvinther/050216/scratch/output/wt")
getFreq("wt_coli.mpileup",
                     CC = c(0,1), #control vs treated
                     CSVout = "wt_coli_out.txt", #sample in the mpileup file
                     minFrac = 0.0001, #minimum fraction of mutatation for analysis
                     minCounts = 1, #minimum number of counts for analysis
                     pvalThres = -1, # test all sites even if non-variable
                     nCores = 50, #
                     nSites = 500, #
                     maxSites = 50000000)

cd /binf-isilon/vintherlab/jvinther/050216/scratch/output/wt
grep -P '^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t\d\d\d' wt_coli_out.txt >wt_coli_100_cut_out.txt


################################
#Parse mpileup with getFreq2000 
R
source("/home/jvinther/scripts/mutationmaps/getFreq2000.R")
setwd("/binf-isilon/vintherlab/jvinther/050216/scratch/output/rsmg")
getFreq("rsmg_coli.mpileup",
                     CC = c(0,1), #control vs treated
                     CSVout = "rsmg_coli_out.txt", #sample in the mpileup file
                     minFrac = 0.0001, #minimum fraction of mutatation for analysis
                     minCounts = 1, #minimum number of counts for analysis
                     pvalThres = -1, # test all sites even if non-variable
                     nCores = 50, #
                     nSites = 500, #
                     maxSites = 50000000)

cd /binf-isilon/vintherlab/jvinther/050216/scratch/output/rsmg
grep -P '^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t\d\d\d' rsmg_coli_out.txt >wt_coli_100_cut_out.txt


#STRICT
#################################
#Parse mpileup with getFreq2000 
R
source("/home/jvinther/scripts/mutationmaps/getFreq2000.R")
setwd("/binf-isilon/vintherlab/jvinther/050216/scratch/soutput/wt")
getFreq("swt_coli.mpileup",
                     CC = c(0,1), #control vs treated
                     CSVout = "swt_coli_out.txt", #sample in the mpileup file
                     minFrac = 0.0001, #minimum fraction of mutatation for analysis
                     minCounts = 1, #minimum number of counts for analysis
                     pvalThres = -1, # test all sites even if non-variable
                     nCores = 50, #
                     nSites = 500, #
                     maxSites = 50000000) #for rRNA analysis




################################
#Parse mpileup with getFreq2000 
R
source("/home/jvinther/scripts/mutationmaps/getFreq2000.R")
setwd("/binf-isilon/vintherlab/jvinther/050216/scratch/soutput/rsmg")
getFreq("srsmg_coli.mpileup",
                     CC = c(0,1), #control vs treated
                     CSVout = "srsmg_coli_out.txt", #sample in the mpileup file
                     minFrac = 0.0001, #minimum fraction of mutatation for analysis
                     minCounts = 1, #minimum number of counts for analysis
                     pvalThres = -1, # test all sites even if non-variable
                     nCores = 50, #
                     nSites = 500, #
                     maxSites = 500000000) #for rRNA analysis


#Use the strictly collapsed data: collapsed on position and barcode





