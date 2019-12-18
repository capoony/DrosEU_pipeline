# DrosEU_pipeline
The bioinformatics pipeline for the generation and analysis of the DrosEU data of 2014.

#  Created by the DrosEU consortium on 13/11/17.
#

############## A) trim, map and re-align around InDels

## 1) Trim raw FASTQ reads for BQ >18 and Minimum length > 75bp with [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

```bash
export PATH=$PATH:scripts/cutadapt-1.8.3/bin

cutadapt \
-q 18 \
--minimum-length 75 \
-o trimmed-read1.fq.gz \
-p trimmed-read2.fq.gz \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
-B CAAGCAGAAGACGGCATACGAGAT \
-O 15 \
-n 3 \
read1.fq.gz \
read2.fq.gz
```

## 2) map trimmed reads with [bwa](https://sourceforge.net/projects/bio-bwa/files/) and filter for intact pairs with MQ > 20 using [samtools](http://samtools.sourceforge.net/) 

```bash
export PATH=$PATH:scripts/samtools-0.1.19
export PATH=$PATH:scripts/bwa-0.7.15

bwa mem \
-M \
-t 24 \
reference.fa.gz \
trimmed-read1.fq.gz \
trimmed-read2.fq.gz \
| samtools view \
-Sbh -q 20 -F 0x100 - > library.bam
```

## 3) sort by position

java \
-Xmx20g \
-Dsnappy.disable=true \
-jar scripts/picard-tools-1.109/SortSam.jar \
I=library.bam \
O=library-sort.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

## 4) de-duplicate

java \
-Xmx20g \
-Dsnappy.disable=true \
-jar scripts/picard-tools-1.109/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
I=library-sort.bam \
O=library-dedup.bam \
M=library-dedup.txt \
VALIDATION_STRINGENCY=SILENT

## 5) add read group to BAM files

java -jar -Xmx10g scripts/picard-tools-1.109/AddOrReplaceReadGroups.jar \
INPUT=librtary-dedup.bam \
OUTPUT=library-dedup_rg.bam \
SORT_ORDER=coordinate \
RGID=library \
RGLB=library \
RGPL=illumina \
RGSM=sample \
RGPU=library \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT

## 6) generate target list of InDel positions

mkdir $out/mapping/realign_list

java -jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R reference.fa \
-I library-dedup_rg.bam \
-o library-dedup_rg.list

## 7) re-align around InDels

java -Xmx20g -jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R reference.fa \
-I library-dedup_rg.bam \
-targetIntervals library-dedup_rg.list \
-o library-dedup_rg_InDel.bam

############## B) Decontamination of libraries with simulans contamination

## 1) obtain simulans genome

wget -O reference/sim_genome.fa http://datadryad.org/bitstream/handle/10255/dryad.62629/dsim-all-chromosome-M252_draft_4-chrnamesok.fa?sequence=1

## 2) add "sim_" to headers

sed 's/>/>sim_/g' reference/sim_genome.fa | gzip -c > reference/sim_genome_prefix.fa.gz

## 3) combine with melanogaster reference

zcat reference/sim_genome_prefix.fa.gz | cat reference.fa - | gzip -c > reference/combined.fa.gz

## 4) extract reads from BAM file

export PATH=$PATH:scripts/bam2fastq-1.1.0
bam2fastq -s -o reads/library# library-dedup_rg_InDel.bam

## 5) competitive mapping

export PATH=$PATH:scripts/bwa-0.7.15
bwa mem -Mt 20 reference/combined.fa.gz reads/library\_1.gz reads/library\_2.gz > library_deSim.sam

## 6) deconvolute the reads

python2.7 scripts/python/FixBAM.py \
--contaminated library-dedup_rg_InDel.bam \
--prefix sim_ \
--detect library_deSim.sam \
--output library_deSim

############## C) Merging BAM files and joint SNP calling

## 1) merge BAM files (in BAMlist.txt) into MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20

export PATH=$PATH:scripts/samtools-0.1.19
samtools mpileup -B \
-f reference.fa \
-b BAMlist.txt \
-q 20 \
-Q 20 \
| gzip > DrosEU.mpileup.gz

## 2) call SNPs with PoolSNP

##  PoolSNP is a heuristic SNP caller, which uses an MPILEUP file and a FASTA references as inputs. NOTE, that the FASTA headers may NOT contain any special characters, such as "/\|,:", or else they will be ignored. Heuristic parameters to be passed to the script are: Minimum coverage across all libraries (min-cov), maximum coverage, which is calculated for every library and chromosomal arm as the percentile of a coverage distribution (max-cov), minimum allele count of the minor allele across all libraries combined (min-count), minimum allele frequency across all libraries combined (min-freq) and missing fraction, whihch is the maximum percentage of libraries that are allowed to NOT full-fill all above criteria (miss-frac).

##  PoolSNP creates a gzipped VCF file v. 4.2 containing  allele counts and frequencies for every position and library, a max-coverage file containing the maximum coverage thressholds for all chromosomal arms and libraries in the mpileup file (separarted by a column) and optionally a "bad-sites" file (by setting the parameter BS=1), which contains a list of (variable and invariable) sites that did not pass the SNP calling criteria. This file can be used to weight windows for the calulation of PopGen statistics (see below).

##  PoolSNP is a python script that is by default single threaded. To process large datasets, I provide a shell script using GNU parallel to utilize multiple threads. Parameters need to be passed to the shell script and all necessary steps will be processed serially (or in parallel, whenever possible). However, the three python scripts being part of the PoolSNP pipeline can also be used as standalone scripts. Please look into the well-documented shellscript (PoolSNP.sh) for more details.

##  PoolSNP has been tested on Mac OSX (10.11) and Linux Ubuntu (16.10). The shellscript only works with a BASH shell and requires Python 2.7 and GNU parallel to be installed in Path

##  to get more help on the different parameters, just execute the shellscript without parameters

bash scripts/PoolSNP/PoolSNP.sh \
mpileup=DrosEU.mpileup.gz \
reference=reference.fa.gz \
names=1_Mauternbach,2_Mauternbach,3_Yesiloz,4_Yesiloz,5_Viltain,7_Viltain,8_Gotheron,9_Sheffield,10_SouthQueensferry,11_Nicosia,12_MarketHarborough,13_Lutterworth,14_Broggingen,15_Broggingen,16_Yalta,18_Yalta,19_Odessa,20_Odessa,21_Odessa,22_Odessa,23_Kyiv,24_Kyiv,25_Varva,26_Piryuatin,27_Drogobych,28_Chernobyl,29_ChernobylYaniv,30_Lund,31_Munich,32_Munich,33_Recarei,34_Lleida,35_Lleida,36_Akaa,37_Akaa,38_Vesanto,39_Karensminde,41_Karensminde,42_ChaletAGobet,43_ChaletAGobet,44_Seeboden,45_Kharkiv,46_Kharkiv,47_ChernobylApple,48_ChernobylPolisske,49_Kyiv,50_Uman,51_Valday,BA_2012_FAT,BA_2012_SPT,FL_sample1,FL_sample2,GA,MA_2012_FAT,MA_2012_SPT,ME_sample1,ME_sample2,NC,NY_2012_FAT,NY_2012_SPT,PA_07_2009,PA_07_2010,PA_07_2011,PA_09_2011,PA_10_2011,PA_11_2009,PA_11_2010,PA_11_2011,PA_2012_FAT,PA_2012_SPT,SC,test,VA_2012_FAT,VA_2012_SPT,VI_2012_FAT,VI_2012_SPT,WI_09_2012,WI_2012_FAT,WI_2012_SPT \
max-cov=0.99 \
min-cov=10 \
min-count=10 \
min-freq=0.001 \
miss-frac=0.2 \
jobs=24 \
BS=1 \
output=SNPs

## 3) identify sites in proximity of InDels

python scripts/DetectIndels.py \
--mpileup DrosEU.mpileup.gz \
--minimum-count 20 \
--mask 5 \
| gzip > InDel-positions_20.txt.gz

## 4) use Repeatmasker to generate GFF with location of TE's

## 5) filter SNPs around InDels and in TE's from original VCF

python2.7 scripts/FilterPosFromVCF.py \
--indel InDel-positions_20.txt.gz \
--te te.gff \
--vcf SNPs.vcf.gz \
| gzip > SNPs_clean.vcf.gz

## 6) annotate with snpEff

java -Xmx4g -jar scripts/snpEff/snpEff.jar \
-ud 2000 \
BDGP6.82 \
-stats  SNPs_clean.html \
SNPs_clean.vcf.gz \
| gzip > SNPs_clean-ann.vcf.gz

############## D) Calculation of unbiased population genetics estimators Pi, Theta and Tajima's D

## 1) convert the VCF to the sync file format

python scripts/VCF2sync.py \
--vcf SNPs_clean-ann.vcf.gz \
| gzip > SNPs.sync.gz

## 2) resample SNPS to a 40x coverage

python scripts/SubsampleSync.py \
--sync SNPs.sync.gz \
--target-cov 40 \
--min-cov 10 \
| gzip > SNPs-40x.sync.gz

## 3) Calculate "true" window-sizes (e.g. for non-overlapping 200kb windows) based on the number of sites that passed the coverage criteria (as calculated from PoolSNP) are not located within TE's and that are not located close to InDels

python TrueWindows.py \
--badcov SNP_BS.txt.gz \
--indel InDel-positions_20.txt.gz \
--te te.gff \
--window 200000 \
--step 200000 \
--output truewindows

## 4) Calculate window-wise Population Genetics parameters Pi, Theta and Tajima's D using Pool-Seq corrections following Kofler et al. 2011

python scripts/PopGen_var.py \
--input SNPs-40x.sync.gz \
--pool-size 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,66,80,80,80,80,80,80,80,80,70,80,80,80 \
--min-count 2 \
--window 200000 \
--step 200000 \
--sitecount truewindows-200000-200000.txt \
--min-sites-frac 0.75 \
--output Popgen

############## E) Inference of Demographic patterns

## 1) isolate SNPs located in Introns < 60 bp length

python scripts/IntronicSnps.py \
--gff dmel-all-filtered-r6.09.gff.gz \
--sync SNPs.sync.gz \
--target-length 60 \
> intron60.sync

## 2) calculate pairwise FST based on the method of Weir and Cockerham 1984

python scripts/FST.py \
--pool-size 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,66,80,80,80,80,80,80,80,80,70,80,80,80 \
--input intron60.sync \
--minimum-count 2 \
--minimum-cov 10 \
| gzip > intron.fst.gz

## 3) average FST across all loci

python scripts/CombineFST.py \
--diff intron.fst.gz \
--stat 0 \
> intron_average.fst

## 4) calculate isolation by distance (IBD)

python scripts/IBD.py \
--fst intron_average.fst \
--meta data/MetaData.txt \
--output IBD_EU

## 5) calculate allele frequencies of the major allele

python scripts/sync2AF.py \
--inp intron60.sync \
--out intron60-af

## 6) calculate PCA in R

echo '''

library(gtools)
library(LEA)

# load data
meta=read.table("data/MetaData.txt",header=T)
freq=read.table("intron60-af_freq.txt",header=F)
rown<-meta[,1]
rownames(freq)<-rown

# calculate PCA
write.lfmm(freq,"test.lfmm")
pc=pca("test.lfmm")
tw=tracy.widom(pc)
a=stars.pval(tw$pvalues)

# identify optimal number of clusters
d=data.frame(pc$eigenvectors[,1:4)
library(mclust)
d_clust=Mclust(as.matrix(d), G=1:4)
m.best <- dim(d_clust$z)[2]

# identify cluster with K-means clustering
comp <- data.frame(pc$eigenvectors[,1:4])
k <- kmeans(comp, 5, nstart=25, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,"Set1"), 0.5))

# plot first two axes
pdf("PCA.pdf",width=14,height=10)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1), heights=c(1.5,1))
par(cex=1.5,mar=c(4,4,2,2))
plot(-1*comp[,1],comp[,2], col=k$clust, pch=16,cex=1.5,xlab="PC1",ylab="PC2")
names=c("Austria","Austria","Turkey","Turkey","France","France","France","UK","UK","Cyprus","UK","UK","Germany","Germany","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Sweden","Germany","Germany","Portugal","Spain","Spain","Finland","Finland","Finland","Denmark","Denmark","Switzerland","Switzerland","Austria","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Ukraine","Russia")

text( -1*comp[,1],comp[,2],names, pos= 3,cex=0.75,pch=19)
barplot(pc$eigenvalues[,1],ylab="Eigenvalues",names.arg=1:nrow(pc$eigenvalues),xlab="Principal components")
abline(h=1,col="red",lwd=2,lty=2)
b=barplot(cumsum(tw$percentage),ylim=c(0,1),names.arg =1:length(tw$percentage),ylab="variance explained",xlab="Principal components")
abline(h=0.8,col="blue",lwd=2,lty=2)
dev.off()

# write PCA scores of first 3 axes to text file
write.table(cbind(k$cluster,comp[,1],comp[,2],comp[,3]),file="PCA-scores.txt",row.names = T, quote=F)

''' > PCA.r

Rscript PCA.r

############## F) Inversion frequencies and correlations with geographical variables

## 1) obtain Marker-SNP positions for full sync dataset

python scripts/OverlapSNPs.py \
--source data/inversion_markers_v6.txt_pos \
--target SNPs.sync.gz \
> inversion_markers.sync

## 2) calculate inversion frequencies based on inversion-specific marker SNPs

python scripts/InvFreq.py \
--sync inversion_markers.sync \
--meta data/MetaData.txt \
--inv data/inversion_markers_v6.txt \
> inversion.freq

## 3) test for correlation of inversion and TE frequencies with geographic variables and account for spatial autocorrelation

python scripts/Test4Correlation.py \
--meta data/MetaData.txt \
> Correlation.test






