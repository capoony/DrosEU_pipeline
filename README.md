# DrosEU bioinformatics pipeline
The bioinformatics pipeline for the generation and analyses of the DrosEU data of 2014. Created by the DrosEU consortium.

## A) trim, map and re-align around InDels

### 1) Trim raw FASTQ reads for BQ >18 and minimum length > 75bp with [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

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

### 2) map trimmed reads with [bwa](https://sourceforge.net/projects/bio-bwa/files/) and filter for propper read pairs with MQ > 20 using [samtools](http://samtools.sourceforge.net/) 

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

### 3) sort BAM by reference position using [picard](https://broadinstitute.github.io/picard/)

```bash
java \
-Xmx20g \
-Dsnappy.disable=true \
-jar scripts/picard-tools-1.109/SortSam.jar \
I=library.bam \
O=library-sort.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT
```

### 4) remove PCR duplicates with [picard](https://broadinstitute.github.io/picard/)

```bash
java \
-Xmx20g \
-Dsnappy.disable=true \
-jar scripts/picard-tools-1.109/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
I=library-sort.bam \
O=library-dedup.bam \
M=library-dedup.txt \
VALIDATION_STRINGENCY=SILENT
```

### 5) add read group tags to BAM files using [picard](https://broadinstitute.github.io/picard/)

```bash
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
```

### 6) generate target list of InDel positions using [GATK](https://software.broadinstitute.org/gatk/)

```bash
mkdir $out/mapping/realign_list

java -jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R reference.fa \
-I library-dedup_rg.bam \
-o library-dedup_rg.list
```

### 7) re-align around InDels using [GATK](https://software.broadinstitute.org/gatk/)

```bash
java -Xmx20g -jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R reference.fa \
-I library-dedup_rg.bam \
-targetIntervals library-dedup_rg.list \
-o library-dedup_rg_InDel.bam
```

## B) Decontamination of libraries with *D. simulans* contamination

### 1) obtain *D. simulans* genome

```bash
wget -O reference/sim_genome.fa http://datadryad.org/bitstream/handle/10255/dryad.62629/dsim-all-chromosome-M252_draft_4-chrnamesok.fa?sequence=1
```


### 2) add *"sim_"* to FASTA headers

```bash
sed 's/>/>sim_/g' reference/sim_genome.fa | gzip -c > reference/sim_genome_prefix.fa.gz
```

### 3) combine with *D. melanogaster* reference

```bash
zcat reference/sim_genome_prefix.fa.gz | cat reference.fa - | gzip -c > reference/combined.fa.gz
```

### 4) extract high confidence mapped reads from BAM file with [bam2fastq](https://github.com/jts/bam2fastq)

```bash
export PATH=$PATH:scripts/bam2fastq-1.1.0

bam2fastq -s -o reads/library# library-dedup_rg_InDel.bam
```

### 5) competitive mapping of extracted reads against combined reference genomes

```bash
export PATH=$PATH:scripts/bwa-0.7.15

bwa mem -Mt 20 reference/combined.fa.gz reads/library\_1.gz reads/library\_2.gz > library_deSim.sam
```

### 6) deconvolute the reads in the original BAM file

```bash
python2.7 scripts/FixBAM.py \
--contaminated library-dedup_rg_InDel.bam \
--prefix sim_ \
--detect library_deSim.sam \
--output library_deSim
```

## C) Merging BAM files and joint SNP calling

### 1) merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20

```bash
export PATH=$PATH:scripts/samtools-0.1.19

samtools mpileup -B \
-f reference.fa \
-b BAMlist.txt \
-q 20 \
-Q 20 \
| gzip > DrosEU.mpileup.gz
```

### 2) call SNPs with [PoolSNP](https://github.com/capoony/PoolSNP)

See [documentation](https://github.com/capoony/PoolSNP/blob/master/README.md) of PoolSNP for further details 

```bash
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
```

### 3) identify sites in proximity of InDels with a minimum count of 20 across all samples pooled and mask sites 5bp up- and downstream of InDel.

```bash
python scripts/DetectIndels.py \
--mpileup DrosEU.mpileup.gz \
--minimum-count 20 \
--mask 5 \
| gzip > InDel-positions_20.txt.gz
```

### 4) use [Repeatmasker](http://www.repeatmasker.org/) to generate a GFF with location of known TE's

#### obtain TE libraries

```bash
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.10_FB2016_02/fasta/dmel-all-transposon-r6.10.fasta.gz
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.10_FB2016_02/fasta/dmel-all-chromosome-r6.10.fasta.gz
```

#### only keep contig name in headers (no spaces)

```bash
python2.7  scripts/adjust-id.py \
dmel-all-transposon-r6.10.fasta \
> dmel-all-transposon-r6.10_fixed-id.fasta
```

#### repeat mask *D. melanogaster* genome using [Repeatmasker](http://www.repeatmasker.org/)

```bash
scripts/RepeatMasker \
-pa 20 \
--lib dmel-all-transposon-r6.10_fixed-id.fasta \
--gff \
--qq \
--no_is \
--nolow \
dmel-all-chromosome-r6.10.fasta
```

### 5) filter SNPs around InDels and in TE's from the original VCF produced with PoolSNP

```bash
python2.7 scripts/FilterPosFromVCF.py \
--indel InDel-positions_20.txt.gz \
--te dmel-all-chromosome-r6.10.fasta.out.gff \
--vcf SNPs.vcf.gz \
| gzip > SNPs_clean.vcf.gz
```

## 6) annotate SNPs with [snpEff](http://snpeff.sourceforge.net/)

```bash
java -Xmx4g -jar scripts/snpEff-4.2/snpEff.jar \
-ud 2000 \
BDGP6.82 \
-stats  SNPs_clean.html \
SNPs_clean.vcf.gz \
| gzip > SNPs_clean-ann.vcf.gz
```

## D) Calculation of unbiased population genetics estimators Tajima's *pi*, Watterson's *Theta* and Tajima's *D*

### 1) convert the VCF to SYNC file format (see Kofler et al. 2011)

```bash
python scripts/VCF2sync.py \
--vcf SNPs_clean-ann.vcf.gz \
| gzip > SNPs.sync.gz
```

### 2) resample SNPS to a 40x coverage

```bash
python scripts/SubsampleSync.py \
--sync SNPs.sync.gz \
--target-cov 40 \
--min-cov 10 \
| gzip > SNPs-40x.sync.gz
```

### 3) Calculate "true" window-sizes (e.g. for non-overlapping 200kb windows) based on the number of sites that passed the coverage criteria (as calculated from [PoolSNP](https://github.com/capoony/PoolSNP)) are not located within TE's and that are not located close to InDels; See Material and Methods in Kapun *et al.* (2018)

```bash
python scripts/TrueWindows.py \
--badcov SNP_BS.txt.gz \
--indel InDel-positions_20.txt.gz \
--te te.gff \
--window 200000 \
--step 200000 \
--output truewindows
```

### 4) Calculate window-wise Population Genetics parameters Tajima's *pi*, Watterson's *Theta* and Tajima's *D* using Pool-Seq corrections following Kofler *et al.* (2011)

```bash
python scripts/PoolGen_var.py \
--input SNPs-40x.sync.gz \
--pool-size 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,66,80,80,80,80,80,80,80,80,70,80,80,80 \
--min-count 2 \
--window 200000 \
--step 200000 \
--sitecount truewindows-200000-200000.txt \
--min-sites-frac 0.75 \
--output Popgen
```

## E) Inference of Demographic patterns

### 1) isolate SNPs located in Introns < 60 bp length using the *D. melanogaster* genome annotation in GFF file format and retain only SNPS with a minimum recombination rate >3 (based on Comeron et al. 2012) and in a minimum Distance of 1 mb to the breakpoints of common chromosomal inversions (based on Corbett-Detig et al. 2014).

```bash
# identify SNPs inside introns of < 60bp length

python scripts/IntronicSnps.py \
--gff dmel-all-filtered-r6.09.gff.gz \
--sync SNPs.sync.gz \
--target-length 60 \
> intron60_all.sync

# remove SNPs within and in 1mb distance to chromosomal inversions and with recombination rates <3

python scripts/FilterByRecomRateNInversion.py \
--inv data/inversions_breakpoints_v5v6.txt \
--RecRates data/DrosEU-SNPs.recomb.gz \
--input intron60_all.sync \
--D 1000000 \
--r 3 \
> intron60.sync
```

### 2) calculate pairwise *F*<sub>ST</sub> based on the method of Weir & Cockerham (1984)

```bash
python scripts/FST.py \
--pool-size 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,66,80,80,80,80,80,80,80,80,70,80,80,80 \
--input intron60.sync \
--minimum-count 2 \
--minimum-cov 10 \
| gzip > intron.fst.gz
```

### 3) average *F*<sub>ST</sub> across all loci

```bash
python scripts/CombineFST.py \
--diff intron.fst.gz \
--stat 0 \
> intron_average.fst
```

### 4) calculate isolation by distance (IBD)

```bash
python scripts/IBD.py \
--fst intron_average.fst \
--meta data/MetaData.txt \
--output IBD_EU
```

### 5) calculate allele frequencies of the major allele

```bash
python scripts/sync2AF.py \
--inp intron60.sync \
--out intron60-af
```

### 6) calculate PCA in *R*

```R
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
```

### analyse population structure and admixture with the *R* package *conStruct*

```R
# install.packages("conStruct")
library(conStruct)

# Load allele frequencies, coordinates of each sampling site and geographic distances in kilometers
Freq=read.table("intron60-af_freq.txt",header=F)
CoordRaw=read.table("data/DrosEU-coord.txt",header=T)
Coord=as.matrix(CoordRaw[,4:3])
colnames(Coord)=c("Lon","Lat")
DistRaw=read.table("data/DrosEU-geo.dist",header=T)
Dist=as.matrix(DistRaw)

# Set working directory
setwd("/Volumes/MartinResearch2/new2014Analyses/analyses/4fold/construct")

# test values of K ranging from 1 to 10 in 8-fold replication with cross-validation
my.xvals <- x.validation(train.prop = 0.9,
    n.reps = 8,
    K = 1:10,
    freqs = as.matrix(Freq),
    data.partitions = NULL,
    geoDist = Dist,
    coords = Coord,
    prefix = "example2",
    n.iter = 1e3,
    make.figs = T,
    save.files = T,
    parallel = TRUE,
    n.nodes = 20)

# load both the results for the spatial and non-spation models
sp.results <- as.matrix(
    read.table("example2_sp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))
nsp.results <- as.matrix(
    read.table("example2_nsp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))

# format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:10 with 8 replicates and visualize results with confidence interval bars
pdf("cross-validate-sp.pdf",width=4,height=12)
plot(rowMeans(sp.results),
pch=19,col="blue",
ylab="predictive accuracy",xlab="values of K",
ylim=range(sp.CIs),
main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
y0 = sp.CIs[1,],
x1 = 1:nrow(sp.results),
y1 = sp.CIs[2,],
col = "blue",lwd=2)
dev.off()

# plot all Admixture plots for values of K ranging from 1:10
for (i in seq(1,10,1)){
      my.run <- conStruct(spatial = TRUE,
            K = i,
            freqs = as.matrix(Freq),
            geoDist = Dist,
            coords = Coord,
            prefix = paste("test_",i,seq=""),
            n.chains = 1,
            n.iter = 1e3,
            make.figs = T,
            save.files = T)

      admix.props <- my.run$chain_1$MAP$admix.proportions
      pdf(paste("STRUCTURE_",i,"_1.pdf"),width=8,height=4)
      make.structure.plot(admix.proportions = admix.props,
            sample.names=CoordRaw$node_label,
            mar = c(6,4,2,2),
            sort.by=1)
      dev.off()

      # plot map with pie-charts showing proportion admixture  
      pdf(paste("AdmPIE_",i,"_1.pdf"),width=9,height=8)
      maps::map(xlim = range(Coord[,1]) + c(-5,5), ylim = range(Coord[,2])+c(-2,2), col="gray")
      make.admix.pie.plot(admix.proportions = admix.props,
            coords = Coord,
            add = TRUE)
      box()
      axis(1)
      axis(2)
      dev.off()
}
```

## F) Correlation with climatic variation using the WorldClim dataset (Hijmans *et al.* 2005) 

### 1) obtain climatic data

```R
# load raster package
library(raster)

# first load WC bio variables at the resolution of 2.5 deg
biod <- getData("worldclim", var="bio", res=2.5)

# read csv file with geographic coordinates
geod<-read.table("data/DrosEU-coord.txt", header=T)

# extact for each coordinate bio clim variables
bio<-extract(biod, geod[,c(4,3)])

# create a full dataset
bio.data<-cbind(geod,bio)

# save into external file
write.table(bio.data,file="data/climate.txt",sep="\t", row.names=FALSE ,quote=FALSE)
```
### 2) calculate PCA 

```R
#read csv file with geographic coordinates
geod<-read.table("data/climate.txt", header=T)

# get variance and mean for PCA
library(FactoMineR)

# prepare dataset
geo<-geod[,5:nrow(geod)]
rownames(geo)<-geod[,2]

# do PCA
pca<-PCA(geo)

# make Figures
pdf("data/climate_scree.pdf",width=6,height=6)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(2,3), heights=c(1.5,1))
plot(pca)
barplot(pca$eig[,1],ylab="Eigenvalues",names.arg=1:nrow(pca$eig),xlab="Principal components")
abline(h=1,col="blue")
barplot(pca$eig[,3],names.arg=1:nrow(pca$eig),ylab="variance explained",xlab="Principal components")
abline(h=80,col="blue")
dev.off()

# export loadings
cat(capture.output(pca$var$coord),file="data/climate.rot",sep="")

# export PC axes
bio.data<-cbind(geod,pca$ind$coord)
write.table(bio.data,file="data/climate_PCA.txt",sep="\t", row.names=FALSE ,quote=FALSE)
```
### 3) convert vcf file to bscan input files (python 2.7, requires PvVCF (https://pyvcf.readthedocs.io/en/latest/))
```bash
scripts/vcf2Bscan.py -filename SNPs.vcf.gz \
	-qual 15 \
	-depth 10 \
	-task bscan \
	-prefix bscan_2R \
	-region 2R
```

### 4) run BayeScEnv
```bash
bayescenv SNPs.in \
	-pilot 1000 \
	-npb 5 \
	-n 2000 \
	-env env_var.txt \
	-od ~/path/to/output/dir/
```

### 5) run GOwinda for top snps
```bash
gowinda --output-file pc1 \
	--annotation-file dmel-all-r6.12.gtf \
	--snp-file AllSNPS.vcf \
	--candidate-snp-file topSNPs.bed \
	--gene-set-file Dmel_funcassociate_go_associations.txt \
	--simulations 1000000 \
	--gene-definition updownstream2000 \
	--threads 3
```

## G) Inversion frequencies and correlations with geographical variables

### 1) obtain Marker-SNP positions in full SYNC dataset (see Kapun *et al.* 2014 for more details)

```bash
python3 scripts/OverlapSNPs.py \
--source data/inversion_markers_v6.txt_pos \
--target SNPs.sync.gz \
> inversion_markers.sync
```

### 2) calculate inversion frequencies based on inversion-specific marker SNPs

```bash
python3 scripts/InvFreq.py \
--sync inversion_markers.sync \
--meta data/MetaData.txt \
--inv data/inversion_markers_v6.txt \
> inversion.freq
```

## 3) test for correlations of inversion and TE frequencies with geographic variables and account for spatial autocorrelation

```bash
python3 scripts/Test4Correlation.py \
--meta data/MetaData.txt \
> Correlation.test
```

## H) Selective sweeps

### 1) Compute pileups from bam files

```bash
samtools mpileup \
      -B \
      -f holo_dmel_6.12.fa \ 
      -q 20 \
      -Q 20 \
      $bamfile \
      > $out/$name.mpileup
```

### 2) Separate pileup files by chromosome

```awk
awk '$1 == "2L" {print $0}' $name.pileup > $name-2L.pileup
awk '$1 == "2R" {print $0}' $name.pileup > $name-2R.pileup
awk '$1 == "3L" {print $0}' $name.pileup > $name-3L.pileup
awk '$1 == "3R" {print $0}' $name.pileup > $name-3R.pileup
awk '$1 == "4" {print $0}' $name.pileup > $name-4.pileup
awk '$1 == "X" {print $0}' $name.pileup > $name-X.pileup
```

### 3) Run Pool-hmm for SFS (requires pool-hmm 1.4.4)

```bash
python pool-hmm.py \
      --prefix $name \ # Name associated to each sample and each chromosome
      -n 80 \ # Depending on the sample size
      --only-spectrum \ 
      --theta 0.005 \
      -r 100
```

### 4) Run Pool-hmm for sweep detection (requires pool-hmm 1.4.4)
--pred: for prediction of selective window

-k 0.000000000000001: per site transition probability between hidden states (more restrictive as lower); this changed depending on the sample used

-s: spectrum file

-e: phred quality (required sanger for new illumina reads)

outputs: $name.post, $name.pred, $name.segemit, $name.stat


```bash
python pool-hmm.py --prefix $name-2L -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
python pool-hmm.py --prefix $name-2R -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
python pool-hmm.py --prefix $name-3L -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
python pool-hmm.py --prefix $name-3R -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
python pool-hmm.py --prefix $name-4 -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
python pool-hmm.py --prefix $name-X -n 80 --pred -k 0.000000000000001 -s $name-all -e sanger
```

### 5) Transform Pool-hmm stat files into bed files (Python 2.7)

```bash 
python read_stat_pool-hmm.py *.stat
```

### 6) Bedtools to get the genes inside our candidate selective sweeps
Concatenate all chromosome bed files for the same sample and run bedtools with Drosophila melanogaster v.6.12 annotation file

```bash 
cat *.bed
bedtools intersect -a dmel-all-r6.12_FlyBase_genes_nopseudo.bed -b sample.bed > sammple_flybase.txt 
```

### 7) Look for genes overlapping among samples and populations
all_sample_genes.txt contains a list of all genes obtained in all populations

sampleXXX_genes.txt: genes found for each sample

Overlap among samples; output is a: FBgn*_samples.txt; containing information of samples where the FBgnXXX gene was found

```bash
while read p; do
      grep "$p" sampleXXX_genes.txt > "$p"_samples.txt
      done <all_sample_genes.txt
```

Overlap among populations (Python 2.7)

```bash
python genes_per_population.py \
      FBgn*_samples.txt
```

### 9) Calculate average Tajima's D in the 30 samples included for the analysis (Python 2.7)

```bash
python average_window_tajimasd.py
```


## References

Boitard S, Kofler R, Françoise P, Robelin D, Schlötterer C, Futschik A (2013) Pool-hmm: a Python program for estimating the allele frequency spectrum and detecting selective sweeps from next generation sequencing of pooled samples. Mol Ecol Resour, 13, 337–340. 

Comeron JM, Ratnappan R, Bailin S. 2012. The Many Landscapes of Recombination in *Drosophila melanogaster*. PLoS Genetics 8:e1002905.

Corbett-Detig RB, Cardeno C, Langley CH. 2012. Sequence-based detection and breakpoint assembly of polymorphic inversions. Genetics 192:131–137.

Hijmans RJ, Cameron SE, Parra JL, Jones PG, Jarvis A. 2005. Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology 25:1965–1978.

Kapun M, van Schalkwyk H, Bryant M, Flatt T, Schlötterer C. 2014. Inference of chromosomal inversion dynamics from Pool-Seq data in natural and laboratory populations of *Drosophila melanogaster*. Molecular Ecology 23:1813–1827.

Kapun M, Barron Aduriz MG, Staubach F, Vieira J, Obbard D, Goubert C, Rota Stabelli O, Kankare M, Haudry A, Wiberg RAW, et al. 2018. Genomic analysis of European *Drosophila* populations reveals longitudinal structure and continent-wide selection. bioRxiv Available from: http://biorxiv.org/lookup/doi/10.1101/313759

Kofler R, Pandey RV, Schlotterer C. 2011. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27:3435–3436.

Weir BS, Cockerham CC. 1984. Estimating *F*-Statistics for the Analysis of Population Structure. Evolution 38:1358.




