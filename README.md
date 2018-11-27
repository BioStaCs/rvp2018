# rvp2018
Companion repo for RVP paper, 2018



## 1 Samples Collection

### 1.1 Criteria for sample collection

The samples for RVP database were collected from multiple projects including our own in-house projects. In total, WES data from **740** male and **519** female with well-documented reproductive history and raw sequencing data (bam/fastq) were collected.



### 1.2 Samples quality control

Samples fullfilled with at least one following QC metrics were excluded from variant calling. Finally, 1109 samples passed QC named high quality bams and used for variant calling.

> - Samples had contamination estimate (FREEMIX) > 0.075. Taged as `CONTAMINATION`
> - Samples had low depth of target region (Average depth < 20X). Taged as `LOW DEPTH`
> - Samples had  low converage of target region (Target converage > 70%). Taged as `LOW CONVERAGE`

|      Tag      | Samples |
| :-----------: | :-----: |
| CONTAMINATION |    4    |
|   LOW DEPTH   |   51    |
| LOW CONVERAGE |   95    |



### 1.3 Consortium of samples in RVP 

The sources of high quality bams

|  Source  | Fertile Samples (Male / Female) | Sequencing Strategy |
| :------: | :-----------------------------: | :-----------------: |
| In house |         184 (118 / 66)          |         WES         |
|   SRA    |          71 (36 / 35)           |         WES         |
|   SSC    |         321 (163 / 158)         |         WES         |
|   PPMI   |         533 (352 / 181)         |         WES         |
|  Total   |        1109 (669 / 440)         |         WES         |

> - `SRA` Sporadic samples collected in SRA database. (https://www.ncbi.nlm.nih.gov/sra)
> - `SSC` Simons Simplex Collection. (https://www.sfari.org/resource/simons-simplex-collection/)
> - `PPMI` Parkinson’s Progression Markers Initiative. (https://www.ppmi-info.org/)
> - `In house` The Pakistan and Chinese samples collection were supported by the National Key Research and Developmental Program of China and the National Natural Science Foundation of China. 





## 2 Data Generation

### 2.1 Alignment and BAM processing



#### 2.1.1 Read slignment and bam sorting

The sequencing reads (fastq files) were aligned to the human genome reference (hg19) using `bwa (v0.7.17)  `and then the bam file were sorted using `samtools (v1.7)`.

For PE data:

> - Program(s)
>
>   - `bwa 0.7.17`
>   - `samtools 1.7`
>
> - Input(s)
>
>   - R1 fastq.gz file: \$fastq1
>   - R2 fastq.gz file: $fastq2
>   - BWA Index: \$bwaIndex
>
> - Output(s)
>
>   - BAM file: \$bam
>   - Sorted BAM file: $bam_sort
>
> - Command(s)
>
>   - ```shell
>     $header='@RG\tID:READGROUP_${sampleID}\tSM:${sampleID}\tPL:ILLUMINA\tLB:LIBRARY'
>     $threads=24
>     ```
>
>
>   - ```shell
>     bwa mem -t $threads -R $header $bwaIndex $fastq1 $fastq2 | samtools view -S -b - > $bam
>     ```
>
>   - ```shell
>     samtools sort -@ $threads -O bam -o $bam_sort $bam
>     ```

For SE data:

> - Program(s)
>
>   - `bwa 0.7.17`
>   - `samtools 1.7`
>
> - Input(s)
>
>   - fastq.gz file: \$fastq
>   - BWA Index: \$bwaIndex
>
> - Output(s)
>
>   - BAM file: \$bam
>   - Sorted BAM file: $bam_sort
>
> - Command(s)
>
>   - ```shell
>     $header='@RG\tID:READGROUP_${sampleID}\tSM:${sampleID}\tPL:ILLUMINA\tLB:LIBRARY'
>     $threads=24
>     ```
>
>   - ```shell
>     bwa mem -t $threads -R $header $bwaIndex $fastq | samtools view -S -b - > $bam
>     ```
>
>   - ```shell
>     samtools sort -@ $threads -O bam -o $bam_sort $bam
>     ```



#### 2.1.2 Read duplicates marking

PCR duplicate reads were then marked for each BAM using `Picard MarkDuplicates`. 

>- Program(s)
>
>   - `MarkDuplicates(Picard 2.10.2)`
>   - `samtools 1.7`
>
>- Input(s)
>
>   - Sorted BAM file: \$bam_sort
>
>- Output(s)
>
>   - Markdup metrics file: \$markdup_metrics
>   - BAM file with duplicates marked: $bam_markdup
>   - BAM index file: ${bam_markdup}.bai
>
>- Command(s)
>
>   - ```shell
>     picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT \
>     INPUT=$bam_sort OUTPUT=$bam_markdup M=$markdup_metrics
>     ```
>   - ```shell
>     samtools index $bam_markdup
>     ```



#### 2.1.3 Bam quality control

Performing quality control on each bam files. A python script `bamQC.py` (https://github.com/BioStaCs/rvp2018/blob/master/src/bamQC.py) produces summary statistics of the following quality metrics used for quality control.

- Alignments
- Mapping rate
- Mapping quality
- Insert size (for PE)
- Duplication rate
- Capture rate
- Depth in target
- Target coverage
- Target>20X

> - Program(s)
>   - `bamQC.py`
>
> - Input(s)
>   - BAM file with duplicates marked: $bam_markdup
>   - Exome calling interval (bed): \$exome_calling_region
>
> - Output(s)
>   - Markdup metrics file: \${prefix}.bamqc.txt
>
> - Command(s)
>
>   - ```shell
>     python bamQC.py -b $bam_markdup -r $exome_calling_region -o ${prefix}.bamqc.txt
>     ```



#### 2.1.4 Insertion/Deletion Realignment

The re-alignment intervals for each BAM was determined using `GATK RealignerTargetCreator` and a list of known indel sites. Using this interval list, local realignment was then performed by `GATK IndelRealigner`.

> - Program(s)
>
>   - `RealignerTargetCreator (GATK 3.7.0)`
>   - `IndelRealigner (GATK 3.7.0)`
>
> - Input(s)
>
>   - Human GRCh37 reference (fasta): \$reference
>   - BAM file with duplicates marked: $bam_markdup
>   - Gold indel interval: $gold_indel
>
> - Output(s)
>
>   - re-alignment intervals file: \${prefix}.IndelRealigner.intervals
>   - Re-aligned BAM file: \$bam_realign
>
> - Command(s)
>
>   - ```shell
>     gatk -T RealignerTargetCreator -R $reference -known $gold_indel \
>     -I $bam_markdup -o ${prefix}.IndelRealigner.intervals
>     ```
>
>   - ```shell
>     gatk -T IndelRealigner --filter_bases_not_stored \
>     -R $reference -known $gold_indel \
>     -targetIntervals ${prefix}.IndelRealigner.intervals \
>     -I $bam_markdup -o $bam_realign
>     ```



#### 2.1.5 Base Quality Recalibration

The base quality scores were then recalibrated using `GATK BaseRecalibrator` and a list of known variant sites. The new base quality scores were then applied using `GATK PrintReads` but retaining the original base quality scores within the BAM.

> - Program(s)
>
>   - `BaseRecalibrator (GATK 3.7.0)`
>   - `PrintReads (GATK 3.7.0)`
>
> - Input(s)
>
>   - Human GRCh37 reference (fasta): \$reference
>   - Re-aligned BAM file: \$bam_realign
>   - Gold indel interval: $gold_indel
>   - Knwon snps in dbsnp: \$known_snp
>
> - Output(s)
>
>   - The parameters file used for BQSR: \${prefix}.recal_data.table
>   - Re-calibrated BAM file: \$bam_bqsr
>
> - Command(s)
>
>   - ```shell
>     gatk -T BaseRecalibrator -R $reference \
>     -knownSites $gold_indel -knownSites $known_snp \
>     -I $bam_realign -o ${prefix}.recal_data.table
>     ```
>
>   - ```shell
>     gatk -T PrintReads -R $reference --BQSR ${prefix}.recal_data.table \
>     -I $bam_realign -o $bam_bqsr
>     ```



#### 2.1.6 Contamination Detection

The sample contimination wes estimated by `VerifyBamID 1.1.3`. Samples with FREEMIX > 0.075 were excluded from variant calling.

> - Program(s)
>
>   - `VerifyBamID 1.1.3`
>
> - Input(s)
>
>   - Common snps from 1000 genomes(phase3_v4_20130502): $common_snp
>   - Re-calibrated BAM file: \$bam_bqsr
>   - Re-calibrated BAM index file: \$bam_bqsr_bai
>
> - Output(s)
>
>   - Contamination estimation results: \$prefix
>
> - Command(s)
>
>   - ```shell
>     verifyBamID --vcf $common_snp --bam $bam_bqsr --bai $bam_bqsr_bai -o $prefix
>     ```





### 2.2 Variant Calling and Recalibration



#### 2.2.1 Single sample variant calling

The `Genome Analysis Toolkit (GATK) v3.7.0` HaplotypeCaller algorithm was used to generate gVCFs for all 1109 BAMs across the [exome interval set defined by ExAC project](ftp://ftp.broadinstitute.org/pub/ExAC_release/current/resources/exome_calling_regions.v1.interval_list) totaling 59,684,884 bp and known sites were annotated with dbSNP138.

>- Program(s)
>
>   - `HaplotypeCaller (gatk 3.7.0)`
>
>- Input(s)
>
>   - Exome calling interval (bed): \$exome_calling_region
>   - Human GRCh37 reference (fasta): \$reference
>   - known variants in dbSNP 138 (vcf): \$known_variant
>   - Re-calibrated BAM file: \$bam_bqsr
>
>- Output(s)
>
>   - gVCF file from each bam: $gvcf_single
>
>- Command(s)
>
>   - ```shell
>     gatk -T HaplotypeCaller --emitRefConfidence GVCF \
>     -R $reference -L $exome_calling_region --dbsnp $known_variant \
>     -I $bam_bqsr -o $gvcf_single
>     ```



#### 2.2.2 Joint genotyping

The sample gVCFs were combined into 33 groups with $\sqrt{samples}$ samples in each group.

And all 33 grouped gVCFs were then used as input for joint genotyping. To save the running time, we performed the genotyping on each chromosome separately.

> - Program(s)
>
>   - `CombineGVCFs(gatk 3.7.0)`
>   - `GenotypeGVCFs(gatk 3.7.0)`
>   - `GatherVcfs (Picard 2.10.2)`
>
> - Input(s)
>
>   - Human GRCh37 reference (fasta): \$reference
>   - known variants in dbSNP 138 (vcf): \$known_variant
>   - Re-calibrated BAM file: \$bam_bqsr
>
>
>   - gVCF file from each bam: $gvcf_single
>
> - Output(s)
>
>   - Grouped gVCF file: \$gvcf_group
>   - Raw VCF (combind VCF of all high-quality BAMs): \$prefix.raw.vcf
>
> - Command(s)
>
>   - ```shell
>     gatk -T CombineGVCFs -R $reference -o $gvcf_group \
>     --variant $gvcf_single1
>     --variant $gvcf_single2
>     ...
>     ```
>
>   - ```shell
>     gatk -T GenotypeGVCFs -R $reference -D $known_variant -L $chr \
>     -o $prefix.$chr.raw.vcf
>     --variant $gvcf_group1
>     --variant $gvcf_group2
>     ...
>     ```
>
>   - ```shell
>     picard GatherVcfs O=$prefix.raw.vcf \
>     I=$prefix.chr1.raw.vcf \
>     I=$prefix.chr2.raw.vcf \
>     ...
>     ```





#### 2.2.3 Variant recalibration

GATK Variant Quality Score Recalibration (VQSR) was used to filter variants. The SNP VQSR model was trained using HapMap3.3 and 1KG Omni 2.5 SNP sites and a 99.5% sensitivity threshold was applied to filter variants, while Mills et. al. 1KG gold standard was used for insertions/deletion sites and a 95.0% sensitivity threshold was used.

>- Program(s)
>   - `VariantRecalibrator (gatk 3.7.0)`
>   - `ApplyRecalibration (gatk 3.7.0)`
>- Input(s)
>   - Exome calling interval (bed): \$exome_calling_region
>   - Human GRCh37 reference (fasta): \$reference
>   - Raw VCF file: \$prefix.raw.vcf
>   - Gold standard genotyped snp in HapMap V3.3: \$hapmap
>   - Gold standard genotyped snp in 1000G (omni 2.5): \$omini
>   - High quality genotyped snp in 1000G (phase1): \$phase1
>   - Recorded snp in dbSNP138: \$dbsnp
>   - Mills and 1000G gold standard indels: \$mills
>   - Raw VCF (combind VCF of all high-quality BAMs): \$prefix.raw.vcf
>
>
>- Output(s)
>
>   - $prefix.snps.VQSR.vcf
>   - $prefix.VQSR.vcf
>
>- Command(s)
>
>   - ```shell
>     gatk -T VariantRecalibrator \
>     --disable_auto_index_creation_and_locking_when_reading_rods \
>     -R $reference -input $prefix.raw.vcf \
>     -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
>     -resource:omini,known=false,training=true,truth=false,prior=12.0 $omini \
>     -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase1 \
>     -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
>     -an QD -an MQ -an MQRankSum -an ReadPosRankSum \
>     -an FS -an SOR -an DP -an InbreedingCoeff \
>     -mode SNP --maxGaussians 6 \
>     -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 \
>     -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \
>     -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
>     -recalFile $prefix.snps.recal \
>     -tranchesFile $prefix.snps.tranches \
>     -rscriptFile $prefix.snps.plots.R
>     ```
>
>   - ```shell
>     gatk -T ApplyRecalibration \
>     -R $reference \
>     -L $exome_calling_region \
>     --ts_filter_level 99.5 \
>     -recalFile $prefix.snps.recal \
>     -tranchesFile $prefix.snps.tranches \
>     -mode SNP \
>     -input $prefix.raw.vcf \
>     -o $prefix.snps.VQSR.vcf
>     ```
>
>   - ```shell
>     gatk -T VariantRecalibrator \
>     -R $reference -input $prefix.snps.VQSR.vcf \
>     -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \
>     -an QD -an MQ -an MQRankSum -an ReadPosRankSum \
>     -an FS -an SOR -an DP -an InbreedingCoeff \
>     -mode INDEL --maxGaussians 6 \
>     -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \
>     -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
>     -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \
>     -tranche 91.0 -tranche 90.0 \
>     -recalFile $prefix.snps.indels.recal \
>     -tranchesFile $prefix.snps.indels.tranches \
>     -rscriptFile $prefix.snps.indels.plots.R
>     ```
>
>   - ```shell
>     gatk -T ApplyRecalibration \
>     -R $reference \
>     -L $exome_calling_region \
>     --ts_filter_level 95.0 \
>     -recalFile $prefix.snps.indels.recal \
>     -tranchesFile $prefix.snps.indels.tranches \
>     -mode INDEL \
>     -input $prefix.snps.VQSR.vcf \
>     -o $prefix.VQSR.vcf
>     ```



### 2.3 Variant Filtering



#### 2.3.1 Assessment of variant filtering

The VQSLOD score is known to have a bias towards common and well behaved variants, thus the 99.5% sensitivity threshold will lose many singleton SNPs for population variant data sets . We performed a series of exploration for the relationship between Variant Quality Score Log Odds (VQSLOD) and other metrics (i.e. titv ratio, singleton transmission ratio and Mendel inconsistancy)  to determine the best VQSLOD threshold for variant filtering as described in ExAC project.

We reduced the stringency moving the VQSLOD threshold from -2.0622 (black dash line) to -6.63 (red dash line) with minimal impact on transition to transversion (TiTv) ratio from 1.92 to 1.84.

The adjustment from standard VQSR filtering was assessed using singleton transmission and Mendel inconsistency. We use 21 trios to assess the transmission of singletons from parents. Reducing the stringency resulted in an improved average transmission rates from 50.7% to 49.9% while the cumulative Mendel inconsistency was with a minor increase from 3.60% to 4.09%.

One additional round of site filtering was performed to filter sites with inbreeding coefiicient (InbreddingCoeff) < -0.2 to remove sites with excess heterozygous individuals and sites with AC (allele count) = 0, which occurs when remaining genotype calls in the release subset does not meet minimum quality threshold of DP >= 10 and GQ >= 20. In contrast, the Indel PASS cut off was not adjusted.

![01_metrics_vs_vqsrank](D:\jarning\ongoing_project\rvp2018\figures\data_generation\01_metrics_vs_vqsrank.png)

#### 2.3.2 Data quality evaluation

The NA12878 exome call set by our pipeline generated from NIST NA12878 WES data on illumina HiSeq platform was evaluated against the NIST v2.18 consensus call set. The  NA12878 exome have 50% ExAC exome v1 intervals with depth >= 20, which is lower than most of our data (black points). In total, there were 255 variants (209 SNPs) that were absent in NIST v2.18. The false discovery rate (FDR) for NA12878 was imporoved from 2.34% to 1.12% after adjusted VQSR and genotype filter. Meanwhile, the sensitivity decreased from 81.35% to 64.87%, which cased by the low depth of raw reads.

![04_NA12878_eval](figures\data_generation\04_NA12878_eval.png)

#### 2.3.3 High quality call set

Variant sites met the following criteria were identified as high-quality: 

(1) they were marked a PASS filter status by adjusted VQSR (see above), 

(2) at least 100 supported chromsomes (which means at least 50 individuals in autosomes) in the dataset had at least depth (DP) >= 10 and genotype quality (GQ) >= 20, 

(3) there was at leasat one individual harboring the alternate allele with DP >= 10 and GQ >= 20. 

The application of this criteria, includes all variaint in chromosome X and Y.



### 2.4 Variant Annotation



#### 2.4.1 Function annotation

Variant annotation was performed using `Variant Effect Predictor` (VEP) version 92 on GRCh37. LoF (loss-of-function) annotation was performed using [`LOFTEE`](https://github.com/konradjk/loftee), a pugin to VEP developed by [MacArthur Lab](http://macarthurlab.org/) to identify LoF variation. VEP is used to determine the following annotations:

- Variant consequence and impact
- cDNA change and amino acid change
- Gene symbol, ensembl ID, transcript ID and whether the transcript is canonical 
- Affected exon and intron numbering
- Affected position on cDNA, CDS or protein level.
- Name of overlapping protein domain
- Allele frequency in 1000 genomes, ESP6500 and genomAD.
- SIFT and Polyphen score
- Existing variant co-located with the variant
- Report Pubmed IDs for publications that cite existing variant.

>- Program(s)
>   - ` ensembl-vep 92.5`
>   - `LOFTEE (vep-plugin)`
>- Input(s)
>   - The final call set with no filtering: $prefix.VQSR.vcf
>
>
>- Output(s)
>
>   - VCF file from VEP+LOFTEE: $prefix.VQSR.vep_loftee.vcf
>
>- Command(s)
>
>   - ```shell
>     vep -i $prefix.VQSR.vcf -o $prefix.VQSR.vep_loftee.vcf \
>     --cache --offline --no_stats --assembly GRCh37 --vcf \
>     --dir_cache $vep_home/cache \
>     --dir_plugin $vep_home/Plugins/loftee \
>     --fasta $vep_home/homo_sapiens/92_GRCh37Homo_sapiens.GRCh37.dna.toplevel.fa.gz \
>     --merged \
>     --allele_number \
>     --af --af_1kg --af_esp --af_gnomad --max_af \
>     --sift b --polyphen b --variant_class \
>     --symbol --total_length --numbers --domains --canonical \
>     --regulatory --biotype --gene_phenotype --pubmed \
>     --plugin LoF,\
>     loftee_path:$vep_home/Plugins/loftee,\
>     human_ancestor_fa:$vep_home/Plugins/loftee/human_ancestor/human_ancestor.fa.rz
>     ```

CADD score is annotated using a python script `04_cadd_annotate_parallel.py`.

To speed up the annotation step, we divided the variant call set into 10Mbp bundles, and annotated them parallelly. 



#### 2.4.2 Allele frequency calculation

A custom python script was used to calculate the allele frequency (AF) between males and females separately in the INFO field of VCF call set.  We also calculated the homozygotes counts (HC) in all samples together with males and females. For variants in not PAR (Pseudoautosomal region) on sex chromosome, we fixed the AF where count the allele number(AN) in males as 1 for each genotype.

| New tag in INFO field | Explanation                  |
| --------------------- | ---------------------------- |
| AC_M                  | male alt allele counts       |
| AC_F                  | female alt allele counts     |
| AN_M                  | male allele number           |
| AN_F                  | female allele number         |
| AC_M                  | male alt allele frequency    |
| AC_F                  | female alt allele frequency  |
| HC                    | homozygote counts            |
| HC_M                  | homozygote counts in males   |
| HC_F                  | homozygote counts in females |



## 3 Protein-truncating variation



## 4 Data Availability

| Call set                                                     | Samples | Filters                                                      | Alt alleles | Annotations                               | Analysis                        |
| ------------------------------------------------------------ | ------- | ------------------------------------------------------------ | :---------- | ----------------------------------------- | ------------------------------- |
| RVP fertile cohort r0.1                                      | 1,109   | None                                                         | 1,027,259   | AF and HC                                 | None                            |
| RVP fertile cohort r0.1 -Variant site filtering              | 1,109   | Adjusted VQSR PASS with imbreeding coff > -0.2               | 916,171     | AF and HC                                 | NA12878 sensitivity/specificity |
| RVP fertile cohort r0.1 -Variant site filtering and Genotype filtering | 1,109   | Adjusted VQSR PASS with imbreeding coff > -0.2 + High quality variants | 836,315     | AF, HC and VEP+CADD functional annotation | All analysis in paper           |

## 5 Acknowledgements



## 6 References

 