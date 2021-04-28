#!/bin/sh

# ---------------------------------------------------------------
## Scripts to merge bams from the same clone across samples
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

# Follows run_bamslice.sh

#JK136

#clone a
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK136_JK202_mergebams/clone_a.bam \
JK136_reg1_tis/clones/outs/subsets/clone_a.bam \
JK136_reg1_org/clones/outs/subsets/clone_a.bam \
JK136_reg2_tis/clones/outs/subsets/clone_a.bam \
JK136_reg2_org/clones/outs/subsets/clone_a.bam \
JK202_reg1_tis/clones/outs/subsets/clone_a.bam \
JK202_reg1_org/clones/outs/subsets/clone_a.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK136_JK202_mergebams/clone_a.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
-o JK136_JK202_mergebams/clone_a.vcf.gz



#clone c
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK136_JK202_mergebams/clone_c.bam \
JK136_reg1_tis/clones/outs/subsets/clone_c.bam \
JK136_reg1_org/clones/outs/subsets/clone_c.bam \
JK136_reg2_tis/clones/outs/subsets/clone_c.bam \
JK136_reg2_org/clones/outs/subsets/clone_c.bam \
JK202_reg1_tis/clones/outs/subsets/clone_c.bam \
JK202_reg1_org/clones/outs/subsets/clone_c.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK136_JK202_mergebams/clone_c.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
-o JK136_JK202_mergebams/clone_c.vcf.gz



#non-malignant cells

#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK136_JK202_mergebams/clone_non-malignant.bam \
JK136_reg1_tis/clones/outs/subsets/clone_non-malignant.bam \
JK136_reg1_org/clones/outs/subsets/clone_non-malignant.bam \
JK136_reg2_tis/clones/outs/subsets/clone_non-malignant.bam \
JK136_reg2_org/clones/outs/subsets/clone_non-malignant.bam \
JK202_reg1_tis/clones/outs/subsets/clone_non-malignant.bam \
JK202_reg1_org/clones/outs/subsets/clone_non-malignant.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK136_JK202_mergebams/clone_non-malignant.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
-o JK136_JK202_mergebams/clone_non-malignant.vcf.gz









#JK142

#clone a
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_a.bam \
JK142_reg1_tis/clones/outs/subsets/clone_a.bam \
JK142_reg1_org/clones/outs/subsets/clone_a.bam \
JK142_reg2_org/clones/outs/subsets/clone_a.bam \
JK196_reg1_tis/clones/outs/subsets/clone_a.bam \
JK196_reg1_org/clones/outs/subsets/clone_a.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_a.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_a.vcf.gz



#clone b
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_b.bam \
JK196_reg1_tis/clones/outs/subsets/clone_b.bam \
JK196_reg1_org/clones/outs/subsets/clone_b.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_b.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_b.vcf.gz



#clone d
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_d.bam \
JK142_reg1_tis/clones/outs/subsets/clone_d.bam \
JK142_reg1_org/clones/outs/subsets/clone_d.bam \
JK142_reg2_org/clones/outs/subsets/clone_d.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_d.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_d.vcf.gz



#clone e
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_e.bam \
JK142_reg1_tis/clones/outs/subsets/clone_e.bam \
JK142_reg1_org/clones/outs/subsets/clone_e.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_e.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_e.vcf.gz



#clone j
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_j.bam \
JK142_reg1_org/clones/outs/subsets/clone_j.bam \
JK142_reg2_org/clones/outs/subsets/clone_j.bam \
JK196_reg1_tis/clones/outs/subsets/clone_j.bam \
JK196_reg1_org/clones/outs/subsets/clone_j.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_j.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_j.vcf.gz



#clone m
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_m.bam \
JK142_reg1_tis/clones/outs/subsets/clone_m.bam \
JK142_reg1_org/clones/outs/subsets/clone_m.bam \
JK142_reg2_org/clones/outs/subsets/clone_m.bam \
JK196_reg1_tis/clones/outs/subsets/clone_m.bam \
JK196_reg1_org/clones/outs/subsets/clone_m.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_m.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_m.vcf.gz



#non-malignant cells

#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 20 -p \
JK142_JK196_mergebams/clone_non-malignant.bam \
JK142_reg1_tis/clones/outs/subsets/clone_non-malignant.bam \
JK142_reg1_org/clones/outs/subsets/clone_non-malignant.bam \
JK142_reg2_tis/clones/outs/subsets/clone_non-malignant.bam \
JK142_reg2_org/clones/outs/subsets/clone_non-malignant.bam \
JK196_reg1_tis/clones/outs/subsets/clone_non-malignant.bam \
JK196_reg1_org/clones/outs/subsets/clone_non-malignant.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK142_JK196_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK142_JK196_mergebams/clone_non-malignant.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK142_JK196_mergebams/clone_non-malignant.vcf.gz






#JK153

#clone b
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK153_mergebams/clone_b.bam \
JK153_reg1_tis/clones/outs/subsets/clone_b.bam \
JK153_reg1_org/clones/outs/subsets/clone_b.bam \
JK153_reg2_tis/clones/outs/subsets/clone_b.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_b.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_b.vcf.gz

/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.het.genotype.output.filtered.g.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_b.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_b.germline.vcf.gz



#clone c
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK153_mergebams/clone_c.bam \
JK153_reg1_tis/clones/outs/subsets/clone_c.bam \
JK153_reg1_org/clones/outs/subsets/clone_c.bam \
JK153_reg2_tis/clones/outs/subsets/clone_c.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_c.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_c.vcf.gz

/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.het.genotype.output.filtered.g.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_c.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_c.germline.vcf.gz



#clone d
#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK153_mergebams/clone_d.bam \
JK153_reg1_tis/clones/outs/subsets/clone_d.bam \
JK153_reg1_org/clones/outs/subsets/clone_d.bam \
JK153_reg2_tis/clones/outs/subsets/clone_d.bam \
JK153_reg2_org/clones/outs/subsets/clone_d.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_d.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_d.vcf.gz

/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.het.genotype.output.filtered.g.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_d.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_d.germline.vcf.gz



#non-malignant cells

#merge bam files
/gsc/software/linux-x86_64-centos6/sambamba-0.6.1/sambamba_v0.6.1 merge \
-t 10 -p \
JK153_mergebams/clone_non-malignant.bam \
JK153_reg1_tis/clones/outs/subsets/clone_non-malignant.bam \
JK153_reg1_org/clones/outs/subsets/clone_non-malignant.bam \
JK153_reg2_tis/clones/outs/subsets/clone_non-malignant.bam \
JK153_reg2_org/clones/outs/subsets/clone_non-malignant.bam

#genotype variants of interest
/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_non-malignant.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_non-malignant.vcf.gz

/projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools mpileup \
-f /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-R /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.het.genotype.output.filtered.g.vcf.gz \
--ignore-RG \
JK153_mergebams/clone_non-malignant.bam \
| /projects/vleblanc_prj/tools/bcftools-1.9/bin/bcftools call \
-m -Oz \
--ploidy GRCh38 \
-o JK153_mergebams/clone_non-malignant.germline.vcf.gz
