#!/bin/sh

#-----------------------------------------------------------
## Scripts to run Mutect2 (from GATK 4.1.4.1) in multisample mode (all samples from each patient)
#-----------------------------------------------------------

gnomad_resource="/projects/vleblanc_prj/genomes/hg38/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.cleaned.noalt.withHeader.vcf.bgz"


#-----------------------------------------------------------
## JK124
#-----------------------------------------------------------

#not including reg2 organoid since low tumour content

#Run mutect2 to get variant calls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-I ./data/exome_seq/bams/JK124_reg1_tissue.bam \
-I ./data/exome_seq/bams/JK124_reg1_organoid.bam \
-I ./data/exome_seq/bams/JK124_reg2_tissue.bam \
-I ./data/exome_seq/bams/JK124_blood.bam \
-normal JK124_blood \
--germline-resource ${gnomad_resource} \
-L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
--native-pair-hmm-threads 32 \
-O ./results/exome/mutect2_multisamp/JK124_noreg2org/somatic.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK124_noreg2org/log

#Run FilterMutectCalls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
-V ./results/exome/mutect2_multisamp/JK124_noreg2org/somatic.vcf.gz \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-O ./results/exome/mutect2_multisamp/JK124_noreg2org/somatic.filter.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK124_noreg2org/filter_log

#Make vcf files with variants passing filter (bgziped and tabix-indexed)
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
./results/exome/mutect2_multisamp/JK124_noreg2org/somatic.filter.vcf.gz \
| bgzip > ./results/exome/mutect2_multisamp/JK124_noreg2org/passed.somatic.filter.vcf.gz

tabix ./results/exome/mutect2_multisamp/JK124_noreg2org/passed.somatic.filter.vcf.gz

sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ./results/exome/mutect2_multisamp/JK124_noreg2org passed.somatic.filter





#-----------------------------------------------------------
## JK136
#-----------------------------------------------------------

#Run mutect2 to get variant calls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-I ./data/exome_seq/bams/JK136_reg1_tissue.bam \
-I ./data/exome_seq/bams/JK136_reg1_organoid.bam \
-I ./data/exome_seq/bams/JK136_reg2_tissue.bam \
-I ./data/exome_seq/bams/JK136_reg2_organoid.bam \
-I ./data/exome_seq/bams/JK202_tissue.bam \
-I ./data/exome_seq/bams/JK202_organoid.bam \
-I ./data/exome_seq/bams/JK136_JK202_blood.bam \
-normal JK136_JK202_blood \
--germline-resource /projects/vleblanc_prj/genomes/hg38/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.cleaned.noalt.withHeader.vcf.bgz \
-L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
--native-pair-hmm-threads 32 \
-O ./results/exome/mutect2_multisamp/JK136/somatic.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK136/log

#Run FilterMutectCalls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
-V ./results/exome/mutect2_multisamp/JK136/somatic.vcf.gz \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-O ./results/exome/mutect2_multisamp/JK136/somatic.filter.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK136/filter_log

#Make vcf files with variants passing filter (bgziped and tabix-indexed)
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
./results/exome/mutect2_multisamp/JK136/somatic.filter.vcf.gz \
| bgzip > ./results/exome/mutect2_multisamp/JK136/passed.somatic.filter.vcf.gz

tabix ./results/exome/mutect2_multisamp/JK136/passed.somatic.filter.vcf.gz

sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ./results/exome/mutect2_multisamp/JK136 passed.somatic.filter






#-----------------------------------------------------------
## JK142
#-----------------------------------------------------------

#not including reg2 tissue since low tumour content

#Run mutect2 to get variant calls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-I ./data/exome_seq/bams/JK142_reg1_tissue.bam \
-I ./data/exome_seq/bams/JK142_reg1_organoid.bam \
-I ./data/exome_seq/bams/JK142_reg2_organoid.bam \
-I ./data/exome_seq/bams/JK196_tissue.bam \
-I ./data/exome_seq/bams/JK196_organoid.bam \
-I ./data/exome_seq/bams/JK142_JK196_blood.bam \
-normal JK142_JK196_blood \
--germline-resource /projects/vleblanc_prj/genomes/hg38/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.cleaned.noalt.withHeader.vcf.bgz \
-L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
--native-pair-hmm-threads 32 \
-O ./results/exome/mutect2_multisamp/JK142_noreg2tis/somatic.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK142_noreg2tis/log

#Run FilterMutectCalls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
-V ./results/exome/mutect2_multisamp/JK142_noreg2tis/somatic.vcf.gz \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-O ./results/exome/mutect2_multisamp/JK142_noreg2tis/somatic.filter.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK142_noreg2tis/filter_log

#Make vcf files with variants passing filter (bgziped and tabix-indexed)
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
./results/exome/mutect2_multisamp/JK142_noreg2tis/somatic.filter.vcf.gz \
| bgzip > ./results/exome/mutect2_multisamp/JK142_noreg2tis/passed.somatic.filter.vcf.gz

tabix ./results/exome/mutect2_multisamp/JK142_noreg2tis/passed.somatic.filter.vcf.gz

sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ./results/exome/mutect2_multisamp/JK142_noreg2tis passed.somatic.filter









#-----------------------------------------------------------
## JK153
#-----------------------------------------------------------

#Run mutect2 to get variant calls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-I ./data/exome_seq/bams/JK153_reg1_tissue.bam \
-I ./data/exome_seq/bams/JK153_reg1_organoid.bam \
-I ./data/exome_seq/bams/JK153_reg2_tissue.bam \
-I ./data/exome_seq/bams/JK153_reg2_organoid.bam \
-I ./data/exome_seq/bams/JK153_blood.bam \
-normal JK153_blood \
--germline-resource /projects/vleblanc_prj/genomes/hg38/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.cleaned.noalt.withHeader.vcf.bgz \
-L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
--native-pair-hmm-threads 32 \
-O ./results/exome/mutect2_multisamp/JK153/somatic.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK153/log

#Run FilterMutectCalls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
-V ./results/exome/mutect2_multisamp/JK153/somatic.vcf.gz \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-O ./results/exome/mutect2_multisamp/JK153/somatic.filter.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK153/filter_log

#Make vcf files with variants passing filter (bgziped and tabix-indexed)
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
./results/exome/mutect2_multisamp/JK153/somatic.filter.vcf.gz \
| bgzip > ./results/exome/mutect2_multisamp/JK153/passed.somatic.filter.vcf.gz

tabix ./results/exome/mutect2_multisamp/JK153/passed.somatic.filter.vcf.gz

sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ./results/exome/mutect2_multisamp/JK153 passed.somatic.filter





#-----------------------------------------------------------
## JK163
#-----------------------------------------------------------

#Run mutect2 to get variant calls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-I ./data/exome_seq/bams/JK163_reg1_tissue.bam \
-I ./data/exome_seq/bams/JK163_reg1_organoid.bam \
-I ./data/exome_seq/bams/JK163_reg2_tissue.bam \
-I ./data/exome_seq/bams/JK163_reg2_organoid.bam \
-I ./data/exome_seq/bams/JK163_blood.bam \
-normal JK163_blood \
--germline-resource /projects/vleblanc_prj/genomes/hg38/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.cleaned.noalt.withHeader.vcf.bgz \
-L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
--native-pair-hmm-threads 32 \
-O ./results/exome/mutect2_multisamp/JK163/somatic.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK163/log

#Run FilterMutectCalls
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
-V ./results/exome/mutect2_multisamp/JK163/somatic.vcf.gz \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-O ./results/exome/mutect2_multisamp/JK163/somatic.filter.vcf.gz \
2>&1 | tee ./results/exome/mutect2_multisamp/JK163/filter_log

#Make vcf files with variants passing filter (bgziped and tabix-indexed)
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
./results/exome/mutect2_multisamp/JK163/somatic.filter.vcf.gz \
| bgzip > ./results/exome/mutect2_multisamp/JK163/passed.somatic.filter.vcf.gz

tabix ./results/exome/mutect2_multisamp/JK163/passed.somatic.filter.vcf.gz

sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ./results/exome/mutect2_multisamp/JK163 passed.somatic.filter
