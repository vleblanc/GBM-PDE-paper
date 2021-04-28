#!/bin/sh

dir=$1
name=$2

#Get effects
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/snpEff.jar \
-v \
-stats ${dir}/${name}.eff.stats.html \
-dataDir /gsc/software/linux-x86/snpEff-4.1/ref_data \
GRCh38.79 \
${dir}/${name}.vcf.gz \
| bgzip > ${dir}/${name}.eff.vcf.gz

/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar annotate \
/gsc/resources/annotation/dbsnp/hg38_no_alt/v149/dbSNP_v149.vcf.gz \
${dir}/${name}.eff.vcf \
| /gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar annotate \
/gsc/QA-bio/bioapps/ref_data/annotation_vcfs/hg38/cosmic_v79_coding.vcf \
| /gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar annotate \
/gsc/resources/annotation/clinvar/hg38_no_alt/20170104/clinvar_20170104.vcf.gz \
| bgzip > ${dir}/${name}.eff.dbSNP_v149.cosmic_v79.clinvar_20170104.annotations.vcf.gz

tabix ${dir}/${name}.eff.dbSNP_v149.cosmic_v79.clinvar_20170104.annotations.vcf.gz
