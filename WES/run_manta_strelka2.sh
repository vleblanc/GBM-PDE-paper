#!/bin/sh

norm_bam=$1
tum_bam=$2
strelka_dir=$3
n_cores=$4
mem=$5

#Run manta
/projects/vleblanc_prj/tools/manta-1.6.0.centos6_x86_64/bin/configManta.py \
--normalBam ${norm_bam} \
--tumorBam ${tum_bam} \
--referenceFasta /projects/alignment_references/9606/hg38_no_alt/genome/fasta/hg38_no_alt.fa \
--exome \
--callRegions /projects/vleblanc_prj/genomes/Exome-IDT-xGen-hg38-targets_merged_sorted.bed.gz \
--runDir ${strelka_dir}/manta

${strelka_dir}/manta/runWorkflow.py -j ${n_cores} -g ${mem}


#Run strelka
/projects/vleblanc_prj/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${norm_bam} \
--tumorBam ${tum_bam} \
--referenceFasta /projects/alignment_references/9606/hg38_no_alt/genome/fasta/hg38_no_alt.fa \
--indelCandidates ${strelka_dir}/manta/results/variants/candidateSmallIndels.vcf.gz \
--exome \
--callRegions /projects/vleblanc_prj/genomes/Exome-IDT-xGen-hg38-targets_merged_sorted.bed.gz \
--runDir ${strelka_dir}

${strelka_dir}/runWorkflow.py -m local -j ${n_cores} -g ${mem}



#Make vcf files with variants passing EVS filter (bgziped and tabix-indexed)
#SNVs
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
${strelka_dir}/results/variants/somatic.snvs.vcf.gz \
| bgzip > ${strelka_dir}/results/variants/passed.somatic.snvs.vcf.gz

tabix ${strelka_dir}/results/variants/passed.somatic.snvs.vcf.gz

#Indels
/gsc/software/linux-x86_64/jre1.7.0_03/bin/java -Xmx16g -jar /gsc/software/linux-x86/snpEff-4.1/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" \
${strelka_dir}/results/variants/somatic.indels.vcf.gz \
| bgzip > ${strelka_dir}/results/variants/passed.somatic.indels.vcf.gz

tabix ${strelka_dir}/results/variants/passed.somatic.indels.vcf.gz



#Annotate vcf files
sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ${strelka_dir}/results/variants somatic.snvs
sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ${strelka_dir}/results/variants passed.somatic.snvs
sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ${strelka_dir}/results/variants somatic.indels
sh /projects/vleblanc_prj/GBM_organoids/code/utils/annotate_vcf.sh ${strelka_dir}/results/variants passed.somatic.indels
