#!/bin/sh

#Get list of normals
norm_list=()
while IFS= read -d $'\0' -r file ; do
  norm_list+=("$file")
done < <(find /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK*_blood/ -type d -name "haplotypecaller" -print0)

#Run HaplotypeCaller on each
for hcaller_dir in "${norm_list[@]}"; do
  #If variant file exists, don't run HaplotypeCaller
  if [ -f "$hcaller_dir/output.g.vcf.gz" ]; then
    printf "$hcaller_dir/output.g.vcf.gz already exists, not running HaplotypeCaller\n"

  else
    samp=`echo "$hcaller_dir" | cut -d/ -f7`

    #Run HaplotypeCaller to get variant calls
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller \
    -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
    -I bams/${samp}.bam \
    -L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
    -ERC GVCF \
    --native-pair-hmm-threads 24 \
    -O "$hcaller_dir/output.g.vcf.gz" \
    2>&1 | tee "$hcaller_dir/log"
  fi
done



#Get list of tumours
tum_list=()
while IFS= read -d $'\0' -r file ; do
  tum_list+=("$file")
done < <(find /projects/vleblanc_prj/GBM_organoids/data/exome_seq/ -maxdepth 5 -path "*/B*_B*/haplotypecaller" -print0)

#Run HaplotypeCaller on each
for hcaller_dir in "${tum_list[@]}"; do
  tum=`echo "$hcaller_dir" | cut -d/ -f7`

  printf "$tum"

  if echo "$tum" | grep JK124; then
    norm="JK124_blood"
  elif echo "$tum" | grep -e JK136 -e JK202; then
    norm="JK136_JK202_blood"
  elif echo "$tum" | grep -e JK142 -e JK196; then
    norm="JK142_JK196_blood"
  elif echo "$tum" | grep JK153; then
    norm="JK153_blood"
  elif echo "$tum" | grep JK163; then
    norm="JK163_blood"
  fi

  #If variant file exists, don't run HaplotypeCaller
  if [ -f "$hcaller_dir/output.g.vcf.gz" ]; then
    printf "$hcaller_dir/output.g.vcf.gz already exists, not running HaplotypeCaller\n"

    #If GVCFs have already been combined, don't run CombineGVCFs
    if [ -f "$hcaller_dir/tum_norm_pair.g.vcf.gz" ]; then
      printf "$hcaller_dir/tum_norm_pair.g.vcf.gz already exists, not running CombineGVCFs\n"

    else
      #Combine tumour and normal GVCF
      /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CombineGVCFs \
      -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
      --variant "$norm/haplotypecaller/output.g.vcf.gz" \
      --variant "$hcaller_dir/output.g.vcf.gz" \
      -O "$hcaller_dir/tum_norm_pair.g.vcf.gz" \
      2>&1 | tee -a "$hcaller_dir/log"
    fi

    #If GVCFs have already been genotyped, don't run GenotypeGVCFs
    if [ -f "$hcaller_dir/output.tum_norm_pair.g.vcf.gz" ]; then
      printf "$hcaller_dir/output.tum_norm_pair.g.vcf.gz already exists, not running GenotypeGVCFs\n"

    else
      #Genotype GVCFs
      /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar GenotypeGVCFs \
      -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
      -V "$hcaller_dir/tum_norm_pair.g.vcf.gz" \
      -O "$hcaller_dir/output.tum_norm_pair.g.vcf.gz"
    fi

  else
    #Run HaplotypeCaller to get variant calls
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller \
    -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
    -I bams/${tum}.bam \
    -L /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
    -ERC GVCF \
    --native-pair-hmm-threads 24 \
    -O "$hcaller_dir/output.g.vcf.gz" \
    2>&1 | tee "$hcaller_dir/log"

    #Combine tumour and normal GVCF
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CombineGVCFs \
    -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
    --variant "$norm/haplotypecaller/output.g.vcf.gz" \
    --variant "$hcaller_dir/output.g.vcf.gz" \
    -O "$hcaller_dir/tum_norm_pair.g.vcf.gz" \
    2>&1 | tee -a "$hcaller_dir/log"

    #Genotype GVCFs
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar GenotypeGVCFs \
    -R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
    -V "$hcaller_dir/tum_norm_pair.g.vcf.gz" \
    -O "$hcaller_dir/output.tum_norm_pair.g.vcf.gz"

    #Separate SNPs and indels for filtering
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
    -V "$hcaller_dir/output.tum_norm_pair.g.vcf.gz" \
    -select-type SNP \
    -O "$hcaller_dir/output.tum_norm_pair.g.snps.vcf.gz"

    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
    -V "$hcaller_dir/output.tum_norm_pair.g.vcf.gz" \
    -select-type INDEL \
    -O "$hcaller_dir/output.tum_norm_pair.g.indels.vcf.gz"

    #Filter SNPs (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216)
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration \
    -V "$hcaller_dir/output.tum_norm_pair.g.snps.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "$hcaller_dir/output.tum_norm_pair.g.snps.filtered.vcf.gz"

    #Filter indels (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216)
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration \
    -V "$hcaller_dir/output.tum_norm_pair.g.indels.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "$hcaller_dir/output.tum_norm_pair.g.indels.filtered.vcf.gz"

    #Cobmine filtered SNPs and indels
    /projects/vleblanc_prj/tools/bcftools-1.9/bcftools concat -a -Oz \
    -o "$hcaller_dir/output.tum_norm_pair.g.filtered.vcf.gz" \
    ${hcaller_dir}/output.tum_norm_pair.g.snps.filtered.vcf.gz \
    ${hcaller_dir}/output.tum_norm_pair.g.indels.filtered.vcf.gz

    tabix "$hcaller_dir/output.tum_norm_pair.g.filtered.vcf.gz"

    #Make file with only passing variants
    /gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
    -V "$hcaller_dir/output.tum_norm_pair.g.filtered.vcf.gz" \
    -select 'vc.isNotFiltered()' \
    -O "$hcaller_dir/passed.output.tum_norm_pair.g.filtered.vcf.gz"

  fi
done







#Genotype GVCFs
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar GenotypeGVCFs \
-R /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt.fa \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/output.g.vcf.gz \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.g.vcf.gz

#Separate SNPs and indels for filtering
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.g.vcf.gz \
-select-type SNP \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.snps.g.vcf.gz

/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.g.vcf.gz \
-select-type INDEL \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.indels.g.vcf.gz

#Filter SNPs (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216)
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.snps.g.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.snps.filtered.g.vcf.gz

#Filter indels (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216)
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.indels.g.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.indels.filtered.g.vcf.gz

#Cobmine filtered SNPs and indels
/projects/vleblanc_prj/tools/bcftools-1.9/bcftools concat -a -Oz \
-o /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.filtered.g.vcf.gz \
/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.snps.filtered.g.vcf.gz \
/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.indels.filtered.g.vcf.gz

tabix /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.filtered.g.vcf.gz

#Make file with only passing variants
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/genotype.output.filtered.g.vcf.gz \
-select 'vc.isNotFiltered()' \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.genotype.output.filtered.g.vcf.gz

#Make file with only heterozygous variants
/gsc/software/linux-x86_64/jre1.8.0_66/bin/java -Xmx128g -jar /projects/vleblanc_prj/tools/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
-V /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.genotype.output.filtered.g.vcf.gz \
-select 'vc.getGenotype("JK153_blood").isHet()' \
-O /projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_blood/haplotypecaller/passed.het.genotype.output.filtered.g.vcf.gz
