#!/bin/sh

tum_list=()
while IFS= read -d $'\0' -r file ; do
  tum_list+=("$file")
done < <(find /projects/vleblanc_prj/GBM_organoids/data/exome_seq/ -maxdepth 4 -path "*/B*_B*" -print0)

for dir in "${tum_list[@]}"; do
  #If variant file exists, don't run 
  if [ -d "$dir/rtgtools" ]; then
    printf "$dir/rtgtools already exists\n"

  else
    tum=`echo "$dir" | cut -d/ -f7`

    /projects/vleblanc_prj/tools/rtg-tools-3.10.1/rtg RTG_MEM=16G vcfeval \
    -b ${dir}/mutect2/somatic.filter.vcf.gz \
    -c ${dir}/strelka2/results/variants/somatic.vcf.gz \
    --sample ${tum},ALT \
    -e /projects/vleblanc_prj/genomes/hg38/exome/Exome-IDT-xGen-hg38-targets_merged_sorted.bed \
    -o ${dir}/rtgtools \
    -t /projects/vleblanc_prj/genomes/hg38/genome/hg38_no_alt_sdf \
    --squash-ploidy \
    --vcf-score-field QUAL \
    -T 10
  fi
done
