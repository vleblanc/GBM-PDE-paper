#!/bin/sh


#JK124
sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK124_blood.bam \
./bams/JK124_reg1_tis.bam \
./JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200/strelka2 \
30 256 2>&1 | tee JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK124_blood.bam \
./bams/JK124_reg1_org.bam \
./JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200/strelka2 \
30 256 2>&1 | tee JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK124_blood.bam \
./bams/JK124_reg2_tis.bam \
./JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200/strelka2 \
30 256 2>&1 | tee JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK124_blood.bam \
./bams/JK124_reg2_org.bam \
./JK124_reg2_organoid/hg38_no_alt/EXC/B58195_B58200/strelka2 \
30 256 2>&1 | tee JK124_reg2_organoid/hg38_no_alt/EXC/B58195_B58200/strelka2/log



#JK136
sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK136_reg1_tis.bam \
./JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201/strelka2 \
30 256 2>&1 | tee JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK136_reg1_org.bam \
./JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201/strelka2 \
30 256 2>&1 | tee JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK136_reg2_tis.bam \
./JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201/strelka2 \
30 256 2>&1 | tee JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK136_reg2_org.bam \
./JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201/strelka2 \
30 256 2>&1 | tee JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK202_tis.bam \
./JK202_tissue/hg38_no_alt/EXC/B58180_B58201/strelka2 \
30 256 2>&1 | tee JK202_tissue/hg38_no_alt/EXC/B58180_B58201/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK136_JK202_blood.bam \
./bams/JK202_org.bam \
./JK202_organoid/hg38_no_alt/EXC/B58181_B58201/strelka2 \
30 256 2>&1 | tee JK202_organoid/hg38_no_alt/EXC/B58181_B58201/strelka2/log



#JK142
sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK142_reg1_tis.bam \
./JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202/strelka2 \
30 256 2>&1 | tee JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK142_reg1_org.bam \
./JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202/strelka2 \
30 256 2>&1 | tee JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK142_reg2_tis.bam \
./JK142_reg2_tissue/hg38_no_alt/EXC/B58184_B58202/strelka2 \
30 256 2>&1 | tee JK142_reg2_tissue/hg38_no_alt/EXC/B58184_B58202/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK142_reg2_org.bam \
./JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202/strelka2 \
30 256 2>&1 | tee JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK196_tis.bam \
./JK196_tissue/hg38_no_alt/EXC/B58186_B58202/strelka2 \
30 256 2>&1 | tee JK196_tissue/hg38_no_alt/EXC/B58186_B58202/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK142_JK196_blood.bam \
./bams/JK196_org.bam \
./JK196_organoid/hg38_no_alt/EXC/B58187_B58202/strelka2 \
30 256 2>&1 | tee JK196_organoid/hg38_no_alt/EXC/B58187_B58202/strelka2/log



#JK153
sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK153_blood.bam \
./bams/JK153_reg1_tis.bam \
./JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203/strelka2 \
30 256 2>&1 | tee JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK153_blood.bam \
./bams/JK153_reg1_org.bam \
./JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203/strelka2 \
30 256 2>&1 | tee JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK153_blood.bam \
./bams/JK153_reg2_tis.bam \
./JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203/strelka2 \
30 256 2>&1 | tee JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK153_blood.bam \
./bams/JK153_reg2_org.bam \
./JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203/strelka2 \
30 256 2>&1 | tee JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203/strelka2/log



#JK124
sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK163_blood.bam \
./bams/JK163_reg1_tis.bam \
./JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204/strelka2 \
30 256 2>&1 | tee JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK163_blood.bam \
./bams/JK163_reg1_org.bam \
./JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204/strelka2 \
30 256 2>&1 | tee JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK163_blood.bam \
./bams/JK163_reg2_tis.bam \
./JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204/strelka2 \
30 256 2>&1 | tee JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204/strelka2/log

sh /projects/vleblanc_prj/GBM_organoids/code/analysis/run_manta_strelka2.sh \
./bams/JK163_blood.bam \
./bams/JK163_reg2_org.bam \
./JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204/strelka2 \
30 256 2>&1 | tee JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204/strelka2/log

