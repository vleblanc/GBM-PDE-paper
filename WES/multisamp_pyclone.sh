#!/bin/sh


#---------------JK124---------------

JK124_reg1_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200"
JK124_reg1_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200"
JK124_reg2_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200"
JK124_reg2_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK124_reg2_organoid/hg38_no_alt/EXC/B58195_B58200"

PyClone setup_analysis --in_files ${JK124_reg1_tis_dir}/pyclone/JK124_reg1_tissue_titan_subclonefil.tsv \
${JK124_reg1_org_dir}/pyclone/JK124_reg1_organoid_titan_subclonefil.tsv \
${JK124_reg2_tis_dir}/pyclone/JK124_reg2_tissue_titan_subclonefil.tsv \
${JK124_reg2_org_dir}/pyclone/JK124_reg2_organoid_titan_subclonefil.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/ \
--tumour_contents 0.98669 0.9445 0.8778 0.2113 \
--samples JK124_reg1_tis JK124_reg1_org JK124_reg2_tis JK124_reg2_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil/tables/loci.tsv \
--table_type loci \
--burnin 1000



#Not including reg2 org (looks like 0 tumour content)
PyClone setup_analysis --in_files ${JK124_reg1_tis_dir}/pyclone/JK124_reg1_tissue_titan_subclonefil_adfil_noreg2org.tsv \
${JK124_reg1_org_dir}/pyclone/JK124_reg1_organoid_titan_subclonefil_adfil_noreg2org.tsv \
${JK124_reg2_tis_dir}/pyclone/JK124_reg2_tissue_titan_subclonefil_adfil_noreg2org.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/ \
--tumour_contents 0.98669 0.9445 0.8778 \
--samples JK124_reg1_tis JK124_reg1_org JK124_reg2_tis \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/loci.tsv \
--table_type loci \
--burnin 1000

#run_citup_iter.py --submit local --maxjobs 10 JK124_subclonefil_adfil_cluster_freqs.txt iter_results.h5





#---------------JK136---------------

JK136_reg1_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201"
JK136_reg1_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201"
JK136_reg2_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201"
JK136_reg2_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201"
JK202_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK202_tissue/hg38_no_alt/EXC/B58180_B58201"
JK202_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK202_organoid/hg38_no_alt/EXC/B58181_B58201"

PyClone setup_analysis --in_files ${JK136_reg1_tis_dir}/pyclone/JK136_reg1_tissue_titan_subclonefil_adfil.tsv \
${JK136_reg1_org_dir}/pyclone/JK136_reg1_organoid_titan_subclonefil_adfil.tsv \
${JK136_reg2_tis_dir}/pyclone/JK136_reg2_tissue_titan_subclonefil_adfil.tsv \
${JK136_reg2_org_dir}/pyclone/JK136_reg2_organoid_titan_subclonefil_adfil.tsv \
${JK202_tis_dir}/pyclone/JK202_tissue_titan_subclonefil_adfil.tsv \
${JK202_org_dir}/pyclone/JK202_organoid_titan_subclonefil_adfil.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/ \
--tumour_contents 0.91299 0.24500 0.89060 0.70010 0.76860 0.86310 \
--samples JK136_reg1_tis JK136_reg1_org JK136_reg2_tis JK136_reg2_org JK202_tis JK202_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK136_subclonefil_adfil/tables/loci.tsv \
--table_type loci \
--burnin 1000

#run_citup_iter.py --submit local --maxjobs 10 JK136_subclonefil_adfil_cluster_freqs.txt iter_results.h5





#---------------JK142---------------

JK142_reg1_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202"
JK142_reg1_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202"
JK142_reg2_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK142_reg2_tissue/hg38_no_alt/EXC/B58184_B58202"
JK142_reg2_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202"
JK196_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK196_tissue/hg38_no_alt/EXC/B58186_B58202"
JK196_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK196_organoid/hg38_no_alt/EXC/B58187_B58202"

PyClone setup_analysis --in_files ${JK142_reg1_tis_dir}/pyclone/JK142_reg1_tissue_titan_subclonefil.tsv \
${JK142_reg1_org_dir}/pyclone/JK142_reg1_organoid_titan_subclonefil.tsv \
${JK142_reg2_tis_dir}/pyclone/JK142_reg2_tissue_titan_subclonefil.tsv \
${JK142_reg2_org_dir}/pyclone/JK142_reg2_organoid_titan_subclonefil.tsv \
${JK196_tis_dir}/pyclone/JK196_tissue_titan_subclonefil.tsv \
${JK196_org_dir}/pyclone/JK196_organoid_titan_subclonefil.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/ \
--tumour_contents 0.5599 0.94911 0.0028 0.96746 0.7601 0.92177 \
--samples JK142_reg1_tis JK142_reg1_org JK142_reg2_tis JK142_reg2_org JK196_tis JK196_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil/tables/loci.tsv \
--table_type loci \
--burnin 1000


#Not including JK142_reg2_tis (very low tumour content)
PyClone setup_analysis --in_files ${JK142_reg1_tis_dir}/pyclone/JK142_reg1_tissue_titan_subclonefil_adfil_noreg2tis.tsv \
${JK142_reg1_org_dir}/pyclone/JK142_reg1_organoid_titan_subclonefil_adfil_noreg2tis.tsv \
${JK142_reg2_org_dir}/pyclone/JK142_reg2_organoid_titan_subclonefil_adfil_noreg2tis.tsv \
${JK196_tis_dir}/pyclone/JK196_tissue_titan_subclonefil_adfil_noreg2tis.tsv \
${JK196_org_dir}/pyclone/JK196_organoid_titan_subclonefil_adfil_noreg2tis.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/ \
--tumour_contents 0.5599 0.94911 0.96746 0.7601 0.92177 \
--samples JK142_reg1_tis JK142_reg1_org JK142_reg2_org JK196_tis JK196_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables/loci.tsv \
--table_type loci \
--burnin 1000

#run_citup_iter.py --max_nodes 13 --submit local --maxjobs 20 JK142_subclonefil_adfil_noreg2tis_cluster_freqs.txt iter_results.h5





#---------------JK153---------------

JK153_reg1_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203"
JK153_reg1_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203"
JK153_reg2_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203"
JK153_reg2_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203"

PyClone setup_analysis --in_files ${JK153_reg1_tis_dir}/pyclone/JK153_reg1_tissue_titan_subclonefil_adfil.tsv \
${JK153_reg1_org_dir}/pyclone/JK153_reg1_organoid_titan_subclonefil_adfil.tsv \
${JK153_reg2_tis_dir}/pyclone/JK153_reg2_tissue_titan_subclonefil_adfil.tsv \
${JK153_reg2_org_dir}/pyclone/JK153_reg2_organoid_titan_subclonefil_adfil.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/ \
--tumour_contents 0.91689 0.94070 0.93361 0.79530 \
--samples JK153_reg1_tis JK153_reg1_org JK153_reg2_tis JK153_reg2_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK153_subclonefil_adfil/tables/loci.tsv \
--table_type loci \
--burnin 1000

#run_citup_iter.py --max_nodes 12 --submit local --maxjobs 20 JK153_subclonefil_adfil_cluster_freqs.txt iter_results.h5





#---------------JK163---------------

JK163_reg1_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204"
JK163_reg1_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204"
JK163_reg2_tis_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204"
JK163_reg2_org_dir="/projects/vleblanc_prj/GBM_organoids/data/exome_seq/JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204"

PyClone setup_analysis --in_files ${JK163_reg1_tis_dir}/pyclone/JK163_reg1_tissue_titan_subclonefil_adfil.tsv \
${JK163_reg1_org_dir}/pyclone/JK163_reg1_organoid_titan_subclonefil_adfil.tsv \
${JK163_reg2_tis_dir}/pyclone/JK163_reg2_tissue_titan_subclonefil_adfil.tsv \
${JK163_reg2_org_dir}/pyclone/JK163_reg2_organoid_titan_subclonefil_adfil.tsv \
--working_dir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/ \
--tumour_contents 0.68440 0.95910 0.85960 0.93565 \
--samples JK163_reg1_tis JK163_reg1_org JK163_reg2_tis JK163_reg2_org \
--init_method connected

PyClone run_analysis --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/config.yaml \
--seed 12345

mkdir /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/tables

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/tables/cluster.tsv \
--table_type cluster \
--burnin 1000

PyClone build_table --config_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/config.yaml \
--out_file /projects/vleblanc_prj/GBM_organoids/results/exome/pyclone/JK163_subclonefil_adfil/tables/loci.tsv \
--table_type loci \
--burnin 1000
