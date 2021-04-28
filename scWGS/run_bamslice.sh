#!/bin/sh

# ---------------------------------------------------------------
## Scripts to run cellranger-dna bamslice on cells from separate clones within each sample
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

#JK136

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK136_reg1_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-A3/outs/possorted_bam.bam \
--localcores=32 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK136_reg1_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-B3/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK136_reg2_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-C3/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK136_reg2_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-D3/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK202_reg1_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-E3/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK202_reg1_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-F3/outs/possorted_bam.bam \
--localcores=24 --localmem=128




#JK142

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK142_reg1_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-A4/outs/possorted_bam.bam \
--localcores=32 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK142_reg1_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-B4/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK142_reg2_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-C4/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK142_reg2_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-D4/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK196_reg1_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-E4/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK196_reg1_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-F4/outs/possorted_bam.bam \
--localcores=24 --localmem=128




#JK153

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK153_reg1_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-G3/outs/possorted_bam.bam \
--localcores=32 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK153_reg1_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-H3/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK153_reg2_tis
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-G4/outs/possorted_bam.bam \
--localcores=24 --localmem=128

cd /projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/JK153_reg2_org
cellranger-dna bamslice \
--id=clones \
--csv=clones.csv \
--bam=./SI-GA-H4/outs/possorted_bam.bam \
--localcores=24 --localmem=128
