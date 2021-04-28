# ---------------------------------------------------------------
## Scripts to load 10X copy number output data into a combined object
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_load_clean.R")

#Load 10X copy number
outs_dir <- function(sample) paste0("./data/cellranger_dna/", sample, "/outs")

all_cnv <- load.10x.cnv(id = "200106_10x_scDNA_full",
                        dirs = c(JK136_reg1_tis = outs_dir("SI-GA-A3"), JK136_reg1_org = outs_dir("SI-GA-B3"),
                                 JK136_reg2_tis = outs_dir("SI-GA-C3"), JK136_reg2_org = outs_dir("SI-GA-D3"),
                                 JK202_reg1_tis = outs_dir("SI-GA-E3"), JK202_reg1_org = outs_dir("SI-GA-F3"),
                                 JK142_reg1_tis = outs_dir("SI-GA-A4"), JK142_reg1_org = outs_dir("SI-GA-B4"),
                                 JK142_reg2_tis = outs_dir("SI-GA-C4"), JK142_reg2_org = outs_dir("SI-GA-D4"),
                                 JK196_reg1_tis = outs_dir("SI-GA-E4"), JK196_reg1_org = outs_dir("SI-GA-F4"),
                                 JK153_reg1_tis = outs_dir("SI-GA-G3"), JK153_reg1_org = outs_dir("SI-GA-H3"),
                                 JK153_reg2_tis = outs_dir("SI-GA-G4"), JK153_reg2_org = outs_dir("SI-GA-H4")),
                        samp_levels = c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org",
                                        "JK202_reg1_tis", "JK202_reg1_org",
                                        "JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_tis", "JK142_reg2_org",
                                        "JK196_reg1_tis", "JK196_reg1_org",
                                        "JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org"),
                        tum_levels = c("JK136", "JK202", "JK142", "JK196", "JK153"),
                        merged_calls = T, filter_cnvs = F,
                        resume = F, n_cores = 20)


