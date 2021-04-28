
#JK124
JK124_vars <- get.var.info.from.table("./results/exome/common_mutations/JK124_noreg2org_strelka_mutect/tp-baseline.table")

make.pyclone.tsv(var_table = JK124_vars$tbl_nosnp[apply(JK124_vars$tbl_nosnp[,grepl(".DP", colnames(JK124_vars$tbl_nosnp))], 1, function(x) x[which.min(x)]) >= 20 &
                                                    apply(JK124_vars$tbl_nosnp[,grepl("_alt_ad", colnames(JK124_vars$tbl_nosnp))], 1, function(x) x[which.max(x)]) >= 10,], 
                 cn_met = "titan",
                 dat_dirs = c("JK124_reg1_tissue" = "./data/exome_seq/JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200/",
                              "JK124_reg1_organoid" = "./data/exome_seq/JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200/",
                              "JK124_reg2_tissue" = "./data/exome_seq/JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200/"),
                 titan_path = "./results/exome/titanCNA/JK124/results/titan/hmm/optimalClusterSolution.txt",
                 sex = "female", filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = "_subclonefil_adfil_noreg2org", return_granges_list = FALSE)

#JK136
JK136_vars <- get.var.info.from.table("./results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.table")

make.pyclone.tsv(var_table = JK136_vars$tbl_nosnp[apply(JK136_vars$tbl_nosnp[,grepl(".DP", colnames(JK136_vars$tbl_nosnp))], 1, function(x) x[which.min(x)]) >= 20 &
                                                    apply(JK136_vars$tbl_nosnp[,grepl("_alt_ad", colnames(JK136_vars$tbl_nosnp))], 1, function(x) x[which.max(x)]) >= 10,], 
                 cn_met = "titan",
                 dat_dirs = c("JK136_reg1_tissue" = "./data/exome_seq/JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201/",
                              "JK136_reg1_organoid" = "./data/exome_seq/JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201/",
                              "JK136_reg2_tissue" = "./data/exome_seq/JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201/",
                              "JK136_reg2_organoid" = "./data/exome_seq/JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201/",
                              "JK202_tissue" = "./data/exome_seq/JK202_tissue/hg38_no_alt/EXC/B58180_B58201/",
                              "JK202_organoid" = "./data/exome_seq/JK202_organoid/hg38_no_alt/EXC/B58181_B58201/"),
                 titan_path = "./results/exome/titanCNA/JK136/results/titan/hmm/optimalClusterSolution.txt",
                 sex = "female", filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = "_subclonefil_adfil", return_granges_list = FALSE)


#JK142
JK142_vars <- get.var.info.from.table("./results/exome/common_mutations/JK142_JK196_noreg2tis_strelka_mutect/tp-baseline.table")

make.pyclone.tsv(var_table = JK142_vars$tbl_nosnp[apply(JK142_vars$tbl_nosnp[,grepl(".DP", colnames(JK142_vars$tbl_nosnp))], 1, function(x) x[which.min(x)]) >= 20 &
                                                    apply(JK142_vars$tbl_nosnp[,grepl("_alt_ad", colnames(JK142_vars$tbl_nosnp))], 1, function(x) x[which.max(x)]) >= 10,], 
                 cn_met = "titan",
                 dat_dirs = c("JK142_reg1_tissue" = "./data/exome_seq/JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202/",
                              "JK142_reg1_organoid" = "./data/exome_seq/JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202/",
                              "JK142_reg2_tissue" = "./data/exome_seq/JK142_reg2_tissue/hg38_no_alt/EXC/B58184_B58202/",
                              "JK142_reg2_organoid" = "./data/exome_seq/JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202/",
                              "JK196_tissue" = "./data/exome_seq/JK196_tissue/hg38_no_alt/EXC/B58186_B58202/",
                              "JK196_organoid" = "./data/exome_seq/JK196_organoid/hg38_no_alt/EXC/B58187_B58202/"),
                 titan_path = "./results/exome/titanCNA/JK142/results/titan/hmm/optimalClusterSolution.txt",
                 sex = "female", filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = "_subclonefil_adfil_noreg2tis", return_granges_list = FALSE)



#JK153
JK153_vars <- get.var.info.from.table("./results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.table")

make.pyclone.tsv(var_table = JK153_vars$tbl_nosnp[apply(JK153_vars$tbl_nosnp[,grepl(".DP", colnames(JK153_vars$tbl_nosnp))], 1, function(x) x[which.min(x)]) >= 20 &
                                                    apply(JK153_vars$tbl_nosnp[,grepl("_alt_ad", colnames(JK153_vars$tbl_nosnp))], 1, function(x) x[which.max(x)]) >= 10,], 
                 cn_met = "titan",
                 dat_dirs = c("JK153_reg1_tissue" = "./data/exome_seq/JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203/",
                              "JK153_reg1_organoid" = "./data/exome_seq/JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203/",
                              "JK153_reg2_tissue" = "./data/exome_seq/JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203/",
                              "JK153_reg2_organoid" = "./data/exome_seq/JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203/"),
                 titan_path = "./results/exome/titanCNA/JK153/results/titan/hmm/optimalClusterSolution.txt",
                 sex = "female", filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = "_subclonefil_adfil", return_granges_list = FALSE)



#JK163
JK163_vars <- get.var.info.from.table("./results/exome/common_mutations/JK163_strelka_mutect/tp-baseline.table")

JK163_tbl <- JK163_vars$tbl_nosnp[apply(JK163_vars$tbl_nosnp[,grepl(".DP", colnames(JK163_vars$tbl_nosnp))], 1, function(x) x[which.min(x)]) >= 20 &
                                    apply(JK163_vars$tbl_nosnp[,grepl("_alt_ad", colnames(JK163_vars$tbl_nosnp))], 1, function(x) x[which.max(x)]) >= 10,]
JK163_tbl <- JK163_tbl[JK163_tbl$JK163_reg2_organoid_alt_ad >= 50,]
make.pyclone.tsv(var_table = JK163_tbl, 
                 cn_met = "titan",
                 dat_dirs = c("JK163_reg1_tissue" = "./data/exome_seq/JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204/",
                              "JK163_reg1_organoid" = "./data/exome_seq/JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204/",
                              "JK163_reg2_tissue" = "./data/exome_seq/JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204/",
                              "JK163_reg2_organoid" = "./data/exome_seq/JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204/"),
                 titan_path = "./results/exome/titanCNA/JK163/results/titan/hmm/optimalClusterSolution.txt",
                 sex = "female", filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = "_subclonefil_adfil50org2", return_granges_list = FALSE)








make.pyclone.tsv <- function(var_dir = NULL, var_table = NULL, cn_met = c("cnvkit", "titan"), dat_dirs, titan_path = NULL, 
                             sex = c("female", "male"), filter_subclonal_cn = TRUE, write_table = TRUE, out_file_suffix = NULL, 
                             return_granges_list = TRUE) {
  #Need to supply either var_dir (will load variants from vcf file) or var_table
  require(vcfR)
  
  #Get variants
  if(!is.null(var_dir)) {
    vars_vcf <- read.vcfR(paste0(var_dir, "tp-baseline.eff.vcf.gz"), verbose = FALSE)
    
    vars_gr <- GRanges(seqnames = vars_vcf@fix[,1], 
                       ranges = IRanges(start = as.numeric(vars_vcf@fix[,2]), 
                                        end = as.numeric(vars_vcf@fix[,2])),
                       strand = "*",
                       gene = sapply(strsplit(sapply(strsplit(vars_vcf@fix[,"INFO"], ";"), function(x) x[grepl("ANN=", x)]), "[|]"), "[", 4),
                       ref = vars_vcf@fix[,4],
                       alt = vars_vcf@fix[,5])
    
    #Get ref/alt allele depths for all relevant samples
    vars_list <- setNames(lapply(names(dat_dirs), function(samp) {
      samp_refd <- sapply(strsplit(sapply(strsplit(vars_vcf@gt[,samp], ":"), "[", 2), ","), "[", 1)
      samp_altd <- sapply(strsplit(sapply(strsplit(vars_vcf@gt[,samp], ":"), "[", 2), ","), "[", 2)
      
      samp_vars_gr <- vars_gr
      mcols(samp_vars_gr) <- cbind(mcols(samp_vars_gr), ref_depth = samp_refd, alt_depth = samp_altd)
      
      return(samp_vars_gr)
    }), names(dat_dirs))
    
  } else {
    vars_gr <- GRanges(seqnames = var_table$CHROM,
                       ranges = IRanges(start = var_table$POS,
                                        end = var_table$POS),
                       strand = "*",
                       gene = var_table$gene,
                       ref = var_table$REF,
                       alt = var_table$ALT,
                       type = var_table$TYPE)
    
    #Get ref/alt allele depths for all relevant samples
    vars_list <- setNames(lapply(names(dat_dirs), function(samp) {
      samp_refd <- var_table[,paste0(samp, "_ref_ad")]
      samp_altd <- var_table[,paste0(samp, "_alt_ad")]
      
      samp_vars_gr <- vars_gr
      mcols(samp_vars_gr) <- cbind(mcols(samp_vars_gr), ref_depth = samp_refd, alt_depth = samp_altd)
      
      return(samp_vars_gr)
    }), names(dat_dirs))
    
  }

  
  #Get copy number tables
  if(cn_met == "cnvkit") {
    cn_list <- setNames(lapply(names(dat_dirs), function(samp) {
      cn <- read.delim(paste0(dat_dirs[samp], "/CNVkit/", samp, ".call.cns"))
      cn_gr <- makeGRangesFromDataFrame(cn[,c("chromosome", "start", "end", "cn1", "cn2")], 
                                        keep.extra.columns = T)
      colnames(mcols(cn_gr)) <- c("major_cn", "minor_cn")
    }), names(dat_dirs))
    
  } else if(cn_met == "titan") {
    titan_optimal <- read.delim(titan_path, stringsAsFactors = FALSE)
    titan_segs <- setNames(titan_optimal$path, sapply(strsplit(titan_optimal$id, "_"), function(x) paste0(x[1:(length(x) - 1)], collapse = "_")))
    
    tum <- strsplit(names(dat_dirs)[1], "_")[[1]][1]
    
    cn_list <- setNames(lapply(names(titan_segs), function(samp) {
      cn <- read.delim(paste0("./results/exome/titanCNA/", tum, "/", titan_segs[samp], ".segs.txt"))[,c("Chromosome", "Start_Position.bp.", 
                                                                                                        "End_Position.bp.", "MinorCN", "MajorCN",
                                                                                                        "Clonal_Cluster", "Cellular_Prevalence")]
      colnames(cn) <- c("chromosome", "start", "end", "minor_cn", "major_cn", "clonal_cluster", "cell_prev")
      cn_gr <- makeGRangesFromDataFrame(cn, keep.extra.columns = T)
    }), names(titan_segs))
    
    #If there is a sample missing from titan, remove variants
    vars_list <- vars_list[names(titan_segs)]
    
  } else {
    stop("cn_met must be one of \"cnvkit\" or \"titan\"")
  }
  
  #Get sample-specific info for each variant
  res_list <- setNames(lapply(names(vars_list), function(samp) {
    res <- vars_list[[samp]]
    
    #Get copy number info
    cn <- cn_list[[samp]]
    
    #Get overlap between variants and copy number objects
    cn_ol <- findOverlaps(cn, res)
    
    #Keep only variants that have an overlapping segment
    res_out <- res[subjectHits(cn_ol)]
    
    #Add copy number info to variants
    mcols(res_out) <- cbind(mcols(res)[subjectHits(cn_ol),], mcols(cn)[queryHits(cn_ol),])
    
    return(res_out)
  }), names(vars_list))
  
  #Make dataframe with required information for tsv output file
  out_list <- lapply(res_list, function(res) {
    if(sex == "female") normal_cn <- ifelse(as.character(seqnames(res)) == "chrY", 0, 2)
    if(sex == "male") normal_cn <- ifelse(as.character(seqnames(res)) %in% c("chrX", "chrY"), 1, 2)
    
    out <- data.frame(mutation_id = paste0(mcols(res)$gene, "_", as.character(seqnames(res)), ":", start(res)),
                      ref_counts = res$ref_depth,
                      var_counts = res$alt_depth,
                      normal_cn = normal_cn,
                      minor_cn = res$minor_cn,
                      major_cn = res$major_cn)
    
    #Add additional metadata columns if applicable
    if(!is.null(var_table)) {
      out <- cbind(out, mcols(res)[,c("type"), drop = FALSE])
    }
    
    if(cn_met == "titan") {
      out <- cbind(out, mcols(res)[,c("clonal_cluster", "cell_prev")])
      
      #Filter variants in regions with subclonal copy number alterations if requested
      if(filter_subclonal_cn) {
        #Only remove if copy number subclone is >25% or <75% (otherwise shouldn't have much of an impact, according to Andrew Roth)
        out <- out[c(which(is.na(out$cell_prev)), which(out$cell_prev < 0.25 | out$cell_prev > 0.75)),]
      }
    }
    
    #Remove variants with no or nonsensical major copy number
    out <- out[!is.na(out$major_cn) & out$major_cn != 0,]
  })
  
  if(write_table) {
    lapply(names(vars_list), function(samp) {
      if(cn_met == "cnvkit") {
        write.table(out_list[[samp]], 
                    file = paste0(dat_dirs[samp], "/pyclone/", samp, "_cnvkit", out_file_suffix, ".tsv"),
                    sep = "\t", quote = F, row.names = F)
        
      } else if(cn_met == "titan") {
        write.table(out_list[[samp]], 
                    file = paste0(dat_dirs[samp], "/pyclone/", samp, "_titan", out_file_suffix, ".tsv"),
                    sep = "\t", quote = F, row.names = F)
        
      }
    })
  }
  
  if(return_granges_list) {
    return(res_list)
  } else {
    return
  }
}
