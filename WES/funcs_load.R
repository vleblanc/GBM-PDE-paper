# ---------------------------------------------------------------
## Functions to load variant outputs
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

library(vcfR)

#------------------------------------------------
##Function to get variant info from VCF table output
#------------------------------------------------

get.var.info.from.table <- function(file) {
  tbl <- read.delim(file, stringsAsFactors = F)
  tbl$multiall <- ifelse(grepl(",", tbl[,grepl(".AF", colnames(tbl))][,1]), TRUE, FALSE)
  
  #Get variant info from annotation column
  tbl$gene <- sapply(strsplit(tbl$ANN, "[|]"), "[", 4)
  tbl$effect <- sapply(strsplit(tbl$ANN, "[|]"), "[", 2)
  tbl$impact <- factor(sapply(strsplit(tbl$ANN, "[|]"), "[", 3), levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
  tbl$HVGS.c <- sapply(strsplit(tbl$ANN, "[|]"), "[", 10)
  tbl$HVGS.p <- sapply(strsplit(tbl$ANN, "[|]"), "[", 11)
  
  #Fix allele-specific values for multiallelic variants (if present)
  if(sum(tbl$multiall > 0)) {
    #Get unique positions
    tbl_multiall_unq <- unique(tbl[tbl$multiall, c("CHROM", "POS")])
    
    for(i in rownames(tbl_multiall_unq)) {
      #Get number of variants
      nvar <- length(strsplit(tbl[i, grepl(".AF", colnames(tbl))][,1], ",")[[1]])
      
      for(j in 1:nvar) {
        row_name <- as.character(as.numeric(i) + (j-1))
        
        #Fix allelic depths
        tbl[row_name, grepl(".AD", colnames(tbl))] <- sapply(strsplit(as.character(tbl[row_name, grepl(".AD", colnames(tbl))]), ","), 
                                                             function(x) paste(x[1], x[j+1], sep = ","))
        
        #Fix allelic fractions
        tbl[row_name, grepl(".AF", colnames(tbl))] <- sapply(strsplit(as.character(tbl[row_name, grepl(".AF", colnames(tbl))]), ","), 
                                                             "[", j)
      }
    }
  }
  
  #Split allele depths into ref and alt
  tbl <- cbind(tbl[, !grepl(".AD", colnames(tbl))],
               do.call(cbind, lapply(colnames(tbl)[grepl(".AD", colnames(tbl))], function(col) {
                 samp <- strsplit(col, "[.]")[[1]][1]
                 out <- data.frame(do.call('rbind', strsplit(as.character(tbl[, col]), ',')), stringsAsFactors = FALSE)
                 colnames(out) <- c(paste0(samp, "_ref_ad"), paste0(samp, "_alt_ad"))
                 return(out)
               })))
  
  #Make object with all SnpEff annotation info
  gene_list <- lapply(strsplit(tbl$ANN, ","), function(var) {
    #Split annotations and keep relevant ones
    ann <- sapply(var, function(ann) strsplit(ann, "[|]"))
    ann <- lapply(ann, function(sngl_ann) sngl_ann[c(2:4,6:8,10,11)])
    
    #Combine into data frame and add informative column names
    ann_df <- do.call(rbind.data.frame, ann)
    colnames(ann_df) <- c("effect", "impact", "gene", "gene_id", "feature_id", "transcript_biotype", "HGVS.c", "HVGS.p")
    
    #Replace empties with NAs (usually HVGS.p if not a coding mutation)
    ann_df[ann_df == ""] <- NA
    return(ann_df)
  })
  
  #Remove ANN column
  tbl <- tbl[,colnames(tbl) != "ANN"]
  
  #Convert relevant columns to numeric
  tbl[,grepl(".AF|.DP|_ref_ad|_alt_ad", colnames(tbl))] <- sapply(tbl[,grepl(".AF|.DP|_ref_ad|_alt_ad", colnames(tbl))], as.numeric)
  
  #Make second table filtering out common SNPs
  tbl_nosnp <- tbl[is.na(tbl$COMMON) | tbl$COMMON == "0",]
  
  #Return table and gene list
  return(list(tbl = tbl, tbl_nosnp = tbl_nosnp, gene_list = gene_list))
}





#------------------------------------------------
##Function to load germline variants from VCF output
#------------------------------------------------

load.germ.cov.vcf <- function(file) {
  germ_cov <- read.vcfR(file, verbose = FALSE)
  germ_cov <- data.frame(CHROM = germ_cov@fix[,"CHROM"],
                         POS = germ_cov@fix[,"POS"],
                         ref_dp = sapply(strsplit(germ_cov@fix[,"INFO"], ";"),
                                         function(x) {
                                           dp <- x[grepl("DP4", x)]
                                           dp_vals <- strsplit(dp, "=")[[1]][2]
                                           dp_vals <- strsplit(dp_vals, ",")[[1]]
                                           sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                         }),
                         alt_dp = sapply(strsplit(germ_cov@fix[,"INFO"], ";"),
                                         function(x) {
                                           dp <- x[grepl("DP4", x)]
                                           dp_vals <- strsplit(dp, "=")[[1]][2]
                                           dp_vals <- strsplit(dp_vals, ",")[[1]]
                                           sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                         }),
                         stringsAsFactors = FALSE)
  germ_cov$ref_af <- germ_cov$ref_dp / (germ_cov$ref_dp + germ_cov$alt_dp)
  germ_cov$alt_af <- germ_cov$alt_dp / (germ_cov$ref_dp + germ_cov$alt_dp)
  
  return(germ_cov)
}
