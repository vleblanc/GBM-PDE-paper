gtex_ref <- as.data.frame(data.table::fread("./data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"))
gtex_annot <- as.data.frame(data.table::fread("./data/GTEx_v7_Annotations_SampleAttributesDS.txt"))[,c("SAMPID", "SMTS", "SMTSD")]
rownames(gtex_annot) <- gtex_annot$SAMPID

#Keep same reference samples as Patel et al (infercnv code)
gtex_annot <- gtex_annot[gtex_annot$SMTSD %in% c("Brain - Cerebellum",
                                                 "Brain - Caudate (basal ganglia)",
                                                 "Brain - Cortex",
                                                 "Brain - Nucleus accumbens (basal ganglia)",
                                                 "Brain - Cerebellar Hemisphere",
                                                 "Brain - Frontal Cortex (BA9)",
                                                 "Brain - Hippocampus"),]
colnames(gtex_annot)[3] <- "cell_type"

gtex_ref$Name <- sapply(strsplit(gtex_ref$Name, "[.]"), "[", 1)
rownames(gtex_ref) <- gtex_ref$Name
gtex_gene_match <- gtex_ref[,1:2]

gtex_ref <- gtex_ref[, colnames(gtex_ref) %in% gtex_annot$SAMPID]

gtex_gene_match$match <- unlist(apply(gtex_gene_match, 1, function(gene){
  if(gene[1] %in% rowData(all_mitofil_orig)$id) {
    rownames(all_mitofil_orig)[rowData(all_mitofil_orig)$id == gene[1]]
  } else if(gene[2] %in% rowData(all_mitofil_orig)$symbol) {
    rownames(all_mitofil_orig)[rowData(all_mitofil_orig)$symbol == gene[2]]
  } else {
    NA
  }
}))
gtex_gene_match <- gtex_gene_match[!is.na(gtex_gene_match$match),]

genes_rem <- c()
for(dup_gene in gtex_gene_match$match[duplicated(gtex_gene_match$match)]) {
  dat <- gtex_gene_match[gtex_gene_match$match == dup_gene, ]
  genes_rem <- c(genes_rem, dat$Name[!dat$Name %in% rowData(all_mitofil_orig)$id])
}

gtex_gene_match <- gtex_gene_match[!rownames(gtex_gene_match) %in% genes_rem,]

gtex_ref <- gtex_ref[gtex_gene_match$Name,]
rownames(gtex_ref) <- gtex_gene_match$match

saveRDS(gtex_ref, file = "./data/gtex_ref_matched.rds")
saveRDS(gtex_annot, file = "./data/gtex_ref_annot.rds")
