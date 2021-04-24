# ---------------------------------------------------------------
## Scripts to calculate signature scores for single cells according to the method and gene sets presented in Neftel et al (Cell 2019)
## Author: Veronique LeBlanc
# ---------------------------------------------------------------


# ---------------------------------------------------------------
## Tissue cells
# ---------------------------------------------------------------

# Load tissue cells object
tis <- readRDS("./data/Robjects/scRNA/clean/tis_mito_dblt_fil_allsamps.rds") #generated in the cell_types.R script

# Get scores
tis_sigs <- get.sig.scores(tis[, !is.na(tis$cell_type)], neftel_sigs, rank_method = "aggregate", n_cores = 4)

#Assign each cell to the signature for which it has the highest score
tis_sigs$assign <- factor(apply(tis_sigs[,c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")], 1, function(cell_scores)
  c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")[which.max(cell_scores)]),
  levels = c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2"))

tis$sig_assign <- tis_sigs[colnames(tis), "assign"]




# ---------------------------------------------------------------
## PDO cells
# ---------------------------------------------------------------

# Load tissue cells object
org <- readRDS("./data/Robjects/scRNA/clean/org_mito_dblt_fil_allsamps.rds") #generated in the cell_types.R script

# Get scores
org_sigs <- get.sig.scores(org[, !is.na(org$cell_type)], neftel_sigs, rank_method = "aggregate", n_cores = 4)

#Assign each cell to the signature for which it has the highest score
org_sigs$assign <- factor(apply(org_sigs[,c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")], 1, function(cell_scores)
  c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")[which.max(cell_scores)]),
  levels = c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2"))

org$sig_assign <- org_sigs[colnames(org), "assign"]




# ---------------------------------------------------------------
## BTICs
# ---------------------------------------------------------------

# Load tissue cells object
btic <- readRDS("./data/Robjects/scRNA/clean/btic_mito_dblt_fil_allsamps.rds") #generated in the cell_types.R script

# Get scores
btic_sigs <- get.sig.scores(btic[, !is.na(btic$cell_type)], neftel_sigs, rank_method = "aggregate", n_cores = 4)

#Assign each cell to the signature for which it has the highest score
btic_sigs$assign <- factor(apply(btic_sigs[,c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")], 1, function(cell_scores)
  c("AC", "MES1", "MES2", "NPC1", "NPC2", "OPC")[which.max(cell_scores)]),
  levels = c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2"))

btic$sig_assign <- btic_sigs[colnames(btic), "assign"]





# ---------------------------------------------------------------
## Figures (Figure 5 + Supplementary Figure 7)
# ---------------------------------------------------------------

# Get 2D representation of signature scores in tissue cells
tis_sigs_2d <- cbind(do.call(rbind.data.frame, apply(tis_sigs[,1:8], 1, function(cell) {
  index <- colnames(tis_sigs)
  
  #y-axis
  D <- max(cell[which(grepl("NPC|OPC", index))]) - max(cell[which(grepl("AC|MES", index))])
  y <- log2(abs(D) + 1)
  if(D < 0) y <- -y
  
  #x-axis
  if(D > 0) {
    diff <- max(cell[which(grepl("OPC", index))]) - max(cell[which(grepl("NPC", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  } else {
    diff <- max(cell[which(grepl("AC", index))]) - max(cell[which(grepl("MES", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  }
  
  return(data.frame(y, x))
})), tis_sigs[,9:ncol(tis_sigs)])



# Get 2D representation of signature scores in PDO cells
org_sigs_2d <- cbind(do.call(rbind.data.frame, apply(org_sigs[,1:8], 1, function(cell) {
  index <- colnames(org_sigs)
  
  #y-axis
  D <- max(cell[which(grepl("NPC|OPC", index))]) - max(cell[which(grepl("AC|MES", index))])
  y <- log2(abs(D) + 1)
  if(D < 0) y <- -y
  
  #x-axis
  if(D > 0) {
    diff <- max(cell[which(grepl("OPC", index))]) - max(cell[which(grepl("NPC", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  } else {
    diff <- max(cell[which(grepl("AC", index))]) - max(cell[which(grepl("MES", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  }
  
  return(data.frame(y, x))
})), org_sigs[,9:ncol(org_sigs)])




# Get 2D representation of signature scores in BTICs
btic_sigs_2d <- cbind(do.call(rbind.data.frame, apply(btic_sigs[,1:8], 1, function(cell) {
  index <- colnames(btic_sigs)
  
  #y-axis
  D <- max(cell[which(grepl("NPC|OPC", index))]) - max(cell[which(grepl("AC|MES", index))])
  y <- log2(abs(D) + 1)
  if(D < 0) y <- -y
  
  #x-axis
  if(D > 0) {
    diff <- max(cell[which(grepl("OPC", index))]) - max(cell[which(grepl("NPC", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  } else {
    diff <- max(cell[which(grepl("AC", index))]) - max(cell[which(grepl("MES", index))])
    x <- log2(abs(diff) + 1)
    if(diff > 0) x <- -x
  }
  
  return(data.frame(y, x))
})), btic_sigs[,9:ncol(btic_sigs)])



# Combine 2d scores from all cells and plot
all_sigs_2d <- rbind(tis_sigs_2d, org_sigs_2d, btic_sigs_2d)

# Plot
ggplot(all_sigs_2d[which(all_sigs_2d$cell_type == "malignant"),], aes(x = x, y = y, colour = source)) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = source_cols) +
  theme_bw() +
  facet_wrap(~patient + reg_stg, ncol = 3)


# Add quadrant info for enrichment analysis
all_sigs_2d$quadrant <- factor(apply(all_sigs_2d[,1:2], 1, function(coord) {
  if(coord[1] > 0) {
    if(coord[2] > 0){
      quad <- "NPC"
    } else {
      quad <- "OPC"
    }
  } else {
    if(coord[2] > 0) {
      quad <- "MES"
    } else {
      quad <- "AC"
    }
  }
  return(quad)
}), levels = c("MES", "AC", "OPC", "NPC"))

# Calculate enrichment between cell types for each quadrant
all_sigs_2d_enrich <- setNames(lapply(levels(all_sigs_2d$patient), function(patient){
  setNames(lapply(levels(all_sigs_2d$reg_stg), function(reg_stg) {
    dat <- all_sigs_2d[all_sigs_2d$cell_type == "malignant" & all_sigs_2d$patient == patient & all_sigs_2d$reg_stg == reg_stg,]
    
    if(nrow(dat) > 0) {
      # Between tissue and PDO cells
      dat_tis_org <- droplevels(dat[dat$source != "btic",])
      
      org_res <- sapply(levels(all_sigs_2d$quadrant), function(quad) {
        org_res <- fisher.test(table(dat_tis_org$source, dat_tis_org$quadrant == quad))
        return(c(odds_ratio = unname(org_res$estimate), p_val = org_res$p.value))
      })
      
      if(sum(dat$source == "btic") > 0) {
        # Between tissue cells and BTICs
        dat_tis_btic <- droplevels(dat[dat$source != "org",])
        
        btic_res <- sapply(levels(all_sigs_2d$quadrant), function(quad) {
          btic_res <- fisher.test(table(dat_tis_btic$source, dat_tis_btic$quadrant == quad))
          return(c(odds_ratio = unname(btic_res$estimate), p_val = btic_res$p.value))
        })
      } else {
        btic_res <- NULL
      }
      
      return(list(tis_org = org_res, tis_btic = btic_res))
    } else {
      return(NULL)
    }
  }), levels(all_sigs_2d$reg_stg))
}), levels(all_sigs_2d$patient))



# Calculate entropy of state assignments
# Classify by cell source and region/recurrence
all_sigs_ent <- data.frame(samp = c(unique(paste(tis$patient, tis$reg_stg, "tis", sep = "_")),
                                    unique(paste(org$patient, org$reg_stg, "org", sep = "_")),
                                    unique(paste(btic$patient, btic$reg_stg, "btic", sep = "_"))),
                           stringsAsFactors = FALSE)

all_sigs_ent$patient <- factor(sapply(strsplit(all_sigs_ent$samp, "_"), "[", 1), levels = levels(tis$patient))
all_sigs_ent$reg_stg <- factor(sapply(strsplit(all_sigs_ent$samp, "_"), "[", 2), levels = levels(tis$reg_stg))
all_sigs_ent$source <- factor(sapply(strsplit(all_sigs_ent$samp, "_"), "[", 3), levels = c("tis", "org", "btic"))
all_sigs_ent$tumour <- factor(ifelse(all_sigs_ent$reg_stg == "rec", 
                                     ifelse(all_sigs_ent$patient == "JK136", "JK202", "JK196"), as.character(all_sigs_ent$patient)),
                              levels = levels(tis$tumour))

all_sigs_ent <- cbind(all_sigs_ent, t(apply(all_sigs_ent, 1, function(samp) {
  # Entropy by quadrant (meta-states)
  quad_ent <- entropy::entropy(table(all_sigs_2d$quadrant[all_sigs_2d$patient == samp[2] & 
                                                            all_sigs_2d$reg_stg == samp[3] & 
                                                            all_sigs_2d$source == samp[4] &
                                                            all_sigs_2d$cell_type == "malignant"]),
                               method = "ML")
  # Entropy by state
  sig_ent <- entropy::entropy(table(all_sigs_2d$assign[all_sigs_2d$patient == samp[2] & 
                                                         all_sigs_2d$reg_stg == samp[3] & 
                                                         all_sigs_2d$source == samp[4] &
                                                         all_sigs_2d$cell_type == "malignant"]),
                              method = "ML")
  
  return(c(quad_ent = quad_ent, sig_ent = sig_ent))
})))

#Remove samples with (almost) no malignant cells
all_sigs_ent <- all_sigs_ent[!all_sigs_ent$samp %in% c("JK124_reg2_org", "JK142_reg2_tis"),]

# Plot
ggplot(all_sigs_ent, aes(x = source, y = entropy)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = patient)) +
  scale_colour_manual(values = tum_cols) +
  labs(x = NULL, y = "Shannon index") +
  theme_bw()