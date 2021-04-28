# ---------------------------------------------------------------
## Scripts to load PyClone results and plot mean cellular prevalence across pairs (Supplementary Figure 1)
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

library(ggplot2)
library(reshape2)

# ---------------------------------------------------------------
## JK124
# ---------------------------------------------------------------

# Load PyClone results
JK124_pyclone_clust <- read.delim("./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/cluster.tsv", stringsAsFactors = FALSE)
JK124_pyclone_clust$cluster_id <- as.factor(JK124_pyclone_clust$cluster_id)
JK124_pyclone_clust$sample_id <- factor(JK124_pyclone_clust$sample_id,
                                        levels = c("JK124_reg1_tis", "JK124_reg1_org", "JK124_reg2_tis"))

# Get mean cellular prevalence across clusters
JK124_pyclone_clust_means <- dcast(JK124_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")
JK124_pyclone_clust_means$cluster_size <- sapply(levels(JK124_pyclone_clust$cluster_id), 
                                                 function(clust) unique(JK124_pyclone_clust$size[JK124_pyclone_clust$cluster_id == clust]))

# Region 1 samples - correlation and plot (Supplementary Figure 1)
cor(JK124_pyclone_clust_means$JK124_reg1_tis, JK124_pyclone_clust_means$JK124_reg1_org)

ggplot(JK124_pyclone_clust_means, aes(x = JK124_reg1_tis, y = JK124_reg1_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK124_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK124 region 1") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()

# Tissue samples - correlation and plot (Supplementary Figure 1)
cor(JK124_pyclone_clust_means$JK124_reg1_tis, JK124_pyclone_clust_means$JK124_reg2_tis)

ggplot(JK124_pyclone_clust_means, aes(x = JK124_reg1_tis, y = JK124_reg2_tis)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK124_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK124 tissues") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()





# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

# Load PyClone results
JK136_pyclone_clust <- read.delim("./results/exome/pyclone/JK136_subclonefil_adfil/tables/cluster.tsv", stringsAsFactors = FALSE)
JK136_pyclone_clust$cluster_id <- as.factor(JK136_pyclone_clust$cluster_id)
JK136_pyclone_clust$sample_id <- factor(JK136_pyclone_clust$sample_id,
                                        levels = c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org",
                                                   "JK202_tis", "JK202_org"))

# Get mean cellular prevalence across clusters
JK136_pyclone_clust_means <- dcast(JK136_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")
JK136_pyclone_clust_means$cluster_size <- sapply(levels(JK136_pyclone_clust$cluster_id), 
                                                 function(clust) unique(JK136_pyclone_clust$size[JK136_pyclone_clust$cluster_id == clust]))

# Region 1 samples - correlation and plot (Supplementary Figure 1)
cor(JK136_pyclone_clust_means$JK136_reg1_tis, JK136_pyclone_clust_means$JK136_reg1_org)

ggplot(JK136_pyclone_clust_means, aes(x = JK136_reg1_tis, y = JK136_reg1_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK136_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK136 region 1") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Region 2 samples - correlation and plot (Supplementary Figure 1)
cor(JK136_pyclone_clust_means$JK136_reg2_tis, JK136_pyclone_clust_means$JK136_reg2_org)

ggplot(JK136_pyclone_clust_means, aes(x = JK136_reg2_tis, y = JK136_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK136_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK136 region 2") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Recurrence samples - correlation and plot (Supplementary Figure 1)
cor(JK136_pyclone_clust_means$JK202_tis, JK136_pyclone_clust_means$JK202_org, use = "complete")

ggplot(JK136_pyclone_clust_means, aes(x = JK202_tis, y = JK202_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK136_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK136 recurrence") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Tissue samples - correlation and plot (Supplementary Figure 1)
cor(JK136_pyclone_clust_means$JK136_reg1_tis, JK136_pyclone_clust_means$JK136_reg2_tis)

ggplot(JK136_pyclone_clust_means, aes(x = JK136_reg1_tis, y = JK136_reg2_tis)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK136_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK136 tissues") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# PDO samples - correlation and plot (Supplementary Figure 1)
cor(JK136_pyclone_clust_means$JK136_reg1_org, JK136_pyclone_clust_means$JK136_reg2_org)

ggplot(JK136_pyclone_clust_means, aes(x = JK136_reg1_org, y = JK136_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK136_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK136 PDOs") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()






# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

# Load PyClone results
JK142_pyclone_clust <- read.delim("./results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables/cluster.tsv", stringsAsFactors = FALSE)
JK142_pyclone_clust$cluster_id <- as.factor(JK142_pyclone_clust$cluster_id)
JK142_pyclone_clust$sample_id <- factor(JK142_pyclone_clust$sample_id,
                                        levels = c("JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_tis", "JK142_reg2_org",
                                                   "JK196_tis", "JK196_org"))

# Get mean cellular prevalence across clusters
JK142_pyclone_clust_means <- dcast(JK142_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")
JK142_pyclone_clust_means$cluster_size <- sapply(levels(JK142_pyclone_clust$cluster_id), 
                                                 function(clust) unique(JK142_pyclone_clust$size[JK142_pyclone_clust$cluster_id == clust]))


# Region 1 samples - correlation and plot (Supplementary Figure 1)
cor(JK142_pyclone_clust_means$JK142_reg1_tis, JK142_pyclone_clust_means$JK142_reg1_org)

ggplot(JK142_pyclone_clust_means, aes(x = JK142_reg1_tis, y = JK142_reg1_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK142_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK142 region 1") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Recurrence samples - correlation and plot (Supplementary Figure 1)
cor(JK142_pyclone_clust_means$JK196_tis, JK142_pyclone_clust_means$JK196_org)

ggplot(JK142_pyclone_clust_means, aes(x = JK196_tis, y = JK196_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK142_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK142 recurrence") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# PDO samples - correlation and plot (Supplementary Figure 1)
cor(JK142_pyclone_clust_means$JK142_reg1_org, JK142_pyclone_clust_means$JK142_reg2_org)

ggplot(JK142_pyclone_clust_means, aes(x = JK142_reg1_org, y = JK142_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK142_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK142 PDOs") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()







# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

# Load PyClone results
JK153_pyclone_clust <- read.delim("./results/exome/pyclone/JK153_subclonefil_adfil/tables/cluster.tsv", stringsAsFactors = FALSE)
JK153_pyclone_clust$cluster_id <- as.factor(JK153_pyclone_clust$cluster_id)
JK153_pyclone_clust$sample_id <- factor(JK153_pyclone_clust$sample_id,
                                        levels = c("JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org"))

# Get mean cellular prevalence across clusters
JK153_pyclone_clust_means <- dcast(JK153_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")
JK153_pyclone_clust_means$cluster_size <- sapply(levels(JK153_pyclone_clust$cluster_id), 
                                                 function(clust) unique(JK153_pyclone_clust$size[JK153_pyclone_clust$cluster_id == clust]))


# Region 1 samples - correlation and plot (Supplementary Figure 1)
cor(JK153_pyclone_clust_means$JK153_reg1_tis, JK153_pyclone_clust_means$JK153_reg1_org)

ggplot(JK153_pyclone_clust_means, aes(x = JK153_reg1_tis, y = JK153_reg1_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK153_pyclone_clust_means)))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK153 region 1") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Region 2 samples - correlation and plot (Supplementary Figure 1)
cor(JK153_pyclone_clust_means$JK153_reg2_tis, JK153_pyclone_clust_means$JK153_reg2_org)

ggplot(JK153_pyclone_clust_means, aes(x = JK153_reg2_tis, y = JK153_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK153_pyclone_clust_means$cluster_id))))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK153 region 2") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Tissue samples - correlation and plot (Supplementary Figure 1)
cor(JK153_pyclone_clust_means$JK153_reg1_tis, JK153_pyclone_clust_means$JK153_reg2_tis)

ggplot(JK153_pyclone_clust_means, aes(x = JK153_reg1_tis, y = JK153_reg2_tis)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK153_pyclone_clust_means$cluster_id))))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK153 tissues") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# PDO samples - correlation and plot (Supplementary Figure 1)
cor(JK153_pyclone_clust_means$JK153_reg1_org, JK153_pyclone_clust_means$JK153_reg2_org)

ggplot(JK153_pyclone_clust_means, aes(x = JK153_reg1_org, y = JK153_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK153_pyclone_clust_means$cluster_id))))) +
  scale_size_area() +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK153 PDOs") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()






# ---------------------------------------------------------------
## JK163
# ---------------------------------------------------------------

# Load PyClone results
JK163_pyclone_clust <- read.delim("./results/exome/pyclone/JK163_subclonefil_adfil50org2/tables/cluster.tsv", stringsAsFactors = FALSE)
JK163_pyclone_clust$cluster_id <- as.factor(JK163_pyclone_clust$cluster_id)
JK163_pyclone_clust$sample_id <- factor(JK163_pyclone_clust$sample_id,
                                        levels = c("JK163_reg1_tis", "JK163_reg1_org", "JK163_reg2_tis", "JK163_reg2_org"))

# Get mean cellular prevalence across clusters
JK163_pyclone_clust_means <- dcast(JK163_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")
JK163_pyclone_clust_means$cluster_size <- sapply(levels(JK163_pyclone_clust$cluster_id), 
                                                 function(clust) unique(JK163_pyclone_clust$size[JK163_pyclone_clust$cluster_id == clust]))


# Region 1 samples - correlation and plot (Supplementary Figure 1)
cor(JK163_pyclone_clust_means$JK163_reg1_tis, JK163_pyclone_clust_means$JK163_reg1_org)

ggplot(JK163_pyclone_clust_means, aes(x = JK163_reg1_tis, y = JK163_reg1_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(nrow(JK163_pyclone_clust_means)))) +
  scale_size_area(limits = c(0, 200)) +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK163 region 1") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Region 2 samples - correlation and plot (Supplementary Figure 1)
cor(JK163_pyclone_clust_means$JK163_reg2_tis, JK163_pyclone_clust_means$JK163_reg2_org)

ggplot(JK163_pyclone_clust_means, aes(x = JK163_reg2_tis, y = JK163_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK163_pyclone_clust_means$cluster_id))))) +
  scale_size_area(limits = c(0, 200)) +
  labs(x = "mean cellular prevalence in tissue", y = "mean cellular prevalence in PDO", title = "JK163 region 2") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# Tissue samples - correlation and plot (Supplementary Figure 1)
cor(JK163_pyclone_clust_means$JK163_reg1_tis, JK163_pyclone_clust_means$JK163_reg2_tis)

ggplot(JK163_pyclone_clust_means, aes(x = JK163_reg1_tis, y = JK163_reg2_tis)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK163_pyclone_clust_means$cluster_id))))) +
  scale_size_area(limits = c(0, 200)) +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK163 tissues") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()


# PDO samples - correlation and plot (Supplementary Figure 1)
cor(JK163_pyclone_clust_means$JK163_reg1_org, JK163_pyclone_clust_means$JK163_reg2_org)

ggplot(JK163_pyclone_clust_means, aes(x = JK163_reg1_org, y = JK163_reg2_org)) +
  geom_point(aes(colour = cluster_id, size = cluster_size)) +
  scale_colour_manual(values = unname(iwanthue(length(levels(JK163_pyclone_clust_means$cluster_id))))) +
  scale_size_area(limits = c(0, 200)) +
  labs(x = "mean cellular prevalence in region 1", y = "mean cellular prevalence in region 2", title = "JK163 PDOs") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()