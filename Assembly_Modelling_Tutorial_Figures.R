#####written September 5, 2025 as a null model tutorial####
rm(list=ls());graphics.off()#clear all
#change line below to directory where your files are
setwd('/Users/grah930/Documents/WSU_fire_paper/assembly_model_code_and_tutorial')
library(reshape2);library(ggplot2);library(vegan)

#assembly processes:
#homog. selection: bNTI < -2
#variable selection: bNTI > 2
#homog. dispersal: for bNTI between -2 and 2, RCbray < -0.95
#dispersal limitation: for bNTI between -2 and 2, RCbray > 0.95
#undominated: for bNTI between -2 and 2, RCbray between -0.95 and 0.95

#####define some functions for use in calculating % assembly, mantel tests with NA values, and plotting######
# Classification function
classify_pairs <- function(bnti_vals, rc_vals, threshold = 2, rc_cut = 0.95) {
  ifelse(bnti_vals < -threshold, "Homog. Selection",
         ifelse(bnti_vals >  threshold, "Variable Selection",
                ifelse(rc_vals >  rc_cut,      "Dispersal Limitation",
                       ifelse(rc_vals < -rc_cut,      "Homog. Dispersal", "Mixed"))))
}

fraction_categories <- function(bnti, rc, samples, threshold = 2, rc_cut = 0.95) {
  # Subset matrices for the samples in this group
  bnti_sub <- bnti[samples, samples, drop = FALSE]
  rc_sub   <- rc[samples, samples, drop = FALSE]
  
  # Extract lower-triangle pairwise values
  vals_bnti <- bnti_sub[lower.tri(bnti_sub)]
  vals_rc   <- rc_sub[lower.tri(rc_sub)]
  
  # Classify each pair
  cats <- classify_pairs(vals_bnti, vals_rc, threshold, rc_cut)
  
  # Force cats to be a factor with the 5 canonical levels
  cats <- factor(cats, levels = c("Homog. Selection","Variable Selection",
                                  "Dispersal Limitation","Homog. Dispersal","Mixed"))
  
  # Compute proportions
  prop <- prop.table(table(cats)) * 100
  
  # Return both
  list(proportions = prop, raw_cats = cats)
}

# Mantel permutation test for two vectors (useful for non-square matrices)
mantel_perm <- function(x, y, nperm = 999) {
  obs_cor <- cor(x, y)
  perm_cor <- replicate(nperm, cor(sample(x), y))
  p_val <- (sum(abs(perm_cor) >= abs(obs_cor)) + 1) / (nperm + 1)
  list(statistic = obs_cor, p_value = p_val)
}

# Extract lower-triangle values from a square matrix
get_lower_tri <- function(mat) {
  mat[lower.tri(mat)]
}

# Filter NA pairs (and optional additional condition)
filter_valid <- function(x, y, cond = NULL) {
  valid <- !is.na(x) & !is.na(y)
  
  if (!is.null(cond)) {
    # Treat NA in cond as FALSE
    cond_clean <- ifelse(is.na(cond), FALSE, cond)
    valid <- valid & cond_clean
  }
  
  list(x = x[valid], y = y[valid])
}

add_threshold_lines <- function(thresholds = c(-2, 2), #default thresholds are -2 and 2
                                labels = NULL,
                                col = "blue", 
                                x_pos = 0.90, 
                                text_pos = 3) {
  abline(h = thresholds, col = col)
  
  if (!is.null(labels)) {
    if (length(labels) != length(thresholds)) stop("Length of labels must match length of thresholds")
    for (i in seq_along(thresholds)) {
      text(
        x = par("usr")[2] * x_pos,
        y = thresholds[i],
        labels = labels[i],
        pos = text_pos,
        col = col
      )
    }
  }
}

# Plot Mantel results with optional regression line
plot_mantel <- function(vals, ylab, xlab, type = c("bnti", "rc")) {
  type <- match.arg(type)
  
  # Set y-axis limits and thresholds based on type
  if (type == "bnti") {
    ylim <- c(min(-2.05, min(vals$y, na.rm = TRUE)), max(2.05,max(vals$y, na.rm = TRUE)))
    thresholds <- c(-2, 2)
    labels <- c("Homog. Selection", "Variable Selection")
  } else if (type == "rc") {
    ylim <- c(-1, 1)
    thresholds <- c(-0.95, 0.95)
    labels <- c("Homog. Dispersal", "Dispersal Limitation")
  }
  
  # Plot points
  plot(vals$x, vals$y,
       xlab = xlab,
       ylab = ylab,
       ylim = ylim)
  
  # Add threshold lines
  add_threshold_lines(thresholds = thresholds, labels = labels)
}

annotate_mantel <- function(mnt, x, y, line_sig = TRUE) {
  text(
    x = max(x, na.rm = TRUE) * 0.7,
    y = max(y, na.rm = TRUE) * 0.8,
    labels = paste0("Mantel r = ", round(mnt$statistic, 2),
                    ", p = ", signif(mnt$p_value, 3))
  )
  
  if (line_sig && mnt$p_value < 0.05) {
    mod <- lm(y ~ x)
    lines(x, predict(mod), col = "red")
  }
}

####load data####
data.set.name = 'Knelman_2019';rare.depth = 10730
rc = read.csv(paste(data.set.name,'_RC_weighted.csv', sep = ""), row.names = 1);
bnti = read.csv(paste(data.set.name,'_bNTI_weighted.csv',sep=""), row.names = 1)
otu = read.csv(paste(data.set.name,'_OTU_rarefied_',rare.depth,'.csv',sep=""), row.names = 1)

#NOTE: NEXT FOLLOWING LINES ARE FOR TUTORIAL ONLY; adding random taxa
#our OTU table does not have taxonomy, so we are added a fake taxonomy column for the last example plot. You should use real data :)
set.seed(123)  # for reproducibility
bacteria_classes <- c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria",
                      "Actinobacteria", "Bacteroidia", "Clostridia", "Firmicutes")
otu <- cbind(otu, taxon = NA)  # initialize with NA
otu[, "taxon"] <- sample(bacteria_classes, nrow(otu), replace = TRUE)
#head(otu[, "taxon"])#check

#####set groups (categorical comparisons), gradients (linear comparisons), and/or taxon of interest######
group = c('bpre','bpost','tpre','tpost')#list all groups
gradient = 'pH' #can switch to any variable of interest
taxon_of_interest = "Alphaproteobacteria"
meta = read.csv('mapping.csv')#file of other data associated with samples for plotting

#reorder metadata to match bnti and rc
sample_ids = sub("^X", "", rownames(bnti))   # row names of bnti without leading "X"
meta_order = meta[match(sample_ids, meta$Sample), ]

# Verify order matches
identical(sample_ids, as.character(meta_order$Sample))
identical(sample_ids, rownames(rc))

# Optional: quick peek if mismatches exist
# setdiff(sample_ids, meta_order$Sample)
# setdiff(sample_ids, rownames(rc))

rm(sample_ids)

###### Precompute reusable objects ######
#make euclidean distance matrix for gradient; Euclidean distance is linear and good for edaphic data
grad.mat <- as.matrix(vegdist(meta_order[[gradient]], method = "euclidean"))
rc.mat   <- as.matrix(rc)
bnti.mat <- as.matrix(bnti)
otu_subset <- otu[otu[, "taxon"] == taxon_of_interest, setdiff(colnames(otu), "taxon")]
otu_dist <- as.matrix(vegdist(t(otu_subset), method = "bray"))

bnti_lt  <- get_lower_tri(bnti.mat)
grad_lt  <- get_lower_tri(grad.mat)
rc_lt    <- get_lower_tri(rc.mat)
otu_lt <- get_lower_tri(otu_dist)

#####calculate % assembly processes in each group of samples######
#assembly processes homog. selection, variable selection, homog. dispersal, dispersal limiation, unknown
#NOTE: samples in bnti and meta must be in the same order
#for rows where group[i] = meta$test; find meta$Sample; 
#use all pairwise comparisons where both samples are in group i; that is the group that i what to find assembly process distributions for using: 

sample_groups <- meta$test
names(sample_groups) <- meta$Sample  # match sample IDs to sample names in bnti

# Compute assembly fractions for each group
assembly <- lapply(group, function(g) {
  group_samples <- names(sample_groups)[sample_groups == g]
  fraction_categories(bnti, rc, group_samples)
})

#put into matrix
# Define all possible assembly process categories
all_categories <- c("Homog. Selection", "Variable Selection",
                    "Dispersal Limitation", "Homog. Dispersal", "Mixed")

# Convert list of named vectors into a full matrix
assembly_mat <- t(sapply(assembly, function(x) {
  # x is a list with elements $proportions and $raw_cats
  prop <- numeric(length(all_categories))
  names(prop) <- all_categories
  
  # Fill in the proportions for existing categories
  prop[names(x$proportions)] <- x$proportions
  prop
}))

# Assign row names as group names
rownames(assembly_mat) <- group
colnames(assembly_mat) <- all_categories

# make boxplot from assembly_mat
# Convert matrix to data frame and melt to long format
assembly_df <- as.data.frame(assembly_mat)
assembly_df$group <- rownames(assembly_df)  # keep groups as a column
long_df <- melt(assembly_df, id.vars = "group", variable.name = "process", value.name = "percent")

# Plot grouped barplot
#pdf(paste0(data.set.name,"_assembly_distrib", ".pdf"))
ggplot(long_df, aes(x = group, y = percent, fill = process)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab("% of Assembly Process") +
  xlab("") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_blank()
  )
#dev.off()

#####boxplot of assembly processes across groups#####
#e.g., what assembly processes operate on communities to change them from pre- to post-burn

# Define all possible group pairs for cross-group comparisons
group_pairs <- t(combn(group, 2))  # each row is a pair: group1 vs group2

# Compute assembly fractions for each cross-group pair
assembly_cross <- lapply(1:nrow(group_pairs), function(i) {
  g1 <- group_pairs[i, 1]
  g2 <- group_pairs[i, 2]
  
  # Samples in each group
  samples_g1 <- names(sample_groups)[sample_groups == g1]
  samples_g2 <- names(sample_groups)[sample_groups == g2]
  
  # Subset numeric blocks
  bnti_sub <- as.matrix(bnti[samples_g1, samples_g2])
  rc_sub   <- as.matrix(rc[samples_g1, samples_g2])
  
  # Include transpose block to get lower triangle values
  bnti_sub_t <- as.matrix(bnti[samples_g2, samples_g1])
  rc_sub_t   <- as.matrix(rc[samples_g2, samples_g1])
  
  # Combine, flatten, and remove NA
  vals_bnti <- c(as.vector(bnti_sub), as.vector(bnti_sub_t))
  vals_rc   <- c(as.vector(rc_sub), as.vector(rc_sub_t))
  
  # Keep only valid (non-NA) pairs
  valid <- !is.na(vals_bnti) & !is.na(vals_rc)
  vals_bnti <- vals_bnti[valid]
  vals_rc   <- vals_rc[valid]
  
  #for debugging:
  # cat("\nPair:", g1, "vs", g2,
  #     "\n  samples_g1:", samples_g1,
  #     "\n  samples_g2:", samples_g2,
  #     "\n  head(vals_bnti):", head(vals_bnti),
  #     "\n  head(vals_rc):", head(vals_rc), "\n")
  
  #classify
  cats <- classify_pairs(vals_bnti, vals_rc)
  
  # Initialize full vector of categories
  prop <- numeric(length(all_categories))
  names(prop) <- all_categories
  
  tab <- prop.table(table(cats)) * 100  # proportions per category
  prop[names(tab)] <- tab               # fill the matching categories
  
  list(pair = paste(g1, g2, sep = "_vs_"), 
       proportions = prop, 
       raw_cats = cats)
})

# Convert list into matrix for plotting
assembly_mat_cross <- t(sapply(assembly_cross, function(x) x$proportions))
rownames(assembly_mat_cross) <- sapply(assembly_cross, function(x) x$pair)
colnames(assembly_mat_cross) <- all_categories

# Convert to long format for ggplot
assembly_df_cross <- as.data.frame(assembly_mat_cross)
assembly_df_cross$pair <- rownames(assembly_df_cross)
long_df_cross <- melt(assembly_df_cross, id.vars = "pair", variable.name = "process", value.name = "percent")

# Plot
ggplot(long_df_cross, aes(x = pair, y = percent, fill = process)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab("% of Assembly Process") +
  xlab("") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

######plot bNTI vs gradient######
#inference: selection changes along gradient (up = divergent; down = convergent)
#can plot line for visualization
#statistical tests needs to be Mantel (or similar) because pairwise comparisons violate independence assumptions of linear regression
#both axes will be pairwise (i.e., amount of change in selection vs. amount of change in gradient)

vals <- filter_valid(grad_lt,bnti_lt)
mnt <- mantel_perm(vals$x, vals$y)

plot_mantel(vals, ylab = "bNTI", xlab = paste("Difference in", gradient), type = "bnti")

######plot bNTI vs gradient######
#inference: selection changes as gradient changes 
#can plot line for visualization
#statistical tests needs to be Mantel (or similar) because pairwise comparisons violate independence assumptions of linear regression
#both axes will be pairwise (i.e., amount of change in disperal vs. amount of change in gradient)
# Filter valid pairs
vals <- filter_valid(otu_lt,bnti_lt)

# Mantel permutation test
mnt <- mantel_perm(vals$x, vals$y)

# Plot with automated thresholds and limits for bNTI
plot_mantel(vals, 
            ylab = "bNTI",
            xlab = paste("Difference in", taxon_of_interest, "(Bray-Curtis, rarefied counts)"),
            type = "bnti")

# Annotate Mantel results and regression if significant
annotate_mantel(mnt, vals$x, vals$y)

######plot RCbray of non-selection-dominated samples vs taxon######
#inference: disperal changes as taxon increase or decrease in relative abundance (up = more limited; down = less limited)
#can plot line for visualization
#statistical tests needs to be Mantel (or similar) because pairwise comparisons violate independence assumptions of linear regression
#both axes will be pairwise (i.e., amount of change in disperal vs. amount of change in taxon)
vals <- filter_valid(
  x = grad_lt,          # x-axis: gradient distances
  y = rc_lt,            # y-axis: RC values
  cond = (bnti_lt > -2 & bnti_lt < 2)
)

# Mantel permutation test
mnt <- mantel_perm(vals$x, vals$y)

# Plot RC vs Gradient with thresholds and y-limits for RC
plot_mantel(vals,
            ylab = "RC",
            xlab = paste("Difference in", gradient),
            type = "rc")  # 'type' will set thresholds and ylim internally

# Annotate Mantel results with regression line if significant
annotate_mantel(mnt, vals$x, vals$y)

######plot RCbray of non-selection-dominated samples vs taxon######
#inference: disperal changes as taxon increase or decrease in relative abundance (up = more limited; down = less limited)
#can plot line for visualization
#statistical tests needs to be Mantel (or similar) because pairwise comparisons violate independence assumptions of linear regression
#both axes will be pairwise (i.e., amount of change in disperal vs. amount of change in taxon)
vals <- filter_valid(
  x = otu_lt,       # x-axis: taxon Bray-Curtis distances
  y = rc_lt,        # y-axis: RC values
  cond = (bnti_lt > -2 & bnti_lt < 2)  # only non-selection-dominated samples
)

# --- Mantel permutation test ---
mnt <- mantel_perm(vals$x, vals$y)

# --- Plot RC vs Taxon with automated thresholds and limits for RC ---
plot_mantel(vals,
            ylab = "RC",
            xlab = paste("Difference in", taxon_of_interest, "(Bray-Curtis, rarefied counts)"),
            type = "rc")  # 'type' sets thresholds (-0.95, 0.95) and ylim (-1, 1)

# --- Annotate Mantel results and add regression line if significant ---
annotate_mantel(mnt, vals$x, vals$y)
