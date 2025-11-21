################################################################################
#                                                                              #
#   PigmentCoordination_MedShrubland                                           #
#   ---------------------------------------------------                        #
#                                                                              #
#   Author:      [Paschalis Chatzopoulos]                                      #
#   Contact:     [paschatzop@gmail.com]                                        #
#   Date:        November 2025                                                 #
#   Paper:       [Pigment disorganization in Salvia rosmarinus under drought   #
#                 and nitrogen addition in a Mediterranean shrubland]          #
#   Data DOI:    [https://doi.org/10.6084/m9.figshare.30636824.v1]             #
#   Manuscript DOI: [10.XXXX/XXXXX]                                            #
#                                                                              #
#   DESCRIPTION:                                                               #
#   This script reproduces the full analytical workflow for the manuscript     #
#   listed above. Specifically, it:                                            #
#    - Imports the curated pigment dataset directly from Figshare via the      #
#      Figshare API.                                                           #
#    - Derives pigment ratios relative to total chlorophyll (a + b).           #
#    - Prepares the data structure (factors, interaction terms) required for   #
#      mixed-effects modelling.                                                #
#    - Fits linear mixed-effects models (ANOVA) for each pigment, summarises   #
#      fixed effects, and exports the resulting ANOVA table (Table 1).         #
#    - Reshapes the data to long format for visualisation and computes means   #
#      and standard errors.                                                    #
#    - Produces Figures 2, 4, 5 and Supplementary Figures S1–S3, including     #
#      bar plots, pigment ratio plots, and coupling-related plots.             #
#    - Quantifies pigment coupling (mean absolute Spearman correlation)        #
#      using EcoCoupleR.                                                      #
#    - Constructs a permutation-based null model (999 permutations) for        #
#      coupling and pigment-level coupling, and derives confidence envelopes.  #
#    - Exports supplementary regression results (Table S1) and network-related #
#      objects for subsequent network visualisation.                           #
#                                                                              #
################################################################################

# Clean environment
rm(list = ls())

# Load packages:
library(tidyverse)
library(Hmisc)
library(EcoCoupleR)
library(purrr)
library(viridis)
library(svglite)
library(nlme)
library(patchwork)
library(httr2)
library(jsonlite)
library(ggpubr)
library(broom)

# Data Import (via Figshare API) --------------------------------------------

# Define Figshare article ID containing the pigment dataset
article_id <- 30636824

# Request metadata for the specified Figshare article
res <- request(paste0("https://api.figshare.com/v2/articles/", article_id)) %>%
  req_perform()

# Parse the JSON response into an R list / data structure
meta <- res %>%
  resp_body_json(simplifyVector = TRUE)

# Select the specific file (first file associated with the record).
file_url <- meta$files$download_url[1] 
data_tidy <- read_csv(file_url, show_col_types = FALSE)

# Data Calculation and Transformation ---------------------------------------

# Calculate pigment ratios relative to total chlorophyll (Chl a + Chl b).
# These ratios express each pigment on a chlorophyll basis, capturing
# shifts in pigment composition independent of absolute chlorophyll content.
pigment_ratios <- data_tidy %>%
  mutate(Neo = (Neoxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Vio = (Violaxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Ant = (Antheraxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Lut = (Lutein / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Zea = (Zeaxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         'β-Car' = (`β-Carotene` / (`Chlorophyll-a` / `Chlorophyll-b`)),
         VAZ = (VAZ / (`Chlorophyll-a` + `Chlorophyll-b`)))

# Prepare data for statistical analyses:
# - Convert N, treatment, and block to factors.
# - Create an interaction factor (Nblock) for specifying the random structure.
# - Reorder columns to place key design variables first.
stats_data <- pigment_ratios %>%
  mutate(
    N = as.factor(N),
    treatment = as.factor(treatment),
    block = as.factor(block),
    # Interaction factor used in the random-effects structure (Nblock / block)
    Nblock = factor(paste(N, block, sep = "."))) %>%
  as.data.frame() %>%
  select(trt_N, Nblock, treatment, N, block, everything()) %>%
  select(-AZ)

# Define colour palettes for N levels and treatments
col_palette_N <- c("#90d743", "#6ece58", "#4ec36b", "#25ab82")
col_palette_treat <- c("Drought" = "#A6611A", "Recovery" = "#018571")

# Define a unified publication-style theme to be reused across figures
theme_pub <- theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.2))

# Define labeller for N levels used in facets (Figures 2 and S2)
n_labels <- c(
  "0"  = "0 kg·N·ha-1·yr-1",
  "10" = "10 kg·N·ha-1·yr-1",
  "20" = "20 kg·N·ha-1·yr-1",
  "50" = "50 kg·N·ha-1·yr-1")

# Statistical Analysis (Mixed Models) ---------------------------------------

# Define the column indices corresponding to pigment variables
# (to be analysed with mixed models).
target_cols <- 6:21 

# Initialize list to store ANOVA tables per pigment
anova_results_list <- list()

# Loop over pigment columns and fit a mixed model for each pigment.
# Fixed effects: N, treatment, and their interaction; response is log-transformed.
# Random effects: intercepts nested as 1 | Nblock / block.
for(i in target_cols) {
  
  # Retrieve the column name (pigment) for labelling results
  col_name <- names(stats_data)[i]
  
  # Assign current pigment to a temporary column used in the model
  stats_data$Y_temp <- stats_data[, i]
  
  # Fit the linear mixed-effects model for the current pigment
  tryCatch({
    model <- nlme::lme(
      fixed = log(Y_temp) ~ N * treatment,
      random = ~ 1 | Nblock / block,
      data = stats_data,
      na.action = na.omit
    )
    
    # Store the ANOVA table in the list using the pigment name as key
    anova_results_list[[col_name]] <- anova(model)
    
  }, error = function(e) {
    message(paste("Error on column index", i, "(", col_name, "):"))
    message(e)
  })
}

# Combine all ANOVA tables into a single wide data frame
table1 <- do.call("cbind", anova_results_list) %>% 
  as.data.frame()

# Identify DF columns:
# - numDF: numerator degrees of freedom
# - denDF: denominator degrees of freedom
idx_num <- grep("\\.numDF$", names(table1))
idx_den <- grep("\\.denDF$", names(table1))

# For DF columns, keep only the first instance of numDF and denDF.
# Remaining duplicated DF columns are removed.
cols_to_remove <- c(idx_num[-1], idx_den[-1])

# Process and tidy the ANOVA table:
table1_final <- table1 %>%
  # Remove redundant DF columns
  select(-all_of(cols_to_remove)) %>%
  # Rename remaining DF columns to generic names
  rename(
    numDF = all_of(idx_num[1]),
    denDF = all_of(idx_den[1])) %>%
  # Remove rows corresponding to the intercept term
  filter(row.names(.) != "(Intercept)") %>%
  # Format numeric values with two decimal places for reporting
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))

print(table1_final)

# Export ANOVA results (Table 1)
write.csv(table1,"exports/tables/Table1.csv")

# Data Preparation ----------------------------------------------------------

# Create the master LONG dataframe¨
# This step reshapes the pigment concentration data from wide to long format,
# keeping treatment, N, and block as grouping variables.
data_long <- data_tidy %>%
  pivot_longer(cols = c(-trt_N, -treatment, -N, -block), 
               names_to = "pigment", values_to = "concentration") %>%
  mutate(block = factor(block),
         N = factor(N),
         treatment = factor(treatment),
         pigment = as.factor(pigment))

# Compute mean concentration and standard error per pigment × treatment × N
# and recode pigment names to shorter labels used in figures.
data_mean <- data_long %>%
  group_by(treatment, N, pigment) %>%
  summarise(mean_conc = mean(concentration, na.rm = TRUE),
            se_conc = sd(concentration, na.rm = TRUE) / sqrt(n())) %>%
  # Recode pigment names
  mutate(pigment = fct_recode(pigment,
                              "Neo" = "Neoxanthin",
                              "Vio" = "Violaxanthin",
                              "Anth" = "Antheraxanthin",
                              "Lut" = "Lutein",
                              "Zea" = "Zeaxanthin",
                              "Chl-b" = "Chlorophyll-b",
                              "Chl-a" = "Chlorophyll-a",
                              "β-Car" = "β-Carotene"))

# Enforce a consistent pigment ordering across all figures
data_mean$pigment <- fct_relevel(data_mean$pigment, 
                                 "Anth", "β-Car", "Chl-a", "Chl-b", "Lut", "Neo", "Vio", "Zea", "AZ", "VAZ", "AZ/VAZ")

# Figure S1: Pigment Concentrations by Nitrogen ------------------------------

# Filter out ratio-based variables (AZ, VAZ) for the main pigment concentration plot
figS1_data <- data_mean %>%
  filter(!pigment %in% c("AZ", "VAZ"))

figureS1 <- ggplot(figS1_data, aes(x = pigment, y = mean_conc, fill = N)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_conc - se_conc, 
                    ymax = mean_conc + se_conc), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = col_palette_N) +
  geom_vline(xintercept = as.numeric(which(levels(figS1_data$pigment) == "Zea")) + 
               0.5, linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(~ treatment, nrow = 2) +
  labs(x = "Pigment",
       fill = expression("kg·N·ha"^-1~yr^-1*"")) +
  scale_y_continuous(name = expression("Concentration ("*mu*"mol·g"^-1*" dry weight)"),
                     sec.axis = sec_axis(
                       ~ . / max(1),
                       breaks = seq(0, 1, by = 0.25),
                       name = "De-epoxidation index"),) +
  theme_pubclean() +
  theme_pub

figureS1

# Export Figure S1 (Supplementary)
ggsave("exports/figures//FigureS1.tiff", plot = figureS1, width = 10, height = 12, dpi = 1200)
ggsave("exports/figures//FigureS1.pdf", plot = figureS1, width = 10, height = 12, dpi = 1200, device = cairo_pdf)

# Supplementary Figure S2: Pigment Ratios by Nitrogen and Treatment ----------

# Prepare data for pigment ratios (relative to total chlorophyll).
figureS2.dat <- pigment_ratios %>%
  select(trt_N, treatment, N, block, Neo, Vio, Ant, Lut, Zea, 'β-Car', VAZ) %>%
  pivot_longer(cols = c(-trt_N, -treatment, -N, -block), 
               names_to = "pigment", values_to = "ratio") %>%
  mutate(block = factor(block),
         N = factor(N),
         treatment = factor(treatment),
         pigment = as.factor(pigment)) %>%
  group_by(treatment, N, pigment) %>%
  summarise(mean_conc = mean(ratio, na.rm = TRUE),
            se_conc = sd(ratio, na.rm = TRUE) / sqrt(n()))

figureS2 <- ggplot(figureS2.dat, aes(x = pigment, y = mean_conc, fill = treatment)) +
  geom_bar(aes(group = treatment), stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_conc - se_conc, 
                    ymax = mean_conc + se_conc), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = col_palette_treat) +
  facet_wrap(~ N, labeller = as_labeller(n_labels), nrow = 2) +
  labs(x = "Pigment",
       fill = expression("kg·N·ha"^-1~yr^-1*"")) +
  scale_y_continuous(name = "Pigment-to-chlorophyll (a + b) ratio") +
  theme_pubclean() +
  theme_pub

figureS2

# Export Figure S2 (Supplementary)
ggsave("exports/figures//FigureS2.tiff", plot = figureS2, width = 10, height = 12, dpi = 1200)
ggsave("exports/figures//FigureS2.pdf", plot = figureS2, width = 10, height = 12, dpi = 1200, device = cairo_pdf)

# Figure 2: --------------------------------------------------------------------
# Data preparation:
# - Exclude AZ and VAZ (ratio variables)
# - Ensure N is a factor
subset_long <- data_long %>%
  filter(pigment != "AZ", pigment != "VAZ") %>%
  mutate(N = factor(N))

# Create an empty data frame to store significance results of paired t-tests
sig_results <- data.frame()

# Identify all unique pigment × N combinations
combinations <- unique(subset_long[, c("pigment", "N")])

# Loop over each pigment × N combination:
# - Reshape data to wide format by treatment (Drought vs Recovery)
# - Perform paired t-tests (within blocks) where both treatments are present
# - Convert p-values to significance codes (*, **, ***)
for (i in 1:nrow(combinations)) {
  pigment_i <- combinations$pigment[i]
  N_i <- combinations$N[i]
  
  # Subset data for the current pigment × N combination
  subset_df <- subset_long %>%
    filter(pigment == pigment_i, N == N_i)
  
  # Pivot to wide format for paired comparison by treatment
  wide_df <- subset_df %>%
    select(block, treatment, concentration) %>%
    pivot_wider(names_from = treatment, values_from = concentration)
  
  # Only perform the test if both treatments are present
  if (all(c("Drought", "Recovery") %in% colnames(wide_df))) {
    test <- t.test(wide_df$Drought, wide_df$Recovery, paired = TRUE)
    pval <- test$p.value
  } else {
    pval <- NA
  }
  
  # Convert p-value to significance asterisks
  sig <- if (is.na(pval)) {
    ""
  } else if (pval <= 0.001) {
    "***"
  } else if (pval <= 0.01) {
    "**"
  } else if (pval <= 0.05) {
    "*"
  } else {
    ""
  }
  
  # Append result to output data frame
  sig_results <- rbind(sig_results, data.frame(
    pigment = pigment_i,
    N = N_i,
    p_value = pval,
    sig = sig
  ))
}

# Determine y-position for the significance labels:
# - For each pigment × N, place the label slightly above the tallest bar.
y_pos_df <- subset_long %>%
  group_by(pigment, N) %>%
  summarise(y_pos = max(concentration, na.rm = TRUE) + 0.03, .groups = "drop") %>%
  # Set manually the y-position for Chlorophyll-b (to improve spacing)
  mutate(y_pos = ifelse(pigment == "Chlorophyll-b", max(y_pos) - 0.4, y_pos))

# Merge significance codes with y-positions and keep only significant comparisons
sig_df2 <- sig_results %>%
  left_join(y_pos_df, by = c("pigment", "N")) %>%
  filter(sig != "")

# Recode pigment names to match plotting labels
sig_df2$pigment <- fct_recode(sig_df2$pigment,
                              "Neo" = "Neoxanthin", "Vio" = "Violaxanthin", "Anth" = "Antheraxanthin",
                              "Lut" = "Lutein", "Zea" = "Zeaxanthin", "Chl-b" = "Chlorophyll-b",
                              "Chl-a" = "Chlorophyll-a", "β-Car" = "β-Carotene")

# Ensure pigment order is consistent with other figures
sig_df2$pigment <- fct_relevel(sig_df2$pigment, 
                               "Anth", "β-Car", "Chl-a", "Chl-b", "Lut", "Neo", "Vio", "Zea", "AZ", "VAZ", "AZ/VAZ")


figure2 <- ggplot(data_mean %>% filter(!pigment %in% c("AZ", "VAZ")),
                  aes(x = pigment, y = mean_conc, fill = treatment)) +
  geom_bar(aes(group = treatment), stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_conc - se_conc, 
                    ymax = mean_conc + se_conc), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = sig_df2, aes(x = pigment, y = y_pos, label = sig),
            position = position_dodge(0.9), vjust = 0, size = 5, inherit.aes = FALSE) +
  scale_fill_manual(values = col_palette_treat) +
  geom_vline(xintercept = as.numeric(which(levels(data_mean$pigment) == "Zea")) + 
               0.5, linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(~ N, labeller = as_labeller(n_labels), nrow = 2) +
  labs(
    x = "Pigment",
    fill = "") +
  scale_y_continuous(name = expression("Concentration ("*mu*"mol·g"^-1*" dry weight)"),
                     sec.axis = sec_axis(
                       ~ . / max(1),
                       breaks = seq(0, 1.2, by = 0.25),
                       name = "De-epoxidation index")) +
  theme_pubclean() +
  theme_pub

figure2

# Export Figure 2
ggsave("exports/figures//Figure2.tiff", plot = figure2, width = 10, height = 8, dpi = 1200)
ggsave("exports/figures//Figure2.pdf", plot = figure2, width = 10, height = 8, dpi = 1200, device = cairo_pdf)

################################################################################
# Calculate coupling ----
################################################################################

# Function to create node table from a data frame of pigment variables
create_nodes <- function(data) {
  data %>%
    # Extract unique column names (pigment names)
    colnames() %>%
    as.data.frame() %>%
    # Rename the first column to "pigment"
    rename(pigment = 1) %>%
    mutate(id = row_number()) %>%
    select(id, pigment)
}

# Data input for coupling estimation:
# - metadata: treatment, N, block, trt_N
# - data: pigment concentration matrix (used to construct correlation networks)
metadata <- data_tidy %>%
  select(trt_N, treatment, N, block)

data <- data_tidy %>%
  select(-treatment, -N, -block)

# Initialize lists to store coupling, network edges, and links per treatment × N
results_list <- list()
edges_list <- list()
links_list <- list()

# Loop through unique trt_N values (treatment × N combinations):
for (i in unique(data$trt_N)) {
  
  # Filter data for the current trt_N and retain only pigment columns
  data_filtered <- data %>% 
    filter(trt_N == i) %>%
    select(-trt_N, -AZ, -VAZ, -`AZ/VAZ`)
  
  nodes <- create_nodes(data_filtered)
  
  # Calculate the Spearman correlation matrix (r and p-values)
  cor_matrix <- rcorr(as.matrix(data_filtered), type = "spearman")
  
  # Calculate ecological coupling (mean absolute correlation) for the full network
  coupling_per_plot <- eco_coupling(cor_matrix$r, data_str = "matrix")
  
  # Extract and process correlation matrices (r, p-values)
  cor_df_raw <- as.data.frame(cor_matrix$r)
  cor_p <- as.data.frame(cor_matrix$P)
  
  # Remove diagonal (self-correlations are not informative)
  diag(cor_df_raw) <- NA  
  
  # Save an untrimmed copy of the correlation matrix to compute pigment-level coupling
  cor_df_full <- cor_df_raw
  
  # Work on a deep copy for connectivity calculations
  cor_df <- cor_df_raw
  
  # Trim lower triangle of the correlation matrix to avoid double counting
  cor_df[lower.tri(cor_df)] <- NA
  
  # Filter significant correlations (p < 0.05)
  cor_df_sig <- cor_df
  cor_df_sig[cor_p > 0.05] <- NA
  
  # Calculate network connectivity as the ratio of significant to total evaluated correlations
  connectivity <- sum(!is.na(cor_df_sig)) / sum(!is.na(cor_df))
  
  #### Pigment coupling:
  # Calculate mean absolute correlation per pigment (row-wise mean of |r|)
  pigment_coupling <- cor_df_full %>%
    mutate(pigment = rownames(cor_df_full)) %>%
    summarise(across(-pigment, ~ mean(abs(.), na.rm = TRUE))) %>%
    bind_rows() %>%
    pivot_longer(cols = everything(), names_to = "pigment", values_to = "pig_coup")
  
  # Create links and edges for every pairwise pigment interaction:
  # - links_df: edge list with weights and significance flag
  # - links_long: reshaped version listing pigments and associated correlations
  links_df <- cor_df %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
    filter(!is.na(weight)) %>%
    # Determine whether each correlation is significant or not
    mutate(sig_cor = if_else(!is.na(cor_df_sig[cbind(from, to)]), "sig", "ns"))
  
  links_long <- links_df %>%
    pivot_longer(
      cols = c(from, to),
      names_to = "direction",    
      values_to = "pigment"      
    ) %>%
    select(-direction) %>%
    distinct()
  
  # Create edge list with numeric node IDs for network visualisation
  edges_df <- links_df %>%
    left_join(nodes, by = c("from" = "pigment")) %>%
    rename(from_id = id) %>%
    left_join(nodes, by = c("to" = "pigment")) %>%
    rename(to_id = id) %>%
    select(from_id, to_id, weight, sig_cor) %>%
    rename(from = from_id, to = to_id)
  
  # Store edges and links in lists keyed by trt_N
  edges_list[[i]] <- edges_df
  links_list[[i]] <- links_df
  
  # Store all results (correlation matrices, coupling metrics, nodes) in a master list
  results_list[[i]] <- list(
    trt_N = i,
    cor_df = cor_df,
    cor_df_sig = cor_df_sig,
    cor_p = cor_p,
    coupling_per_plot = coupling_per_plot,
    connectivity = connectivity,
    pigment_coupling = pigment_coupling,
    nodes = nodes
  )
}

# Combine coupling and connectivity results into a single data frame:
coupling <- bind_rows(lapply(results_list, function(x) {
  data.frame(trt_N = x$trt_N,
             connectivity = x$connectivity,
             coupling = x$coupling_per_plot$coupling)
}))

# Combine pigment-level coupling results:
coupling_pig_df <- bind_rows(lapply(results_list, function(x) {
  data.frame(trt_N = x$trt_N, 
             pigment = x$pigment_coupling$pigment, 
             pig_coup = x$pigment_coupling$pig_coup)
}))

################################################################################
# NULL MODEL ----
################################################################################

# Number of permutations for the null model
num_permutations <- 999

# Create permuted datasets:
# - For each permutation, pigment values are shuffled within each trt_N
#   (preserving marginal distributions and sample size).
permutations <- map(1:num_permutations, ~ data %>%
                      group_by(trt_N) %>%
                      mutate_at(vars(-group_cols()), ~ sample(., replace = FALSE)) %>%
                      ungroup())

# Lists to store results for permutations
perm_results_list <- vector("list", length = num_permutations)
coupling_pig_df_perm <- data.frame()

# Loop through each permuted dataset and recompute coupling metrics
for (perm_idx in seq_along(permutations)) {
  perm_data <- permutations[[perm_idx]]
  
  # Loop through each unique trt_N within the permuted dataset
  for (i in unique(perm_data$trt_N)) {
    data_filtered_perm <- perm_data %>%
      filter(trt_N == i) %>%
      select(-trt_N, -AZ, -VAZ, -`AZ/VAZ`)
    
    # Calculate Spearman correlation matrix in the permuted data
    cor_matrix_perm <- rcorr(as.matrix(data_filtered_perm), type = "spearman")
    
    coupling_perm <- eco_coupling(cor_matrix_perm$r, data_str = "matrix")
    
    # Extract correlation matrix and p-values
    cor_df_raw_perm <- as.data.frame(cor_matrix_perm$r)
    cor_p_perm <- as.data.frame(cor_matrix_perm$P)
    
    # Remove diagonal
    diag(cor_df_raw_perm) <- NA
    
    # Save a full copy for pigment-level coupling
    cor_df_full_perm <- cor_df_raw_perm
    
    # Trim lower triangle
    cor_df_raw_perm[lower.tri(cor_df_raw_perm)] <- NA
    
    # Filter significant correlations
    cor_df_sig_perm <- cor_df_raw_perm
    cor_df_sig_perm[cor_p_perm > 0.05] <- NA
    
    # Calculate connectivity (ratio of significant correlations) in permuted data
    connectivity_perm <- sum(!is.na(cor_df_sig_perm)) / sum(!is.na(cor_df_raw_perm))
    
    #### Pigment coupling (permutation)
    coupling_per_pigment_perm <- cor_df_full_perm %>%
      mutate(pigment = rownames(cor_df_full_perm)) %>%
      summarise(across(-pigment, ~ mean(abs(.), na.rm = TRUE))) %>%
      bind_rows() %>%
      pivot_longer(cols = everything(), names_to = "pigment", values_to = "pig_coup")
    
    # Store coupling and connectivity results for this permutation and trt_N
    perm_results_list[[perm_idx]][[i]] <- list(
      trt_N = i,  
      coupling = coupling_perm,
      connectivity_perm = connectivity_perm,
      coupling_per_pigment_perm = coupling_per_pigment_perm
    )
  }
}

# Store coupling from permutations:
# - Combine permutation-level coupling and connectivity into a single data frame
final_result_df <- bind_rows(lapply(seq_along(perm_results_list), function(perm_idx) {
  bind_rows(lapply(perm_results_list[[perm_idx]], function(x) {
    data.frame(
      permutation = perm_idx,
      trt_N = x$trt_N,  
      coupling = x$coupling,
      connectivity = x$connectivity_perm
    )
  }))
}))

# Calculate null-model quantiles for coupling (per trt_N)
final_result_df <- final_result_df %>%
  group_by(trt_N) %>%
  summarise(
    quan_coup_5.0 = quantile(coupling, probs = 0.050),
    quan_coup_50.0 = quantile(coupling, probs = 0.500),
    quan_coup_95.0 = quantile(coupling, probs = 0.950))

# Merge observed coupling values with null-model quantiles and
# separate trt_N into treatment and N.
coupling_df <- final_result_df %>%
  left_join(coupling, by = "trt_N")%>% 
  separate(trt_N, into = c("treatment", "N"), sep = "_")

# Store pigment coupling from permutations:
# - Combine pigment-level coupling over all permutations and trt_N
pigment_final_df <- bind_rows(lapply(seq_along(perm_results_list), function(perm_idx) {
  bind_rows(lapply(perm_results_list[[perm_idx]], function(x) {
    data.frame(
      permutation = perm_idx,
      trt_N = x$trt_N,  
      pigment = x$coupling_per_pigment_perm$pigment,
      pig_coup = x$coupling_per_pigment_perm$pig_coup
    )
  }))
}))

# Calculate null-model quantiles for pigment-level coupling
pigment_final_df.2 <- pigment_final_df %>%
  group_by(pigment, trt_N) %>%
  summarise(
    quan_coup_5.0 = quantile(pig_coup, probs = 0.050),
    quan_coup_50.0 = quantile(pig_coup, probs = 0.500),
    quan_coup_95.0 = quantile(pig_coup, probs = 0.950)) 

# Obtain observed pigment coupling values and compare with null envelopes:
# - Combine observed pig_coup with permutation quantiles
# - Split trt_N into treatment and N
# - Derive a significance flag (coup_sig) indicating whether observed coupling
#   lies outside the 5–95% null interval.
pigment_df <- pigment_final_df.2 %>%
  left_join(coupling_pig_df, by = c("trt_N", "pigment")) %>%
  separate(trt_N, into = c("treatment", "N"), sep = "_") %>%
  arrange(treatment, N, pigment) %>%
  mutate(coup_sig = if_else(
    pig_coup < quan_coup_5.0 | pig_coup > quan_coup_95.0,
    "sig",
    "ns"
  ))

# Figure 4 ----
# Visualisation of pigment coupling (observed values and null
# envelopes) across N levels for drought and recovery treatments.
figure4 <- ggplot(coupling_df, aes(x = N, y = coupling, color = treatment, group = treatment, fill = treatment)) +
  geom_line(size = 0.6, show.legend = FALSE, linewidth = 0.6, alpha = 1) +
  geom_ribbon(aes(ymin = quan_coup_5.0, ymax = quan_coup_95.0), alpha = 0.3, linewidth = 0.8,
              show.legend = FALSE, linetype = 0) +
  geom_point(aes(x = N, y = coupling), 
             size = 2.5, shape = 20, alpha = 1,
             show.legend = FALSE) +
  scale_color_manual(values = col_palette_treat) +
  scale_fill_manual(values = col_palette_treat) +
  scale_x_discrete(expand = c(0.03, 0.03)) +
  facet_wrap(~ treatment, ncol = 2) +
  labs(x = expression("kg·N·ha"^-1~yr^-1*""),
       y = "Pigment coupling (|Spearman correlation|)") +
  ylim(0, 0.8) +
  theme_pubclean() +
  theme_pub +
  theme(panel.spacing = unit(0.5, "lines")) +
  # Text annotations to help interpret coupling gradients
  annotate("text", x = 2.25, y = 0.5, color = "antiquewhite4",
           label = "Coupled", size = 2, hjust = 0) +
  annotate("text", x = 2.2, y = 0.375, color = "antiquewhite4",
           label = "Decoupled", size = 2, hjust = 0) +
  annotate("text", x = 2.18, y = 0.25, color = "antiquewhite4",
           label = "Anticoupled", size = 2, hjust = 0) +
  theme(aspect.ratio = 1,
        axis.title.y = element_text(colour = "black", size = 8),
        axis.title.x = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 8),
        strip.text = element_text(colour = "black", size = 8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank())
figure4

# Export Figure 4
ggsave("exports/figures//Figure4.tiff", figure4, width = 16, height = 8, unit = "cm", dpi = 1200)
ggsave("exports/figures//Figure4.pdf", figure4, width = 16, height = 8, unit = "cm", dpi = 1200, device = cairo_pdf)

# Supplementary figure 3 ----
# Define a custom theme for pigment-level coupling plots
custom_theme <- theme(
  axis.text.y = element_text(colour = "black", size = 14),
  axis.text.x = element_text(colour = "black", size = 14),
  strip.text = element_text(size = 12),
  legend.text = element_text(colour = "black", size =12),
  legend.position = "bottom",
  axis.title.y = element_text(size = 14, colour = "black"),
  axis.title.x = element_text(size = 14, colour = "black"),
  panel.background = element_blank(),
  legend.key = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_line(size = 0.4, linewidth = 0.5),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank())

# Create facet labels combining pigment and treatment, later simplified in the labeller
pigment_df$facet_group <- paste(pigment_df$pigment, pigment_df$treatment, sep = "_")

# Visualize single-pigment coupling across N levels and treatments,
# including null-model ribbons (S3).
figureS3 <- ggplot(pigment_df, aes(x = N, y = pig_coup, 
                                   group = treatment, 
                                   color = treatment, 
                                   fill = treatment)) +
  geom_line(size = 1, linewidth = 1) +
  geom_ribbon(aes(ymin = quan_coup_5.0, ymax = quan_coup_95.0), alpha = 0.2,
              show.legend = FALSE, linetype = 0) +
  geom_point(aes(x = N, y = pig_coup), 
             size = 3, shape = 20, alpha = 1, show.legend = FALSE) +
  labs(x = expression("kg N ha"^-1~yr^-1*""),
       y = "Single pigment coupling (|Spearman correlation|)",
       color = "") +
  ylim(0, 1) +
  scale_x_discrete(expand = c(0.01, 0.05)) + 
  scale_color_manual(values = col_palette_treat) +
  scale_fill_manual(values = col_palette_treat) +
  facet_wrap(~ facet_group, ncol = 4, 
             labeller = labeller(facet_group = function(x) {
               # Regex to keep only the text BEFORE the "_" (the pigment name)
               sub("_.*", "", x) 
             })) +
  theme_minimal() +
  custom_theme +
  theme(aspect.ratio = 1,
        axis.title.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        legend.text = element_text(colour = "black", size =10),
        strip.text = element_text(
          margin = margin(t = 0, b = 2, r = 0, l = 0),
          vjust = 0, size = 10
        ))

figureS3

# Export Supplementary Figure S3
ggsave("exports/figures//FigureS3.tiff", figureS3, width = 30, height = 24, unit = "cm", dpi = 300)
ggsave("exports/figures//FigureS3.pdf", figureS3, width = 30, height = 24, unit = "cm", dpi = 300, device = cairo_pdf)

# Reshape pigment coupling data for pairwise comparison of N levels (Figure 5)
# Here, N = 0 (control) is used as the x-axis reference and compared against
# higher N loads (10, 20, 50 kg N ha-1 yr-1).
pigment_df_long <- pigment_df %>%
  filter(N %in% c("0", "10", "20", "50")) %>%
  select(pigment, treatment, N, pig_coup) %>%
  # Pivot into wide format with columns N_0, N_10, N_20, N_50
  pivot_wider(names_from = N, values_from = pig_coup, names_prefix = "N ") %>%
  # Pivot into long format keeping N 0 as the control axis
  pivot_longer(cols = c(`N 10`, `N 20`, `N 50`), 
               names_to = "comparison", 
               values_to = "pig_coup")

n_labels_fig_4 <- c(
  "N 10" = "10~kg~N~ha^{-1}~yr^{-1}",
  "N 20" = "20~kg~N~ha^{-1}~yr^{-1}",
  "N 50" = "50~kg~N~ha^{-1}~yr^{-1}")


# Figure 5 ----
# Scatterplots comparing pigment coupling at N = 0 (control) vs each N load.
# - 1:1 line: no change in coupling
# - Linear regression line summarises directional shifts in coupling.
figure5 <- ggplot(pigment_df_long, aes(x = `N 0`, y = pig_coup)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  geom_point(aes(color = pigment), size = 1) +
  facet_grid(treatment ~ comparison, 
             scales = "fixed", 
             labeller = labeller(comparison = as_labeller(n_labels_fig_4, label_parsed))) +
  labs(x = "Single pigment coupling (control)", 
       y = "Single pigment coupling (N load)", color = "") +
  ylim(0.2, 0.9) +
  xlim(0.2, 0.9) +
  theme_minimal() +
  coord_fixed() +
  custom_theme +
  theme(
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.text.x = element_text(colour = "black", size = 7),
    strip.text = element_text(colour = "black", size = 7),
    legend.text = element_text(colour = "black", size =7),
    legend.key.size = unit(0.2, "cm"))

print(figure5)

# Export Figure 5
ggsave("exports/figures//Figure5.tiff", figure5, width = 20, height = 10, unit = "cm", dpi = 1200)
ggsave("exports/figures//Figure5.pdf", figure5, width = 20, height = 10, unit = "cm", dpi = 1200, device = cairo_pdf)

# Run the models
# For each treatment × N-load comparison (column "comparison"), fit a linear
# model: pig_coup (N-load) ~ pig_coup (N = 0). Extract slope, SE, p-value,
# and R² to summarise how pigment coupling scales with N loading.
lm_results <- pigment_df_long %>%
  group_by(treatment, comparison) %>%
  nest() %>%
  mutate(model = map(data, ~lm(pig_coup ~ `N 0`, data = .x))) %>%
  mutate(tidied = map(model, tidy),
         glanced = map(model, glance))

# Create a clean summary table (Table S1) with slope estimates and model fit
tableS1 <- lm_results %>%
  unnest(tidied) %>%
  filter(term == "`N 0`") %>% 
  select(treatment, comparison, estimate, std.error, p.value) %>%
  left_join(
    lm_results %>% unnest(glanced) %>% select(treatment, comparison, r.squared, adj.r.squared),
    by = c("treatment", "comparison")) %>%
  rename(N_load = comparison) %>%
  mutate(across(c(estimate, std.error, p.value, r.squared, adj.r.squared), ~ round(., 3)))

# View the regression summary (optional)
print(tableS1)

# Export regression summary as Supplementary Table S1
write.csv(tableS1, "exports/tables/TableS1.csv", row.names = FALSE)

# Export network-related objects for downstream network visualisation
save(edges_list, pigment_df, results_list,
     file = "data/network_data.rdata")
