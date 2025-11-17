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
#   This script reproduces the statistical analysis (ANOVA/mixed models) and   #
#   figures for the manuscript listed above.                                   #
#                                                                              #
#   DATA SOURCE:                                                               #
#   Data is retrieved automatically via the Figshare API (ID: 30636824).       #
#   No manual download is required.                                            #
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

# Data Import (via Figshare API) --------------------------------------------

# Define Article ID in Figshare:
article_id <- 30636824

# Fetch metadata
res <- request(paste0("https://api.figshare.com/v2/articles/", article_id)) %>%
  req_perform()
meta <- res %>%
  resp_body_json(simplifyVector = TRUE)

# Select the specific file (First file in the record)
# Note: If the repo has multiple files, verify index [1] is always the dataset
file_url <- meta$files$download_url[1] 
data_tidy <- read_csv(file_url, show_col_types = FALSE)

# Data Calculation and Transformation ---------------------------------------
# Calculate Pigment Ratios relative to Total Chlorophyll (a + b)
pigment_ratios <- data_tidy %>%
  mutate(Neo = (Neoxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Vio = (Violaxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Ant = (Antheraxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Lut = (Lutein / (`Chlorophyll-a` + `Chlorophyll-b`)),
         Zea = (Zeaxanthin / (`Chlorophyll-a` + `Chlorophyll-b`)),
         'β-Car' = (`β-Carotene` / (`Chlorophyll-a` / `Chlorophyll-b`)),
         VAZ = (VAZ / (`Chlorophyll-a` + `Chlorophyll-b`)))

# Prepare data for Statistics (Factors and Interactions)
stats_data <- pigment_ratios %>%
  mutate(
    N = as.factor(N),
    treatment = as.factor(treatment),
    block = as.factor(block),
    # Interaction factor for random effects
    Nblock = factor(paste(N, block, sep = "."))) %>%
  as.data.frame() %>%
  select(trt_N, Nblock, treatment, N, block, everything()) %>%
  select(-AZ)

# Define Colors
col_palette_N <- c("#90d743", "#6ece58", "#4ec36b", "#25ab82")
col_palette_treat <- c("Drought" = "#A6611A", "Recovery" = "#018571")

# Define a Unified Theme
theme_pub <- theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Define Labeller for Facets
n_labels <- c(
  "0"  = "0 kg·N·ha-1·yr-1",
  "10" = "10 kg·N·ha-1·yr-1",
  "20" = "20 kg·N·ha-1·yr-1",
  "50" = "50 kg·N·ha-1·yr-1")

# Statistical Analysis (Mixed Models) ---------------------------------------

# Define the column indices for the pigments
target_cols <- 6:21 

# Initialize list
anova_results_list <- list()

# Loop through indices
for(i in target_cols) {
  
  # Get the actual name of the column (for saving results)
  col_name <- names(stats_data)[i]
  
  # Assign data to a temporary column
  stats_data$Y_temp <- stats_data[, i]
  
  # Run Model
  # We model 'Y_temp', but save the result under 'col_name'
  tryCatch({
    model <- nlme::lme(
      fixed = log(Y_temp) ~ N * treatment,
      random = ~ 1 | Nblock / block,
      data = stats_data,
      na.action = na.omit
    )
    
    # Save result
    anova_results_list[[col_name]] <- anova(model)
    
  }, error = function(e) {
    message(paste("Error on column index", i, "(", col_name, "):"))
    message(e)
  })
}

# Bind list to wide dataframe
table1 <- do.call("cbind", anova_results_list) %>% 
  as.data.frame()

# Identify DF Columns
# Find indices of all columns ending in .numDF or .denDF
idx_num <- grep("\\.numDF$", names(table1))
idx_den <- grep("\\.denDF$", names(table1))

# We keep the first instance of each, and list the rest for removal
cols_to_remove <- c(idx_num[-1], idx_den[-1])

# Process Table
table1_final <- table1 %>%
  # Remove the redundant DF columns 
  select(-all_of(cols_to_remove)) %>%
  # Rename the remaining first DF columns to generic names
  rename(
    numDF = all_of(idx_num[1]),
    denDF = all_of(idx_den[1])) %>%
  # Remove Intercept rows
  filter(row.names(.) != "(Intercept)") %>%
  # Round numeric values:
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))

# Export Table 1:
write.csv(table1,"exports/tables/Table1.csv")

# Data Preparation ----------------------------------------------------------
# Create the master LONG dataframe¨
data_long <- data_tidy %>%
  pivot_longer(cols = c(-trt_N, -treatment, -N, -block), 
               names_to = "pigment", values_to = "concentration") %>%
  mutate(block = factor(block),
         N = factor(N),
         treatment = factor(treatment),
         pigment = as.factor(pigment))

data_mean <- data_long %>%
  group_by(treatment, N, pigment) %>%
  summarise(mean_conc = mean(concentration, na.rm = TRUE),
            se_conc = sd(concentration, na.rm = TRUE) / sqrt(n())) %>%
  # recode pigment names
  mutate(pigment = fct_recode(pigment,
                              "Neo" = "Neoxanthin",
                              "Vio" = "Violaxanthin",
                              "Anth" = "Antheraxanthin",
                              "Lut" = "Lutein",
                              "Zea" = "Zeaxanthin",
                              "Chl-b" = "Chlorophyll-b",
                              "Chl-a" = "Chlorophyll-a",
                              "β-Car" = "β-Carotene"))

# Consistent pigment order
data_mean$pigment <- fct_relevel(data_mean$pigment, 
                               "Anth", "β-Car", "Chl-a", "Chl-b", "Lut", "Neo", "Vio", "Zea", "AZ", "VAZ", "AZ/VAZ")

# Figure 1: Pigment Concentrations by Nitrogen ------------------------------
# Filter out ratio metrics for the main plot
fig1_data <- data_mean %>%
  filter(!pigment %in% c("AZ", "VAZ"))

figureS1 <- ggplot(fig1_data, aes(x = pigment, y = mean_conc, fill = N)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = mean_conc - se_conc, 
                    ymax = mean_conc + se_conc), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = col_palette_N) +
  geom_vline(xintercept = as.numeric(which(levels(fig1_data$pigment) == "Zea")) + 
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

# Export Figure S1:
ggsave("exports/figures//figureS1.jpg", plot = figureS1, width = 10, height = 12, dpi = 300)

# Supplementary Figure S2: Pigment Ratios by Nitrogen and Treatment ------------
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

# Export Figure S2:
ggsave("exports/figures//figureS2.jpg", plot = figureS2, width = 10, height = 12, dpi = 300)

# Figure 1: --------------------------------------------------------------------
# Data preparation:
subset_long <- data_long %>%
  filter(pigment != "AZ", pigment != "VAZ") %>%
  mutate(N = factor(N))

# Create an empty data frame to store results
sig_results <- data.frame()

# Get all unique combinations of pigment and N
combinations <- unique(subset_long[, c("pigment", "N")])

# Loop over each combination
for (i in 1:nrow(combinations)) {
  pigment_i <- combinations$pigment[i]
  N_i <- combinations$N[i]
  
  # Subset data for that pigment × N
  subset_df <- subset_long %>%
    filter(pigment == pigment_i, N == N_i)
  
  # Pivot to wide for paired comparison
  wide_df <- subset_df %>%
    select(block, treatment, concentration) %>%
    pivot_wider(names_from = treatment, values_from = concentration)
  
  # Check if both treatments are present
  if (all(c("Drought", "Recovery") %in% colnames(wide_df))) {
    test <- t.test(wide_df$Drought, wide_df$Recovery, paired = TRUE)
    pval <- test$p.value
  } else {
    pval <- NA
  }
  
  # Convert p-value to asterisk
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
  
  # Add result to output data frame
  sig_results <- rbind(sig_results, data.frame(
    pigment = pigment_i,
    N = N_i,
    p_value = pval,
    sig = sig
  ))
}

# Calculate y-position above the tallest bar for each pigment × N
y_pos_df <- subset_long %>%
  group_by(pigment, N) %>%
  summarise(y_pos = max(concentration, na.rm = TRUE) + 0.03, .groups = "drop") %>%
  # set manually the y-position for AZ/VAZ
  mutate(y_pos = ifelse(pigment == "Chlorophyll-b", max(y_pos) - 0.4, y_pos))

# Merge significance and y-position
sig_df2 <- sig_results %>%
  left_join(y_pos_df, by = c("pigment", "N")) %>%
  filter(sig != "")

# Re-code pigment names to match plot
sig_df2$pigment <- fct_recode(sig_df2$pigment,
                             "Neo" = "Neoxanthin", "Vio" = "Violaxanthin", "Anth" = "Antheraxanthin",
                             "Lut" = "Lutein", "Zea" = "Zeaxanthin", "Chl-b" = "Chlorophyll-b",
                             "Chl-a" = "Chlorophyll-a", "β-Car" = "β-Carotene")

# Ensure pigment order is consistent
sig_df2$pigment <- fct_relevel(sig_df2$pigment, 
                              "Anth", "β-Car", "Chl-a", "Chl-b", "Lut", "Neo", "Vio", "Zea", "AZ", "VAZ", "AZ/VAZ")


figure1 <- ggplot(data_mean %>% filter(!pigment %in% c("AZ", "VAZ")),
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

figure1

ggsave("exports/figures//Figure1.jpg", plot = figure1, width = 10, height = 8, dpi = 300)

################################################################################
# Calculate coupling ----
################################################################################
# Function to create nodes from data
create_nodes <- function(data) {
  data %>%
    # Extract unique col names (pigment names)
    colnames() %>%
    as.data.frame() %>%
    # Rename the first column to pigment
    rename(pigment = 1) %>%
    mutate(id = row_number()) %>%
    select(id, pigment)
}

# Data input for coupling estimation:
metadata <- data_tidy %>%
  select(trt_N, treatment, N, block)

data <- data_tidy %>%
  select(-treatment, -N, -block)

# Initialize lists to store results:
results_list <- list()
edges_list <- list()
links_list <- list()

# Loop through unique trt_N values:
for (i in unique(data$trt_N)) {
  
  # Filter data for the current trt_N and select only pigments:
  data_filtered <- data %>% 
    filter(trt_N == i) %>%
    select(-trt_N, -AZ, -VAZ, -`AZ/VAZ`)
  
  nodes <- create_nodes(data_filtered)
  
  # Calculate the Spearman correlation matrix:
  cor_matrix <- rcorr(as.matrix(data_filtered), type = "spearman")
  
  # Calculate ecological coupling for all correlations:
  coupling_per_plot <- eco_coupling(cor_matrix$r, data_str = "matrix")
  
  # Extract and process correlation matrices (r, p-values, n):
  cor_df_raw <- as.data.frame(cor_matrix$r)
  cor_p <- as.data.frame(cor_matrix$P)
  
  # Remove diagonal
  diag(cor_df_raw) <- NA  
  
  # Save a copy of untrimmed correlation matrix to calculate pigment coupling later:
  cor_df_full <- cor_df_raw
  
  # Make a deep copy of the raw correlation matrix:
  cor_df <- cor_df_raw
  
  # Trim lower triangle of the correlation matrix:
  cor_df[lower.tri(cor_df)] <- NA
  
  # Filter significant correlations (p < 0.05):
  cor_df_sig <- cor_df
  cor_df_sig[cor_p > 0.05] <- NA
  
  # Calculate network connectivity (ratio of significant correlations):
  connectivity <- sum(!is.na(cor_df_sig)) / sum(!is.na(cor_df))
  
  #### Pigment coupling:
  # Calculate mean absolute correlation per pigment:
  pigment_coupling <- cor_df_full %>%
    # Add a pigment column (rownames will be added as a new column)
    mutate(pigment = rownames(cor_df_full)) %>%
    # Calculate the mean of the absolute correlations across columns, omitting NA values:
    summarise(across(-pigment, ~ mean(abs(.), na.rm = TRUE))) %>%
    # Ensure that the output has the correct structure for pivot_longer
    bind_rows() %>%
    # Pivot longer to reshape the data frame
    pivot_longer(cols = everything(), names_to = "pigment", values_to = "pig_coup")
  
  # Create links and edges for every interaction:
  links_df <- cor_df %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
    filter(!is.na(weight)) %>%
    # take the significance from cor_df_sig:
    mutate(sig_cor = if_else(!is.na(cor_df_sig[cbind(from, to)]), "sig", "ns"))
  
  # Reshape links_all by combining 'from' and 'to' into one 'sp_code' column
  links_long <- links_df %>%
    pivot_longer(
      cols = c(from, to),
      names_to = "direction",    # Temporary column to hold 'from'/'to' labels
      values_to = "pigment"      # Final column for pigment names
    ) %>%
    # Optionally remove the 'direction' column if not needed
    select(-direction) %>%
    # Remove duplicate rows if 'weight' and 'sig' were the same in 'from' and 'to' pairs
    distinct()
  
  # Create edges by merging with nodes and using the significance column
  edges_df <- links_df %>%
    left_join(nodes, by = c("from" = "pigment")) %>%
    rename(from_id = id) %>%
    left_join(nodes, by = c("to" = "pigment")) %>%
    rename(to_id = id) %>%
    select(from_id, to_id, weight, sig_cor) %>%
    rename(from = from_id, to = to_id)
  
  # Store edges and links in the list
  edges_list[[i]] <- edges_df
  links_list[[i]] <- links_df
  
  # Store results in the list:
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

# Combine coupling and connectivity results:
coupling <- bind_rows(lapply(results_list, function(x) {
  data.frame(trt_N = x$trt_N,
             connectivity = x$connectivity,
             coupling = x$coupling_per_plot$coupling)
}))

# Combine pigment coupling results:
coupling_pig_df <- bind_rows(lapply(results_list, function(x) {
  data.frame(trt_N = x$trt_N, 
             pigment = x$pigment_coupling$pigment, 
             pig_coup = x$pigment_coupling$pig_coup)
}))

################################################################################
# NULL MODEL ----
################################################################################

num_permutations <- 999

# Create permutations using purrr::map
permutations <- map(1:num_permutations, ~ data %>%
                      group_by(trt_N) %>%
                      mutate_at(vars(-group_cols()), ~ sample(., replace = FALSE)) %>%
                      ungroup())

# Lists to store results for permutations
perm_results_list <- vector("list", length = num_permutations)
coupling_pig_df_perm <- data.frame()

# Loop through each permuted dataset
for (perm_idx in seq_along(permutations)) {
  perm_data <- permutations[[perm_idx]]
  
  # Loop through each unique trt_N
  for (i in unique(perm_data$trt_N)) {
    data_filtered_perm <- perm_data %>%
      filter(trt_N == i) %>%
      select(-trt_N, -AZ, -VAZ, -`AZ/VAZ`)
    
    # Calculate Spearman correlation matrix
    cor_matrix_perm <- rcorr(as.matrix(data_filtered_perm), type = "spearman")
    
    coupling_perm <- eco_coupling(cor_matrix_perm$r, data_str = "matrix")
    
    # Extract correlation matrix and p-values
    cor_df_raw_perm <- as.data.frame(cor_matrix_perm$r)
    cor_p_perm <- as.data.frame(cor_matrix_perm$P)
    
    # Remove diagonal
    diag(cor_df_raw_perm) <- NA
    
    # Save a copy for pigment coupling
    cor_df_full_perm <- cor_df_raw_perm
    
    # Trim lower triangle of the correlation matrix
    cor_df_raw_perm[lower.tri(cor_df_raw_perm)] <- NA
    
    # Filter significant correlations (p < 0.05)
    cor_df_sig_perm <- cor_df_raw_perm
    cor_df_sig_perm[cor_p_perm > 0.05] <- NA
    
    # Calculate connectivity (ratio of significant correlations)
    connectivity_perm <- sum(!is.na(cor_df_sig_perm)) / sum(!is.na(cor_df_raw_perm))
    
    #### Pigment coupling
    coupling_per_pigment_perm <- cor_df_full_perm %>%
      mutate(pigment = rownames(cor_df_full_perm)) %>%
      summarise(across(-pigment, ~ mean(abs(.), na.rm = TRUE))) %>%
      bind_rows() %>%
      pivot_longer(cols = everything(), names_to = "pigment", values_to = "pig_coup")
    
    # Store results, including trt_N
    perm_results_list[[perm_idx]][[i]] <- list(
      trt_N = i,  # Store trt_N here
      coupling = coupling_perm,
      connectivity_perm = connectivity_perm,
      coupling_per_pigment_perm = coupling_per_pigment_perm
    )
  }
}

# Store coupling from permutations: 
# Combine results into data frames
final_result_df <- bind_rows(lapply(seq_along(perm_results_list), function(perm_idx) {
  bind_rows(lapply(perm_results_list[[perm_idx]], function(x) {
    data.frame(
      permutation = perm_idx,
      trt_N = x$trt_N,  # Include trt_N
      coupling = x$coupling,
      connectivity = x$connectivity_perm
    )
  }))
}))

# Calculate quantiles for coupling
final_result_df <- final_result_df %>%
  group_by(trt_N) %>%
  summarise(
    quan_coup_5.0 = quantile(coupling, probs = 0.050),
    quan_coup_50.0 = quantile(coupling, probs = 0.500),
    quan_coup_95.0 = quantile(coupling, probs = 0.950))

# Obtain original coupling values from coupling df:
coupling_df <- final_result_df %>%
  left_join(coupling, by = "trt_N")%>% 
  separate(trt_N, into = c("treatment", "N"), sep = "_")

# Store pigment coupling from permutations:
# Combine pigment coupling results into a data frame, including trt_N
pigment_final_df <- bind_rows(lapply(seq_along(perm_results_list), function(perm_idx) {
  bind_rows(lapply(perm_results_list[[perm_idx]], function(x) {
    data.frame(
      permutation = perm_idx,
      trt_N = x$trt_N,  # Include trt_N
      pigment = x$coupling_per_pigment_perm$pigment,
      pig_coup = x$coupling_per_pigment_perm$pig_coup
    )
  }))
}))

# Calculate quantiles for pigment coupling
pigment_final_df <- pigment_final_df %>%
  group_by(pigment, trt_N) %>%
  summarise(
    quan_coup_5.0 = quantile(pig_coup, probs = 0.050),
    quan_coup_50.0 = quantile(pig_coup, probs = 0.500),
    quan_coup_95.0 = quantile(pig_coup, probs = 0.950)) 

# Obtain original pigment coupling values from coupling df:
pigment_df <- pigment_final_df %>%
  left_join(coupling_pig_df, by = c("trt_N", "pigment")) %>%
  # Split trt_N into treatment and N
  separate(trt_N, into = c("treatment", "N"), sep = "_") %>%
  # order
  arrange(treatment, N, pigment) %>%
  # create a significance variable
  mutate(coup_sig = if_else(
    pig_coup < quan_coup_5.0 | pig_coup > quan_coup_95.0,
    "sig",
    "ns"
  ))

# Figure 3 ----
# Visualize:
figure3 <- ggplot(coupling_df, aes(x = N, y = coupling, color = treatment, group = treatment, fill = treatment)) +
  geom_line(size = 1, show.legend = FALSE, linewidth = 1, alpha = 1) +
  geom_ribbon(aes(ymin = quan_coup_5.0, ymax = quan_coup_95.0), alpha = 0.4,
              show.legend = FALSE) +
  geom_point(aes(x = N, y = coupling), 
             size = 3, shape = 20, alpha = 1,
             show.legend = FALSE) +
  scale_color_manual(values = col_palette_treat) +
  scale_fill_manual(values = col_palette_treat) +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  facet_wrap(~ treatment, ncol = 2) +
  labs(x = expression("kg·N·ha"^-1~yr^-1*""),
       y = "Pigment coupling (|Spearman correlation|)") +
  ylim(0, 1) +
  theme_pubclean() + 
  theme_pub +
  theme(panel.spacing = unit(0.5, "lines")) +
  # Add text annotations
  annotate("text", x = 2.26, y = 0.6, color = "antiquewhite4",
           label = "Coupled", size = 4, hjust = 0) +
  annotate("text", x = 2.2, y = 0.375, color = "antiquewhite4",
           label = "Decoupled", size = 4, hjust = 0) +
  annotate("text", x = 2.18, y = 0.15, color = "antiquewhite4",
           label = "Anticoupled", size = 4, hjust = 0)

figure3

# Export
ggsave("exports/figures//figure3.jpeg", figure3, width = 20, height = 12, unit = "cm", dpi = 300)

# Supplementary figure 1 ----
# Set theme for plot
custom_theme <- theme(
  axis.text.y = element_text(colour = "black", size = 14),
  axis.text.x = element_text(colour = "black", size = 14),
  strip.text = element_text(size = 12),
  legend.text = element_text(colour = "black", size =12),
  legend.position = "bottom",
  axis.title.y = element_text(size = 14),
  axis.title.x = element_text(size = 14, colour = "black"),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, size = 1),
  legend.key = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_line(size = 0.4),
  axis.line = element_line(size = 0.5)
)

# Visualize pigment coupling:
figureS3 <- ggplot(pigment_df, aes(x = N, y = pig_coup, 
                                  group = treatment, 
                                  color = treatment, 
                                  fill = treatment)) +
  geom_line(size = 1, linewidth = 1) +
  geom_ribbon(aes(ymin = quan_coup_5.0, ymax = quan_coup_95.0), alpha = 0.2,
              show.legend = FALSE) +
  geom_point(aes(x = N, y = pig_coup), 
             size = 3, shape = 20, alpha = 1, show.legend = FALSE) +
  labs(x = expression("kg N ha"^-1~yr^-1*""),
       y = "Single pigment coupling (|Spearman correlation|)",
       color = "") +
  ylim(0, 1) +
  scale_x_discrete(expand = c(0.01, 0.05)) + 
  scale_color_manual(values = col_palette_treat) +
  scale_fill_manual(values = col_palette_treat) +
  facet_wrap(~ pigment + treatment, ncol = 4, labeller = labeller(treatment = function(x) "")) +
  theme_minimal() +
  custom_theme

figureS3

# Export:
ggsave("exports/figures//figureS3.png", figureS3, width = 30, height = 24, unit = "cm", dpi = 300)

# Reshape the data into long format
pigment_df_long <- pigment_df %>%
  filter(N %in% c("0", "10", "20", "50")) %>%
  select(pigment, treatment, N, pig_coup) %>%
  # Pivot into wide format and use N_X as columns
  pivot_wider(names_from = N, values_from = pig_coup, names_prefix = "N ") %>%
  # Pivot into long but keeping N_0 for comparisons
  pivot_longer(cols = c(`N 10`, `N 20`, `N 50`), 
               names_to = "comparison", 
               values_to = "pig_coup")

# Figure 5 ----
# Plot all comparisons with N=0 in one plot, ensuring the same scale
figure4 <- ggplot(pigment_df_long, aes(x = `N 0`, y = pig_coup)) +
  geom_point(aes(color = pigment), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_grid(treatment ~ comparison, scales = "fixed") +
  labs(x = "Single pigment coupling (control)", 
       y = "Single pigment coupling (N load)", color = "") +
  ylim(0.2, 0.9) +
  xlim(0.2, 0.9) +
  # add 1:1 abline
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  custom_theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.text = element_text(colour = "black", size =10),
    legend.key.size = unit(0.2, "cm"))

print(figure4)

ggsave("exports/figures//figure4.png", figure4, width = 20, height = 12, unit = "cm", dpi = 300)

# Visualize networks:
library(igraph)
library(ggraph)
library(tidygraph)

pos <- "#2ECC71"
neg <- "#FF4040"

muestreo_un <- unique(paste(pigment_df$treatment, pigment_df$N, sep = "_"))

net_tidy_list <- list()

# Loop through each link and edge df
for (trt_load in muestreo_un) {
  
  # Extract metadata
  muestreo <- strsplit(trt_load, "_")[[1]]
  treatment <- factor(muestreo[1])
  N <- factor(muestreo[2])
  
  # Retrieve nodes and edges from the results
  nodes <- results_list[[trt_load]]$nodes %>%
    as.data.frame() %>%
    mutate(treatment = as.character(treatment), 
           N = as.character(N)) %>%
    left_join(pigment_df, by = c("treatment", "N", "pigment"))

  edges_df <- edges_list[[trt_load]] 
  
  # Debugging check for edge_alpha
  print(edges_df)
  
  # Create tbl_graph for significant edges
  net_tidy <- tbl_graph(nodes = nodes, edges = edges_df, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(N = N, treatment = treatment,
           node_fill_color = ifelse(coup_sig == "sig", "#1f78b4", "gray")) %>% # Node color
    activate(edges) %>%
      mutate(sign_cor = if_else(
        weight > 0, "pos", "neg"),
        edge_color = ifelse(sign_cor == "pos", pos, neg),
        edge_alpha = if_else(sig_cor == "sig", 0.8, 0.3),  # Transparency based on significance
        weight = abs(weight))
  
  # Save in the list
  net_tidy_list[[trt_load]] <- net_tidy
}

network_plots <- list()

# Visualize networks
for (trt_load in names(net_tidy_list)) {
  
  # Extract the tidygraph object for the current treatment and load
  net_tidy <- net_tidy_list[[trt_load]]
  
  # Extract metadata
  muestreo <- strsplit(trt_load, "_")[[1]]
  Muestreo <- factor(muestreo[1])
  N <- factor(muestreo[2])
  
  # Visualize the network using ggraph
  network_plot <- ggraph(net_tidy, layout = "circle") +  
    geom_edge_arc(
      aes(edge_color = sign_cor, edge_alpha = I(edge_alpha), edge_width = weight),
      strength = 0.2) +
    scale_edge_width_continuous(
      range = c(0.2, 2.5), limits = c(0, 1), name = "Correlation Weight") +
    scale_edge_color_manual(
      values = c("pos" = pos, "neg" = neg),
      name = "Correlation Sign",
      labels = c("Negative", "Positive")) +
    geom_node_point(
      aes(fill = I(node_fill_color)), shape = 21, color = "white", stroke = 1.2, size = 2.5) +
    geom_node_text(
      aes(label = pigment), repel = TRUE, size = 3, fontface = "bold", color = "darkgrey") +
    theme_void() +
    labs(title = paste(Muestreo, N)) +
    coord_fixed() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.title.position = "plot")
  
  # Store the plot in the list
  network_plots[[trt_load]] <- network_plot
}

# Combine all plots into a grid layout
combined_network_grid <- wrap_plots(network_plots, ncol = 4) +
  plot_annotation(theme = theme(legend.position = "right")) +
  plot_layout(guides = "collect")

# Display the grid
combined_network_grid

# Export
ggsave(
  "exports/figures//figure2.jpeg", 
  plot = combined_network_grid, 
  width = 16,
  height = 12,
  dpi = 300
)
