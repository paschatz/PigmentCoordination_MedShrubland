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
  # Add text annotations
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
        strip.text = element_text(colour = "black", size = 8))
figure3

# Export
ggsave("exports/figures//figure3.jpeg", figure3, width = 16, height = 8, unit = "cm", dpi = 300)

# Supplementary figure 3 ----
# Set theme for plot
custom_theme <- theme(
  axis.text.y = element_text(colour = "black", size = 14),
  axis.text.x = element_text(colour = "black", size = 14),
  strip.text = element_text(size = 12),
  legend.text = element_text(colour = "black", size =12),
  legend.position = "bottom",
  axis.title.y = element_text(size = 14, colour = "black"),
  axis.title.x = element_text(size = 14, colour = "black"),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, size = 1),
  legend.key = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_line(size = 0.4),
  axis.line = element_line(size = 0)
)

pigment_df$facet_group <- paste(pigment_df$pigment, pigment_df$treatment, sep = "_")

# Visualize pigment coupling:
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
               # Regex to keep only the text BEFORE the "_" (The Pigment Name)
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

n_labels_fig_4 <- c(
  "N 10" = "10~kg~N~ha^{-1}~yr^{-1}",
  "N 20" = "20~kg~N~ha^{-1}~yr^{-1}",
  "N 50" = "50~kg~N~ha^{-1}~yr^{-1}")


# Figure 4 ----
# Plot all comparisons with N=0 in one plot, ensuring the same scale
figure4 <- ggplot(pigment_df_long, aes(x = `N 0`, y = pig_coup)) +
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
  # add 1:1 abline
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

print(figure4)

ggsave("exports/figures//figure4.tiff", figure4, width = 20, height = 10, unit = "cm", dpi = 300)

library(broom) # Essential for cleaning model outputs

# Run the models
lm_results <- pigment_df_long %>%
  # Group by the factors you want to split the models by
  group_by(treatment, comparison) %>%
  nest() %>%
  # Run the linear model on each nested group
  # Formula: y ~ x  (Treatment Coupling ~ Control Coupling)
  mutate(model = map(data, ~lm(pig_coup ~ `N 0`, data = .x))) %>%
  # Extract Coefficients (Slope & Intercept)
  mutate(tidied = map(model, tidy),
         glanced = map(model, glance))

# Create a Clean Summary Table
tableS1 <- lm_results %>%
  unnest(tidied) %>%
  filter(term == "`N 0`") %>% #
  select(treatment, comparison, estimate, std.error, p.value) %>%
  # Join with R-squared values
  left_join(
    lm_results %>% unnest(glanced) %>% select(treatment, comparison, r.squared, adj.r.squared),
    by = c("treatment", "comparison")) %>%
  # rename comparison to N load:
  rename(N_load = comparison) %>%
  # round numeric collumns to 0.001
  mutate(across(c(estimate, std.error, p.value, r.squared, adj.r.squared), ~ round(., 3)))

# View the results
print(tableS1)

# Export as supplementary tables S1
write.csv(tableS1, "exports/tables/TableS1.csv", row.names = FALSE)

# Supplementary Figure 2: Network Plots ----

library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(tidyverse)
library(cowplot)

# Setup Colors
pos <- "#2ECC71"
neg <- "#FF4040"

# Pigment shotcuts:
pigment.shorts <- c(
  "Neoxanthin" = "Neo",
  "Violaxanthin" = "Vio",
  "Antheraxanthin" = "Ant",
  "Lutein" = "Lut",
  "Zeaxanthin" = "Zea",
  "Chlorophyll-b" = "Chl-b",
  "Chlorophyll-a" = "Chl-a",
  "β-Carotene" = "β-Car")

theme.net <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 12, hjust = 0.5),
  legend.position = "none",
  axis.text = element_blank(),
  axis.title = element_blank(),
  plot.margin = margin(10, 10, 10, 10),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
)

# Define a fixed order for pigments so they appear in the same spot in the circle every time
fixed_pigment_order <- sort(unique(pigment_df$pigment))

# Process network data:
net_tidy_list <- list()
muestreo_un <- unique(paste(pigment_df$treatment, pigment_df$N, sep = "_"))

for (trt_load in muestreo_un) {
  
  # Metadata
  parts <- strsplit(trt_load, "_")[[1]]
  curr_treatment <- parts[1]
  curr_N <- parts[2]
  
  # --- STEP A: Prepare Nodes & Edges ---
  current_meta <- pigment_df %>%
    filter(
      as.character(treatment) == as.character(curr_treatment),
      as.character(N) == as.character(curr_N)
    ) %>%
    select(pigment, coup_sig, treatment, N) %>%  # Keep only the pigment name and the significance column
    distinct() # ensure one row per pigment
  
  # PREPARE NODES
  nodes_raw <- results_list[[trt_load]]$nodes %>%
    as.data.frame() 
  
  # JOIN
  # Now we only join by "pigment", because we already filtered treatment/N manually
  nodes_annotated <- nodes_raw %>%
    left_join(current_meta, by = "pigment") %>%
    mutate(
      # Assign Loop Variables manually
      treatment = curr_treatment,
      N = curr_N,
      
      # Create Labels
      pigment_label = pigment.shorts[as.character(pigment)],
      
      # Create Color Logic
      node_fill_color = case_when(
        coup_sig == "sig" ~ "#1f78b4",       # Significant -> Blue
        TRUE ~ "gray90"                      # Not Significant -> Gray
      )
    )
  
  # Get edges
  edges_df <- edges_list[[trt_load]]
  
  # Create Graph Object to calculate topology metrics
  temp_graph <- tbl_graph(nodes = nodes_annotated, edges = edges_df, directed = FALSE)
  
  # --- STEP B: Enhance Graph Data ---
  net_tidy <- temp_graph %>%
    activate(nodes) %>%
    mutate(
      # 1. Keep the original full name
      pigment_full = pigment, 
      # 2. New label
      pigment_label = pigment.shorts[as.character(pigment)],
      # Color logic
      treatment = curr_treatment,
      N = curr_N,
      node_fill_color = ifelse(coup_sig == "sig", "#1f78b4", "gray90")
    ) %>%
    arrange(pigment) %>% # Sort nodes so they plot in order
    activate(edges) %>%
    mutate(
      sign_cor = if_else(weight > 0, "pos", "neg"),
      edge_color = ifelse(sign_cor == "pos", pos, neg),
      # Make non-significant edges faint, significant ones solid
      edge_alpha_val = if_else(sig_cor == "sig", 1, 0.5), 
      edge_linetype = if_else(sig_cor == "sig", "solid", "dotted"),
      weight = abs(weight)
    )
  
  net_tidy_list[[trt_load]] <- net_tidy
}

# Define the Layout (Linear Layout used to make a Circle)
# We use a linear layout sorted by our factor, then wrap it into a circle in coords

network_plots <- list()

for (trt_load in names(net_tidy_list)) {
  
  net_tidy <- net_tidy_list[[trt_load]]
  
  # Extract Title Info
  parts <- strsplit(trt_load, "_")[[1]]
  plot_title <- paste(parts[1], "-", parts[2], "kg N")
  
  # --- NETWORK ---
  p <- ggraph(net_tidy, layout = 'linear', circular = TRUE) + 
    
    # --- EDGES ---
    geom_edge_arc(
      aes(color = sign_cor,               # Set color based on the sign
          alpha = I(edge_alpha_val),      # Set transparency based on significance
          width = weight,                 # Edge width by absolute correlation
          linetype = I(edge_linetype))) + 
    
    scale_edge_width_continuous(
      range = c(0.2, 2), 
      limits = c(0, 1), 
      breaks = c(0.25, 0.5, 0.75, 1),
      name = "Correlation"
    ) +
    
    scale_edge_color_manual(
      values = c("pos" = pos, "neg" = neg), 
      labels = c("Negative", "Positive"),
      name = "Sign", guide = "none"
    ) +
    
    # --- NODES ---
    geom_node_point(
      aes(fill = I(node_fill_color), size = 0.5), # Fill based on pigment coupling sig.
      shape = 21, color = "white", stroke = 1) +
    
    # --- LABELS ---
    # Push labels slightly outward
    geom_node_text(
      aes(
        label = pigment_label,
        # Dynamic Nudge: Pushes the "ideal" label position outward
        nudge_x = x * 1.2, 
        nudge_y = y * 1.2
      ),
      
      # Keep Repel ON to fix overlaps
      repel = TRUE,
      
      # keep it clean
      segment.color = NA,
      box.padding = 0.1,
      point.padding = 0.1, 
      force = 1,
      
      size = 3, 
      color = "black"
    ) +
    
    # LEGEND
    guides(
      edge_color = "none",
      edge_width = "none",
      edge_alpha = "none",
      fill = "none",
      size = "none",
      color = "none"
    ) +
    
    # --- THEME ---
    theme_minimal() +
    coord_fixed(clip = "off", 
                xlim = c(-1.1, 1.1),
                ylim = c(-1.1, 1.1)) +
    labs(title = plot_title) +
    theme.net
  
  network_plots[[trt_load]] <- p

}


# Arrange Logically
desired_order <- c(
  "Drought_0", "Drought_10", "Drought_20", "Drought_50",
  "Recovery_0", "Recovery_10", "Recovery_20", "Recovery_50"
)

# Filter/Sort 
ordered_plots <- network_plots[desired_order[desired_order %in% names(network_plots)]]


# Define Column Titles
col_titles <- c("0 kg·N·ha-1·yr-1",
                "10 kg·N·ha-1·yr-1",
                "20 kg·N·ha-1·yr-1",
                "50 kg·N·ha-1·yr-1")

# Separate list into Top (Drought) and Bottom (Recovery) rows
row1_plots <- ordered_plots[1:4]
row2_plots <- ordered_plots[5:8]

# Loop through TOP ROW to add titles
for(i in 1:4) {
  row1_plots[[i]] <- row1_plots[[i]] + 
    ggtitle(col_titles[i]) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b=10)),
      plot.margin = margin(5, 5, 5, 5) 
    )
}

# Loop through BOTTOM ROW to ensure NO titles
for(i in 1:4) {
  row2_plots[[i]] <- row2_plots[[i]] + 
    labs(title = NULL) + # Remove title if it exists
    theme(plot.margin = margin(5, 5, 5, 5))
}

# Function to create a vertical label
create_row_label <- function(label_text) {
  ggplot() + 
    annotate("text", 
             x = 1, 
             y = 0.5, 
             label = label_text, 
             angle = 90, 
             size = 6, 
             color = "black",
             hjust = 0.5) +
    theme_void() +
    # Fix the limits so x=1 is exactly the edge
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + 
    # Negative Right Margin pulls the next plot closer
    theme(plot.margin = margin(t = 0, r = 14, b = 0, l = 0))
}

label_drought <- create_row_label("Drought")
label_recov   <- create_row_label("Recovery")

# Build Row 1 
# widths = c(0.5, 15) makes the label column much narrower
row1_combined <- label_drought + wrap_plots(row1_plots, nrow = 1) + 
  plot_layout(widths = c(0.5, 30)) 

# Build Row 2 
row2_combined <- label_recov + wrap_plots(row2_plots, nrow = 1) + 
  plot_layout(widths = c(0.5, 30))

figure2 <- (row1_combined / row2_combined) +
  plot_layout(guides = "collect") 

# View and Save
print(figure2)

ggsave("exports/figures/figure2.tiff", figure2, width = 16, height = 10)
 
