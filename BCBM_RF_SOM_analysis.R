# ====
# BCBM RADIOMIC PHENOTYPING VIA SELF-ORGANIZING MAPS & RANDOM FOREST
# ====
# Dataset : BCBM-RadioGenomics (Taha et al., Scientific Data 2025,
#           DOI: 10.1038/s41597-025-06066-6)
# Platform: Posit Cloud (free tier)
#
# HOW TO RUN:
#   1. Save this file as BCBM_pipeline.R on your computer
#   2. Posit Cloud: Session -> Restart R
#   3. Files pane -> Upload -> upload BCBM_pipeline.R + BOTH .xlsx files
#   4. Console: source("BCBM_pipeline.R")
#   5. Outputs (PNG 600 dpi + CSV) -> BCBM_SOM_Outputs3/
#
# ANALYSIS (identical to the reported results):
#   seed 2025 | ntree = 1000 | mtry = sqrt(p) | SOM 5x5 hexagonal, rlen = 500
#   Ward's D2 on codebook -> k = 4 Radiomic Phenotypes
#
# REVISION 3 (what changed vs "after revision 2"):
#   - Sections 0-9: main results are UNCHANGED (same seeds, same RNG order).
#   - Fig 2  : larger, full matrix, labels outside the matrix (no overlap).
#   - Fig 9  : cut line placed BETWEEN merge 21 and merge 22 (true k = 4);
#              label is centred in the image and shows the exact value.
#   - Fig 8  : neighbour search uses a numeric tolerance (hexagonal grid).
#   - Fig 16 : left margin added so the first axis label is not clipped.
#   - Table S1: p-values of 0 are written as "<1e-300".
#   - Extra printouts: patients before/after filtering, NZV removed, NAs imputed.
#   - NEW Section 10: extra analyses for the remaining reviewer issues
#       (patients per phenotype, field strength vs biology, stratified tests,
#        cluster-robust tests, RF grouped CV and stability, SOM stability,
#        k sensitivity, decorrelated / rank-transformed / Laplacian-without-
#        Coarseness / ComBat-then-select pipelines, lesion-size check).
#     Every Section 10 block is wrapped in safe_run(), so one failure cannot
#     stop the rest. Outputs: S_*.csv and Fig22_*.png in BCBM_SOM_Outputs3/.
# ====


# ====
# SECTION 0: ENVIRONMENT SETUP
# ====

set.seed(2025)

options(error = function() {
  cat("\n[TRACEBACK - copy everything below when reporting an error]\n")
  traceback()
  cat("[END TRACEBACK]\n")
})

WORKING_DIR    <- getwd()
RADIOMICS_FILE <- file.path(WORKING_DIR,
                            "BCBM-RadioGenomic_Radiomics_Data_v01_20250110.xlsx")
CLINICAL_FILE  <- file.path(WORKING_DIR,
                            "BCBM-RadioGenomics-Clinical-data_v01_20250321v2.xlsx")
OUTPUT_DIR     <- file.path(WORKING_DIR, "BCBM_SOM_Outputs3")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Fig 9 option: FALSE = correct cut for k = 4 (recommended).
# TRUE = draw the old line at the 20th merge height (4.042); the console
# will print how many clusters that height really gives.
FIG9_USE_LEGACY_CUT <- FALSE

pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "tibble",
          "ggplot2", "caret", "randomForest", "kohonen",
          "corrplot", "viridis", "ggdendro", "gridExtra", "scales")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

for (f in c(RADIOMICS_FILE, CLINICAL_FILE)) {
  if (!file.exists(f)) stop("File not found: ", f,
                            " - upload it via the Files pane.")
}
message("[OK] Both input files found.")

strip_special_classes <- function(df) {
  df[] <- lapply(df, function(x) {
    if (inherits(x, c("pillar_num", "vctrs_vctr"))) as.numeric(x) else x
  })
  tibble::as_tibble(df)
}

theme_pub <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      axis.title        = element_text(face = "bold", size = base_size + 2),
      axis.text         = element_text(size = base_size - 1, color = "grey10"),
      legend.title      = element_text(face = "bold", size = base_size),
      legend.text       = element_text(size = base_size - 1),
      legend.key.size   = unit(0.9, "lines"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA)
    )
}

# PNG only, 600 dpi
save_plot <- function(plot_obj, filename, width = 7, height = 6, dpi = 600) {
  ggsave(file.path(OUTPUT_DIR, paste0(filename, ".png")),
         plot_obj, width = width, height = height, dpi = dpi,
         units = "in", bg = "white")
  message("[PLOT SAVED] ", filename, ".png (600 dpi)")
}

# PNG only, 600 dpi, base-R plots
save_base_plot <- function(plot_fn, filename, width = 7, height = 6.5,
                           pointsize = 15) {
  png(file.path(OUTPUT_DIR, paste0(filename, ".png")),
      width = width, height = height, units = "in", res = 600,
      pointsize = pointsize, bg = "white")
  plot_fn()
  dev.off()
  message("[PLOT SAVED] ", filename, ".png (600 dpi)")
}

# Safe chi-square helpers
# NOTE: Monte Carlo p-values are used when the smallest OBSERVED cell count
# is below 5 (min(ct) < 5). Describe it this way in the Methods text.
drop_empty <- function(ct) ct[rowSums(ct) > 0, colSums(ct) > 0, drop = FALSE]

cramer_v <- function(ct) {
  ct <- drop_empty(ct)
  if (nrow(ct) < 2L || ncol(ct) < 2L) return(NA_real_)
  tt <- suppressWarnings(chisq.test(ct))
  n  <- sum(ct)
  k  <- min(nrow(ct) - 1L, ncol(ct) - 1L)
  as.numeric(sqrt(tt$statistic / (n * k)))
}

safe_chisq <- function(ct, B = 2000) {
  ct <- drop_empty(ct)
  if (nrow(ct) < 2L || ncol(ct) < 2L)
    return(list(statistic = NA_real_, parameter = NA_integer_, p.value = NA_real_))
  chisq.test(ct, simulate.p.value = (min(ct) < 5), B = B)
}

chi_test <- function(ct, B = 2000) {
  ct <- drop_empty(ct)
  if (nrow(ct) < 2L || ncol(ct) < 2L)
    return(list(chi2 = NA_real_, df = NA_integer_, p = NA_real_, V = NA_real_))
  tt <- suppressWarnings(chisq.test(ct, simulate.p.value = (min(ct) < 5), B = B))
  list(chi2 = unname(tt$statistic),
       df   = (nrow(ct) - 1L) * (ncol(ct) - 1L),
       p    = tt$p.value,
       V    = cramer_v(ct))
}


# ====
# SECTION 1: DATA INGESTION, MERGING, COHORT FLOW
# ====

message("\n=== SECTION 1: DATA INGESTION & MERGING ===\n")

radiomics_raw <- strip_special_classes(read_excel(RADIOMICS_FILE))
radiomics_raw <- radiomics_raw %>%
  rename_with(~ str_trim(str_squish(.x))) %>%
  mutate(across(where(is.character), str_trim))

clinical_raw <- strip_special_classes(
  read_excel(CLINICAL_FILE, sheet = "Clinical+Genetics"))
clinical_raw <- clinical_raw %>%
  rename_with(~ str_trim(str_squish(.x))) %>%
  mutate(across(where(is.character), str_trim))

radiomics_raw <- radiomics_raw %>%
  mutate(FilenamePrefix = str_trim(FilenamePrefix))
clinical_raw  <- clinical_raw %>% rename(FilenamePrefix = ID)

merged_raw <- inner_join(radiomics_raw, clinical_raw, by = "FilenamePrefix") %>%
  rename(MagneticFieldStrength = `Magnetic Field Strength ID`)

message(sprintf("[INFO] After inner join: %d lesion-level rows", nrow(merged_raw)))

# --- 1.1 Missing-data QC (Figure 1) and cohort-flow diagnostics (Issue A) ---
miss_flag <- function(x) is.na(x) | as.character(x) == "NA"
qc_df <- data.frame(
  Variable = c("ER", "PR", "HER2", "Age", "MagneticFieldStrength"),
  Missing  = c(sum(miss_flag(merged_raw$ER)),
               sum(miss_flag(merged_raw$PR)),
               sum(miss_flag(merged_raw$HER2)),
               sum(miss_flag(merged_raw$Age)),
               sum(miss_flag(merged_raw$MagneticFieldStrength)))
) %>% mutate(Percent = round(Missing / nrow(merged_raw) * 100, 1))
print(qc_df)

p1 <- ggplot(qc_df, aes(x = reorder(Variable, Missing), y = Missing)) +
  geom_col(fill = "#546E7A", width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", Missing, Percent)),
            hjust = -0.1, size = 5, fontface = "bold") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Lesions with missing value (n)") +
  theme_pub(16)
save_plot(p1, "Fig01_Missing_Data_QC", width = 8, height = 5)

cohort_flow <- data.frame(
  Step = c("Lesions after inner join",
           "Missing ER", "Missing PR", "Missing HER2",
           "Lesions with any of ER/PR/HER2 missing",
           "Analytical lesions (complete ER/PR/HER2)"),
  N = c(nrow(merged_raw),
        sum(miss_flag(merged_raw$ER)),
        sum(miss_flag(merged_raw$PR)),
        sum(miss_flag(merged_raw$HER2)),
        sum(miss_flag(merged_raw$ER) | miss_flag(merged_raw$PR) |
              miss_flag(merged_raw$HER2)),
        nrow(merged_raw) -
          sum(miss_flag(merged_raw$ER) | miss_flag(merged_raw$PR) |
                miss_flag(merged_raw$HER2)))
)
cat("\n--- Cohort flow (Issue A) ---\n")
print(cohort_flow)
write.csv(cohort_flow, file.path(OUTPUT_DIR, "CohortFlow.csv"), row.names = FALSE)

# --- 1.2 Filter and derive variables ---
merged_filtered <- merged_raw %>%
  filter(!miss_flag(ER), !miss_flag(PR), !miss_flag(HER2)) %>%
  mutate(
    ER_status   = factor(ifelse(ER   == "+", "ER+",   "ER-")),
    PR_status   = factor(ifelse(PR   == "+", "PR+",   "PR-")),
    HER2_status = factor(ifelse(HER2 == "+", "HER2+", "HER2-")),
    Manufacturer_clean = factor(str_trim(as.character(Manufacturer)))
  )

mfs_num <- suppressWarnings(as.numeric(
  str_extract(as.character(merged_filtered$MagneticFieldStrength),
              "[0-9]+\\.?[0-9]*")))
merged_filtered$MagField_clean <- factor(
  ifelse(!is.na(mfs_num) & abs(mfs_num - 1.5) < 0.26, "1.5T",
         ifelse(!is.na(mfs_num) & abs(mfs_num - 3)   < 0.51, "3T", NA_character_)),
  levels = c("1.5T", "3T"))

merged_filtered <- merged_filtered %>%
  mutate(
    MolecularSubtype = case_when(
      (ER == "+" | PR == "+") & HER2 == "-" ~ "Luminal A",
      (ER == "+" | PR == "+") & HER2 == "+" ~ "Luminal B",
      ER == "-" & PR == "-"   & HER2 == "+" ~ "HER2-Enriched",
      ER == "-" & PR == "-"   & HER2 == "-" ~ "Triple Negative",
      TRUE ~ NA_character_
    ),
    MolecularSubtype = factor(MolecularSubtype,
                              levels = c("Luminal A", "Luminal B",
                                         "HER2-Enriched", "Triple Negative"))
  ) %>%
  filter(!is.na(MolecularSubtype)) %>%
  mutate(MolecularSubtype = droplevels(MolecularSubtype))

message("[INFO] Molecular subtype distribution:")
print(table(merged_filtered$MolecularSubtype))
message(sprintf("[INFO] Final analytical dataset: %d lesions", nrow(merged_filtered)))

# --- 1.2b Patients and sessions before / after filtering (NEW) ---
pid_fun <- function(x) sub("-[^-]*$", "", as.character(x))
patient_flow <- data.frame(
  Step     = c("After inner join", "Analytical dataset"),
  Lesions  = c(nrow(merged_raw), nrow(merged_filtered)),
  Sessions = c(length(unique(merged_raw$FilenamePrefix)),
               length(unique(merged_filtered$FilenamePrefix))),
  Patients = c(length(unique(pid_fun(merged_raw$FilenamePrefix))),
               length(unique(pid_fun(merged_filtered$FilenamePrefix))))
)
cat("\n--- Patients / sessions before and after filtering ---\n")
print(patient_flow)
write.csv(patient_flow, file.path(OUTPUT_DIR, "CohortFlow_Patients.csv"),
          row.names = FALSE)

# --- 1.3 HER2 coding check (Issue C) ---
cat("\n--- HER2 raw coding (Issue C) ---\n")
print(table(merged_raw$HER2, useNA = "always"))
cat("\n--- HER2 raw coding by molecular subtype ---\n")
print(table(merged_filtered$HER2, merged_filtered$MolecularSubtype))


# ====
# SECTION 2: PRE-PROCESSING & RANDOM FOREST FEATURE SELECTION
# ====

message("\n=== SECTION 2: PRE-PROCESSING & FEATURE REDUCTION ===\n")

feature_cols <- names(radiomics_raw)[3:109]
feature_cols <- intersect(feature_cols, names(merged_filtered))
message(sprintf("[INFO] %d radiomic feature columns identified", length(feature_cols)))

X_features_raw <- merged_filtered[, feature_cols]
Y_subtype      <- merged_filtered$MolecularSubtype
meta_data      <- merged_filtered %>%
  select(FilenamePrefix, Segmentation_Name, MolecularSubtype,
         ER_status, PR_status, HER2_status,
         Manufacturer_clean, MagField_clean, Age)

nzv_idx <- nearZeroVar(X_features_raw, saveMetrics = TRUE)
features_to_remove_nzv <- rownames(nzv_idx)[nzv_idx$nzv == TRUE]
feature_cols_clean <- setdiff(feature_cols, features_to_remove_nzv)
X_features_clean   <- X_features_raw[, feature_cols_clean]
message(sprintf("[INFO] Near-zero-variance features removed: %d",
                length(features_to_remove_nzv)))

n_na_imputed <- sum(is.na(X_features_clean))
message(sprintf("[INFO] Missing feature values imputed with the median: %d",
                n_na_imputed))
X_features_clean <- X_features_clean %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

preproc_obj <- preProcess(X_features_clean, method = c("center", "scale"))
X_scaled    <- predict(preproc_obj, X_features_clean)
message(sprintf("[INFO] Features after NZV removal and scaling: %d", ncol(X_scaled)))

message("[INFO] Training Random Forest (ntree = 1000)...")
rf_model <- randomForest(
  x = X_scaled, y = Y_subtype, ntree = 1000,
  mtry = floor(sqrt(ncol(X_scaled))), importance = TRUE, keep.forest = TRUE
)

# MDA = permutation importance scaled by its standard error (randomForest
# default, scale = TRUE). MDA_unscaled = raw mean decrease in accuracy.
rf_importance <- importance(rf_model, type = 1)
rf_importance_df <- data.frame(
  Feature      = rownames(rf_importance),
  MDA          = rf_importance[, "MeanDecreaseAccuracy"],
  MDA_unscaled = importance(rf_model, type = 1, scale = FALSE)[, "MeanDecreaseAccuracy"],
  MDG          = importance(rf_model, type = 2)[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>% arrange(desc(MDA))
write.csv(rf_importance_df, file.path(OUTPUT_DIR, "RF_Importance_All_Features.csv"),
          row.names = FALSE)

oob_err <- round(rf_model$err.rate[nrow(rf_model$err.rate), "OOB"] * 100, 2)
message("[INFO] Random Forest OOB Error Rate: ", oob_err, "%")

TOP_N_FEATURES <- 15
top_features   <- rf_importance_df$Feature[1:TOP_N_FEATURES]
X_top          <- X_scaled[, top_features]

cat("\n--- Top 15 features by MDA ---\n")
top15 <- rf_importance_df[1:TOP_N_FEATURES, c("Feature", "MDA")]
top15$Feature <- sub("^original_", "", top15$Feature)
print(top15, row.names = FALSE)


# ====
# SECTION 3: SOM TRAINING & CLUSTERING
# ====

message("\n=== SECTION 3: SOM TRAINING & HIERARCHICAL CLUSTERING ===\n")

SOM_XDIM <- 5
SOM_YDIM <- 5
som_grid <- somgrid(xdim = SOM_XDIM, ydim = SOM_YDIM,
                    topo = "hexagonal", toroidal = FALSE)

X_som_matrix <- as.matrix(X_top)
som_model <- som(
  X = X_som_matrix, grid = som_grid, rlen = 500,
  alpha = c(0.05, 0.01), radius = c(SOM_XDIM, 1), keep.data = TRUE
)

node_assignments <- som_model$unit.classif
codebook_vectors <- getCodes(som_model)

N_PHENOTYPES   <- 4
codebook_dist  <- dist(codebook_vectors, method = "euclidean")
hclust_obj     <- hclust(codebook_dist, method = "ward.D2")
node_phenotype <- cutree(hclust_obj, k = N_PHENOTYPES)
obs_phenotype  <- node_phenotype[node_assignments]

meta_data$SOM_Node      <- node_assignments
meta_data$RadiomicPheno <- factor(paste0("RP-", obs_phenotype),
                                  levels = paste0("RP-", 1:N_PHENOTYPES))
message("[INFO] Radiomic Phenotype distribution:")
print(table(meta_data$RadiomicPheno))

test_vars <- c("MolecularSubtype", "ER_status", "PR_status",
               "HER2_status", "MagField_clean")

chi_results <- list()
for (var in test_vars) {
  ct    <- table(meta_data$RadiomicPheno, meta_data[[var]])
  ctest <- safe_chisq(ct)
  chi_results[[var]] <- list(
    table   = ct,
    chi2    = round(ctest$statistic, 3),
    df      = ctest$parameter,
    p_value = round(ctest$p.value, 4),
    p_adj   = NA
  )
}
p_vals_raw <- sapply(chi_results, function(x) x$p_value)
p_vals_adj <- p.adjust(p_vals_raw, method = "BH")
for (var in test_vars) chi_results[[var]]$p_adj <- round(p_vals_adj[var], 4)

stat_summary <- data.frame(
  Variable       = test_vars,
  Chi2_Statistic = sapply(chi_results, function(x) x$chi2),
  DF             = sapply(chi_results, function(x) x$df),
  P_Value        = sapply(chi_results, function(x) x$p_value),
  P_Adj_BH       = sapply(chi_results, function(x) x$p_adj),
  Significant    = sapply(chi_results, function(x) ifelse(x$p_value < 0.05, "*", "ns"))
)
print(stat_summary)
write.csv(stat_summary, file.path(OUTPUT_DIR, "Statistical_Summary_ChiSquare.csv"),
          row.names = FALSE)


# ====
# SECTION 4: MAIN-TEXT FIGURES (Fig02 - Fig15)
# ====

message("\n=== SECTION 4: FIGURES ===\n")

subtype_colors <- c("Luminal A" = "#2196F3", "Luminal B" = "#4CAF50",
                    "HER2-Enriched" = "#FF9800", "Triple Negative" = "#F44336")
phenotype_colors <- c("RP-1" = "#7B2D8B", "RP-2" = "#1976D2",
                      "RP-3" = "#388E3C", "RP-4" = "#D32F2F")

# ── Fig 02: correlation matrix of selected features (REVISED: clearer) ──
# Full matrix, feature names outside the matrix (left = horizontal,
# top = vertical), larger figure, smaller coefficient text.
message("[PLOT 2] Correlation matrix...")
cor_matrix <- cor(X_top, use = "pairwise.complete.obs", method = "pearson")
colnames(cor_matrix) <- gsub("^original_", "", top_features)
rownames(cor_matrix) <- colnames(cor_matrix)

save_base_plot(function() {
  par(mar = c(0, 0, 1, 0))
  corrplot(cor_matrix, method = "color", type = "full", order = "hclust",
           hclust.method = "ward.D2",
           col = colorRampPalette(c("#053061", "#2166AC", "#F7F7F7",
                                    "#D6604D", "#67001F"))(200),
           addCoef.col = "black", number.cex = 0.8, number.digits = 2,
           tl.pos = "lt", tl.col = "black", tl.cex = 1.0, tl.srt = 90,
           cl.cex = 1.1, cl.ratio = 0.12, cl.align.text = "l",
           diag = FALSE)
}, "Fig02_Correlation_Heatmap", width = 15, height = 14, pointsize = 14)

# ── Fig 03: RF feature importance ──
message("[PLOT 3] RF feature importance...")
rf_top20 <- rf_importance_df[1:min(20, nrow(rf_importance_df)), ] %>%
  mutate(
    Feature_short = gsub("^original_", "", Feature),
    Feature_short = factor(Feature_short, levels = rev(Feature_short)),
    Top15 = ifelse(Feature %in% top_features, "Selected Top 15", "Other")
  )

p3 <- ggplot(rf_top20, aes(x = MDA, y = Feature_short, fill = Top15)) +
  geom_col(width = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Selected Top 15" = "#1565C0", "Other" = "#90A4AE"),
                    name = "Feature selection") +
  labs(x = "Mean decrease in accuracy (MDA)", y = NULL) +
  theme_pub(16) +
  theme(legend.position = "bottom")
save_plot(p3, "Fig03_RF_Feature_Importance", width = 9, height = 7.5)

# ── Fig 04: SOM training convergence ──
message("[PLOT 4] SOM training progress...")
training_df <- data.frame(Iteration = seq_along(som_model$changes),
                          MeanDistance = som_model$changes)

p4 <- ggplot(training_df, aes(x = Iteration, y = MeanDistance)) +
  geom_line(color = "#1565C0", linewidth = 1.1) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.15, se = FALSE,
              color = "#B71C1C", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = max(training_df$Iteration) * 0.65,
           y = max(training_df$MeanDistance) * 0.9,
           label = sprintf("Final MDD: %.4f", tail(training_df$MeanDistance, 1)),
           color = "#B71C1C", fontface = "bold", size = 5.5) +
  labs(x = "Training iteration", y = "Mean distance to BMU") +
  theme_pub(16)
save_plot(p4, "Fig04_SOM_Training_Progress", width = 8, height = 5.5)

# ── Fig 05: SOM node counts ──
message("[PLOT 5] SOM counts...")
grid_pos <- expand.grid(col = 1:SOM_XDIM, row = 1:SOM_YDIM)
grid_pos$node <- 1:(SOM_XDIM * SOM_YDIM)
grid_pos$n <- as.integer(table(factor(node_assignments,
                                      levels = grid_pos$node)))

p5 <- ggplot(grid_pos, aes(x = col, y = row, fill = n)) +
  geom_tile(color = "white", linewidth = 1, width = 0.92, height = 0.92) +
  geom_text(aes(label = n), color = "white", size = 5.5, fontface = "bold") +
  scale_fill_viridis_c(option = "mako", name = "Lesions\n(n)", direction = -1) +
  labs(x = "SOM column", y = "SOM row") +
  coord_fixed() +
  theme_pub(16) +
  theme(panel.grid = element_blank())
save_plot(p5, "Fig05_SOM_Counts", width = 7.5, height = 7)

# ── Fig 06: SOM quality (mean quantisation error per node) ──
message("[PLOT 6] SOM quality...")
qe_node <- tapply(som_model$distances,
                  factor(node_assignments, levels = grid_pos$node), mean)
grid_pos$qe <- ifelse(is.na(qe_node), NA, as.numeric(qe_node))

p6 <- ggplot(grid_pos, aes(x = col, y = row, fill = qe)) +
  geom_tile(color = "white", linewidth = 1, width = 0.92, height = 0.92) +
  geom_text(aes(label = ifelse(is.na(qe), "-", sprintf("%.2f", qe))),
            color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1", na.value = "grey85",
                      name = "Mean\ndistance\nto BMU") +
  labs(x = "SOM column", y = "SOM row") +
  coord_fixed() +
  theme_pub(16) +
  theme(panel.grid = element_blank())
save_plot(p6, "Fig06_SOM_Quality", width = 7.5, height = 7)

# ── Fig 07: SOM codebook fan plot (drawn directly; version-independent) ──
draw_fan_map <- function() {
  # kohonen >= 3.0 no longer stores $pos (it uses $pts), so node coordinates
  # are built here and checked against kohonen's own hexagonal adjacency.
  # Node i -> column (i-1) %% xdim + 1, row (i-1) %/% xdim + 1; row 1 at the bottom.
  n_nodes_fan <- SOM_XDIM * SOM_YDIM
  col_i <- (seq_len(n_nodes_fan) - 1) %% SOM_XDIM + 1
  row_i <- (seq_len(n_nodes_fan) - 1) %/% SOM_XDIM + 1
  make_pos <- function(shift_even) {
    shift <- if (shift_even) (row_i %% 2 == 0) else (row_i %% 2 == 1)
    cbind(col_i + 0.5 * shift, row_i * sqrt(3) / 2)
  }
  ud_ref  <- as.matrix(unit.distances(som_grid))
  adj_ref <- abs(ud_ref - 1) < 1e-6
  pos <- make_pos(TRUE)
  mism_even <- sum(xor(abs(as.matrix(dist(make_pos(TRUE))) - 1) < 1e-6, adj_ref))
  mism_odd  <- sum(xor(abs(as.matrix(dist(make_pos(FALSE))) - 1) < 1e-6, adj_ref))
  if (mism_odd < mism_even) pos <- make_pos(FALSE)
  message(sprintf("[CHECK] Fan-plot layout vs kohonen adjacency: mismatches = %d",
                  min(mism_even, mism_odd)))
  codes <- pmax(pmin(codebook_vectors, 2), -2)
  ramp  <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  nf    <- ncol(codes)
  par(mar = c(1, 1, 1, 1))
  plot.new()
  plot.window(xlim = range(pos[, 1]) + c(-0.8, 0.8),
              ylim = range(pos[, 2]) + c(-0.8, 0.8), asp = 1)
  for (i in seq_len(nrow(pos))) {
    cx <- pos[i, 1]; cy <- pos[i, 2]
    for (j in seq_len(nf)) {
      th  <- seq((j - 1) * 2 * pi / nf, j * 2 * pi / nf, length.out = 8)
      idx <- round((codes[i, j] + 2) / 4 * 100) + 1
      polygon(c(cx, cx + 0.5 * cos(th)), c(cy, cy + 0.5 * sin(th)),
              col = ramp[idx], border = NA)
    }
    hx <- seq(pi / 6, 2 * pi + pi / 6, length.out = 7)
    polygon(cx + 0.577 * cos(hx), cy + 0.577 * sin(hx),
            border = phenotype_colors[node_phenotype[i]], lwd = 5, col = NA)
  }
  legend("topright", legend = paste0("RP-", 1:N_PHENOTYPES),
         lwd = 5, col = phenotype_colors, title = "Phenotype",
         bty = "n", cex = 1.1)
}
message("[PLOT 7] SOM codebook fan plot...")
save_base_plot(draw_fan_map, "Fig07_SOM_Codes_Fan", width = 9, height = 8.5,
               pointsize = 16)

# ── Fig 08: SOM U-matrix (mean codebook distance to hexagonal neighbours) ──
# REVISED: neighbours are found with a numeric tolerance, because distances on
# a hexagonal grid can differ from 1 by floating-point error.
message("[PLOT 8] U-matrix...")
unit_d   <- as.matrix(unit.distances(som_grid))
cb_d     <- as.matrix(dist(codebook_vectors))
n_nb     <- integer(nrow(cb_d))
umat     <- sapply(seq_len(nrow(cb_d)), function(i) {
  nb <- which(abs(unit_d[i, ] - 1) < 1e-6)
  n_nb[i] <<- length(nb)
  mean(cb_d[i, nb])
})
message(sprintf("[CHECK] Neighbours per node (hexagonal): min %d, max %d",
                min(n_nb), max(n_nb)))
grid_pos$umat <- umat

p8 <- ggplot(grid_pos, aes(x = col, y = row, fill = umat)) +
  geom_tile(color = "white", linewidth = 1, width = 0.92, height = 0.92) +
  geom_text(aes(label = sprintf("%.2f", umat)), size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#4E342E",
                      name = "Mean\nneighbour\ndistance") +
  labs(x = "SOM column", y = "SOM row") +
  coord_fixed() +
  theme_pub(16) +
  theme(panel.grid = element_blank())
save_plot(p8, "Fig08_SOM_UMatrix", width = 7.5, height = 7)

# ── Fig 09: dendrogram (REVISED) ──
# The cut for k clusters must lie BETWEEN the (n-k+1)-th and (n-k+2)-th merge
# heights. The old code used the (n-k)-th merge height, which gives k + 1
# clusters. The label is now centred in the image, with the exact value.
message("[PLOT 9] Dendrogram...")
h_sorted      <- hclust_obj$height
n_h           <- length(h_sorted)
legacy_cut    <- h_sorted[n_h - N_PHENOTYPES]
k_at_legacy   <- length(unique(cutree(hclust_obj, h = legacy_cut)))
cut_lo        <- h_sorted[n_h - N_PHENOTYPES + 1]
cut_hi        <- h_sorted[n_h - N_PHENOTYPES + 2]
cutoff_height <- if (FIG9_USE_LEGACY_CUT) legacy_cut else (cut_lo + cut_hi) / 2
message(sprintf("[CHECK] Old cut height %.3f gives %d clusters (target = %d).",
                legacy_cut, k_at_legacy, N_PHENOTYPES))
message(sprintf("[CHECK] Correct k = %d range: %.3f to %.3f; line drawn at %.3f.",
                N_PHENOTYPES, cut_lo, cut_hi, cutoff_height))
if (!FIG9_USE_LEGACY_CUT) {
  same_partition <- length(unique(paste(cutree(hclust_obj, h = cutoff_height),
                                        node_phenotype))) == N_PHENOTYPES
  stopifnot(length(unique(cutree(hclust_obj, h = cutoff_height))) == N_PHENOTYPES,
            same_partition)
}

save_base_plot(function() {
  par(mar = c(5, 5, 2, 2), cex = 1.2)
  plot(as.dendrogram(hclust_obj), leaflab = "none",
       ylab = "Euclidean distance (Ward's D2)",
       xlab = "SOM nodes (leaves unlabelled)")
  abline(h = cutoff_height, lty = 2, col = "#B71C1C", lwd = 2)
  usr   <- par("usr")
  lab   <- sprintf("Cut height = %.3f (k = %d)", cutoff_height, N_PHENOTYPES)
  x_mid <- mean(usr[1:2])
  y_lab <- cutoff_height + 0.06 * (usr[4] - usr[3])
  bw    <- strwidth(lab) * 1.08
  bh    <- strheight(lab) * 2.0
  rect(x_mid - bw / 2, y_lab - bh / 2, x_mid + bw / 2, y_lab + bh / 2,
       col = "white", border = "#B71C1C", lwd = 1.5)
  text(x_mid, y_lab, lab, col = "#B71C1C", font = 2)
}, "Fig09_Dendrogram_Codebook", width = 8, height = 6, pointsize = 14)

# ── Helper: SOM grid property maps ──
build_som_property_df <- function(meta_col, node_vec = node_assignments,
                                  xdim = SOM_XDIM, ydim = SOM_YDIM) {
  n_nodes <- xdim * ydim
  df <- data.frame(node = node_vec, property = as.character(meta_col))
  node_summary <- df %>%
    filter(!is.na(property)) %>%
    group_by(node) %>%
    summarise(
      dominant = names(which.max(table(property))),
      prop     = max(table(property)) / n(),
      n        = n(),
      .groups  = "drop"
    )
  grid_pos2 <- expand.grid(col = 1:xdim, row = 1:ydim)
  grid_pos2$node <- 1:n_nodes
  left_join(grid_pos2, node_summary, by = "node") %>%
    mutate(dominant = ifelse(is.na(dominant), "Empty", dominant),
           prop     = ifelse(is.na(prop), 0, prop),
           n        = ifelse(is.na(n), 0, n))
}

plot_som_grid <- function(grid_df, fill_values, fill_name = "Dominant\ngroup") {
  ggplot(grid_df, aes(x = col, y = row, fill = dominant, alpha = prop)) +
    geom_tile(color = "white", linewidth = 1.0, width = 0.92, height = 0.92) +
    geom_text(aes(label = n), color = "black", size = 5.5, fontface = "bold") +
    scale_alpha_continuous(range = c(0.35, 1.0), name = "Proportion\n(dominant)") +
    scale_fill_manual(values = fill_values, name = fill_name, na.value = "grey85") +
    labs(x = "SOM column", y = "SOM row") +
    coord_fixed() +
    theme_pub(16) +
    theme(panel.grid = element_blank())
}

# ── Fig 10: HER2 map ──
message("[PLOT 10] HER2 status map...")
p10 <- plot_som_grid(build_som_property_df(meta_data$HER2_status),
                     c("HER2+" = "#E65100", "HER2-" = "#37474F", "Empty" = "grey85"))
save_plot(p10, "Fig10_SOM_HER2_Status", width = 8, height = 7.5)

# ── Fig 11: ER/PR map ──
message("[PLOT 11] ER/PR status map...")
meta_data$ER_PR <- paste(meta_data$ER_status, meta_data$PR_status, sep = "/")
p11 <- plot_som_grid(build_som_property_df(meta_data$ER_PR),
                     c("ER+/PR+" = "#1B5E20", "ER+/PR-" = "#66BB6A",
                       "ER-/PR+" = "#AB47BC", "ER-/PR-" = "#37474F",
                       "Empty" = "grey85"))
save_plot(p11, "Fig11_SOM_ERPR_Status", width = 8, height = 7.5)

# ── Fig 12: TNBC map ──
message("[PLOT 12] TNBC distribution...")
meta_data$TNBC <- ifelse(meta_data$MolecularSubtype == "Triple Negative",
                         "TNBC", "Non-TNBC")
p12 <- plot_som_grid(build_som_property_df(meta_data$TNBC),
                     c("TNBC" = "#C62828", "Non-TNBC" = "#90A4AE",
                       "Empty" = "grey85"))
save_plot(p12, "Fig12_SOM_TNBC_Distribution", width = 8, height = 7.5)

# ── Fig 14: magnetic field strength map ──
message("[PLOT 14] Field strength...")
p14 <- plot_som_grid(build_som_property_df(meta_data$MagField_clean),
                     c("1.5T" = "#0D47A1", "3T" = "#880E4F", "Empty" = "grey85"),
                     fill_name = "Dominant\nfield strength")
save_plot(p14, "Fig14_SOM_FieldStrength", width = 8, height = 7.5)

# ── Fig 15: consensus map + subtype enrichment ──
message("[PLOT 15] Consensus map + enrichment...")
pheno_grid <- build_som_property_df(meta_data$RadiomicPheno)

p15a <- ggplot(pheno_grid, aes(x = col, y = row, fill = dominant)) +
  geom_tile(color = "white", linewidth = 1.2, width = 0.92, height = 0.92) +
  geom_text(aes(label = ifelse(n > 0, paste0(dominant, "\nn=", n), "")),
            color = "white", size = 5, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c(phenotype_colors, "Empty" = "grey85"),
                    name = "Radiomic\nphenotype") +
  labs(x = "SOM column", y = "SOM row") +
  coord_fixed() +
  theme_pub(20) +
  theme(panel.grid = element_blank())

enrichment_df <- meta_data %>%
  group_by(RadiomicPheno, MolecularSubtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(RadiomicPheno) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p15b <- ggplot(enrichment_df, aes(x = RadiomicPheno, y = pct, fill = MolecularSubtype)) +
  geom_col(position = "stack", color = "white", width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 6, fontface = "bold") +
  scale_fill_manual(values = subtype_colors, name = "Molecular\nsubtype") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0, 0)) +
  labs(x = "Radiomic phenotype", y = "Lesions (%)") +
  theme_pub(20) +
  theme(legend.position = "right")

p15_combined <- gridExtra::arrangeGrob(p15a, p15b, ncol = 2, widths = c(1, 1.15))
save_plot(p15_combined, "Fig15_Consensus_Map_Enrichment", width = 15, height = 7.5)


# ====
# SECTION 5: STATISTICAL SUMMARY (Tables 2-3 source data)
# ====

message("\n=== SECTION 5: STATISTICAL SUMMARY ===\n")

for (var in test_vars) {
  cat(sprintf("\n--- %s ---\n", var))
  print(chi_results[[var]]$table)
  cat(sprintf("  Chi-Square: %s | df: %s | p: %s | p-adj (BH): %s\n",
              chi_results[[var]]$chi2, chi_results[[var]]$df,
              chi_results[[var]]$p_value, chi_results[[var]]$p_adj))
}

pheno_stats <- meta_data %>%
  group_by(RadiomicPheno) %>%
  summarise(
    N_lesions      = n(),
    N_LuminalA     = sum(MolecularSubtype == "Luminal A"),
    N_LuminalB     = sum(MolecularSubtype == "Luminal B"),
    N_HER2Enriched = sum(MolecularSubtype == "HER2-Enriched"),
    N_TripleNeg    = sum(MolecularSubtype == "Triple Negative"),
    Pct_TripleNeg  = round(mean(MolecularSubtype == "Triple Negative") * 100, 1),
    Pct_HER2pos    = round(mean(HER2_status == "HER2+") * 100, 1),
    Pct_ERpos      = round(mean(ER_status == "ER+") * 100, 1),
    Pct_PRpos      = round(mean(PR_status == "PR+") * 100, 1),
    Mean_Age       = round(mean(as.numeric(Age), na.rm = TRUE), 1),
    .groups = "drop"
  )
print(as.data.frame(pheno_stats))
write.csv(pheno_stats, file.path(OUTPUT_DIR, "Descriptive_Stats_per_Phenotype.csv"),
          row.names = FALSE)

cat(sprintf("\nOOB error rate: %.2f%%\n", oob_err))
cat(sprintf("Phenotype x subtype (%%):\n"))
print(round(prop.table(table(meta_data$RadiomicPheno,
                             meta_data$MolecularSubtype), margin = 1) * 100, 1))

# free memory no longer needed
rm(X_features_raw, nzv_idx, preproc_obj, merged_raw, radiomics_raw, clinical_raw)
invisible(gc())

message("\n[COMPLETE] Main-text pipeline (Sections 0-5) finished.")


# ====
# SECTION 6: PHENOTYPE FEATURE CHARACTERIZATION (R1-Results-2)
# Outputs: Fig16 + TableS1 (Kruskal-Wallis) + TableS2 (feature means)
# ====

message("\n=== SECTION 6: PHENOTYPE FEATURE CHARACTERIZATION ===\n")

feat_by_pheno <- X_top %>%
  mutate(RadiomicPheno = meta_data$RadiomicPheno) %>%
  group_by(RadiomicPheno) %>%
  summarise(across(everything(), mean), .groups = "drop")

feat_long <- feat_by_pheno %>%
  pivot_longer(-RadiomicPheno, names_to = "Feature", values_to = "MeanZ") %>%
  mutate(Feature = gsub("^original_", "", Feature))

kw_mat <- sapply(top_features, function(f) {
  kt <- suppressWarnings(kruskal.test(X_top[[f]], meta_data$RadiomicPheno))
  c(chi2 = unname(kt$statistic), p = kt$p.value)
})
fmt_p <- function(p) ifelse(p < 1e-300, "<1e-300", formatC(p, format = "e", digits = 2))
kw_p_raw <- unname(kw_mat["p", ])
kw_p_bh  <- p.adjust(kw_p_raw, method = "BH")
kw_df <- data.frame(
  Feature   = gsub("^original_", "", top_features),
  KW_chi2   = round(kw_mat["chi2", ], 2),
  P_raw     = signif(kw_p_raw, 3),
  P_BH      = signif(kw_p_bh, 3),
  P_raw_fmt = fmt_p(kw_p_raw),
  P_BH_fmt  = fmt_p(kw_p_bh),
  row.names = NULL
) %>% arrange(desc(KW_chi2))
print(kw_df)

write.csv(kw_df, file.path(OUTPUT_DIR, "TableS1_KruskalWallis_Features_by_Phenotype.csv"),
          row.names = FALSE)
write.csv(as.data.frame(feat_by_pheno),
          file.path(OUTPUT_DIR, "TableS2_Phenotype_Feature_Means.csv"), row.names = FALSE)

# REVISED: extra left margin so the first x-axis label is not clipped
p16 <- ggplot(feat_long, aes(x = Feature, y = RadiomicPheno, fill = MeanZ)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", MeanZ)), size = 4, color = "grey10") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Mean\nz-score") +
  labs(x = NULL, y = NULL) +
  theme_pub(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = ggplot2::margin(10, 15, 10, 70))
save_plot(p16, "Fig16_Phenotype_Feature_Heatmap", width = 12, height = 5.5)


# ====
# SECTION 7: ACQUISITION AUDIT & COMBAT HARMONIZATION (R1-M2 ; R2-1)
# ====

message("\n=== SECTION 7: ACQUISITION AUDIT & COMBAT ===\n")

# ComBat, parametric empirical Bayes (Johnson, Li & Rabinovic 2007),
# base-R implementation (no sva/edgeR installation needed).
combat_eb <- function(dat, batch, mod = NULL, eps = 1e-8) {
  stopifnot(is.matrix(dat), length(batch) == ncol(dat))
  if (anyNA(batch)) stop("[ComBat] NA in batch vector.")
  batch <- factor(batch)
  if (nlevels(batch) < 2L) stop("[ComBat] batch must have >= 2 levels.")
  
  n.array   <- ncol(dat)
  n.batches <- as.vector(table(batch))
  n.batch   <- length(n.batches)
  
  batchmod <- model.matrix(~ -1 + batch)
  design   <- if (is.null(mod)) batchmod else cbind(batchmod, mod)
  design   <- design[, !apply(design, 2, function(x) all(x == 1)), drop = FALSE]
  if (qr(design)$rank < ncol(design))
    stop("[ComBat] design not full rank (batch confounded with covariates).")
  batch_cols <- seq_len(ncol(batchmod))
  
  B.hat <- solve(crossprod(design), crossprod(design, t(dat)))
  grand.mean <- as.vector((n.batches / n.array) %*% B.hat[batch_cols, , drop = FALSE])
  
  stand.mean <- matrix(grand.mean, nrow = nrow(dat), ncol = n.array)
  nb <- setdiff(seq_len(ncol(design)), batch_cols)
  if (length(nb) > 0)
    stand.mean <- stand.mean + t(design[, nb, drop = FALSE] %*% B.hat[nb, , drop = FALSE])
  
  var.pooled <- pmax(as.vector(((dat - t(design %*% B.hat))^2) %*%
                                 rep(1 / n.array, n.array)), eps)
  s.data <- (dat - stand.mean) / sqrt(var.pooled)
  
  gamma.hat <- rowsum(t(s.data), batch) / n.batches
  delta.hat <- pmax((rowsum(t(s.data^2), batch) -
                       sweep(gamma.hat^2, 1, n.batches, `*`)) / (n.batches - 1), eps)
  
  gamma.bar <- rowMeans(gamma.hat)
  t2        <- pmax(apply(gamma.hat, 1, stats::var), eps)
  m_d       <- rowMeans(delta.hat)
  v_d       <- pmax(apply(delta.hat, 1, stats::var), eps)
  a.prior   <- m_d^2 / v_d + 2
  b.prior   <- m_d * (a.prior - 1)
  
  gamma.star <- (sweep(gamma.hat, 1, n.batches * t2, `*`) +
                   sweep(delta.hat, 1, gamma.bar, `*`)) /
    sweep(delta.hat, 1, n.batches * t2, `+`)
  
  delta.star <- delta.hat
  for (i in seq_len(n.batch)) {
    idx <- which(as.integer(batch) == i)
    dev <- s.data[, idx, drop = FALSE] - gamma.hat[i, ]
    delta.star[i, ] <- (b.prior[i] + 0.5 * rowSums(dev^2)) /
      (n.batches[i] / 2 + a.prior[i] - 1)
  }
  delta.star <- pmax(delta.star, eps)
  
  adj <- s.data
  for (i in seq_len(n.batch)) {
    idx <- which(as.integer(batch) == i)
    adj[, idx] <- (s.data[, idx, drop = FALSE] - gamma.star[i, ]) /
      sqrt(delta.star[i, ])
  }
  out <- adj * sqrt(var.pooled) + stand.mean
  dimnames(out) <- dimnames(dat)
  out
}

# --- 7.1 Acquisition audit ---
mk_counts <- function(x, label) {
  tb <- table(x, useNA = "ifany")
  data.frame(Parameter = label,
             Level     = names(tb),
             N         = as.integer(tb),
             Percent   = round(as.integer(tb) / sum(tb) * 100, 1))
}
acq_audit <- rbind(mk_counts(meta_data$Manufacturer_clean, "Manufacturer"),
                   mk_counts(meta_data$MagField_clean, "Magnetic field strength"))
acq_audit$Level[is.na(acq_audit$Level)] <- "Missing"
print(acq_audit)
write.csv(acq_audit, file.path(OUTPUT_DIR, "Acquisition_Parameter_Audit.csv"),
          row.names = FALSE)

# --- 7.2 Choose batch variable ---
batch_candidates <- list(
  "field strength x manufacturer" = interaction(meta_data$MagField_clean,
                                                meta_data$Manufacturer_clean, drop = TRUE),
  "magnetic field strength"       = meta_data$MagField_clean,
  "manufacturer"                  = meta_data$Manufacturer_clean
)

batch_var <- NULL; batch_label <- ""; batch_tab <- NULL
for (nm in names(batch_candidates)) {
  tb <- table(batch_candidates[[nm]])
  tb <- tb[tb > 0]
  if (length(tb) >= 2L && min(tb) >= 10L) {
    batch_var <- batch_candidates[[nm]]; batch_label <- nm; batch_tab <- tb
    break
  }
}

chi_before_vec <- sapply(chi_results, function(x) x$chi2)
v_before_vec   <- round(sapply(chi_results, function(x) cramer_v(x$table)), 3)

if (is.null(batch_var)) {
  # BRANCH A: single batch -> ComBat not applicable
  message("[FINDING] Single acquisition batch: ComBat NOT APPLICABLE.")
  tab3 <- data.frame(Variable = test_vars, Chi2_Before = chi_before_vec,
                     Chi2_After = NA_real_, Pct_Change = NA_real_,
                     V_Before = v_before_vec, V_After = NA_real_,
                     P_After = NA_real_, P_After_BH = NA_real_, N_After = NA_integer_,
                     Note = "ComBat not applicable: single acquisition batch",
                     row.names = NULL)
  write.csv(tab3, file.path(OUTPUT_DIR, "Table3_ChiSquare_Before_After_ComBat.csv"),
            row.names = FALSE)
  
  p17 <- ggplot(acq_audit, aes(x = Level, y = N, fill = Parameter)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, Percent)),
              vjust = -0.3, size = 5, fontface = "bold") +
    facet_wrap(~ Parameter, scales = "free_x") +
    scale_fill_manual(values = c("Manufacturer" = "#3949AB",
                                 "Magnetic field strength" = "#880E4F"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Lesions (n)") +
    theme_pub(16)
  save_plot(p17, "Fig17_Acquisition_Batch_Audit", width = 9, height = 5.5)
  
} else {
  # BRANCH B: ComBat feasible
  message(sprintf("[INFO] ComBat batch variable: %s", batch_label))
  ok <- !is.na(batch_var)
  batch_sub <- droplevels(factor(batch_var[ok]))
  meta_sub  <- meta_data[ok, ]
  
  X_unscaled_top <- as.data.frame(X_features_clean[ok, top_features])
  mod_bio <- model.matrix(~ MolecularSubtype, data = meta_sub)
  
  combat_mat <- combat_eb(dat = t(as.matrix(X_unscaled_top)),
                          batch = batch_sub, mod = mod_bio)
  X_combat <- as.data.frame(t(combat_mat))
  colnames(X_combat) <- top_features
  message("[OK] ComBat harmonization complete.")
  
  preproc_cb      <- preProcess(X_combat, method = c("center", "scale"))
  X_combat_scaled <- predict(preproc_cb, X_combat)
  
  set.seed(2026)
  som_combat <- som(X = as.matrix(X_combat_scaled), grid = som_grid, rlen = 500,
                    alpha = c(0.05, 0.01), radius = c(SOM_XDIM, 1), keep.data = TRUE)
  cb_nodes    <- som_combat$unit.classif
  cb_codebook <- getCodes(som_combat)
  cb_hclust   <- hclust(dist(cb_codebook, method = "euclidean"), method = "ward.D2")
  cb_node_ph  <- cutree(cb_hclust, k = N_PHENOTYPES)
  
  RPc <- rep(NA_character_, nrow(meta_data))
  RPc[ok] <- paste0("RPc-", cb_node_ph[cb_nodes])
  meta_data$RadiomicPheno_ComBat <- factor(RPc)
  print(table(meta_data$RadiomicPheno_ComBat))
  
  chi_after <- lapply(test_vars, function(v)
    chi_test(table(meta_data$RadiomicPheno_ComBat[ok], meta_sub[[v]])))
  names(chi_after) <- test_vars
  p_after_bh <- p.adjust(sapply(chi_after, function(x) x$p), method = "BH")
  
  compare_combat <- data.frame(
    Variable    = test_vars,
    Chi2_Before = chi_before_vec,
    Chi2_After  = round(sapply(chi_after, function(x) x$chi2), 3),
    Pct_Change  = round((sapply(chi_after, function(x) x$chi2) - chi_before_vec) /
                          chi_before_vec * 100, 1),
    V_Before    = v_before_vec,
    V_After     = round(sapply(chi_after, function(x) x$V), 3),
    P_After     = signif(sapply(chi_after, function(x) x$p), 3),
    P_After_BH  = signif(p_after_bh, 3),
    N_After     = sum(ok),
    row.names   = NULL
  )
  print(compare_combat)
  write.csv(compare_combat, file.path(OUTPUT_DIR, "Table3_ChiSquare_Before_After_ComBat.csv"),
            row.names = FALSE)
  
  comp_long <- compare_combat %>%
    select(Variable, Chi2_Before, Chi2_After) %>%
    pivot_longer(-Variable, names_to = "Stage", values_to = "Chi2") %>%
    filter(!is.na(Chi2)) %>%
    mutate(Stage = factor(Stage, levels = c("Chi2_Before", "Chi2_After"),
                          labels = c("Before ComBat", "After ComBat")))
  
  p17 <- ggplot(comp_long, aes(x = Variable, y = Chi2, fill = Stage)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    geom_text(aes(label = round(Chi2, 1)),
              position = position_dodge(width = 0.72),
              vjust = -0.4, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Before ComBat" = "#B71C1C",
                                 "After ComBat" = "#1565C0"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = expression(chi^2 ~ "statistic")) +
    theme_pub(16) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))
  save_plot(p17, "Fig17_ChiSquare_BeforeAfter_ComBat", width = 10, height = 6)
  
  pheno_colors_cb <- setNames(unname(phenotype_colors), paste0("RPc-", 1:N_PHENOTYPES))
  mfs_grid_cb   <- build_som_property_df(meta_sub$MagField_clean, node_vec = cb_nodes)
  p18a <- plot_som_grid(mfs_grid_cb,
                        c("1.5T" = "#0D47A1", "3T" = "#880E4F", "Empty" = "grey85"),
                        fill_name = "Dominant\nfield strength")
  
  pheno_grid_cb <- build_som_property_df(meta_data$RadiomicPheno_ComBat[ok],
                                         node_vec = cb_nodes)
  p18b <- ggplot(pheno_grid_cb, aes(x = col, y = row, fill = dominant)) +
    geom_tile(color = "white", linewidth = 1.2, width = 0.92, height = 0.92) +
    geom_text(aes(label = ifelse(n > 0, paste0(dominant, "\nn=", n), "")),
              color = "white", size = 4.5, fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = c(pheno_colors_cb, "Empty" = "grey85"),
                      name = "Phenotype\n(post-ComBat)") +
    labs(x = "SOM column", y = "SOM row") +
    coord_fixed() +
    theme_pub(16) +
    theme(panel.grid = element_blank())
  
  p18 <- gridExtra::arrangeGrob(p18a, p18b, ncol = 2, widths = c(1, 1))
  save_plot(p18, "Fig18_PostComBat_FieldStrength_and_Phenotypes", width = 15, height = 7.5)
  
  rm(combat_mat, X_combat, X_combat_scaled, som_combat, cb_nodes,
     cb_codebook, cb_hclust, cb_node_ph, X_unscaled_top)
  invisible(gc())
}


# ====
# SECTION 8: LABEL-FREE FEATURE-SELECTION CONTROLS (R1-M1 ; R2-2)
# Outputs: Fig19, Fig20, Table4, Laplacian scores CSV
# ====

message("\n=== SECTION 8: LABEL-FREE FEATURE-SELECTION CONTROLS ===\n")

laplacian_scores <- function(X, k = 5) {
  X <- as.matrix(X)
  n <- nrow(X)
  dmat <- as.matrix(dist(X))
  t_par <- median(dmat[dmat > 0])^2
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nn <- order(dmat[i, ])[2:(k + 1)]
    W[i, nn] <- exp(-dmat[i, nn]^2 / t_par)
  }
  W <- pmax(W, t(W))
  rm(dmat); gc()
  dvec <- rowSums(W)
  dtot <- sum(dvec)
  apply(X, 2, function(fr) {
    f1  <- fr - sum(fr * dvec) / dtot
    Lf  <- dvec * f1 - as.vector(W %*% f1)
    num <- sum(f1 * Lf)
    den <- sum(dvec * f1 * f1)
    if (den == 0) return(Inf)
    num / den
  })
}

message("[INFO] Computing Laplacian scores (~1-3 min)...")
lap_scores <- laplacian_scores(X_scaled, k = 5)
lap_df <- data.frame(Feature = colnames(X_scaled),
                     LaplacianScore = round(lap_scores, 6),
                     row.names = NULL) %>%
  arrange(LaplacianScore)
top_lap <- lap_df$Feature[1:TOP_N_FEATURES]
print(lap_df[1:TOP_N_FEATURES, ])
write.csv(lap_df, file.path(OUTPUT_DIR, "Laplacian_Scores_All_Features.csv"),
          row.names = FALSE)

var_df <- data.frame(
  Feature  = colnames(X_features_clean),
  Variance = sapply(X_features_clean, function(col) var(col, na.rm = TRUE)),
  row.names = NULL
) %>% arrange(desc(Variance))
top_var <- var_df$Feature[1:TOP_N_FEATURES]
print(var_df[1:TOP_N_FEATURES, ])

run_som_pipeline <- function(feat_set, Xmat_scaled, seed, tag) {
  message(sprintf("[INFO] Control SOM run: %s", tag))
  set.seed(seed)
  m  <- som(as.matrix(Xmat_scaled[, feat_set]), grid = som_grid, rlen = 500,
            alpha = c(0.05, 0.01), radius = c(SOM_XDIM, 1), keep.data = TRUE)
  nd <- m$unit.classif
  cb <- getCodes(m)
  hc <- hclust(dist(cb, method = "euclidean"), method = "ward.D2")
  ph <- factor(paste0("RP-", cutree(hc, k = N_PHENOTYPES)[nd]))
  res_list <- list()
  for (v in test_vars) {
    r <- chi_test(table(ph, meta_data[[v]]))
    res_list[[v]] <- c(chi2 = r$chi2, p = r$p, V = r$V)
  }
  list(model = m, pheno = ph, res = as.data.frame(do.call(rbind, res_list)))
}

lap_res <- run_som_pipeline(top_lap, X_scaled, seed = 2028, tag = "Laplacian score")
var_res <- run_som_pipeline(top_var, X_scaled, seed = 2029, tag = "Maximum variance")

fs_compare <- data.frame(
  Variable         = test_vars,
  Chi2_RF_guided   = sapply(chi_results, function(x) x$chi2),
  Chi2_Laplacian   = round(lap_res$res$chi2, 3),
  Chi2_MaxVariance = round(var_res$res$chi2, 3),
  P_Laplacian      = signif(lap_res$res$p, 3),
  P_MaxVariance    = signif(var_res$res$p, 3),
  V_RF_guided      = round(sapply(chi_results, function(x) cramer_v(x$table)), 3),
  V_Laplacian      = round(lap_res$res$V, 3),
  V_MaxVariance    = round(var_res$res$V, 3),
  row.names = NULL
)
print(fs_compare)
write.csv(fs_compare, file.path(OUTPUT_DIR, "Table4_FeatureSelection_Control_ChiSquare.csv"),
          row.names = FALSE)

fs_long <- bind_rows(
  data.frame(Variable = test_vars, Method = "RF-guided (main)",
             V = sapply(chi_results, function(x) cramer_v(x$table))),
  data.frame(Variable = test_vars, Method = "Laplacian score", V = lap_res$res$V),
  data.frame(Variable = test_vars, Method = "Maximum variance", V = var_res$res$V)
) %>%
  filter(!is.na(V)) %>%
  mutate(Method = factor(Method, levels = c("RF-guided (main)",
                                            "Laplacian score", "Maximum variance")))

p19 <- ggplot(fs_long, aes(x = Variable, y = V, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.2f", V)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c("RF-guided (main)" = "#B71C1C",
                               "Laplacian score" = "#1565C0",
                               "Maximum variance" = "#388E3C"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Cramer's V (association strength)") +
  theme_pub(16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p19, "Fig19_FeatureSelection_Control_Comparison", width = 10, height = 6)

lap_top20 <- lap_df[1:min(20, nrow(lap_df)), ] %>%
  mutate(
    Feature_short = gsub("^original_", "", Feature),
    Feature_short = factor(Feature_short, levels = rev(Feature_short)),
    Selected = ifelse(Feature %in% top_lap, "Selected Top 15", "Other")
  )

p20 <- ggplot(lap_top20, aes(x = LaplacianScore, y = Feature_short, fill = Selected)) +
  geom_col(width = 0.7, color = "white") +
  scale_fill_manual(values = c("Selected Top 15" = "#1565C0", "Other" = "#90A4AE"),
                    name = "Feature selection") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Laplacian score (lower = better)", y = NULL) +
  theme_pub(16) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        plot.margin = ggplot2::margin(10, 20, 10, 10))
save_plot(p20, "Fig20_Laplacian_Score_Ranking", width = 11, height = 8)


# ====
# SECTION 9: PATIENT-LEVEL SENSITIVITY (R2-3)
# Outputs: Fig21, Table5, Table5b (excluding dominant patient)
# ====

message("\n=== SECTION 9: PATIENT-LEVEL SENSITIVITY ===\n")

meta_data$PatientID <- sub("-[^-]*$", "", as.character(meta_data$FilenamePrefix))
n_patients <- length(unique(meta_data$PatientID))
lesions_per_patient <- table(meta_data$PatientID)

cat(sprintf("\nUnique patients: %d | lesions/patient median %.0f, IQR %.0f-%.0f, max %.0f\n",
            n_patients, median(lesions_per_patient),
            quantile(lesions_per_patient, 0.25), quantile(lesions_per_patient, 0.75),
            max(lesions_per_patient)))

mode_fn <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

patient_level <- meta_data %>%
  group_by(PatientID) %>%
  summarise(
    N_lesions         = n(),
    N_phenotypes      = length(unique(RadiomicPheno)),
    Purity            = max(table(RadiomicPheno)) / n(),
    DominantPhenotype = mode_fn(as.character(RadiomicPheno)),
    MolecularSubtype  = mode_fn(as.character(MolecularSubtype)),
    ER_status         = mode_fn(as.character(ER_status)),
    PR_status         = mode_fn(as.character(PR_status)),
    HER2_status       = mode_fn(as.character(HER2_status)),
    MagField_clean    = mode_fn(as.character(MagField_clean)),
    .groups = "drop"
  ) %>%
  mutate(DominantPhenotype = factor(DominantPhenotype,
                                    levels = paste0("RP-", 1:N_PHENOTYPES)))

pct_multi <- round(mean(patient_level$N_phenotypes > 1) * 100, 1)
cat(sprintf("Patients with >1 phenotype: %.1f%% | mean purity: %.2f\n",
            pct_multi, mean(patient_level$Purity)))
write.csv(as.data.frame(patient_level),
          file.path(OUTPUT_DIR, "PatientLevel_DominantPhenotype.csv"), row.names = FALSE)

pat_tests <- lapply(test_vars, function(v)
  chi_test(table(patient_level$DominantPhenotype, patient_level[[v]])))
names(pat_tests) <- test_vars
p_pat_bh <- p.adjust(sapply(pat_tests, function(x) x$p), method = "BH")

patient_compare <- data.frame(
  Variable     = test_vars,
  Chi2_Lesion  = sapply(chi_results, function(x) x$chi2),
  Chi2_Patient = round(sapply(pat_tests, function(x) x$chi2), 3),
  DF_Patient   = sapply(pat_tests, function(x) x$df),
  P_Patient    = signif(sapply(pat_tests, function(x) x$p), 3),
  P_Patient_BH = signif(p_pat_bh, 3),
  V_Lesion     = round(sapply(chi_results, function(x) cramer_v(x$table)), 3),
  V_Patient    = round(sapply(pat_tests, function(x) x$V), 3),
  N_Patients   = n_patients,
  row.names = NULL
)
print(patient_compare)
write.csv(patient_compare, file.path(OUTPUT_DIR, "Table5_PatientLevel_Sensitivity.csv"),
          row.names = FALSE)

sens_long <- patient_compare %>%
  select(Variable, V_Lesion, V_Patient) %>%
  pivot_longer(-Variable, names_to = "Level", values_to = "V") %>%
  filter(!is.na(V)) %>%
  mutate(Level = factor(Level, levels = c("V_Lesion", "V_Patient"),
                        labels = c(sprintf("Lesion level (n=%d)", nrow(meta_data)),
                                   sprintf("Patient level (n=%d)", n_patients))))

p21 <- ggplot(sens_long, aes(x = Variable, y = V, fill = Level)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", V)),
            position = position_dodge(width = 0.72),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("#EF6C00", "#00695C"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Cramer's V (association strength)") +
  theme_pub(16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p21, "Fig21_PatientLevel_vs_LesionLevel", width = 10, height = 6)

# --- 9.1 Sensitivity excluding the dominant patient (Issue B) ---
big_pt <- names(which.max(table(meta_data$PatientID)))
keep   <- meta_data$PatientID != big_pt
cat(sprintf("\nExcluding dominant patient %s: %d lesions removed, %d remain\n",
            big_pt, sum(!keep), sum(keep)))

excl_res <- lapply(test_vars, function(v)
  chi_test(table(meta_data$RadiomicPheno[keep], meta_data[[v]][keep])))
names(excl_res) <- test_vars

table5b <- data.frame(
  Variable        = test_vars,
  Chi2_Full       = chi_before_vec,
  Chi2_ExclDom    = round(sapply(excl_res, function(x) x$chi2), 3),
  DF_ExclDom      = sapply(excl_res, function(x) x$df),
  P_ExclDom       = signif(sapply(excl_res, function(x) x$p), 3),
  V_Full          = v_before_vec,
  V_ExclDom       = round(sapply(excl_res, function(x) x$V), 3),
  N_Lesions_Used  = sum(keep),
  row.names = NULL
)
print(table5b)
write.csv(table5b, file.path(OUTPUT_DIR, "Table5b_Sensitivity_ExclDominantPatient.csv"),
          row.names = FALSE)


# ====
# SECTION 10 (NEW): ADDITIONAL ANALYSES FOR THE REMAINING REVIEWER ISSUES
# All blocks are independent and wrapped in safe_run().
# Approximate run time on the free tier: 15-30 minutes (mostly 10.6 and 10.7).
# To shorten: set N_STAB and N_SOM_SEEDS to 5.
# ====

message("\n=== SECTION 10: ADDITIONAL ANALYSES ===\n")

safe_run <- function(label, expr) {
  tryCatch(expr, error = function(e)
    message(sprintf("[WARN] Block '%s' failed: %s", label, conditionMessage(e))))
}
save_csv <- function(df, name) {
  write.csv(df, file.path(OUTPUT_DIR, name), row.names = FALSE)
  message("[CSV SAVED] ", name)
}
ari <- function(a, b) {   # adjusted Rand index (no extra packages)
  tb  <- table(a, b); n <- sum(tb)
  sc  <- function(x) sum(choose(x, 2))
  si  <- sc(rowSums(tb)); sj <- sc(colSums(tb)); sij <- sc(tb)
  ex  <- si * sj / choose(n, 2)
  (sij - ex) / (0.5 * (si + sj) - ex)
}
chi_stat <- function(x, y) {   # fast chi-square statistic only
  tb <- table(x, y)
  tb <- tb[rowSums(tb) > 0, colSums(tb) > 0, drop = FALSE]
  if (nrow(tb) < 2L || ncol(tb) < 2L) return(NA_real_)
  as.numeric(suppressWarnings(chisq.test(tb)$statistic))
}
bio_vars  <- c("MolecularSubtype", "ER_status", "PR_status", "HER2_status")
mtry_main <- floor(sqrt(ncol(X_scaled)))
pat_names <- patient_level$PatientID
pat_idx   <- match(meta_data$PatientID, pat_names)

# ---- 10.1 Cohort descriptives, RF baseline and per-class performance ----
safe_run("10.1 cohort descriptives", {
  sub_tab <- table(meta_data$MolecularSubtype)
  sub_counts <- data.frame(Subtype = names(sub_tab), N = as.integer(sub_tab),
                           Percent = round(as.integer(sub_tab) / sum(sub_tab) * 100, 1))
  cat("\n--- Subtype counts (add to Methods / Supplement) ---\n"); print(sub_counts)
  save_csv(sub_counts, "S_Subtype_Counts.csv")
  
  lpp <- as.integer(table(meta_data$PatientID))
  lpp_sum <- data.frame(
    N_patients = length(lpp), Median = median(lpp),
    Q1 = unname(quantile(lpp, 0.25)), Q3 = unname(quantile(lpp, 0.75)),
    Min = min(lpp), Max = max(lpp),
    Max_share_pct = round(max(lpp) / sum(lpp) * 100, 1))
  cat("\n--- Lesions per patient ---\n"); print(lpp_sum)
  save_csv(lpp_sum, "S_Lesions_Per_Patient.csv")
  
  cm <- rf_model$confusion[, seq_len(nlevels(Y_subtype)), drop = FALSE]
  recall <- diag(cm) / rowSums(cm)
  base_acc <- max(prop.table(table(Y_subtype)))
  rf_perf <- data.frame(
    Metric = c("RF OOB accuracy (lesion level)", "RF OOB balanced accuracy",
               "Majority-class baseline accuracy"),
    Value  = round(c(sum(diag(cm)) / sum(cm), mean(recall), base_acc), 4))
  cat("\n--- RF performance vs. baseline ---\n"); print(rf_perf)
  save_csv(rf_perf, "S_RF_OOB_vs_Baseline.csv")
  save_csv(data.frame(Class = rownames(cm), N = rowSums(cm),
                      OOB_Recall = round(recall, 3)), "S_RF_OOB_PerClass.csv")
})

# ---- 10.2 Patients per phenotype and the dominant patient ----
safe_run("10.2 patients per phenotype", {
  pheno_pat <- meta_data %>%
    group_by(RadiomicPheno) %>%
    summarise(
      N_lesions = n(),
      N_patients = n_distinct(PatientID),
      Top_patient_lesions = max(table(PatientID)),
      Top_patient_share_pct = round(max(table(PatientID)) / n() * 100, 1),
      Effective_N_patients = round(1 / sum((as.numeric(table(PatientID)) / n())^2), 1),
      .groups = "drop")
  cat("\n--- Patients per phenotype ---\n"); print(as.data.frame(pheno_pat))
  save_csv(as.data.frame(pheno_pat), "S_Patients_Per_Phenotype.csv")
  
  is_dom <- meta_data$PatientID == big_pt
  dom_df <- rbind(
    data.frame(Variable = "Phenotype", as.data.frame(table(meta_data$RadiomicPheno[is_dom]),
                                                     stringsAsFactors = FALSE)),
    data.frame(Variable = "Subtype", as.data.frame(table(meta_data$MolecularSubtype[is_dom]),
                                                   stringsAsFactors = FALSE)),
    data.frame(Variable = "FieldStrength", as.data.frame(table(meta_data$MagField_clean[is_dom]),
                                                         stringsAsFactors = FALSE)))
  names(dom_df) <- c("Variable", "Level", "N_lesions")
  cat(sprintf("\n--- Dominant patient %s ---\n", big_pt)); print(dom_df)
  save_csv(dom_df, "S_Dominant_Patient_Profile.csv")
  
  top5 <- names(sort(table(meta_data$PatientID), decreasing = TRUE))[1:5]
  sel  <- meta_data$PatientID %in% top5
  top5_tab <- as.data.frame.matrix(table(meta_data$PatientID[sel],
                                         meta_data$RadiomicPheno[sel]))
  top5_tab <- cbind(Patient = rownames(top5_tab), top5_tab)
  print(top5_tab)
  save_csv(top5_tab, "S_Top5_Patients_Phenotype_Composition.csv")
  
  fig22_df <- meta_data %>%
    mutate(Source = ifelse(PatientID == big_pt, "Dominant patient", "Other patients")) %>%
    dplyr::count(RadiomicPheno, Source)
  p22 <- ggplot(fig22_df, aes(x = RadiomicPheno, y = n, fill = Source)) +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 5) +
    scale_fill_manual(values = c("Dominant patient" = "#B71C1C",
                                 "Other patients" = "#90A4AE"), name = NULL) +
    labs(x = "Radiomic phenotype", y = "Lesions (n)") +
    theme_pub(16) + theme(legend.position = "bottom")
  save_plot(p22, "Fig22_Phenotype_Composition_by_DominantPatient", width = 8, height = 6)
})

# ---- 10.3 Is field strength itself associated with biology? (confounding) ----
safe_run("10.3 field strength vs biology", {
  rows <- lapply(bio_vars, function(v) {
    r1 <- chi_test(table(meta_data$MagField_clean, meta_data[[v]]))
    r2 <- chi_test(table(patient_level$MagField_clean, patient_level[[v]]))
    rbind(
      data.frame(Variable = v, Level = "Lesion", N = nrow(meta_data),
                 Chi2 = round(r1$chi2, 2), DF = r1$df, P = signif(r1$p, 3), V = round(r1$V, 3)),
      data.frame(Variable = v, Level = "Patient", N = nrow(patient_level),
                 Chi2 = round(r2$chi2, 2), DF = r2$df, P = signif(r2$p, 3), V = round(r2$V, 3)))
  })
  fs_bio <- do.call(rbind, rows)
  cat("\n--- Field strength vs biology (confounding path) ---\n"); print(fs_bio)
  save_csv(fs_bio, "S_FieldStrength_vs_Biology.csv")
})

# ---- 10.4 Phenotype vs biology WITHIN each field strength (+ CMH test) ----
safe_run("10.4 stratified analysis", {
  rows <- list()
  for (fs in levels(meta_data$MagField_clean)) {
    for (scope in c("All lesions", "Excluding dominant patient")) {
      idx <- which(meta_data$MagField_clean == fs)
      if (scope != "All lesions") idx <- idx[meta_data$PatientID[idx] != big_pt]
      for (v in bio_vars) {
        r <- chi_test(table(meta_data$RadiomicPheno[idx], meta_data[[v]][idx]))
        rows[[length(rows) + 1]] <- data.frame(
          FieldStrength = fs, Scope = scope, Variable = v, N = length(idx),
          Chi2 = round(r$chi2, 2), DF = r$df, P = signif(r$p, 3), V = round(r$V, 3))
      }
    }
  }
  strat <- do.call(rbind, rows)
  cat("\n--- Phenotype vs biology within each field strength ---\n"); print(strat)
  save_csv(strat, "S_Stratified_by_FieldStrength.csv")
  
  cmh <- do.call(rbind, lapply(bio_vars, function(v) {
    arr <- table(meta_data$RadiomicPheno, meta_data[[v]], meta_data$MagField_clean)
    r <- tryCatch(mantelhaen.test(arr), error = function(e) NULL)
    data.frame(Variable = v,
               CMH_chi2 = if (is.null(r)) NA else round(unname(r$statistic), 2),
               DF = if (is.null(r)) NA else unname(r$parameter),
               P = if (is.null(r)) NA else signif(r$p.value, 3))
  }))
  cat("\n--- Cochran-Mantel-Haenszel (phenotype x biology | field strength) ---\n")
  print(cmh)
  save_csv(cmh, "S_CMH_Phenotype_Biology_given_FieldStrength.csv")
})

# ---- 10.5 Cluster-robust inference (patients, not lesions, are the units) ----
safe_run("10.5 cluster-robust tests", {
  B_PERM <- 2000; B_BOOT <- 1000
  set.seed(2040)
  lesion_idx <- split(seq_len(nrow(meta_data)), meta_data$PatientID)
  tt_pat <- table(meta_data$PatientID)
  w_les  <- 1 / as.numeric(tt_pat[meta_data$PatientID])
  rows <- lapply(test_vars, function(v) {
    lab     <- as.character(patient_level[[v]])
    lab_les <- lab[pat_idx]
    obs  <- chi_stat(meta_data$RadiomicPheno, lab_les)
    null <- replicate(B_PERM, chi_stat(meta_data$RadiomicPheno, sample(lab)[pat_idx]))
    p_perm <- (1 + sum(null >= obs, na.rm = TRUE)) / (B_PERM + 1)
    vb <- replicate(B_BOOT, {
      ii <- unlist(lesion_idx[sample(names(lesion_idx), replace = TRUE)], use.names = FALSE)
      cramer_v(table(meta_data$RadiomicPheno[ii], meta_data[[v]][ii]))
    })
    wdf <- data.frame(w = w_les, Pheno = meta_data$RadiomicPheno, X = meta_data[[v]])
    wt  <- xtabs(w ~ Pheno + X, data = wdf)
    wt  <- wt * (nrow(patient_level) / sum(wt))
    data.frame(
      Variable = v,
      V_lesion = round(cramer_v(table(meta_data$RadiomicPheno, meta_data[[v]])), 3),
      V_boot_lo = round(quantile(vb, 0.025, na.rm = TRUE), 3),
      V_boot_hi = round(quantile(vb, 0.975, na.rm = TRUE), 3),
      Chi2_obs_patient_labels = round(obs, 2),
      P_patient_permutation = signif(p_perm, 3),
      V_patient_weighted = round(suppressWarnings(cramer_v(wt)), 3))
  })
  cr <- do.call(rbind, rows)
  cr$P_perm_BH <- signif(p.adjust(cr$P_patient_permutation, "BH"), 3)
  cat("\n--- Cluster-robust tests (patient-level permutation, cluster bootstrap CI) ---\n")
  print(cr)
  save_csv(cr, "S_ClusterRobust_Tests.csv")
})

# ---- 10.6 RF robustness: patient-grouped CV and feature-selection stability ----
safe_run("10.6a grouped CV", {
  set.seed(2031)
  pts     <- unique(meta_data$PatientID)
  fold_pt <- setNames(sample(rep(1:5, length.out = length(pts))), pts)
  fold    <- as.integer(fold_pt[meta_data$PatientID])
  pred_cv <- factor(rep(NA_character_, nrow(X_scaled)), levels = levels(Y_subtype))
  for (k in 1:5) {
    tr <- fold != k
    m  <- randomForest(x = X_scaled[tr, ], y = Y_subtype[tr], ntree = 500, mtry = mtry_main)
    pred_cv[!tr] <- predict(m, X_scaled[!tr, ])
  }
  cm_cv <- table(True = Y_subtype, Pred = pred_cv)
  cv_perf <- data.frame(
    Metric = c("Patient-grouped 5-fold CV accuracy",
               "Patient-grouped 5-fold CV balanced accuracy",
               "Majority-class baseline accuracy"),
    Value  = round(c(sum(diag(cm_cv)) / sum(cm_cv),
                     mean(diag(cm_cv) / rowSums(cm_cv)),
                     max(prop.table(table(Y_subtype)))), 4))
  cat("\n--- RF with patient-grouped cross-validation ---\n"); print(cv_perf); print(cm_cv)
  save_csv(cv_perf, "S_RF_PatientGrouped_CV.csv")
})

safe_run("10.6b selection stability", {
  N_STAB <- 10
  rank_mat <- sapply(seq_len(N_STAB), function(i) {
    set.seed(2100 + i)
    m   <- randomForest(x = X_scaled, y = Y_subtype, ntree = 500,
                        mtry = mtry_main, importance = TRUE)
    imp <- importance(m, type = 1)[, "MeanDecreaseAccuracy"]
    rank(-imp[colnames(X_scaled)], ties.method = "min")
  })
  stab <- data.frame(
    Feature = rownames(rank_mat),
    Mean_rank = round(rowMeans(rank_mat), 1),
    Selection_frequency = round(rowMeans(rank_mat <= TOP_N_FEATURES), 2),
    In_main_top15 = rownames(rank_mat) %in% top_features,
    row.names = NULL) %>% arrange(Mean_rank)
  cat(sprintf("\n--- RF selection stability over %d seeds (top 25) ---\n", N_STAB))
  print(head(stab, 25))
  message(sprintf("[INFO] Main top-15 features selected in >= 80%% of runs: %d of %d",
                  sum(stab$In_main_top15 & stab$Selection_frequency >= 0.8), TOP_N_FEATURES))
  save_csv(stab, "S_RF_Selection_Stability.csv")
})

# ---- 10.7 SOM robustness: seeds, k, quantization / topographic error, outliers ----
safe_run("10.7a SOM seed stability", {
  N_SOM_SEEDS <- 10
  rows <- lapply(seq_len(N_SOM_SEEDS), function(i) {
    set.seed(3000 + i)
    m  <- som(X_som_matrix, grid = som_grid, rlen = 500, alpha = c(0.05, 0.01),
              radius = c(SOM_XDIM, 1), keep.data = TRUE)
    cb <- getCodes(m)
    ph <- cutree(hclust(dist(cb), "ward.D2"), k = N_PHENOTYPES)[m$unit.classif]
    data.frame(Seed = 3000 + i,
               ARI_vs_main = round(ari(ph, obs_phenotype), 3),
               V_subtype = round(cramer_v(table(ph, meta_data$MolecularSubtype)), 3),
               V_fieldstrength = round(cramer_v(table(ph, meta_data$MagField_clean)), 3))
  })
  seed_tab <- do.call(rbind, rows)
  cat("\n--- SOM seed stability (ARI vs the main run) ---\n"); print(seed_tab)
  save_csv(seed_tab, "S_SOM_Seed_Stability.csv")
})

safe_run("10.7b k sensitivity", {
  have_cl <- requireNamespace("cluster", quietly = TRUE)
  rows <- lapply(3:6, function(k) {
    cl  <- cutree(hclust_obj, k = k)
    lab <- cl[node_assignments]
    data.frame(k = k,
               Mean_silhouette_codebook = if (have_cl)
                 round(mean(cluster::silhouette(cl, codebook_dist)[, 3]), 3) else NA,
               ARI_vs_k4 = round(ari(lab, obs_phenotype), 3),
               V_subtype = round(cramer_v(table(lab, meta_data$MolecularSubtype)), 3),
               V_HER2 = round(cramer_v(table(lab, meta_data$HER2_status)), 3),
               V_fieldstrength = round(cramer_v(table(lab, meta_data$MagField_clean)), 3))
  })
  ktab <- do.call(rbind, rows)
  cat("\n--- Sensitivity to the number of phenotypes (k) ---\n"); print(ktab)
  save_csv(ktab, "S_k_Sensitivity.csv")
})

safe_run("10.7c SOM quality and outliers", {
  bmu_d <- sqrt(rowSums((X_som_matrix - codebook_vectors[node_assignments, , drop = FALSE])^2))
  sq    <- outer(rowSums(X_som_matrix^2), rowSums(codebook_vectors^2), "+") -
    2 * X_som_matrix %*% t(codebook_vectors)
  two   <- t(apply(sq, 1, function(r) order(r)[1:2]))
  adj   <- abs(unit_d[cbind(two[, 1], two[, 2])] - 1) < 1e-6
  maxabs <- apply(abs(X_som_matrix), 1, max)
  node_hi <- grid_pos$node[which.max(grid_pos$qe)]
  in_hi   <- node_assignments == node_hi
  qual <- data.frame(
    Metric = c("Quantization error (mean Euclidean distance to BMU, z-units)",
               "Topographic error (2nd BMU not adjacent)",
               "Share of lesions with any |z| > 4 (all)",
               sprintf("Node %d (highest mean distance): lesions", node_hi),
               sprintf("Node %d: patients", node_hi),
               sprintf("Node %d: share from the top patient (%%)", node_hi),
               sprintf("Node %d: share of lesions with any |z| > 4 (%%)", node_hi)),
    Value = round(c(mean(bmu_d), mean(!adj), mean(maxabs > 4),
                    sum(in_hi), length(unique(meta_data$PatientID[in_hi])),
                    100 * max(table(meta_data$PatientID[in_hi])) / sum(in_hi),
                    100 * mean(maxabs[in_hi] > 4)), 4))
  cat("\n--- SOM quality diagnostics ---\n"); print(qual)
  save_csv(qual, "S_SOM_Quality.csv")
  
  skew <- sapply(as.data.frame(X_top), function(x) mean((x - mean(x))^3) / sd(x)^3)
  tail_share <- sapply(as.data.frame(X_top), function(x) mean(abs(x) > 4))
  dist_diag <- data.frame(Feature = sub("^original_", "", names(skew)),
                          Skewness = round(skew, 2),
                          Share_abs_z_gt_4 = round(tail_share, 4), row.names = NULL)
  print(dist_diag)
  save_csv(dist_diag, "S_Feature_Distribution_Diagnostics.csv")
})

# ---- 10.8 Alternative pipelines (does the finding survive?) ----
safe_run("10.8 sensitivity pipelines", {
  collect <- function(res, tag, n_feat) {
    data.frame(Pipeline = tag, Variable = test_vars, N_features = n_feat,
               Chi2 = round(res$res$chi2, 2), P = signif(res$res$p, 3),
               V = round(res$res$V, 3),
               ARI_vs_main = round(ari(as.character(res$pheno), obs_phenotype), 3),
               row.names = NULL)
  }
  out <- list(data.frame(Pipeline = "Main (RF top 15)", Variable = test_vars,
                         N_features = TOP_N_FEATURES,
                         Chi2 = round(as.numeric(chi_before_vec), 2),
                         P = signif(p_vals_raw, 3), V = as.numeric(v_before_vec),
                         ARI_vs_main = 1, row.names = NULL))
  
  # (a) decorrelated: drop features with |r| > 0.9 to an already kept (higher-MDA) feature
  cmat <- abs(cor(X_top)); kept <- character(0)
  for (f in top_features)
    if (length(kept) == 0 || max(cmat[f, kept]) <= 0.9) kept <- c(kept, f)
  message(sprintf("[INFO] Decorrelated set keeps %d of %d features", length(kept), TOP_N_FEATURES))
  out[[length(out) + 1]] <- collect(
    run_som_pipeline(kept, X_scaled, 2033, "Decorrelated (|r| <= 0.9)"),
    "Decorrelated |r|<=0.9", length(kept))
  
  # (b) rank-based inverse-normal transform of the 15 features (heavy tails)
  X_rint <- as.data.frame(lapply(as.data.frame(X_top), function(x)
    qnorm((rank(x, ties.method = "average") - 0.5) / length(x))))
  out[[length(out) + 1]] <- collect(
    run_som_pipeline(colnames(X_rint), X_rint, 2032, "Rank-INT transformed"),
    "Rank-based inverse normal", ncol(X_rint))
  
  # (c) Laplacian control without ngtdm_Coarseness
  top_lap2 <- setdiff(lap_df$Feature, "original_ngtdm_Coarseness")[1:TOP_N_FEATURES]
  out[[length(out) + 1]] <- collect(
    run_som_pipeline(top_lap2, X_scaled, 2034, "Laplacian without Coarseness"),
    "Laplacian without Coarseness", TOP_N_FEATURES)
  
  # (d) ComBat on ALL 107 features, then RF selection, then SOM
  if (exists("batch_sub") && exists("mod_bio") && length(batch_sub) == nrow(meta_data)) {
    combat_all <- combat_eb(dat = t(as.matrix(X_features_clean)),
                            batch = batch_sub, mod = mod_bio)
    X_h <- as.data.frame(t(combat_all))
    X_h_scaled <- predict(preProcess(X_h, method = c("center", "scale")), X_h)
    set.seed(2035)
    rf_h <- randomForest(x = X_h_scaled, y = Y_subtype, ntree = 1000,
                         mtry = floor(sqrt(ncol(X_h_scaled))), importance = TRUE)
    imp_h <- importance(rf_h, type = 1)[, "MeanDecreaseAccuracy"]
    top_h <- names(sort(imp_h, decreasing = TRUE))[1:TOP_N_FEATURES]
    message(sprintf("[INFO] ComBat-then-select: %d of %d features overlap with the main top 15",
                    length(intersect(top_h, top_features)), TOP_N_FEATURES))
    out[[length(out) + 1]] <- collect(
      run_som_pipeline(top_h, X_h_scaled, 2036, "ComBat all features, RF, SOM"),
      "ComBat(107) -> RF -> SOM", TOP_N_FEATURES)
  }
  
  sens <- do.call(rbind, out)
  cat("\n--- Sensitivity pipelines ---\n"); print(sens)
  save_csv(sens, "S_Sensitivity_Pipelines.csv")
})

# ---- 10.9 Lesion size: are phenotypes just lesion volume? ----
safe_run("10.9 lesion size", {
  vol <- X_features_clean[["original_shape_VoxelVolume"]]
  size_by_pheno <- data.frame(Phenotype = meta_data$RadiomicPheno, Volume = vol) %>%
    group_by(Phenotype) %>%
    summarise(N = n(), Median = median(Volume),
              Q1 = quantile(Volume, 0.25), Q3 = quantile(Volume, 0.75), .groups = "drop")
  kw_vol <- kruskal.test(vol, meta_data$RadiomicPheno)
  cat("\n--- Lesion volume by phenotype ---\n"); print(as.data.frame(size_by_pheno))
  cat(sprintf("Kruskal-Wallis chi2 = %.1f, p = %.3g\n", kw_vol$statistic, kw_vol$p.value))
  save_csv(as.data.frame(size_by_pheno), "S_LesionVolume_by_Phenotype.csv")
  
  rho <- sapply(top_features, function(f)
    suppressWarnings(cor(X_features_clean[[f]], vol, method = "spearman")))
  rho_df <- data.frame(Feature = sub("^original_", "", top_features),
                       Spearman_rho_with_volume = round(rho, 3), row.names = NULL)
  print(rho_df)
  save_csv(rho_df, "S_Feature_vs_Volume.csv")
  
  tert <- cut(vol, quantile(vol, 0:3 / 3), include.lowest = TRUE,
              labels = c("Small", "Medium", "Large"))
  rows <- lapply(levels(tert), function(tl) {
    idx <- which(tert == tl)
    do.call(rbind, lapply(c("MolecularSubtype", "HER2_status", "MagField_clean"), function(v) {
      r <- chi_test(table(meta_data$RadiomicPheno[idx], meta_data[[v]][idx]))
      data.frame(VolumeTertile = tl, Variable = v, N = length(idx),
                 Chi2 = round(r$chi2, 2), P = signif(r$p, 3), V = round(r$V, 3))
    }))
  })
  tert_tab <- do.call(rbind, rows)
  cat("\n--- Phenotype associations within lesion-volume tertiles ---\n"); print(tert_tab)
  save_csv(tert_tab, "S_Stratified_by_LesionVolume.csv")
})


# ====
# FINAL WRAP-UP
# ====

all_files <- list.files(OUTPUT_DIR, pattern = "\\.png$|\\.csv$", full.names = FALSE)
cat("\n", strrep("=", 70), "\n")
cat("  OUTPUT FILES (", length(all_files), "):\n", sep = "")
cat(strrep("=", 70), "\n")
for (f in sort(all_files)) cat(sprintf("  %s\n", f))
cat(sprintf("\n  Saved to: %s\n", OUTPUT_DIR))

message("\n[COMPLETE] Full pipeline finished.")