# Gene-Level Burden Analysis — Bayesian (Jeffreys prior) with fractional control counts
# Modified: log10 OR plotting, safe p-floor, capped OR, limited labeling, CSV/PNG names adjusted
rm(list = ls())

## -------------- CONFIG --------------
SEED <- 12345
set.seed(SEED)

NSIM_GLOBAL <- 1e6
CHUNK_SIZE <- 5e4
PRIOR_ALPHA <- 0.5
PRIOR_BETA  <- 0.5
PRIOR_ALPHA_CTRL <- 0.5
PRIOR_BETA_CTRL  <- 0.5
CLIP_EPS <- 1e-12
VERBOSE <- TRUE

# plotting / labeling params
TOP_LABEL_N <- 22          # number of FDR-significant genes to label (top by -log10 FDR)
LOG10OR_CAP <- 6           # cap for log10(OR) in plot (±LOG10OR_CAP)
NEGLOG10FDR_CAP <- 300     # cap for -log10(FDR) (very large values)
MIN_NONZERO_P <- 1e-300    # safe floor for p-values equal to 0

# output directory constant
RESULTS_DIR <- "/path/to/output/directory"

required_packages <- c("parallel", "data.table", "dplyr", "readr", "ggplot2", "ggrepel", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

## -------------- Safe CSV reader --------------
safe_read_csv <- function(file_path) {
  tryCatch({
    if (exists("VERBOSE") && VERBOSE) cat("Reading:", basename(file_path), "\n")
    df <- readr::read_csv(file_path, show_col_types = FALSE, col_types = cols(.default = col_character()))
    if (!("Gene Name" %in% colnames(df))) {
      alt <- grep("gene", colnames(df), ignore.case = TRUE, value = TRUE)
      if (length(alt) >= 1) names(df)[which(colnames(df) == alt[1])] <- "Gene Name"
    }
    if (!("Variant" %in% colnames(df))) df$Variant <- paste0("var_", seq_len(nrow(df)))
    if (!("gnomAD_AF" %in% colnames(df))) df$gnomAD_AF <- NA
    if (!("PopFreqMax" %in% colnames(df))) df$PopFreqMax <- NA
    if (!("Gene Name" %in% colnames(df))) df$`Gene Name` <- NA

    df$`Gene Name` <- trimws(df$`Gene Name`)
    df <- df[!is.na(df$`Gene Name`) & df$`Gene Name` != "", , drop = FALSE]

    na_like <- function(x) x %in% c(NA, "NA", "N/A", "n/a", "")
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF[na_like(df$gnomAD_AF)] <- "0"
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax[na_like(df$PopFreqMax)] <- "0"
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF <- suppressWarnings(as.numeric(df$gnomAD_AF))
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax <- suppressWarnings(as.numeric(df$PopFreqMax))
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF[is.na(df$gnomAD_AF)] <- 0
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax[is.na(df$PopFreqMax)] <- 0

    if (exists("VERBOSE") && VERBOSE) cat("Successfully read", nrow(df), "variants from", basename(file_path), "\n")
    return(df)
  }, error = function(e) {
    if (exists("VERBOSE") && VERBOSE) cat("Error reading file:", basename(file_path), "-", e$message, "\n")
    return(NULL)
  })
}

## -------------- sample-level burden --------------
calculate_sample_burden <- function(sample_data, sample_name = "Unknown") {
  if (is.null(sample_data) || nrow(sample_data) == 0) {
    return(data.frame(Gene = character(0), HasBurden = logical(0)))
  }
  gene_burden <- sample_data %>%
    distinct(`Gene Name`) %>%
    mutate(HasBurden = TRUE) %>%
    rename(Gene = `Gene Name`)
  if (exists("VERBOSE") && VERBOSE) cat("Sample", sample_name, "has burden variants in", nrow(gene_burden), "genes\n")
  return(gene_burden)
}

## -------------- population burden estimate --------------
estimate_population_burden <- function(all_variants_data, sample_size = 62L) {
  if (exists("VERBOSE") && VERBOSE) cat("Estimating population burden for", length(unique(all_variants_data$`Gene Name`)), "genes...\n")
  all_variants_data <- all_variants_data %>%
    mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))

  gene_pop_freq <- all_variants_data %>%
    group_by(`Gene Name`) %>%
    summarise(variant_afs = list(variant_AF), variant_count = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
      log_non = sum(log1p(-unlist(carrier_probs))),
      individual_burden_freq = 1 - exp(log_non),
      expected_burden_count = sample_size * individual_burden_freq
    ) %>%
    ungroup() %>%
    select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)
  if (exists("VERBOSE") && VERBOSE) cat("Population burden estimated for", nrow(gene_pop_freq), "genes\n")
  return(gene_pop_freq)
}

## -------------- Bayesian per-gene test (fractional controls) --------------
perform_gene_bayesian_test <- function(gene_name,
                                       cases_burden, controls_burden_est,
                                       total_cases = 62L, total_controls = 62L,
                                       nsim = NSIM_GLOBAL,
                                       prior_alpha = PRIOR_ALPHA, prior_beta = PRIOR_BETA,
                                       prior_alpha_ctrl = PRIOR_ALPHA_CTRL, prior_beta_ctrl = PRIOR_BETA_CTRL,
                                       chunk_size = CHUNK_SIZE) {
  cases_burden <- as.integer(cases_burden)
  cases_no_burden <- as.integer(total_cases - cases_burden)

  controls_burden_est <- as.numeric(controls_burden_est)
  if (is.na(controls_burden_est)) controls_burden_est <- 0.0
  controls_no_burden_est <- as.numeric(total_controls - controls_burden_est)

  a_post <- prior_alpha + cases_burden
  b_post <- prior_beta  + cases_no_burden
  a_post_ctrl <- prior_alpha_ctrl + controls_burden_est
  b_post_ctrl <- prior_beta_ctrl  + controls_no_burden_est

  if (a_post <= 0 || b_post <= 0 || a_post_ctrl <= 0 || b_post_ctrl <= 0) {
    return(data.frame(
      Gene = gene_name, Cases_Burden = cases_burden, Cases_No_Burden = cases_no_burden,
      Controls_Burden = controls_burden_est, Controls_No_Burden = controls_no_burden_est,
      Posterior_Median_OR = NA_real_, OR_CI_Lower = NA_real_, OR_CI_Upper = NA_real_,
      Posterior_Prob_OR_gt_1 = NA_real_, Bayes_p = NA_real_, stringsAsFactors = FALSE
    ))
  }

  nsim <- as.integer(nsim)
  if (nsim <= 0) stop("nsim must be > 0")

  or_samples <- numeric(nsim)
  filled <- 0L
  while (filled < nsim) {
    this_draw <- min(as.integer(chunk_size), nsim - filled)
    p_case_samps <- rbeta(this_draw, shape1 = a_post, shape2 = b_post)
    p_ctrl_samps <- rbeta(this_draw, shape1 = a_post_ctrl, shape2 = b_post_ctrl)

    p_case_samps <- pmin(pmax(p_case_samps, CLIP_EPS), 1 - CLIP_EPS)
    p_ctrl_samps <- pmin(pmax(p_ctrl_samps, CLIP_EPS), 1 - CLIP_EPS)

    odds_case <- p_case_samps / (1 - p_case_samps)
    odds_ctrl <- p_ctrl_samps / (1 - p_ctrl_samps)
    or_chunk <- odds_case / odds_ctrl

    idx0 <- filled + seq_len(this_draw)
    or_samples[idx0] <- or_chunk
    filled <- filled + this_draw
  }

  posterior_prob_or_gt1 <- mean(or_samples > 1, na.rm = TRUE)
  bayes_p_two_sided <- 2 * min(posterior_prob_or_gt1, 1 - posterior_prob_or_gt1)
  median_or <- median(or_samples, na.rm = TRUE)
  or_ci <- as.numeric(quantile(or_samples, probs = c(0.025, 0.975), na.rm = TRUE))

  res_df <- data.frame(
    Gene = gene_name,
    Cases_Burden = as.integer(cases_burden),
    Cases_No_Burden = as.integer(cases_no_burden),
    Controls_Burden = controls_burden_est,
    Controls_No_Burden = controls_no_burden_est,
    Posterior_Median_OR = median_or,
    OR_CI_Lower = or_ci[1],
    OR_CI_Upper = or_ci[2],
    Posterior_Prob_OR_gt_1 = posterior_prob_or_gt1,
    Bayes_p = bayes_p_two_sided,
    stringsAsFactors = FALSE
  )

  rm(or_samples); gc()
  return(res_df)
}

## -------------- Parallel cluster setup --------------
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
if (exists("VERBOSE") && VERBOSE) cat("Setting up parallel processing with", num_cores, "cores...\n")
cl <- NULL
tryCatch({
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    clusterEvalQ(cl, {
      library(data.table); library(dplyr); library(readr)
      TRUE
    })
    # Export functions & globals (VERY IMPORTANT: include VERBOSE and SEED so workers can use them)
    clusterExport(cl, varlist = c(
      "safe_read_csv", "calculate_sample_burden", "estimate_population_burden", "perform_gene_bayesian_test",
      "NSIM_GLOBAL", "CHUNK_SIZE", "PRIOR_ALPHA", "PRIOR_BETA",
      "PRIOR_ALPHA_CTRL", "PRIOR_BETA_CTRL", "CLIP_EPS", "VERBOSE", "SEED"
    ), envir = environment())
    parallel::clusterSetRNGStream(cl, iseed = SEED)
    if (exists("VERBOSE") && VERBOSE) cat("Parallel cluster initialized and functions/globals exported.\n")
  } else {
    if (exists("VERBOSE") && VERBOSE) cat("Single-core mode (no cluster created)\n")
  }
}, error = function(e) {
  if (exists("VERBOSE") && VERBOSE) cat("Warning: Failed to create parallel cluster. Using sequential processing.\n")
  cl <<- NULL
})

## -------------- Main pipeline --------------
main_burden_analysis <- function() {
  if (exists("VERBOSE") && VERBOSE) cat("=== Starting Bayesian Gene-Level Burden Analysis ===\n")

  base_path <- "/path/to/results/directory"
  srr_names_file <- "/path/to/srr_names.txt"
  results_dir <- RESULTS_DIR

  if (!file.exists(srr_names_file)) stop("SRR names file not found: ", srr_names_file)
  srr_names <- readLines(srr_names_file, warn = FALSE)
  srr_names <- trimws(srr_names[srr_names != ""])
  if (exists("VERBOSE") && VERBOSE) cat("Found", length(srr_names), "SRR samples to process\n")

  file_paths <- file.path(base_path, srr_names, "ANNOTATION", "Franklin_Results_Without_Benign_LikelyBenign_0.001.csv")
  existing_files <- file.exists(file_paths)
  if (sum(!existing_files) > 0 && exists("VERBOSE") && VERBOSE) cat("Warning:", sum(!existing_files), "files missing\n")
  valid_paths <- file_paths[existing_files]; valid_srrs <- srr_names[existing_files]
  if (exists("VERBOSE") && VERBOSE) cat("Proceeding with", length(valid_paths), "existing files\n")

  if (!is.null(cl)) {
    if (exists("VERBOSE") && VERBOSE) cat("Reading files in parallel...\n")
    sample_data_list <- parLapply(cl, valid_paths, safe_read_csv)
  } else {
    if (exists("VERBOSE") && VERBOSE) cat("Reading files sequentially...\n")
    sample_data_list <- lapply(valid_paths, safe_read_csv)
  }
  names(sample_data_list) <- valid_srrs
  sample_data_list <- sample_data_list[!sapply(sample_data_list, is.null)]
  if (length(sample_data_list) == 0) stop("No valid sample data loaded.")

  sample_burdens <- list()
  for (i in seq_along(sample_data_list)) {
    sample_name <- names(sample_data_list)[i]
    sample_burdens[[sample_name]] <- calculate_sample_burden(sample_data_list[[i]], sample_name)
  }

  all_variants <- do.call(rbind, sample_data_list)
  if (exists("VERBOSE") && VERBOSE) cat("Combined data contains", nrow(all_variants), "variant observations\n")
  all_variants_unique <- all_variants %>% distinct(Variant, `Gene Name`, .keep_all = TRUE)
  if (exists("VERBOSE") && VERBOSE) cat("Unique variants:", nrow(all_variants_unique), "\n")

  all_genes <- unique(all_variants_unique$`Gene Name`)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  if (exists("VERBOSE") && VERBOSE) cat("Total genes to analyze:", length(all_genes), "\n")

  library(data.table)
  present_long <- rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]; if (nrow(sb) == 0) return(NULL); data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_burden_counts_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_burden_counts <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
    left_join(as.data.frame(gene_burden_counts_dt), by = "Gene") %>%
    mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))
  if (exists("VERBOSE") && VERBOSE) cat("Burden counts calculated for", nrow(gene_burden_counts), "genes\n")

  population_burden <- estimate_population_burden(all_variants_unique, sample_size = length(sample_burdens))

  gene_analysis_data <- gene_burden_counts %>%
    left_join(population_burden, by = "Gene") %>%
    mutate(expected_burden_count = ifelse(is.na(expected_burden_count), 0.0, expected_burden_count),
           Controls_Burden_Fisher = pmax(1L, round(expected_burden_count)),
           Controls_Burden_Bayesian = expected_burden_count)

  if (exists("VERBOSE") && VERBOSE) cat("Analysis data prepared for", nrow(gene_analysis_data), "genes\n")

  total_genes <- nrow(gene_analysis_data)
  if (exists("VERBOSE") && VERBOSE) cat("Will perform Bayesian sampling for", total_genes, "genes (nsim per gene =", NSIM_GLOBAL, ")\n")

  run_bayes_for_row <- function(i_row) {
    row <- gene_analysis_data[i_row, , drop = FALSE]
    perform_gene_bayesian_test(
      gene_name = row$Gene,
      cases_burden = as.integer(row$Cases_Burden),
      controls_burden_est = as.numeric(row$Controls_Burden_Bayesian),
      total_cases = length(sample_burdens),
      total_controls = length(sample_burdens),
      nsim = NSIM_GLOBAL,
      prior_alpha = PRIOR_ALPHA, prior_beta = PRIOR_BETA,
      prior_alpha_ctrl = PRIOR_ALPHA_CTRL, prior_beta_ctrl = PRIOR_BETA_CTRL,
      chunk_size = CHUNK_SIZE
    )
  }

  bayes_results_list <- vector("list", total_genes)
  if (!is.null(cl) && num_cores > 1) {
    if (exists("VERBOSE") && VERBOSE) cat("Running Bayesian sampling in parallel across genes...\n")
    block_size <- max(1, floor(total_genes / num_cores))
    indices <- seq_len(total_genes)
    idx_chunks <- split(indices, ceiling(seq_along(indices)/block_size))
    k <- 1
    for (chunk in idx_chunks) {
      if (exists("VERBOSE") && VERBOSE) cat("Processing gene chunk", k, "of", length(idx_chunks), " (", length(chunk), "genes )\n")
      res_chunk <- parLapply(cl, chunk, run_bayes_for_row)
      for (j in seq_along(chunk)) bayes_results_list[[ chunk[j] ]] <- res_chunk[[j]]
      k <- k + 1
    }
  } else {
    if (exists("VERBOSE") && VERBOSE) cat("Running Bayesian sampling sequentially across genes...\n")
    for (i in seq_len(total_genes)) {
      bayes_results_list[[i]] <- run_bayes_for_row(i)
      if (i %% 50 == 0 && exists("VERBOSE") && VERBOSE) cat("Completed Bayesian sampling for", i, "genes\n")
    }
  }

  bayes_results <- data.table::rbindlist(bayes_results_list, fill = TRUE)
  bayes_results <- as.data.frame(bayes_results)

  valid_results <- bayes_results[!is.na(bayes_results$Bayes_p), , drop = FALSE]
  if (nrow(valid_results) > 0) {
    valid_results$P_Value_Bonferroni <- p.adjust(valid_results$Bayes_p, method = "bonferroni")
    valid_results$P_Value_FDR <- p.adjust(valid_results$Bayes_p, method = "fdr")
    valid_results$Significant_Bonferroni <- valid_results$P_Value_Bonferroni < 0.05
    valid_results$Significant_FDR <- valid_results$P_Value_FDR < 0.05
    valid_results <- valid_results[order(valid_results$Bayes_p), ]
    if (exists("VERBOSE") && VERBOSE) cat("Multiple testing correction applied to", nrow(valid_results), "genes\n")
  } else stop("No valid Bayes_p values were calculated. Please check your data.")

  # save CSVs with requested names
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  readr::write_csv(valid_results, file.path(results_dir, "Gene_Burden_Analysis_Results.csv"))
  readr::write_csv(valid_results[valid_results$Significant_FDR == TRUE, , drop = FALSE], file.path(results_dir, "Significant_Genes_FDR_0.05.csv"))
  readr::write_csv(valid_results[valid_results$Bayes_p < 0.05, , drop = FALSE], file.path(results_dir, "Significant_Genes_Raw_P_0.05.csv"))
  if (exists("VERBOSE") && VERBOSE) cat("Saved CSV outputs to", results_dir, "\n")

  return(valid_results)
}

## -------------- Plotting (modified: log10 OR, safe p-floor, cap, limited labels) --------------
## -------------- Plotting (label exactly 6 genes) --------------
create_summary_plots <- function(results, top_n = TOP_LABEL_N) {
  if (is.null(results) || nrow(results) == 0) { 
    cat("[ERROR] 'results' is NULL or empty — nothing to plot.\n"); 
    return(invisible(NULL)) 
  }
  results_dir <- RESULTS_DIR
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # numeric safety
  results$Posterior_Median_OR <- suppressWarnings(as.numeric(results$Posterior_Median_OR))
  if (!"P_Value_FDR" %in% colnames(results)) {
    if ("Bayes_p" %in% colnames(results)) results$P_Value_FDR <- results$Bayes_p
    else stop("No Bayes_p or P_Value_FDR in results — cannot plot.")
  }
  results$P_Value_FDR <- suppressWarnings(as.numeric(results$P_Value_FDR))

  # safe p-floor and -log10(FDR)
  nonzero_min <- min(results$P_Value_FDR[results$P_Value_FDR > 0 & !is.na(results$P_Value_FDR)], na.rm = TRUE)
  if (is.infinite(nonzero_min) || is.na(nonzero_min)) nonzero_min <- MIN_NONZERO_P
  floor_p <- max(MIN_NONZERO_P, nonzero_min * 0.1)
  results$P_Value_FDR_safe <- ifelse(is.na(results$P_Value_FDR), NA_real_, ifelse(results$P_Value_FDR <= 0, floor_p, results$P_Value_FDR))
  results$neglog10FDR <- -log10(pmax(results$P_Value_FDR_safe, MIN_NONZERO_P))
  results$neglog10FDR <- pmin(results$neglog10FDR, NEGLOG10FDR_CAP)

  # compute log10 OR and cap
  results$log10OR <- ifelse(!is.na(results$Posterior_Median_OR) & results$Posterior_Median_OR > 0,
                            log10(results$Posterior_Median_OR), NA_real_)
  results$log10OR_cap <- pmin(pmax(results$log10OR, -LOG10OR_CAP), LOG10OR_CAP)

  # flags and cases
  results$SignificantFDR <- ifelse(!is.na(results$P_Value_FDR) & results$P_Value_FDR < 0.05, TRUE, FALSE)
  results$Cases_Burden <- ifelse(is.na(results$Cases_Burden), 0L, as.integer(results$Cases_Burden))

  # --- HERE: choose EXACT genes to label ---
  genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
  # filter only those present in our plot_data (and keep order by neglog10FDR so the most significant among them appear first)
  plot_data <- results %>% filter(!is.na(log10OR_cap) & !is.na(neglog10FDR))
  label_data <- NULL
  if (nrow(plot_data) > 0) {
    label_pool <- plot_data %>% filter(Gene %in% genes_to_label)
    if (nrow(label_pool) > 0) {
      # sort label_pool by significance and keep at most the genes in genes_to_label (preserves those exact genes)
      label_data <- label_pool %>% arrange(desc(neglog10FDR))

    }

    # prepare theme
    base_size <- 14
    theme_journal <- theme_bw(base_size = base_size) + theme(
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
    bar_fill <- "#9fd7f5"; bar_border <- "black"; sig_color <- "#e15759"; non_sig_color <- "gray60"

    p_volcano <- ggplot(plot_data, aes(x = log10OR_cap, y = neglog10FDR)) +
      geom_point(aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)), alpha = 0.85) +
      scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = non_sig_color), guide = FALSE) +
      scale_size_continuous(range = c(1.2, 5), guide = guide_legend(title = "Cases with burden")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      labs(title = "Gene Burden Analysis — Volcano Plot (Bayesian)",
           subtitle = paste0("Plotted genes: ", nrow(plot_data), "  |  FDR-significant: ", sum(plot_data$SignificantFDR, na.rm = TRUE)),
           x = paste0("log10(Posterior Median OR) [capped ±", LOG10OR_CAP, "]"),
           y = "-log10(P_Value_FDR)") +
      theme_journal + coord_cartesian(xlim = c(-LOG10OR_CAP, LOG10OR_CAP), expand = TRUE)

    if (!is.null(label_data) && nrow(label_data) > 0) {
      label_data_plot <- label_data %>% mutate(x = log10OR_cap, y = pmin(neglog10FDR, NEGLOG10FDR_CAP))
      p_volcano <- p_volcano + ggrepel::geom_text_repel(
        data = label_data_plot,
        aes(x = x, y = y, label = Gene),
        size = 3.4,
        max.overlaps = 20,
        box.padding = 0.6,
        point.padding = 0.3,
        segment.size = 0.4,
        segment.color = "gray40",
        arrow = grid::arrow(length = grid::unit(0.015, "npc"), type = "closed")
      )
    }

    ggsave(filename = file.path(results_dir, "Gene_Burden_Volcano_Plot_annotated.png"),
           plot = p_volcano, width = 11, height = 8, dpi = 300)
    message("Saved volcano:", file.path(results_dir, "Gene_Burden_Volcano_Plot_annotated.png"))
  } else {
    message("No valid points to plot (volcano)")
  }

  # (rest unchanged) top barplot and distribution follow same logic as before
  plot_genes <- results %>% filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) plot_genes <- results %>% filter(Bayes_p < 0.05 & !is.na(Bayes_p))
  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% arrange(P_Value_FDR) %>% head(50) %>% mutate(score = -log10(pmax(Bayes_p, MIN_NONZERO_P)))
    p_bar <- ggplot(plot_genes2, aes(x = reorder(Gene, score), y = score)) +
      geom_col(fill = bar_fill, color = bar_border) + coord_flip() +
      labs(title = "Top Genes by Significance (Bayesian)", x = "Gene", y = "-log10(Bayes_p)") +
      theme_journal + theme(axis.text.y = element_text(size = base_size - 1))
    height_needed <- max(6, 0.35 * nrow(plot_genes2) + 2)
    ggsave(filename = file.path(results_dir, "Top_Significant_Genes_Styled.png"), plot = p_bar, width = 12, height = height_needed, dpi = 300)
    message("Saved barplot:", file.path(results_dir, "Top_Significant_Genes_Styled.png"))
  } else {
    message("No significant genes to plot (barplot)")
  }

  if ("Cases_Burden" %in% colnames(results)) {
    p_dist <- ggplot(results, aes(x = Cases_Burden)) +
      geom_histogram(bins = 30, fill = bar_fill, color = bar_border) +
      labs(title = "Distribution of Gene Burden Counts (Cases)",
           subtitle = paste("Across", nrow(results), "genes"),
           x = "Number of cases with burden", y = "Number of genes") +
      theme_journal
    ggsave(filename = file.path(results_dir, "Burden_Distribution_Styled.png"), plot = p_dist, width = 10, height = 6, dpi = 300)
    message("Saved distribution:", file.path(results_dir, "Burden_Distribution_Styled.png"))
  } else {
    message("Cases_Burden not found -> skipping distribution")
  }

  return(invisible(NULL))
}

## -------------- MAIN EXECUTION --------------
results <- NULL
tryCatch({
  if (exists("VERBOSE") && VERBOSE) cat("Starting Bayesian analysis (nsim =", NSIM_GLOBAL, ") ...\n")
  results <- main_burden_analysis()
  if (exists("VERBOSE") && VERBOSE) cat("\nCreating summary plots...\n")
  create_summary_plots(results, top_n = TOP_LABEL_N)
  if (exists("VERBOSE") && VERBOSE) cat("\nANALYSIS COMPLETED SUCCESSFULLY!\n")
}, error = function(e) {
  cat("ERROR IN ANALYSIS:\n", e$message, "\n")
}, finally = {
  if (exists("cl") && !is.null(cl)) {
    tryCatch({ stopCluster(cl); if (exists("VERBOSE") && VERBOSE) cat("Parallel processing cluster stopped.\n") }, error = function(e) { cat("Warning: Error stopping cluster:", e$message, "\n") })
  }
})

cat("\n=== SESSION INFO ===\n"); print(sessionInfo()); cat("\nScript execution completed.\n")
