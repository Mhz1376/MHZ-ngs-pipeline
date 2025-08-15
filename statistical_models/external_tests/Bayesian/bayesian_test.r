# Gene-Level Burden Analysis — Bayesian (Jeffreys prior) with fractional control counts
# Modified: log10 OR plotting, safe p-floor, capped OR, limited labeling, CSV/PNG names adjusted
# Author: Updated for user's pipeline (Aug 2025)
rm(list = ls())

# ---------------------------
# CONFIG
# ---------------------------
SEED <- 12345
set.seed(SEED)

NSIM_GLOBAL <- 1e6      # <- note: very large; consider reducing for speed/memory
CHUNK_SIZE <- 5e4
PRIOR_ALPHA <- 0.5
PRIOR_BETA  <- 0.5
PRIOR_ALPHA_CTRL <- 0.5
PRIOR_BETA_CTRL  <- 0.5
CLIP_EPS <- 1e-12
VERBOSE <- TRUE

# plotting / labeling params
TOP_LABEL_N <- 22
LOG10OR_CAP <- 6
NEGLOG10FDR_CAP <- 300
MIN_NONZERO_P <- 1e-300

RESULTS_DIR <- "/path/to/output/directory"

# ---------------------------
# Required packages
# ---------------------------
required_packages <- c(
  "parallel","data.table","dplyr","readr",
  "ggplot2","ggrepel","scales","showtext","sysfonts"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# optional high-quality raster backend
use_ragg <- requireNamespace("ragg", quietly = TRUE)
if (!use_ragg && VERBOSE) message("Note: 'ragg' not installed — TIFF will be produced via ggsave; install.packages('ragg') to improve TIFF rendering.")

# ---------------------------
# Font registration + showtext DPI
# ---------------------------
candidate_paths <- c(
  "C:/Windows/Fonts/arial.ttf",
  "C:/Windows/Fonts/Arial.ttf",
  "/Library/Fonts/Arial.ttf",
  "/System/Library/Fonts/Supplemental/Arial.ttf",
  "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
  "/usr/share/fonts/truetype/freefont/FreeSans.ttf"
)
font_path_found <- ""
for (p in candidate_paths) if (file.exists(p)) { font_path_found <- p; break }
if (nzchar(font_path_found)) {
  font_family <- "CustomSans"
  sysfonts::font_add(family = font_family, regular = font_path_found)
  if (VERBOSE) message("Registered font: ", font_path_found)
} else {
  font_family <- "sans"
  if (VERBOSE) message("No TTF found — using family = 'sans'")
}
showtext::showtext_auto(TRUE)
showtext::showtext_opts(dpi = 1200)

# ---------------------------
# Utility: cm -> in + save_highres (PDF/EPS + TIFF@1200)
# ---------------------------
cm_to_in <- function(x) x / 2.54

save_highres <- function(plot_obj, filename_base,
                         width_cm = 17.6, height_cm = 12.0,
                         save_vector = TRUE, tiff_res = 1200, family = font_family) {
  out_dir <- dirname(filename_base)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  w_in <- cm_to_in(width_cm)
  h_in <- cm_to_in(height_cm)

  # vector
  if (save_vector) {
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".pdf"), plot = plot_obj, device = cairo_pdf,
                      width = w_in, height = h_in, units = "in", family = family)
      if (VERBOSE) message("Saved vector: ", paste0(filename_base, ".pdf"))
    }, error = function(e) warning("PDF save failed: ", e$message))
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".eps"), plot = plot_obj, device = cairo_ps,
                      width = w_in, height = h_in, units = "in", family = family)
      if (VERBOSE) message("Saved vector: ", paste0(filename_base, ".eps"))
    }, error = function(e) warning("EPS save failed: ", e$message))
  }

  # raster TIFF @ specified DPI
  out_tiff <- paste0(filename_base, ".tiff")
  if (use_ragg) {
    tryCatch({
      ragg::agg_tiff(filename = out_tiff, width = w_in, height = h_in, units = "in", res = tiff_res, background = "white")
      print(plot_obj)
      grDevices::dev.off()
      if (VERBOSE) message("Saved raster (ragg): ", out_tiff, " (DPI=", tiff_res, ")")
    }, error = function(e) warning("ragg TIFF failed: ", e$message))
  } else {
    tryCatch({
      ggplot2::ggsave(out_tiff, plot = plot_obj, device = "tiff",
                      width = w_in, height = h_in, units = "in", dpi = tiff_res,
                      limitsize = FALSE, compression = "lzw")
      if (VERBOSE) message("Saved raster (ggsave): ", out_tiff, " (DPI=", tiff_res, ")")
    }, error = function(e) warning("TIFF save failed: ", e$message))
  }

  invisible(TRUE)
}

# ---------------------------
# Safe CSV reader
# ---------------------------
safe_read_csv <- function(file_path) {
  tryCatch({
    if (VERBOSE) cat("Reading:", basename(file_path), "\n")
    df <- readr::read_csv(file_path, show_col_types = FALSE,
                         col_types = readr::cols(.default = readr::col_character()))
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

    if (VERBOSE) cat(" -> read", nrow(df), "variants from", basename(file_path), "\n")
    return(df)
  }, error = function(e) {
    if (VERBOSE) cat("Error reading file:", basename(file_path), "-", e$message, "\n")
    return(NULL)
  })
}

# ---------------------------
# Sample burden, population burden, Bayesian test
# ---------------------------
calculate_sample_burden <- function(sample_data, sample_name = "Unknown") {
  if (is.null(sample_data) || nrow(sample_data) == 0) {
    return(data.frame(Gene = character(0), HasBurden = logical(0), stringsAsFactors = FALSE))
  }
  gene_burden <- sample_data %>%
    dplyr::distinct(`Gene Name`) %>%
    dplyr::mutate(HasBurden = TRUE) %>%
    dplyr::rename(Gene = `Gene Name`)
  if (VERBOSE) cat("Sample", sample_name, "has burden variants in", nrow(gene_burden), "genes\n")
  return(gene_burden)
}

estimate_population_burden <- function(all_variants_data, sample_size = 62L) {
  if (VERBOSE) cat("Estimating population burden for", length(unique(all_variants_data$`Gene Name`)), "genes...\n")
  all_variants_data <- all_variants_data %>%
    dplyr::mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    dplyr::mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))

  gene_pop_freq <- all_variants_data %>%
    dplyr::group_by(`Gene Name`) %>%
    dplyr::summarise(variant_afs = list(variant_AF), variant_count = dplyr::n(), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
      log_non = sum(log1p(-unlist(carrier_probs))),
      individual_burden_freq = 1 - exp(log_non),
      expected_burden_count = sample_size * individual_burden_freq
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)
  if (VERBOSE) cat("Population burden estimated for", nrow(gene_pop_freq), "genes\n")
  return(gene_pop_freq)
}

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

# ---------------------------
# Parallel cluster setup
# ---------------------------
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
if (VERBOSE) cat("Setting up parallel processing with", num_cores, "cores...\n")
cl <- NULL
tryCatch({
  if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterEvalQ(cl, { library(data.table); library(dplyr); library(readr); TRUE })
    parallel::clusterExport(cl, varlist = c(
      "safe_read_csv","calculate_sample_burden","estimate_population_burden","perform_gene_bayesian_test",
      "NSIM_GLOBAL","CHUNK_SIZE","PRIOR_ALPHA","PRIOR_BETA","PRIOR_ALPHA_CTRL","PRIOR_BETA_CTRL",
      "CLIP_EPS","VERBOSE","SEED"
    ), envir = environment())
    parallel::clusterSetRNGStream(cl, iseed = SEED)
    if (VERBOSE) cat("Cluster initialized and functions/globals exported.\n")
  } else {
    if (VERBOSE) cat("Single-core mode (no cluster created)\n")
  }
}, error = function(e) {
  if (VERBOSE) cat("Warning: Failed to create parallel cluster — using sequential processing.\n")
  cl <<- NULL
})

# ---------------------------
# Main Bayesian pipeline
# ---------------------------
main_burden_analysis <- function() {
  if (VERBOSE) cat("=== Starting Bayesian Gene-Level Burden Analysis ===\n")

  srr_file <- "/path/to/srr_names.txt"
  base_path <- "/path/to/results/directory"
  results_dir <- RESULTS_DIR

  if (!file.exists(srr_file)) stop("SRR names file not found: ", srr_file)
  srr_names <- readLines(srr_file, warn = FALSE) %>% trimws()
  srr_names <- srr_names[srr_names != ""]
  if (VERBOSE) cat("Found", length(srr_names), "SRR samples\n")

  file_paths <- file.path(base_path, srr_names, "ANNOTATION", "Franklin_Results_Without_Benign_LikelyBenign_0.001.csv")
  existing_files <- file.exists(file_paths)
  if (sum(!existing_files) > 0 && VERBOSE) cat("Warning:", sum(!existing_files), "files missing\n")
  valid_paths <- file_paths[existing_files]; valid_srrs <- srr_names[existing_files]
  if (VERBOSE) cat("Proceeding with", length(valid_paths), "files\n")

  if (!is.null(cl)) {
    sample_data_list <- parallel::parLapply(cl, valid_paths, safe_read_csv)
  } else {
    sample_data_list <- lapply(valid_paths, safe_read_csv)
  }
  names(sample_data_list) <- valid_srrs
  sample_data_list <- sample_data_list[!sapply(sample_data_list, is.null)]
  if (length(sample_data_list) == 0) stop("No sample data loaded")

  sample_burdens <- list()
  for (i in seq_along(sample_data_list)) {
    name_i <- names(sample_data_list)[i]
    sample_burdens[[name_i]] <- calculate_sample_burden(sample_data_list[[i]], sample_name = name_i)
  }

  all_variants <- do.call(rbind, sample_data_list)
  if (VERBOSE) cat("Combined rows:", nrow(all_variants), "\n")
  all_variants_unique <- all_variants %>% dplyr::distinct(Variant, `Gene Name`, .keep_all = TRUE)

  all_genes <- unique(all_variants_unique$`Gene Name`)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  if (VERBOSE) cat("Genes to analyze:", length(all_genes), "\n")

  present_long <- data.table::rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]; if (nrow(sb) == 0) return(NULL); data.table::data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_cases_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_cases <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(as.data.frame(gene_cases_dt), by = "Gene") %>%
    dplyr::mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))

  population_burden <- estimate_population_burden(all_variants_unique, sample_size = length(sample_burdens))

  gene_analysis_data <- gene_cases %>%
    dplyr::left_join(population_burden, by = "Gene") %>%
    dplyr::mutate(expected_burden_count = ifelse(is.na(expected_burden_count), 0.0, expected_burden_count),
                  Controls_Burden_Fisher = pmax(1L, round(expected_burden_count)),
                  Controls_Burden_Bayesian = expected_burden_count)

  if (VERBOSE) cat("Prepared analysis data for", nrow(gene_analysis_data), "genes\n")

  total_genes <- nrow(gene_analysis_data)
  if (VERBOSE) cat("Will run Bayesian sampling for", total_genes, "genes (nsim=", NSIM_GLOBAL, ")\n")

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
    # split into reasonable chunks to avoid overloading the cluster
    idx <- seq_len(total_genes)
    block_size <- max(1, floor(total_genes / num_cores))
    idx_chunks <- split(idx, ceiling(seq_along(idx)/block_size))
    k <- 1
    for (chunk in idx_chunks) {
      if (VERBOSE) cat("Processing chunk", k, "of", length(idx_chunks), " (", length(chunk), "genes)\n")
      res_chunk <- parallel::parLapply(cl, chunk, run_bayes_for_row)
      for (j in seq_along(chunk)) bayes_results_list[[ chunk[j] ]] <- res_chunk[[j]]
      k <- k + 1
    }
  } else {
    for (i in seq_len(total_genes)) {
      bayes_results_list[[i]] <- run_bayes_for_row(i)
      if (i %% 50 == 0 && VERBOSE) cat("Completed", i, "genes\n")
    }
  }

  bayes_results <- data.table::rbindlist(bayes_results_list, fill = TRUE) %>% as.data.frame()

  valid_results <- bayes_results[!is.na(bayes_results$Bayes_p), , drop = FALSE]
  if (nrow(valid_results) > 0) {
    valid_results$P_Value_Bonferroni <- p.adjust(valid_results$Bayes_p, method = "bonferroni")
    valid_results$P_Value_FDR <- p.adjust(valid_results$Bayes_p, method = "fdr")
    valid_results$Significant_Bonferroni <- valid_results$P_Value_Bonferroni < 0.05
    valid_results$Significant_FDR <- valid_results$P_Value_FDR < 0.05
    valid_results <- valid_results[order(valid_results$Bayes_p), ]
    if (VERBOSE) cat("Applied multiple testing correction for", nrow(valid_results), "genes\n")
  } else stop("No valid Bayes_p values; check data or sampling parameters.")

  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  readr::write_csv(valid_results, file.path(results_dir, "Gene_Burden_Analysis_Results.csv"))
  readr::write_csv(valid_results[valid_results$Significant_FDR == TRUE, , drop = FALSE], file.path(results_dir, "Significant_Genes_FDR_0.05.csv"))
  readr::write_csv(valid_results[valid_results$Bayes_p < 0.05, , drop = FALSE], file.path(results_dir, "Significant_Genes_Raw_P_0.05.csv"))
  if (VERBOSE) cat("Saved CSVs to", results_dir, "\n")

  return(valid_results)
}

# ---------------------------
# Plotting (log10 OR, safe p-floor, cap, only 6 forced labels)
# ---------------------------
create_summary_plots <- function(results, top_n = TOP_LABEL_N) {
  if (is.null(results) || nrow(results) == 0) { message("No results to plot"); return(invisible(NULL)) }
  results_dir <- RESULTS_DIR
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  results$Posterior_Median_OR <- suppressWarnings(as.numeric(results$Posterior_Median_OR))
  if (!"P_Value_FDR" %in% colnames(results)) {
    if ("Bayes_p" %in% colnames(results)) results$P_Value_FDR <- results$Bayes_p
    else stop("No Bayes_p/P_Value_FDR found")
  }
  results$P_Value_FDR <- suppressWarnings(as.numeric(results$P_Value_FDR))

  # safe p-floor
  nonzero_min <- tryCatch(min(results$P_Value_FDR[results$P_Value_FDR > 0 & !is.na(results$P_Value_FDR)], na.rm = TRUE), error = function(e) NA_real_)
  if (is.infinite(nonzero_min) || is.na(nonzero_min)) nonzero_min <- MIN_NONZERO_P
  floor_p <- max(MIN_NONZERO_P, nonzero_min * 0.1)
  results$P_Value_FDR_safe <- ifelse(is.na(results$P_Value_FDR), NA_real_, ifelse(results$P_Value_FDR <= 0, floor_p, results$P_Value_FDR))
  results$neglog10FDR <- -log10(pmax(results$P_Value_FDR_safe, MIN_NONZERO_P))
  results$neglog10FDR <- pmin(results$neglog10FDR, NEGLOG10FDR_CAP)

  # compute + cap log10 OR
  results$log10OR <- ifelse(!is.na(results$Posterior_Median_OR) & results$Posterior_Median_OR > 0,
                            log10(results$Posterior_Median_OR), NA_real_)
  results$log10OR_cap <- pmin(pmax(results$log10OR, -LOG10OR_CAP), LOG10OR_CAP)

  results$SignificantFDR <- ifelse(!is.na(results$P_Value_FDR) & results$P_Value_FDR < 0.05, TRUE, FALSE)
  results$Cases_Burden <- ifelse(is.na(results$Cases_Burden), 0L, as.integer(results$Cases_Burden))

  plot_data <- results %>% dplyr::filter(!is.na(log10OR_cap) & !is.na(neglog10FDR))
  message("Volcano points:", nrow(plot_data), " FDR-significant:", sum(plot_data$SignificantFDR, na.rm = TRUE))

  # EXACT genes to label:
  genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
  label_data <- NULL
  if (nrow(plot_data) > 0) {
    label_data <- plot_data %>% dplyr::filter(Gene %in% genes_to_label) %>% dplyr::arrange(desc(neglog10FDR))
  }

  base_size <- 10
  theme_journal <- ggplot2::theme_bw(base_size = base_size) + ggplot2::theme(
    text = ggplot2::element_text(family = font_family),
    axis.title = ggplot2::element_text(size = base_size + 1),
    axis.text = ggplot2::element_text(size = base_size),
    legend.position = "right",
    panel.grid.major = ggplot2::element_line(color = "gray90"),
    panel.grid.minor = ggplot2::element_blank()
  )

  # Volcano: NO title or subtitle (journal requires captions in manuscript)
  if (nrow(plot_data) > 0) {
    p_volcano <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log10OR_cap, y = neglog10FDR)) +
      ggplot2::geom_point(ggplot2::aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)), alpha = 0.9) +
      ggplot2::scale_color_manual(values = c("TRUE" = "#e15759", "FALSE" = "gray60"), guide = FALSE) +
      ggplot2::scale_size_continuous(range = c(1.2, 5), guide = ggplot2::guide_legend(title = "Cases with burden")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.4) +
      ggplot2::xlab(paste0("log10(Posterior median OR) (capped ±", LOG10OR_CAP, ")")) +
      ggplot2::ylab("-log10(FDR)") +
      theme_journal +
      ggplot2::coord_cartesian(xlim = c(-LOG10OR_CAP, LOG10OR_CAP), expand = TRUE)

    if (!is.null(label_data) && nrow(label_data) > 0) {
      label_plot <- label_data %>% dplyr::mutate(x = log10OR_cap, y = pmin(neglog10FDR, NEGLOG10FDR_CAP))
      p_volcano <- p_volcano + ggrepel::geom_text_repel(
        data = label_plot,
        mapping = ggplot2::aes(x = x, y = y, label = Gene),
        size = 3.0,
        box.padding = 0.6,
        point.padding = 0.3,
        segment.size = 0.35,
        segment.color = "gray40",
        arrow = grid::arrow(length = grid::unit(0.02, "npc"), type = "closed"),
        max.overlaps = 30
      )
    }

    volcano_base <- file.path(results_dir, "Fig1_Gene_Burden_Volcano")
    save_highres(p_volcano, volcano_base, width_cm = 17.6, height_cm = 12.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
    message("Saved volcano (PDF/EPS/TIFF@1200): ", volcano_base)
  } else {
    message("No valid points to plot (volcano)")
  }

  # Top genes barplot (no title inside figure)
  plot_genes <- results %>% dplyr::filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) plot_genes <- results %>% dplyr::filter(Bayes_p < 0.05 & !is.na(Bayes_p))
  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% dplyr::arrange(P_Value_FDR) %>% head(50) %>% dplyr::mutate(score = -log10(pmax(Bayes_p, MIN_NONZERO_P)))
    p_bar <- ggplot2::ggplot(plot_genes2, ggplot2::aes(x = reorder(Gene, score), y = score)) +
      ggplot2::geom_col(fill = "#9fd7f5", color = "black", width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Gene") + ggplot2::ylab("-log10(Bayes_p)") +
      theme_journal
    bar_base <- file.path(results_dir, "Fig2_Top_Significant_Genes")
    height_cm_bar <- max(6, 0.6 * nrow(plot_genes2) + 2)
    save_highres(p_bar, bar_base, width_cm = 8.8, height_cm = height_cm_bar, save_vector = TRUE, tiff_res = 1200, family = font_family)
    message("Saved barplot (PDF/EPS/TIFF@1200): ", bar_base)
  } else {
    message("No significant genes to plot (barplot)")
  }

  # Distribution (no title/caption)
  if ("Cases_Burden" %in% colnames(results)) {
    p_dist <- ggplot2::ggplot(results, ggplot2::aes(x = Cases_Burden)) +
      ggplot2::geom_histogram(bins = 30, fill = "#9fd7f5", color = "black") +
      ggplot2::xlab("Number of cases with burden") + ggplot2::ylab("Number of genes") +
      theme_journal
    dist_base <- file.path(results_dir, "Fig3_Burden_Distribution")
    save_highres(p_dist, dist_base, width_cm = 8.8, height_cm = 6.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
    message("Saved distribution (PDF/EPS/TIFF@1200): ", dist_base)
  } else {
    message("Cases_Burden missing -> skip distribution")
  }

  invisible(NULL)
}

# ---------------------------
# MAIN EXECUTION
# ---------------------------
results <- NULL
tryCatch({
  if (VERBOSE) cat("Starting Bayesian analysis (nsim =", NSIM_GLOBAL, ") ...\n")
  results <- main_burden_analysis()
  if (VERBOSE) cat("Creating publication-ready figures (no internal titles/captions)...\n")
  create_summary_plots(results, top_n = TOP_LABEL_N)
  if (VERBOSE) cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
}, error = function(e) {
  cat("ERROR IN ANALYSIS:\n", e$message, "\n")
}, finally = {
  if (exists("cl") && !is.null(cl)) {
    tryCatch({ parallel::stopCluster(cl); if (VERBOSE) cat("Cluster stopped.\n") }, error = function(e) cat("Warning stopping cluster:", e$message, "\n"))
  }
  cat("\n=== SESSION INFO ===\n"); print(sessionInfo()); cat("\nScript finished.\n")
})
