# Gene-Level Burden Analysis using Barnard's test (DescTools::BarnardTest)
# Label only top FDR-significant genes (by -log10 FDR), axis x fixed to [0,15]
# Outputs CSV: Gene_Burden_Analysis_Results.csv, Significant_Genes_FDR_0.05.csv, Significant_Genes_Raw_P_0.05.csv
# Outputs PNG: Burden_Distribution_Styled.png, Gene_Burden_Volcano_Plot_annotated.png, Top_Significant_Genes_Styled.png
# Author: Updated for user's pipeline (Aug 2025)

rm(list = ls())

# -----------------------
# CONFIG
# -----------------------
RESULTS_DIR <- "/path/to/output/directory"
BASE_PATH   <- "/path/to/results/directory"
SRR_FILE    <- "/path/to/srr_names.txt"
VERBOSE     <- TRUE

HALDANE_PC  <- 0.6    # Haldane continuity for OR reporting
LOG2OR_MIN  <- 0
LOG2OR_MAX  <- 15
TOP_LABEL_N <- 17     # not strictly used for the 6 forced labels, kept for API compatibility

# -----------------------
# REQUIRED PACKAGES
# -----------------------
required_packages <- c(
  "parallel","data.table","dplyr","readr",
  "ggplot2","ggrepel","scales","DescTools",
  "showtext","sysfonts"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# optional high-quality raster device
use_ragg <- requireNamespace("ragg", quietly = TRUE)
if (!use_ragg && VERBOSE) message("Note: 'ragg' not installed — TIFF will be produced via ggsave (showtext handles font rasterization). To improve TIFF rendering, install.packages('ragg').")

# -----------------------
# FONT REGISTRATION + showtext DPI
# -----------------------
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
  if (VERBOSE) message("No TTF found in common paths — using family = 'sans'.")
}
showtext::showtext_auto(TRUE)
showtext::showtext_opts(dpi = 1200)

# -----------------------
# UTIL: cm->inch + save_highres (PDF/EPS + TIFF@1200)
# -----------------------
cm_to_in <- function(x) x / 2.54

save_highres <- function(plot_obj, filename_base,
                         width_cm = 17.6, height_cm = 12.0,
                         save_vector = TRUE, tiff_res = 1200, family = font_family) {
  out_dir <- dirname(filename_base)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  w_in <- cm_to_in(width_cm)
  h_in <- cm_to_in(height_cm)

  # Vector: PDF, EPS (cairo)
  if (save_vector) {
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".pdf"), plot = plot_obj,
                      device = cairo_pdf, width = w_in, height = h_in,
                      units = "in", family = family)
    }, error = function(e) warning("PDF save failed: ", e$message))
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".eps"), plot = plot_obj,
                      device = cairo_ps, width = w_in, height = h_in,
                      units = "in", family = family)
    }, error = function(e) warning("EPS save failed: ", e$message))
  }

  # Raster TIFF at requested DPI
  out_tiff <- paste0(filename_base, ".tiff")
  if (use_ragg) {
    tryCatch({
      ragg::agg_tiff(filename = out_tiff, width = w_in, height = h_in, units = "in", res = tiff_res, background = "white")
      print(plot_obj)
      grDevices::dev.off()
    }, error = function(e) warning("ragg TIFF failed: ", e$message))
  } else {
    tryCatch({
      ggplot2::ggsave(out_tiff, plot = plot_obj,
                      device = "tiff", width = w_in, height = h_in,
                      units = "in", dpi = tiff_res, limitsize = FALSE,
                      compression = "lzw")
    }, error = function(e) warning("ggsave TIFF failed: ", e$message))
  }
  invisible(TRUE)
}

# -----------------------
# PARALLEL SETUP
# -----------------------
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
if (VERBOSE) message("Using ", num_cores, " cores")
cl <- NULL
if (num_cores > 1) {
  cl <- parallel::makeCluster(num_cores)
  parallel::clusterEvalQ(cl, { library(data.table); library(dplyr); library(readr); library(DescTools); TRUE })
  parallel::clusterExport(cl, varlist = c("VERBOSE"), envir = environment())
}

# -----------------------
# HELPERS: safe read + population burden
# -----------------------
safe_read_csv <- function(path) {
  tryCatch({
    if (VERBOSE) message("Reading: ", basename(path))
    df <- readr::read_csv(path, show_col_types = FALSE,
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

    na_like <- function(x) x %in% c(NA,"NA","N/A","n/a","")
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF[na_like(df$gnomAD_AF)] <- "0"
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax[na_like(df$PopFreqMax)] <- "0"
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF <- suppressWarnings(as.numeric(df$gnomAD_AF))
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax <- suppressWarnings(as.numeric(df$PopFreqMax))
    if ("gnomAD_AF" %in% colnames(df)) df$gnomAD_AF[is.na(df$gnomAD_AF)] <- 0
    if ("PopFreqMax" %in% colnames(df)) df$PopFreqMax[is.na(df$PopFreqMax)] <- 0

    if (VERBOSE) message(" -> read ", nrow(df), " variants")
    df
  }, error = function(e) {
    if (VERBOSE) message("Read error: ", e$message)
    NULL
  })
}

estimate_population_burden <- function(all_variants, sample_size = 62L) {
  if (VERBOSE) message("Estimating population burden...")
  all_variants <- all_variants %>%
    dplyr::mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    dplyr::mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))
  gp <- all_variants %>%
    dplyr::group_by(`Gene Name`) %>%
    dplyr::summarise(variant_afs = list(variant_AF), variant_count = dplyr::n(), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
                  log_non = sum(log1p(-unlist(carrier_probs))),
                  individual_burden_freq = 1 - exp(log_non),
                  expected_burden_count = sample_size * individual_burden_freq) %>%
    dplyr::ungroup() %>%
    dplyr::select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)
  if (VERBOSE) message(" -> estimated for ", nrow(gp), " genes")
  gp
}

# -----------------------
# Barnard wrapper (with Fisher fallback)
# -----------------------
perform_gene_barnard_test <- function(gene_name, cases_burden, controls_burden_frac,
                                      total_cases = 62L, total_controls = 62L,
                                      haldane_pc = HALDANE_PC) {
  cases_burden <- as.integer(cases_burden)
  cases_no_burden <- as.integer(total_cases - cases_burden)
  controls_frac <- as.numeric(controls_burden_frac)
  controls_int <- as.integer(round(controls_frac))
  controls_int <- pmax(0L, controls_int)
  controls_no <- as.integer(total_controls - controls_int)
  if (controls_no < 0L) controls_no <- 0L

  if ((cases_burden == 0 && controls_int == 0) || (cases_burden == total_cases && controls_int == total_controls)) {
    return(data.frame(Gene = gene_name,
                      Cases_Burden = cases_burden, Cases_No_Burden = cases_no_burden,
                      Controls_Burden = controls_int, Controls_No_Burden = controls_no,
                      Controls_Burden_Fraction = controls_frac,
                      P_Value = 1.0, Odds_Ratio = NA_real_, CI_Lower = NA_real_, CI_Upper = NA_real_,
                      stringsAsFactors = FALSE))
  }

  tab <- matrix(c(cases_burden, controls_int, cases_no_burden, controls_no), nrow = 2, byrow = TRUE)
  pval <- NA_real_; ok <- FALSE
  try({
    bt <- DescTools::BarnardTest(tab, method = "z-pooled")
    if (!is.null(bt$p.value)) { pval <- as.numeric(bt$p.value); ok <- TRUE }
  }, silent = TRUE)

  if (!ok) {
    # fallback to Fisher with Haldane continuity
    a <- cases_burden + haldane_pc; b <- controls_int + haldane_pc
    c_ <- cases_no_burden + haldane_pc; d <- controls_no + haldane_pc
    ft <- tryCatch(stats::fisher.test(matrix(c(a,b,c_,d), nrow = 2, byrow = TRUE)), error = function(e) NULL)
    pval <- if (!is.null(ft)) ft$p.value else NA_real_
  }

  # Haldane OR and fisher-based CI for reporting
  a_h <- cases_burden + haldane_pc; b_h <- controls_int + haldane_pc
  c_h <- cases_no_burden + haldane_pc; d_h <- controls_no + haldane_pc
  or_est <- (a_h / c_h) / (b_h / d_h)
  ci_l <- NA_real_; ci_u <- NA_real_
  ft2 <- tryCatch(stats::fisher.test(matrix(c(a_h,b_h,c_h,d_h), nrow=2, byrow=TRUE)), error = function(e) NULL)
  if (!is.null(ft2) && !is.null(ft2$conf.int)) { ci_l <- ft2$conf.int[1]; ci_u <- ft2$conf.int[2] }

  data.frame(Gene = gene_name,
             Cases_Burden = cases_burden, Cases_No_Burden = cases_no_burden,
             Controls_Burden = controls_int, Controls_No_Burden = controls_no,
             Controls_Burden_Fraction = controls_frac,
             P_Value = pval, Odds_Ratio = or_est, CI_Lower = ci_l, CI_Upper = ci_u,
             stringsAsFactors = FALSE)
}

# -----------------------
# MAIN pipeline
# -----------------------
main_burden_analysis <- function() {
  if (!file.exists(SRR_FILE)) stop("SRR names file not found: ", SRR_FILE)
  srr <- readLines(SRR_FILE, warn = FALSE) %>% trimws()
  srr <- srr[srr != ""]
  message("Found ", length(srr), " samples")

  paths <- file.path(BASE_PATH, srr, "ANNOTATION", "Franklin_Results_Without_Benign_LikelyBenign_0.001.csv")
  exist <- file.exists(paths)
  if (any(!exist)) warning(sum(!exist), " files missing")
  paths <- paths[exist]; srr <- srr[exist]
  if (length(paths) == 0) stop("No files to read")

  if (!is.null(cl)) parallel::clusterExport(cl, varlist = c("safe_read_csv","VERBOSE"), envir = environment())

  if (!is.null(cl)) {
    sample_list <- parallel::parLapply(cl, paths, safe_read_csv)
  } else sample_list <- lapply(paths, safe_read_csv)
  names(sample_list) <- srr
  sample_list <- sample_list[!sapply(sample_list, is.null)]
  if (length(sample_list) == 0) stop("No sample data loaded")

  # per-sample burden genes
  sample_burdens <- lapply(names(sample_list), function(nm) {
    df <- sample_list[[nm]]
    df %>% dplyr::distinct(`Gene Name`) %>% dplyr::rename(Gene = `Gene Name`)
  })
  names(sample_burdens) <- names(sample_list)

  # combine and unique
  combined <- do.call(rbind, sample_list)
  combined_unique <- combined %>% dplyr::distinct(Variant, `Gene Name`, .keep_all = TRUE)
  genes <- unique(combined_unique$`Gene Name`) %>% na.omit()
  genes <- genes[genes != ""]

  # cases per gene
  present_long <- data.table::rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]; if (nrow(sb) == 0) return(NULL); data.table::data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_cases_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_cases <- data.frame(Gene = genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(as.data.frame(gene_cases_dt), by = "Gene") %>%
    dplyr::mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))

  pop_burden <- estimate_population_burden(combined_unique, sample_size = length(sample_burdens))

  gene_df <- gene_cases %>%
    dplyr::left_join(pop_burden, by = "Gene") %>%
    dplyr::mutate(expected_burden_count = ifelse(is.na(expected_burden_count), 0.0, expected_burden_count),
                  Controls_Burden_Fraction = expected_burden_count)

  message("Prepared ", nrow(gene_df), " genes")

  if (!is.null(cl)) parallel::clusterExport(cl, varlist = c("perform_gene_barnard_test", "gene_df", "HALDANE_PC"), envir = environment())

  n <- nrow(gene_df); idx <- seq_len(n)
  if (!is.null(cl) && num_cores > 1) {
    res_list <- parallel::parLapply(cl, idx, function(i) {
      row <- gene_df[i, , drop = FALSE]
      perform_gene_barnard_test(row$Gene, as.integer(row$Cases_Burden), as.numeric(row$Controls_Burden_Fraction),
                                total_cases = length(sample_burdens), total_controls = length(sample_burdens))
    })
  } else {
    res_list <- vector("list", n)
    for (i in idx) {
      row <- gene_df[i, , drop = FALSE]
      res_list[[i]] <- perform_gene_barnard_test(row$Gene, as.integer(row$Cases_Burden), as.numeric(row$Controls_Burden_Fraction),
                                                 total_cases = length(sample_burdens), total_controls = length(sample_burdens))
      if (i %% 200 == 0 && VERBOSE) message("Processed ", i, " genes")
    }
  }

  results <- data.table::rbindlist(res_list, fill = TRUE) %>% as.data.frame()

  # multiple testing correction
  valid <- results[!is.na(results$P_Value), , drop = FALSE]
  if (nrow(valid) == 0) stop("No valid p-values")
  valid$P_Value_Bonferroni <- p.adjust(valid$P_Value, method = "bonferroni")
  valid$P_Value_FDR <- p.adjust(valid$P_Value, method = "fdr")
  valid$Significant_Bonferroni <- valid$P_Value_Bonferroni < 0.05
  valid$Significant_FDR <- valid$P_Value_FDR < 0.05
  valid <- valid[order(valid$P_Value), ]

  if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
  readr::write_csv(valid, file.path(RESULTS_DIR, "Gene_Burden_Analysis_Results.csv"))
  readr::write_csv(valid[valid$Significant_FDR == TRUE, , drop = FALSE], file.path(RESULTS_DIR, "Significant_Genes_FDR_0.05.csv"))
  readr::write_csv(valid[valid$P_Value < 0.05, , drop = FALSE], file.path(RESULTS_DIR, "Significant_Genes_Raw_P_0.05.csv"))
  message("Saved CSV outputs to ", RESULTS_DIR)

  valid
}

# -----------------------
# PLOTTING (no titles/captions in figures; axis labels allowed)
# -----------------------
create_summary_plots <- function(results, top_n = TOP_LABEL_N) {
  if (is.null(results) || nrow(results) == 0) { message("No results to plot"); return(invisible(NULL)) }

  if (!("Controls_Burden_Fraction" %in% colnames(results))) results$Controls_Burden_Fraction <- results$Controls_Burden
  if (!("Cases_No_Burden" %in% colnames(results))) results$Cases_No_Burden <- NA_real_

  per_row_total <- ifelse(!is.na(results$Cases_No_Burden), (results$Cases_Burden + results$Cases_No_Burden), NA_real_)
  fallback_total <- ifelse(all(is.na(per_row_total)), NA_real_, max(per_row_total, na.rm = TRUE))
  total_controls_row <- ifelse(!is.na(per_row_total), per_row_total, fallback_total)
  if (is.na(total_controls_row[1])) total_controls_row <- rep(max(results$Cases_Burden, na.rm = TRUE), nrow(results))
  results$Controls_No_Burden_Fraction <- total_controls_row - results$Controls_Burden_Fraction
  results$Controls_No_Burden_Fraction <- pmax(results$Controls_No_Burden_Fraction, 0)

  results <- results %>% dplyr::mutate(
    OR_plot = ((Cases_Burden + HALDANE_PC) / (pmax(Cases_No_Burden, 0) + HALDANE_PC)) /
              ((Controls_Burden_Fraction + HALDANE_PC) / (Controls_No_Burden_Fraction + HALDANE_PC)),
    OR_plot = ifelse(is.infinite(OR_plot) | is.nan(OR_plot), NA_real_, OR_plot),
    log2OR = ifelse(!is.na(OR_plot) & OR_plot > 0, log2(OR_plot), NA_real_)
  )

  results$log2OR_plot <- results$log2OR
  results$log2OR_plot[is.na(results$log2OR_plot)] <- NA_real_
  results$log2OR_plot <- pmin(pmax(results$log2OR_plot, LOG2OR_MIN), LOG2OR_MAX)

  if (!"P_Value_FDR" %in% colnames(results)) results$P_Value_FDR <- p.adjust(results$P_Value, method = "fdr")
  results$neglog10FDR <- -log10(pmax(results$P_Value_FDR, 1e-300))
  results$SignificantFDR <- ifelse(!is.na(results$P_Value_FDR) & results$P_Value_FDR < 0.05, TRUE, FALSE)
  results$Cases_Burden <- ifelse(is.na(results$Cases_Burden), 0L, as.integer(results$Cases_Burden))

  plot_data <- results %>% dplyr::filter(!is.na(log2OR_plot) & !is.na(neglog10FDR))
  message("Volcano will plot ", nrow(plot_data), " genes; FDR-significant: ", sum(plot_data$SignificantFDR, na.rm = TRUE))

  # label exactly these 6 genes if present
  genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
  label_data <- NULL
  if (nrow(plot_data) > 0) label_data <- plot_data %>% dplyr::filter(Gene %in% genes_to_label) %>% dplyr::arrange(desc(neglog10FDR))

  base_size <- 10
  theme_journal <- ggplot2::theme_bw(base_size = base_size) + ggplot2::theme(
    text = ggplot2::element_text(family = font_family),
    axis.title = ggplot2::element_text(size = base_size + 1),
    axis.text = ggplot2::element_text(size = base_size),
    legend.position = "right",
    panel.grid.major = ggplot2::element_line(color = "gray90"),
    panel.grid.minor = ggplot2::element_blank()
  )

  # Volcano (no title/subtitle inside figure)
  if (nrow(plot_data) > 0) {
    p_volcano <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log2OR_plot, y = neglog10FDR)) +
      ggplot2::geom_point(ggplot2::aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)), alpha = 0.9) +
      ggplot2::scale_color_manual(values = c("TRUE" = "#e15759", "FALSE" = "gray60"), guide = FALSE) +
      ggplot2::scale_size_continuous(range = c(1.2, 5), guide = ggplot2::guide_legend(title = "Cases with burden")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.4) +
      ggplot2::geom_vline(xintercept = c(1,2), linetype = "dashed", size = 0.4) +
      ggplot2::xlab("log2(Odds Ratio) (clipped)") + ggplot2::ylab("-log10(FDR)") +
      theme_journal +
      ggplot2::coord_cartesian(xlim = c(LOG2OR_MIN, LOG2OR_MAX), expand = TRUE)

    if (!is.null(label_data) && nrow(label_data) > 0) {
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = label_data,
          mapping = ggplot2::aes(x = log2OR_plot, y = neglog10FDR, label = Gene),
          size = 3.0,
          max.overlaps = 30,
          box.padding = 0.6,
          point.padding = 0.3,
          segment.size = 0.35,
          segment.color = "gray40",
          arrow = grid::arrow(length = grid::unit(0.02, "npc"), type = "closed")
        )
    }

    volcano_base <- file.path(RESULTS_DIR, "Fig1_Gene_Burden_Volcano")
    save_highres(p_volcano, volcano_base, width_cm = 17.6, height_cm = 12.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
    message("Saved volcano (PDF/EPS/TIFF@1200): ", volcano_base)
  } else {
    message("No volcano points to plot")
  }

  # Top barplot (no title/caption)
  plot_genes <- results %>% dplyr::filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) plot_genes <- results %>% dplyr::filter(P_Value < 0.05 & !is.na(P_Value))
  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% dplyr::arrange(P_Value) %>% head(50) %>% dplyr::mutate(score = -log10(pmax(P_Value, 1e-300)))
    p_bar <- ggplot2::ggplot(plot_genes2, ggplot2::aes(x = reorder(Gene, score), y = score)) +
      ggplot2::geom_col(fill = "#9fd7f5", color = "black", width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Gene") + ggplot2::ylab("-log10(P)") +
      theme_journal

    height_cm_bar <- max(6, 0.5 * nrow(plot_genes2) + 2)
    bar_base <- file.path(RESULTS_DIR, "Fig2_Top_Significant_Genes")
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

    dist_base <- file.path(RESULTS_DIR, "Fig3_Burden_Distribution")
    save_highres(p_dist, dist_base, width_cm = 8.8, height_cm = 6.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
    message("Saved distribution (PDF/EPS/TIFF@1200): ", dist_base)
  } else {
    message("Cases_Burden missing -> skip distribution")
  }

  invisible(TRUE)
}

# -----------------------
# RUN
# -----------------------
res <- NULL
tryCatch({
  res <- main_burden_analysis()
  create_summary_plots(res, top_n = TOP_LABEL_N)
  message("Analysis completed successfully.")
}, error = function(e) {
  message("ERROR: ", e$message)
}, finally = {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    message("Cluster stopped.")
  }
  message("Session info:"); print(sessionInfo())
})
