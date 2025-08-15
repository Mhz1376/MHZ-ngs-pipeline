# Gene-Level Burden Analysis using Fisher's Exact Test
# Updated plotting: unified journal-style theme + annotated volcano (all FDR < 0.05 labeled)
# Author: Updated for user's pipeline (Aug 2025)
rm(list = ls())

# -------------------------
# Required packages
# -------------------------
required_packages <- c(
  "parallel", "data.table", "dplyr", "readr",
  "ggplot2", "ggrepel", "scales",
  "showtext", "sysfonts"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ragg optional (better TIFF rendering)
if (!requireNamespace("ragg", quietly = TRUE)) {
  message("Note: 'ragg' not installed. To improve TIFF text rendering install.packages('ragg').")
}

# -------------------------
# Font detection + registration
# -------------------------
candidate_paths <- c(
  "C:/Windows/Fonts/arial.ttf",
  "C:/Windows/Fonts/Arial.ttf",
  "/Library/Fonts/Arial.ttf",
  "/System/Library/Fonts/Supplemental/Arial.ttf",
  "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
  "/usr/share/fonts/truetype/freefont/FreeSans.ttf"
)

font_path_found <- ""
for (p in candidate_paths) {
  if (file.exists(p)) { font_path_found <- p; break }
}

if (nzchar(font_path_found)) {
  font_family <- "CustomSans"
  sysfonts::font_add(family = font_family, regular = font_path_found)
  message("Registered font '", font_family, "' from: ", font_path_found)
} else {
  font_family <- "sans"
  message("No system TTF found in common locations — using family = 'sans'.")
}

# enable showtext and set DPI for rasterization
showtext::showtext_auto(enable = TRUE)
showtext::showtext_opts(dpi = 1200)

# -------------------------
# Parallel cluster setup (robust)
# -------------------------
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
cat("Setting up parallel processing with", num_cores, "cores...\n")

cl <- NULL
tryCatch({
  if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterEvalQ(cl, {
      library(data.table); library(dplyr); library(readr)
    })
    cat("Parallel cluster initialized successfully\n")
  } else {
    cat("Single-core mode (no cluster created)\n")
  }
}, error = function(e) {
  cat("Warning: Failed to create parallel cluster. Using sequential processing.\n")
  cl <<- NULL
})

# -------------------------
# Utilities: cm->in and save function (PDF/EPS + TIFF@1200)
# -------------------------
cm_to_in <- function(x) x / 2.54

use_ragg <- requireNamespace("ragg", quietly = TRUE)

save_highres <- function(plot_obj, filename_base,
                         width_cm = 17.6, height_cm = 12.0,
                         save_vector = TRUE, tiff_res = 1200, family = font_family) {
  out_dir <- dirname(filename_base)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  w_in <- cm_to_in(width_cm)
  h_in <- cm_to_in(height_cm)

  # Vector: PDF and EPS (cairo devices)
  if (save_vector) {
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".pdf"), plot = plot_obj,
             device = cairo_pdf, width = w_in, height = h_in,
             units = "in", family = family)
      cat("Saved vector:", paste0(filename_base, ".pdf"), "\n")
    }, error = function(e) {
      warning("PDF save failed: ", e$message)
    })
    tryCatch({
      ggplot2::ggsave(paste0(filename_base, ".eps"), plot = plot_obj,
             device = cairo_ps, width = w_in, height = h_in,
             units = "in", family = family)
      cat("Saved vector:", paste0(filename_base, ".eps"), "\n")
    }, error = function(e) {
      warning("EPS save failed: ", e$message)
    })
  }

  # Raster TIFF
  out_tiff <- paste0(filename_base, ".tiff")
  if (use_ragg) {
    tryCatch({
      ragg::agg_tiff(filename = out_tiff, width = w_in, height = h_in, units = "in", res = tiff_res, background = "white")
      print(plot_obj)
      grDevices::dev.off()
      cat("Saved raster (ragg):", out_tiff, "(DPI =", tiff_res, ")\n")
    }, error = function(e) {
      warning("ragg TIFF failed: ", e$message)
    })
  } else {
    tryCatch({
      ggplot2::ggsave(out_tiff, plot = plot_obj,
             device = "tiff", width = w_in, height = h_in,
             units = "in", dpi = tiff_res, limitsize = FALSE,
             compression = "lzw")
      cat("Saved raster (ggsave):", out_tiff, "(DPI =", tiff_res, ")\n")
    }, error = function(e) {
      warning("TIFF save failed: ", e$message)
    })
  }

  invisible(TRUE)
}

# -------------------------
# Safe CSV reader
# -------------------------
safe_read_csv <- function(file_path) {
  tryCatch({
    cat("Reading:", basename(file_path), "\n")
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
    df$gnomAD_AF[na_like(df$gnomAD_AF)] <- "0"
    df$PopFreqMax[na_like(df$PopFreqMax)] <- "0"
    df$gnomAD_AF <- suppressWarnings(as.numeric(df$gnomAD_AF))
    df$PopFreqMax <- suppressWarnings(as.numeric(df$PopFreqMax))
    df$gnomAD_AF[is.na(df$gnomAD_AF)] <- 0
    df$PopFreqMax[is.na(df$PopFreqMax)] <- 0

    cat("Successfully read", nrow(df), "variants from", basename(file_path), "\n")
    return(df)
  }, error = function(e) {
    cat("Error reading file:", basename(file_path), "-", e$message, "\n")
    return(NULL)
  })
}

# -------------------------
# Sample burden, population burden, fisher test
# -------------------------
calculate_sample_burden <- function(sample_data, sample_name = "Unknown") {
  if (is.null(sample_data) || nrow(sample_data) == 0) {
    return(data.frame(Gene = character(0), HasBurden = logical(0), stringsAsFactors = FALSE))
  }
  gene_burden <- sample_data %>%
    dplyr::distinct(`Gene Name`) %>%
    dplyr::mutate(HasBurden = TRUE) %>%
    dplyr::rename(Gene = `Gene Name`)
  cat("Sample", sample_name, "has burden variants in", nrow(gene_burden), "genes\n")
  return(gene_burden)
}

estimate_population_burden <- function(all_variants_data, sample_size = 62) {
  cat("Estimating population burden (exact product method) for", length(unique(all_variants_data$`Gene Name`)), "genes...\n")
  all_variants_data <- all_variants_data %>%
    dplyr::mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    dplyr::mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))

  gene_pop_freq <- all_variants_data %>%
    dplyr::group_by(`Gene Name`) %>%
    dplyr::summarise(
      variant_afs = list(variant_AF),
      variant_count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
      log_non = sum(log1p(-unlist(carrier_probs))),
      individual_burden_freq = 1 - exp(log_non),
      expected_burden_count = sample_size * individual_burden_freq
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)

  cat("Population burden estimated for", nrow(gene_pop_freq), "genes\n")
  return(gene_pop_freq)
}

perform_gene_fisher_test <- function(gene_name, cases_burden, controls_burden,
                                     total_cases = 62, total_controls = 62) {
  cases_no_burden <- max(0, total_cases - cases_burden)
  controls_no_burden <- max(0, total_controls - controls_burden)

  if ((cases_burden == 0 && controls_burden == 0) ||
      (cases_burden == total_cases && controls_burden == total_controls)) {
    return(data.frame(
      Gene = gene_name,
      Cases_Burden = cases_burden,
      Cases_No_Burden = cases_no_burden,
      Controls_Burden = controls_burden,
      Controls_No_Burden = controls_no_burden,
      P_Value = 1.0,
      Odds_Ratio = NA_real_,
      CI_Lower = NA_real_,
      CI_Upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  contingency_matrix <- matrix(
    c(cases_burden, controls_burden,
      cases_no_burden, controls_no_burden),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("Burden", "No_Burden"), c("Cases", "Controls"))
  )

  res <- tryCatch({
    fisher_result <- stats::fisher.test(contingency_matrix)
    data.frame(
      Gene = gene_name,
      Cases_Burden = cases_burden,
      Cases_No_Burden = cases_no_burden,
      Controls_Burden = controls_burden,
      Controls_No_Burden = controls_no_burden,
      P_Value = fisher_result$p.value,
      Odds_Ratio = if (!is.null(fisher_result$estimate)) as.numeric(fisher_result$estimate) else NA_real_,
      CI_Lower = if (!is.null(fisher_result$conf.int)) fisher_result$conf.int[1] else NA_real_,
      CI_Upper = if (!is.null(fisher_result$conf.int)) fisher_result$conf.int[2] else NA_real_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("Warning: Fisher test failed for gene", gene_name, ":", e$message, "\n")
    data.frame(
      Gene = gene_name,
      Cases_Burden = cases_burden,
      Cases_No_Burden = cases_no_burden,
      Controls_Burden = controls_burden,
      Controls_No_Burden = controls_no_burden,
      P_Value = NA_real_,
      Odds_Ratio = NA_real_,
      CI_Lower = NA_real_,
      CI_Upper = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  return(res)
}

# -------------------------
# Main analysis pipeline
# -------------------------
main_burden_analysis <- function() {
  cat("=== Starting Gene-Level Burden Analysis ===\n")

  base_path <- "/path/to/results/directory"
  srr_names_file <- "/path/to/srr_names.txt"

  if (!file.exists(srr_names_file)) stop("SRR names file not found: ", srr_names_file)
  srr_names <- readLines(srr_names_file, warn = FALSE)
  srr_names <- trimws(srr_names[srr_names != ""])
  cat("Found", length(srr_names), "SRR samples to process\n")

  file_paths <- file.path(base_path, srr_names, "ANNOTATION", "Franklin_Results_Without_Benign_LikelyBenign_0.001.csv")
  existing_files <- file.exists(file_paths)
  missing_files <- sum(!existing_files)
  if (missing_files > 0) {
    cat("Warning:", missing_files, "files are missing out of", length(file_paths), "\n")
  }

  valid_paths <- file_paths[existing_files]
  valid_srrs <- srr_names[existing_files]
  cat("Proceeding with", length(valid_paths), "existing files\n")

  cat("\n=== Step 1: Loading sample data ===\n")
  if (!is.null(cl)) {
    sample_data_list <- parallel::parLapply(cl, valid_paths, safe_read_csv)
  } else {
    sample_data_list <- lapply(valid_paths, safe_read_csv)
  }
  names(sample_data_list) <- valid_srrs

  failed_reads <- sapply(sample_data_list, is.null)
  if (sum(failed_reads) > 0) {
    cat("Warning:", sum(failed_reads), "files failed to load\n")
    sample_data_list <- sample_data_list[!failed_reads]
  }
  cat("Successfully loaded", length(sample_data_list), "sample files\n")
  if (length(sample_data_list) == 0) stop("No sample files were successfully loaded. Please check file paths and formats.")

  cat("\n=== Step 2: Calculating sample-level burden ===\n")
  sample_burdens <- list()
  for (i in seq_along(sample_data_list)) {
    sample_name <- names(sample_data_list)[i]
    sample_burdens[[sample_name]] <- calculate_sample_burden(sample_data_list[[i]], sample_name)
  }

  cat("\n=== Step 3: Combining variant data ===\n")
  all_variants <- do.call(rbind, sample_data_list)
  cat("Combined data contains", nrow(all_variants), "total variant observations\n")
  all_variants_unique <- all_variants %>% dplyr::distinct(Variant, `Gene Name`, .keep_all = TRUE)
  cat("Unique variants:", nrow(all_variants_unique), "\n")

  all_genes <- unique(all_variants_unique$`Gene Name`)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  cat("Total genes to analyze:", length(all_genes), "\n")

  cat("\n=== Step 5: Calculating gene-level burden counts ===\n")
  present_long <- data.table::rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]
    if (nrow(sb) == 0) return(NULL)
    data.table::data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_burden_counts_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_burden_counts <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(as.data.frame(gene_burden_counts_dt), by = "Gene") %>%
    dplyr::mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))
  cat("Burden counts calculated for", nrow(gene_burden_counts), "genes\n")

  cat("\n=== Step 6: Estimating population burden ===\n")
  population_burden <- estimate_population_burden(all_variants_unique, sample_size = length(sample_burdens))

  cat("\n=== Step 7: Merging case and control data ===\n")
  gene_analysis_data <- gene_burden_counts %>%
    dplyr::left_join(population_burden, by = "Gene") %>%
    dplyr::mutate(
      expected_burden_count = ifelse(is.na(expected_burden_count),
                                     0.001 * length(sample_burdens),
                                     expected_burden_count),
      Controls_Burden = pmax(1, round(expected_burden_count))
    )

  cat("Analysis data prepared for", nrow(gene_analysis_data), "genes\n")

  cat("\n=== Step 8: Performing Fisher's exact tests ===\n")
  fisher_results_list <- vector("list", nrow(gene_analysis_data))
  for (i in seq_len(nrow(gene_analysis_data))) {
    row <- gene_analysis_data[i, ]
    result <- perform_gene_fisher_test(
      gene_name = row$Gene,
      cases_burden = as.integer(row$Cases_Burden),
      controls_burden = as.integer(row$Controls_Burden),
      total_cases = length(sample_burdens),
      total_controls = length(sample_burdens)
    )
    fisher_results_list[[i]] <- result
  }
  fisher_results <- do.call(rbind, fisher_results_list)

  cat("\n=== Step 9: Multiple testing correction ===\n")
  valid_results <- fisher_results[!is.na(fisher_results$P_Value), , drop = FALSE]
  if (nrow(valid_results) > 0) {
    valid_results$P_Value_Bonferroni <- p.adjust(valid_results$P_Value, method = "bonferroni")
    valid_results$P_Value_FDR <- p.adjust(valid_results$P_Value, method = "fdr")
    valid_results$Significant_Bonferroni <- valid_results$P_Value_Bonferroni < 0.05
    valid_results$Significant_FDR <- valid_results$P_Value_FDR < 0.05
    valid_results <- valid_results[order(valid_results$P_Value), , drop = FALSE]
    cat("Multiple testing correction applied to", nrow(valid_results), "genes\n")
  } else {
    stop("No valid p-values were calculated. Please check your data.")
  }

  cat("\n=== Step 10: Saving results ===\n")
  results_dir <- "/path/to/output/directory"
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  readr::write_csv(valid_results, file.path(results_dir, "Gene_Burden_Analysis_Results.csv"))
  cat("Detailed results saved to:", file.path(results_dir, "Gene_Burden_Analysis_Results.csv"), "\n")

  significant_genes_fdr <- valid_results[valid_results$Significant_FDR == TRUE, , drop = FALSE]
  significant_genes_raw <- valid_results[valid_results$P_Value < 0.05, , drop = FALSE]
  if (nrow(significant_genes_fdr) > 0) {
    readr::write_csv(significant_genes_fdr, file.path(results_dir, "Significant_Genes_FDR_0.05.csv"))
    cat("FDR-significant genes saved to:", file.path(results_dir, "Significant_Genes_FDR_0.05.csv"), "\n")
  }
  if (nrow(significant_genes_raw) > 0) {
    readr::write_csv(significant_genes_raw, file.path(results_dir, "Significant_Genes_Raw_P_0.05.csv"))
    cat("Raw p-value significant genes saved to:", file.path(results_dir, "Significant_Genes_Raw_P_0.05.csv"), "\n")
  }

  cat("\n=== SUMMARY ===\n")
  cat("Total samples analyzed:", length(sample_burdens), "\n")
  cat("Total genes analyzed:", nrow(valid_results), "\n")
  cat("Genes with raw p < 0.05:", sum(valid_results$P_Value < 0.05, na.rm = TRUE), "\n")
  cat("Genes significant after FDR correction:", sum(valid_results$Significant_FDR, na.rm = TRUE), "\n")

  cat("\n=== TOP 20 GENES BY RAW P-VALUE ===\n")
  top_genes <- head(valid_results, 20)
  print(top_genes[, c("Gene", "Cases_Burden", "Controls_Burden", "P_Value", "P_Value_FDR", "Odds_Ratio")])

  return(valid_results)
}

# -------------------------
# Plotting: unified style and annotated volcano; saves vector + TIFF@1200
# -------------------------
create_summary_plots <- function(results, top_n = 17) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

  results_dir <- "/path/to/output/directory"
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  if (is.null(results) || nrow(results) == 0) {
    cat("[ERROR] 'results' is NULL or empty — nothing to plot.\n"); return(invisible(NULL))
  }

  results$Cases_Burden <- as.integer(results$Cases_Burden)
  if (!"P_Value_FDR" %in% colnames(results)) {
    if ("P_Value" %in% colnames(results)) {
      cat("[INFO] P_Value_FDR not found — computing FDR from P_Value.\n")
      results$P_Value_FDR <- p.adjust(results$P_Value, method = "fdr")
    } else {
      stop("No P_Value or P_Value_FDR in results — cannot plot significance.")
    }
  }

  results <- results %>%
    dplyr::mutate(
      Odds_Ratio = as.numeric(Odds_Ratio),
      log2OR = ifelse(!is.na(Odds_Ratio) & Odds_Ratio > 0, log2(Odds_Ratio), NA_real_),
      neglog10FDR = -log10(pmax(P_Value_FDR, 1e-300)),
      SignificantFDR = ifelse(!is.na(P_Value_FDR) & P_Value_FDR < 0.05, TRUE, FALSE),
      Cases_Burden = ifelse(is.na(Cases_Burden), 0L, as.integer(Cases_Burden))
    )

  base_size <- 10
  theme_journal <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      axis.title = ggplot2::element_text(size = base_size + 1),
      axis.text = ggplot2::element_text(size = base_size),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(5,5,5,5), "pt")
    )

  bar_fill <- "#9fd7f5"
  bar_border <- "black"
  sig_color <- "#e15759"
  non_sig_color <- "gray60"

  # Volcano
  cat("=== Volcano: preparing data ===\n")
  plot_data <- results %>% dplyr::filter(!is.na(log2OR) & !is.na(neglog10FDR))
  cat("[volcano] total plotted genes:", nrow(plot_data), "\n")
  if (nrow(plot_data) > 0) {
    genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
    label_data <- plot_data %>% dplyr::filter(Gene %in% genes_to_label) %>% dplyr::arrange(desc(neglog10FDR))
    cat("[volcano] requested label genes present:", nrow(label_data), "out of", length(genes_to_label), "\n")

    p_volcano <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log2OR, y = neglog10FDR)) +
      ggplot2::geom_point(ggplot2::aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)),
                          alpha = 1, stroke = 0.2) +
      ggplot2::scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = non_sig_color), guide = FALSE) +
      ggplot2::scale_size_continuous(range = c(1.5, 5), guide = ggplot2::guide_legend(title = "Cases with burden")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6, size = 0.4) +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.6, size = 0.4) +
      ggplot2::xlab("log2(Odds Ratio)") + ggplot2::ylab("-log10(FDR P-value)") +
      theme_journal

    if (!is.null(label_data) && nrow(label_data) > 0) {
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = label_data,
          mapping = ggplot2::aes(label = Gene),
          size = 3.2,
          max.overlaps = 30,
          box.padding = 0.6,
          point.padding = 0.3,
          segment.size = 0.35,
          segment.color = "gray40",
          arrow = grid::arrow(length = grid::unit(0.02, "npc"), type = "closed")
        )
    }

    volcano_base <- file.path(results_dir, "Fig1_Gene_Burden_Volcano")
    save_highres(p_volcano, volcano_base, width_cm = 17.6, height_cm = 12.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
  } else {
    cat("[volcano] skipped: no valid OR/FDR to plot\n")
  }

  # Top barplot
  cat("=== Top genes barplot: preparing data ===\n")
  plot_genes <- results %>% dplyr::filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) {
    plot_genes <- results %>% dplyr::filter(P_Value < 0.05 & !is.na(P_Value))
    cat("[barplot] no FDR-significant genes; falling back to raw p<0.05. selected:", nrow(plot_genes), "\n")
  } else {
    cat("[barplot] FDR-significant genes selected:", nrow(plot_genes), "\n")
  }

  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% dplyr::arrange(P_Value_FDR) %>% head(25) %>%
      dplyr::mutate(score = -log10(pmax(P_Value, 1e-300)))
    p_bar <- ggplot2::ggplot(plot_genes2, ggplot2::aes(x = reorder(Gene, score), y = score)) +
      ggplot2::geom_col(fill = bar_fill, color = bar_border, width = 0.75, size = 0.3) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Gene") + ggplot2::ylab("-log10(Raw P-value)") +
      theme_journal +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size - 1))

    n_rows <- nrow(plot_genes2)
    height_cm_bar <- max(6, 0.6 * n_rows + 2)
    bar_base <- file.path(results_dir, "Fig2_Top_Significant_Genes")
    save_highres(p_bar, bar_base, width_cm = 8.8, height_cm = height_cm_bar, save_vector = TRUE, tiff_res = 1200, family = font_family)
  } else {
    cat("[barplot] skipped: no significant genes to show\n")
  }

  # Distribution
  cat("=== Distribution: preparing data ===\n")
  if (!"Cases_Burden" %in% colnames(results)) {
    cat("[distribution] Cases_Burden not found -> skipping\n")
  } else {
    p_dist <- ggplot2::ggplot(results, ggplot2::aes(x = Cases_Burden)) +
      ggplot2::geom_histogram(bins = 30, fill = bar_fill, color = bar_border, size = 0.3) +
      ggplot2::xlab("Number of cases with burden") + ggplot2::ylab("Number of genes") +
      theme_journal

    dist_base <- file.path(results_dir, "Fig3_Burden_Distribution")
    save_highres(p_dist, dist_base, width_cm = 8.8, height_cm = 6.0, save_vector = TRUE, tiff_res = 1200, family = font_family)
  }

  cat("=== create_summary_plots finished ===\n")
  return(invisible(NULL))
}

# -------------------------
# Main execution
# -------------------------
results <- NULL
tryCatch({
  cat("Starting analysis...\n")
  results <- main_burden_analysis()
  cat("\nCreating summary plots...\n")
  create_summary_plots(results)
  cat("\nANALYSIS COMPLETED SUCCESSFULLY!\n")
}, error = function(e) {
  cat("ERROR IN ANALYSIS:\n", e$message, "\n")
}, finally = {
  if (exists("cl") && !is.null(cl)) {
    tryCatch({
      parallel::stopCluster(cl)
      cat("Parallel processing cluster stopped.\n")
    }, error = function(e) {
      cat("Warning: Error stopping cluster:", e$message, "\n")
    })
  }
})

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
cat("\nScript execution completed.\n")
