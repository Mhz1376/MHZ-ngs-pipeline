# Gene-Level Burden Analysis using Fisher's Exact Test
# Updated plotting: unified journal-style theme + annotated volcano (all FDR < 0.05 labeled)
# Author: Updated for user's pipeline (Aug 2025)

rm(list = ls())

required_packages <- c("parallel", "data.table", "dplyr", "readr", "ggplot2", "ggrepel", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Setup parallel (robust)
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
cat("Setting up parallel processing with", num_cores, "cores...\n")

cl <- NULL
tryCatch({
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    clusterEvalQ(cl, {
      library(data.table)
      library(dplyr)
      library(readr)
    })
    cat("Parallel cluster initialized successfully\n")
  } else {
    cat("Single-core mode (no cluster created)\n")
  }
}, error = function(e) {
  cat("Warning: Failed to create parallel cluster. Using sequential processing.\n")
  cl <<- NULL
})

# Safe CSV reader (no classification filtering assumed; input pre-filtered)
safe_read_csv <- function(file_path) {
  tryCatch({
    cat("Reading:", basename(file_path), "\n")
    df <- read_csv(file_path, show_col_types = FALSE,
                   col_types = cols(.default = col_character()))
    if (!("Gene Name" %in% colnames(df))) {
      alt <- grep("gene", colnames(df), ignore.case = TRUE, value = TRUE)
      if (length(alt) >= 1) names(df)[which(colnames(df) == alt[1])] <- "Gene Name"
    }
    if (!("Variant" %in% colnames(df))) df$Variant <- paste0("var_", seq_len(nrow(df)))
    if (!("gnomAD_AF" %in% colnames(df))) df$gnomAD_AF <- NA
    if (!("PopFreqMax" %in% colnames(df))) df$PopFreqMax <- NA
    if (!("Gene Name" %in% colnames(df))) df$`Gene Name` <- NA

    df$`Gene Name` <- trimws(df$`Gene Name`)
    df <- df[!is.na(df$`Gene Name`) & df$`Gene Name` != "", ]

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

# Calculate sample-level gene burden (list of genes that have >=1 burden variant)
calculate_sample_burden <- function(sample_data, sample_name = "Unknown") {
  if (is.null(sample_data) || nrow(sample_data) == 0) {
    return(data.frame(Gene = character(0), HasBurden = logical(0)))
  }
  gene_burden <- sample_data %>%
    distinct(`Gene Name`) %>%
    mutate(HasBurden = TRUE) %>%
    rename(Gene = `Gene Name`)
  cat("Sample", sample_name, "has burden variants in", nrow(gene_burden), "genes\n")
  return(gene_burden)
}

# Exact population burden using carrier probs and product formula
estimate_population_burden <- function(all_variants_data, sample_size = 62) {
  cat("Estimating population burden (exact product method) for", length(unique(all_variants_data$`Gene Name`)), "genes...\n")

  all_variants_data <- all_variants_data %>%
    mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))

  gene_pop_freq <- all_variants_data %>%
    group_by(`Gene Name`) %>%
    summarise(
      variant_afs = list(variant_AF),
      variant_count = n(),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
      log_non = sum(log1p(-unlist(carrier_probs))),
      individual_burden_freq = 1 - exp(log_non),
      expected_burden_count = sample_size * individual_burden_freq
    ) %>%
    ungroup() %>%
    select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)

  cat("Population burden estimated for", nrow(gene_pop_freq), "genes\n")
  return(gene_pop_freq)
}

# Corrected Fisher test
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
    fisher_result <- fisher.test(contingency_matrix)
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

# -----------------------
# Main analysis pipeline
# -----------------------
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
    sample_data_list <- parLapply(cl, valid_paths, safe_read_csv)
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
  all_variants_unique <- all_variants %>% distinct(Variant, `Gene Name`, .keep_all = TRUE)
  cat("Unique variants:", nrow(all_variants_unique), "\n")

  all_genes <- unique(all_variants_unique$`Gene Name`)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  cat("Total genes to analyze:", length(all_genes), "\n")

  cat("\n=== Step 5: Calculating gene-level burden counts ===\n")
  library(data.table)
  present_long <- rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]
    if (nrow(sb) == 0) return(NULL)
    data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_burden_counts_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_burden_counts <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
    left_join(as.data.frame(gene_burden_counts_dt), by = "Gene") %>%
    mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))
  cat("Burden counts calculated for", nrow(gene_burden_counts), "genes\n")

  cat("\n=== Step 6: Estimating population burden ===\n")
  population_burden <- estimate_population_burden(all_variants_unique, sample_size = length(sample_burdens))

  cat("\n=== Step 7: Merging case and control data ===\n")
  gene_analysis_data <- gene_burden_counts %>%
    left_join(population_burden, by = "Gene") %>%
    mutate(
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
  valid_results <- fisher_results[!is.na(fisher_results$P_Value), ]
  if (nrow(valid_results) > 0) {
    valid_results$P_Value_Bonferroni <- p.adjust(valid_results$P_Value, method = "bonferroni")
    valid_results$P_Value_FDR <- p.adjust(valid_results$P_Value, method = "fdr")
    valid_results$Significant_Bonferroni <- valid_results$P_Value_Bonferroni < 0.05
    valid_results$Significant_FDR <- valid_results$P_Value_FDR < 0.05
    valid_results <- valid_results[order(valid_results$P_Value), ]
    cat("Multiple testing correction applied to", nrow(valid_results), "genes\n")
  } else {
    stop("No valid p-values were calculated. Please check your data.")
  }

  cat("\n=== Step 10: Saving results ===\n")
  results_dir <- "/path/to/output/directory"
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  write_csv(valid_results, file.path(results_dir, "Gene_Burden_Analysis_Results.csv"))
  cat("Detailed results saved to:", file.path(results_dir, "Gene_Burden_Analysis_Results.csv"), "\n")

  significant_genes_fdr <- valid_results[valid_results$Significant_FDR == TRUE, ]
  significant_genes_raw <- valid_results[valid_results$P_Value < 0.05, ]
  if (nrow(significant_genes_fdr) > 0) {
    write_csv(significant_genes_fdr, file.path(results_dir, "Significant_Genes_FDR_0.05.csv"))
    cat("FDR-significant genes saved to:", file.path(results_dir, "Significant_Genes_FDR_0.05.csv"), "\n")
  }
  if (nrow(significant_genes_raw) > 0) {
    write_csv(significant_genes_raw, file.path(results_dir, "Significant_Genes_Raw_P_0.05.csv"))
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

# Plotting: unified style and annotated volcano (label all FDR < 0.05)
# ----------------- VOLCANO (modified: label only 6 specific genes) -----------------
create_summary_plots <- function(results, top_n = 17) {
  # dependencies
  if (!require(ggplot2, quietly = TRUE)) stop("ggplot2 required")
  if (!require(ggrepel, quietly = TRUE)) { install.packages("ggrepel"); library(ggrepel) }
  if (!require(scales, quietly = TRUE)) { install.packages("scales"); library(scales) }

  results_dir <- "/path/to/output/directory"
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # sanity checks / prepare columns
  if (is.null(results) || nrow(results) == 0) {
    cat("[ERROR] 'results' is NULL or empty — nothing to plot.\n"); return(invisible(NULL))
  }

  # ensure numeric types (same as original)
  results$Cases_Burden <- as.integer(results$Cases_Burden)
  if (!"P_Value_FDR" %in% colnames(results)) {
    if ("P_Value" %in% colnames(results)) {
      cat("[INFO] P_Value_FDR not found — computing FDR from P_Value.\n")
      results$P_Value_FDR <- p.adjust(results$P_Value, method = "fdr")
    } else {
      stop("No P_Value or P_Value_FDR in results — cannot plot significance.")
    }
  }

  # compute log2OR and -log10(FDR) safely (same as original)
  results <- results %>%
    mutate(
      Odds_Ratio = as.numeric(Odds_Ratio),
      log2OR = ifelse(!is.na(Odds_Ratio) & Odds_Ratio > 0, log2(Odds_Ratio), NA_real_),
      neglog10FDR = -log10(pmax(P_Value_FDR, 1e-300)),
      SignificantFDR = ifelse(!is.na(P_Value_FDR) & P_Value_FDR < 0.05, TRUE, FALSE),
      Cases_Burden = ifelse(is.na(Cases_Burden), 0L, as.integer(Cases_Burden))
    )

  # plotting theme and colors (unchanged)
  base_size <- 14
  theme_journal <- theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  bar_fill <- "#9fd7f5"
  bar_border <- "black"
  sig_color <- "#e15759"
  non_sig_color <- "gray60"

  # ----------------- VOLCANO -----------------
  cat("=== Volcano: preparing data ===\n")
  plot_data <- results %>% filter(!is.na(log2OR) & !is.na(neglog10FDR))
  cat("[volcano] total plotted genes:", nrow(plot_data), "\n")
  if (nrow(plot_data) > 0) {
    # ---- MODIFIED PART: choose exactly these 6 genes to label (if present) ----
    genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")

    # restrict label candidates to the intersection with plot_data
    label_data <- plot_data %>%
      filter(Gene %in% genes_to_label) %>%
      arrange(desc(neglog10FDR))

    cat("[volcano] requested label genes present:", nrow(label_data), "out of", length(genes_to_label), "\n")

    # build volcano plot (same style as before)
    p_volcano <- ggplot(plot_data, aes(x = log2OR, y = neglog10FDR)) +
      geom_point(aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)), alpha = 0.85) +
      scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = non_sig_color), guide = FALSE) +
      scale_size_continuous(range = c(1.2, 4.5), guide = guide_legend(title = "Cases with burden")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.6) +
      labs(title = "Gene Burden Analysis — Volcano Plot",
           subtitle = paste0("Plotted genes: ", nrow(plot_data), "  |  FDR-significant: ", sum(plot_data$SignificantFDR, na.rm = TRUE)),
           x = "log2(Odds Ratio)",
           y = "-log10(FDR P-value)") +
      theme_journal

    # add only the selected labels (if any)
    if (!is.null(label_data) && nrow(label_data) > 0) {
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = Gene),
          size = 3.2,
          max.overlaps = 30,
          box.padding = 1,
          point.padding = 0.3,
          segment.size = 0.45,
          segment.color = "gray40",
          arrow = grid::arrow(length = grid::unit(0.015, "npc"), type = "closed")
        )
    }

    ggsave(filename = file.path(results_dir, "Gene_Burden_Volcano_Plot_annotated.png"),
           plot = p_volcano, width = 11, height = 8, dpi = 300)
    cat("Saved:", file.path(results_dir, "Gene_Burden_Volcano_Plot_annotated.png"), "\n")
  } else {
    cat("[volcano] skipped: no valid OR/FDR to plot\n")
  }

  # ----------------- TOP SIGNIFICANT BARPLOT (unchanged) -----------------
  cat("=== Top genes barplot: preparing data ===\n")
  plot_genes <- results %>% filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) {
    plot_genes <- results %>% filter(P_Value < 0.05 & !is.na(P_Value))
    cat("[barplot] no FDR-significant genes; falling back to raw p<0.05. selected:", nrow(plot_genes), "\n")
  } else {
    cat("[barplot] FDR-significant genes selected:", nrow(plot_genes), "\n")
  }

  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% arrange(P_Value_FDR) %>% head(25) %>%
      mutate(score = -log10(pmax(P_Value, 1e-300)))
    p_bar <- ggplot(plot_genes2, aes(x = reorder(Gene, score), y = score)) +
      geom_col(fill = bar_fill, color = bar_border, width = 0.75) +
      coord_flip() +
      labs(title = "Top Genes by Significance (raw/FDR)",
           subtitle = paste0(nrow(plot_genes2), " top genes"),
           x = "Gene", y = "-log10(Raw P-value)") +
      theme_journal +
      theme(axis.text.y = element_text(size = base_size - 1))

    height_needed <- max(6, 0.4 * nrow(plot_genes2))
    ggsave(filename = file.path(results_dir, "Top_Significant_Genes_Styled.png"),
           plot = p_bar, width = 12, height = height_needed, dpi = 300)
    cat("Saved:", file.path(results_dir, "Top_Significant_Genes_Styled.png"), "\n")
  } else {
    cat("[barplot] skipped: no significant genes to show\n")
  }

  # ----------------- DISTRIBUTION (unchanged) -----------------
  cat("=== Distribution: preparing data ===\n")
  if (!"Cases_Burden" %in% colnames(results)) {
    cat("[distribution] Cases_Burden not found -> skipping\n")
  } else {
    p_dist <- ggplot(results, aes(x = Cases_Burden)) +
      geom_histogram(bins = 30, fill = bar_fill, color = bar_border) +
      labs(title = "Distribution of Gene Burden Counts (Cases)",
           subtitle = paste("Across", nrow(results), "genes"),
           x = "Number of cases with burden", y = "Number of genes") +
      theme_journal

    ggsave(filename = file.path(results_dir, "Burden_Distribution_Styled.png"),
           plot = p_dist, width = 10, height = 6, dpi = 300)
    cat("Saved:", file.path(results_dir, "Burden_Distribution_Styled.png"), "\n")
  }

  cat("=== create_summary_plots finished ===\n")
  return(invisible(NULL))
}

# ============================
# MAIN EXECUTION
# ============================
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
      stopCluster(cl)
      cat("Parallel processing cluster stopped.\n")
    }, error = function(e) {
      cat("Warning: Error stopping cluster:", e$message, "\n")
    })
  }
})

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
cat("\nScript execution completed.\n")
