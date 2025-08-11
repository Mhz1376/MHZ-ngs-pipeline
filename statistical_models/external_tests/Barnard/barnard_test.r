# Gene-Level Burden Analysis using Barnard's test (DescTools::BarnardTest)
# Label only top FDR-significant genes (by -log10 FDR), axis x fixed to [0,15]
# Outputs CSV: Gene_Burden_Analysis_Results.csv, Significant_Genes_FDR_0.05.csv, Significant_Genes_Raw_P_0.05.csv
# Outputs PNG: Burden_Distribution_Styled.png, Gene_Burden_Volcano_Plot_annotated.png, Top_Significant_Genes_Styled.png
# Author: adapted Aug 2025

rm(list = ls())

## ---------- CONFIG ----------
RESULTS_DIR <- "/path/to/output/directory"
BASE_PATH   <- "/path/to/results/directory"
SRR_FILE    <- "/path/to/srr_names.txt"  # file with SRR names, one per line
VERBOSE     <- TRUE

HALDANE_PC  <- 0.5    # Haldane continuity for OR reporting
LOG2OR_MIN  <- 0
LOG2OR_MAX  <- 15
TOP_LABEL_N <- 17     # number of FDR-significant genes to label (top by -log10 FDR)

## ---------- PACKAGES ----------
req <- c("parallel","data.table","dplyr","readr","ggplot2","ggrepel","scales","DescTools")
for (p in req) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
}

## ---------- PARALLEL ----------
available_cores <- parallel::detectCores(logical = FALSE)
num_cores <- min(12, ifelse(is.na(available_cores), 1, max(1, available_cores - 1)))
if (VERBOSE) message("Using ", num_cores, " cores")
cl <- NULL
if (num_cores > 1) {
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, { library(data.table); library(dplyr); library(readr); library(DescTools); TRUE })
  clusterExport(cl, varlist = c("VERBOSE"), envir = environment())
}

## ---------- HELPERS ----------
safe_read_csv <- function(path) {
  tryCatch({
    if (VERBOSE) message("Reading: ", basename(path))
    df <- readr::read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
    # find gene column
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

    if (VERBOSE) message(" -> read ", nrow(df), " variants from ", basename(path))
    df
  }, error = function(e) {
    if (VERBOSE) message("Read error for ", basename(path), ": ", e$message)
    NULL
  })
}

estimate_population_burden <- function(all_variants, sample_size = 62L) {
  if (VERBOSE) message("Estimating population burden...")
  all_variants <- all_variants %>%
    mutate(variant_AF = pmax(gnomAD_AF, PopFreqMax, na.rm = TRUE)) %>%
    mutate(variant_AF = pmin(pmax(variant_AF, 0), 1))
  gp <- all_variants %>%
    group_by(`Gene Name`) %>%
    summarise(variant_afs = list(variant_AF), variant_count = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(carrier_probs = list(pmin(2 * unlist(variant_afs), 1)),
           log_non = sum(log1p(-unlist(carrier_probs))),
           individual_burden_freq = 1 - exp(log_non),
           expected_burden_count = sample_size * individual_burden_freq) %>%
    ungroup() %>%
    select(Gene = `Gene Name`, variant_count, individual_burden_freq, expected_burden_count)
  if (VERBOSE) message(" -> estimated for ", nrow(gp), " genes")
  gp
}

# Barnard wrapper
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

  # Edge cases
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
    # fallback to Fisher with Haldane correction
    a <- cases_burden + haldane_pc; b <- controls_int + haldane_pc
    c_ <- cases_no_burden + haldane_pc; d <- controls_no + haldane_pc
    ft <- tryCatch(stats::fisher.test(matrix(c(a,b,c_,d), nrow = 2, byrow = TRUE)), error = function(e) NULL)
    pval <- if (!is.null(ft)) ft$p.value else NA_real_
  }

  # Haldane OR for reporting (use rounded control count for OR stability)
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

## ---------- MAIN ----------
main_burden_analysis <- function() {
  message("Starting Barnard pipeline...")
  if (!file.exists(SRR_FILE)) stop("SRR names file not found: ", SRR_FILE)
  srr <- readLines(SRR_FILE, warn = FALSE) %>% trimws()
  srr <- srr[srr != ""]
  message("Found ", length(srr), " samples")

  paths <- file.path(BASE_PATH, srr, "ANNOTATION", "Franklin_Results_Without_Benign_LikelyBenign_0.001.csv")
  exist <- file.exists(paths)
  if (any(!exist)) warning(sum(!exist), " files missing")
  paths <- paths[exist]; srr <- srr[exist]
  if (length(paths) == 0) stop("No files to read")

  # ensure workers know safe_read_csv if using cluster
  if (!is.null(cl)) clusterExport(cl, varlist = c("safe_read_csv","VERBOSE"), envir = environment())

  # read samples
  if (!is.null(cl)) {
    sample_list <- parLapply(cl, paths, safe_read_csv)
  } else sample_list <- lapply(paths, safe_read_csv)
  names(sample_list) <- srr
  sample_list <- sample_list[!sapply(sample_list, is.null)]
  if (length(sample_list) == 0) stop("No sample data loaded")

  # sample burdens (per-sample genes with at least one variant)
  sample_burdens <- lapply(names(sample_list), function(nm) {
    df <- sample_list[[nm]]
    df %>% distinct(`Gene Name`) %>% rename(Gene = `Gene Name`)
  })
  names(sample_burdens) <- names(sample_list)

  # combine variants and unique set
  combined <- do.call(rbind, sample_list)
  combined_unique <- combined %>% distinct(Variant, `Gene Name`, .keep_all = TRUE)
  genes <- unique(combined_unique$`Gene Name`) %>% na.omit()
  genes <- genes[genes != ""]

  # cases per gene
  present_long <- data.table::rbindlist(lapply(names(sample_burdens), function(s) {
    sb <- sample_burdens[[s]]; if (nrow(sb) == 0) return(NULL); data.table(Gene = sb$Gene, Sample = s)
  }), use.names = TRUE, fill = TRUE)
  gene_cases_dt <- present_long[, .(Cases_Burden = uniqueN(Sample)), by = Gene]
  gene_cases <- data.frame(Gene = genes, stringsAsFactors = FALSE) %>%
    left_join(as.data.frame(gene_cases_dt), by = "Gene") %>%
    mutate(Cases_Burden = ifelse(is.na(Cases_Burden), 0L, Cases_Burden))

  # estimate population burden (fractional expected controls)
  pop_burden <- estimate_population_burden(combined_unique, sample_size = length(sample_burdens))

  # merge - keep fractional expected for plotting; Barnard uses rounded internally
  gene_df <- gene_cases %>%
    left_join(pop_burden, by = "Gene") %>%
    mutate(expected_burden_count = ifelse(is.na(expected_burden_count), 0.0, expected_burden_count),
           Controls_Burden_Fraction = expected_burden_count)

  message("Prepared ", nrow(gene_df), " genes")

  # export to workers (function + data)
  if (!is.null(cl)) {
    clusterExport(cl, varlist = c("perform_gene_barnard_test", "gene_df", "HALDANE_PC"), envir = environment())
  }

  # run tests (parallel over indices)
  n <- nrow(gene_df)
  idx <- seq_len(n)
  if (!is.null(cl) && num_cores > 1) {
    res_list <- parLapply(cl, idx, function(i) {
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
      if (i %% 200 == 0 && VERBOSE) message("Done ", i, " genes")
    }
  }

  results <- data.table::rbindlist(res_list, fill = TRUE) %>% as.data.frame()

  # multiple testing corrections
  valid <- results[!is.na(results$P_Value), , drop = FALSE]
  if (nrow(valid) == 0) stop("No valid p-values")
  valid$P_Value_Bonferroni <- p.adjust(valid$P_Value, method = "bonferroni")
  valid$P_Value_FDR <- p.adjust(valid$P_Value, method = "fdr")
  valid$Significant_Bonferroni <- valid$P_Value_Bonferroni < 0.05
  valid$Significant_FDR <- valid$P_Value_FDR < 0.05
  valid <- valid[order(valid$P_Value), ]

  # save CSVs with requested names
  if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
  readr::write_csv(valid, file.path(RESULTS_DIR, "Gene_Burden_Analysis_Results.csv"))
  readr::write_csv(valid[valid$Significant_FDR == TRUE, , drop = FALSE], file.path(RESULTS_DIR, "Significant_Genes_FDR_0.05.csv"))
  readr::write_csv(valid[valid$P_Value < 0.05, , drop = FALSE], file.path(RESULTS_DIR, "Significant_Genes_Raw_P_0.05.csv"))
  message("Saved CSV outputs to ", RESULTS_DIR)

  valid
}

## ---------- PLOTTING (label like asked) ----------
## ---------- PLOTTING (label like asked — only 6 specific genes) ----------
create_summary_plots <- function(results, top_n = TOP_LABEL_N) {
  if (is.null(results) || nrow(results) == 0) { message("No results to plot"); return(invisible(NULL)) }
  # prepare plotting OR using fractional controls & Haldane
  if (!("Controls_Burden_Fraction" %in% colnames(results))) results$Controls_Burden_Fraction <- results$Controls_Burden

  # Ensure Cases_No_Burden present (should be from function); if missing, infer from total = Cases_Burden + Cases_No_Burden in rows where available.
  if (!("Cases_No_Burden" %in% colnames(results))) {
    # try to infer total from max per-row if possible; fallback to NA
    results$Cases_No_Burden <- NA_real_
  }

  # For each row, define total_controls as row-specific total (Cases_Burden + Cases_No_Burden) if available, otherwise use max observed total
  per_row_total <- ifelse(!is.na(results$Cases_No_Burden), (results$Cases_Burden + results$Cases_No_Burden), NA_real_)
  fallback_total <- ifelse(all(is.na(per_row_total)), NA_real_, max(per_row_total, na.rm = TRUE))
  total_controls_row <- ifelse(!is.na(per_row_total), per_row_total, fallback_total)
  # if still NA, use length of samples inferred from Cases_Burden max (safe fallback)
  if (is.na(total_controls_row[1])) {
    total_controls_row <- rep(max(results$Cases_Burden, na.rm = TRUE), nrow(results))
  }
  results$Controls_No_Burden_Fraction <- total_controls_row - results$Controls_Burden_Fraction
  # avoid negatives
  results$Controls_No_Burden_Fraction <- pmax(results$Controls_No_Burden_Fraction, 0)

  results <- results %>% mutate(
    OR_plot = ((Cases_Burden + HALDANE_PC) / (pmax(Cases_No_Burden, 0) + HALDANE_PC)) /
              ((Controls_Burden_Fraction + HALDANE_PC) / (Controls_No_Burden_Fraction + HALDANE_PC)),
    OR_plot = ifelse(is.infinite(OR_plot) | is.nan(OR_plot), NA_real_, OR_plot),
    log2OR = ifelse(!is.na(OR_plot) & OR_plot > 0, log2(OR_plot), NA_real_)
  )

  # clip to [LOG2OR_MIN, LOG2OR_MAX] and map negative to LOG2OR_MIN
  results$log2OR_plot <- results$log2OR
  results$log2OR_plot[is.na(results$log2OR_plot)] <- NA_real_
  results$log2OR_plot <- pmin(pmax(results$log2OR_plot, LOG2OR_MIN), LOG2OR_MAX)

  # FDR and -log10
  if (!"P_Value_FDR" %in% colnames(results)) results$P_Value_FDR <- p.adjust(results$P_Value, method = "fdr")
  results$neglog10FDR <- -log10(pmax(results$P_Value_FDR, 1e-300))
  results$SignificantFDR <- ifelse(!is.na(results$P_Value_FDR) & results$P_Value_FDR < 0.05, TRUE, FALSE)
  results$Cases_Burden <- ifelse(is.na(results$Cases_Burden), 0L, as.integer(results$Cases_Burden))

  # Volcano data
  plot_data <- results %>% filter(!is.na(log2OR_plot) & !is.na(neglog10FDR))
  message("Volcano will plot ", nrow(plot_data), " genes; FDR-significant: ", sum(plot_data$SignificantFDR, na.rm = TRUE))

  # ----------------- MODIFIED LABELING: only these 6 genes (if present) -----------------
  genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")

  # select label rows from plot_data (only if present there)
  label_data <- NULL
  if (nrow(plot_data) > 0) {
    label_data <- plot_data %>% filter(Gene %in% genes_to_label) %>% arrange(desc(neglog10FDR))
  }

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

  if (nrow(plot_data) > 0) {
    p_volcano <- ggplot(plot_data, aes(x = log2OR_plot, y = neglog10FDR)) +
      geom_point(aes(color = SignificantFDR, size = pmax(Cases_Burden, 1)), alpha = 0.85) +
      scale_color_manual(values = c("TRUE" = "#e15759", "FALSE" = "gray60"), guide = FALSE) +
      scale_size_continuous(range = c(1.2, 5), guide = guide_legend(title = "Cases with burden")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(1,2), linetype = "dashed") +
      labs(title = "Gene Burden Analysis — Volcano Plot (Barnard)",
           subtitle = paste0("Plotted genes: ", nrow(plot_data), "  |  FDR-significant: ", sum(plot_data$SignificantFDR, na.rm = TRUE)),
           x = paste0("log2(Odds Ratio) [clipped to ", LOG2OR_MIN, "..", LOG2OR_MAX, "]"),
           y = "-log10(FDR P-value)") +
      theme_journal + coord_cartesian(xlim = c(LOG2OR_MIN, LOG2OR_MAX), expand = TRUE)

    if (!is.null(label_data) && nrow(label_data) > 0) {
      # use the clipped x (log2OR_plot) for label x position
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = label_data,
          aes(x = log2OR_plot, y = neglog10FDR, label = Gene),
          size = 3.2,
          max.overlaps = 30,
          box.padding = 0.8,
          point.padding = 0.3,
          segment.size = 0.45,
          segment.color = "gray40",
          arrow = grid::arrow(length = grid::unit(0.015, "npc"), type = "closed")
        )
    }

    ggsave(filename = file.path(RESULTS_DIR, "Gene_Burden_Volcano_Plot_annotated.png"),
           plot = p_volcano, width = 11, height = 8, dpi = 300)
    message("Saved volcano: ", file.path(RESULTS_DIR, "Gene_Burden_Volcano_Plot_annotated.png"))
  } else message("No valid points to plot (volcano)")

  # Top barplot (unchanged)
  plot_genes <- results %>% filter(SignificantFDR == TRUE)
  if (nrow(plot_genes) == 0) plot_genes <- results %>% filter(P_Value < 0.05 & !is.na(P_Value))
  if (nrow(plot_genes) > 0) {
    plot_genes2 <- plot_genes %>% arrange(P_Value) %>% head(50) %>% mutate(score = -log10(pmax(P_Value, 1e-300)))
    p_bar <- ggplot(plot_genes2, aes(x = reorder(Gene, score), y = score)) +
      geom_col(fill = "#9fd7f5", color = "black") + coord_flip() +
      labs(title = "Top Genes by Significance (Barnard)", x = "Gene", y = "-log10(P)") + theme_journal
    height_needed <- max(6, 0.35 * nrow(plot_genes2) + 2)
    ggsave(filename = file.path(RESULTS_DIR, "Top_Significant_Genes_Styled.png"), plot = p_bar, width = 12, height = height_needed, dpi = 300)
    message("Saved barplot: ", file.path(RESULTS_DIR, "Top_Significant_Genes_Styled.png"))
  } else message("No significant genes to plot (barplot)")

  # Distribution (unchanged)
  if ("Cases_Burden" %in% colnames(results)) {
    p_dist <- ggplot(results, aes(x = Cases_Burden)) +
      geom_histogram(bins = 30, fill = "#9fd7f5", color = "black") +
      labs(title = "Distribution of Gene Burden Counts (Cases)", subtitle = paste("Across", nrow(results), "genes"),
           x = "Number of cases with burden", y = "Number of genes") +
      theme_journal
    ggsave(filename = file.path(RESULTS_DIR, "Burden_Distribution_Styled.png"), plot = p_dist, width = 10, height = 6, dpi = 300)
    message("Saved distribution: ", file.path(RESULTS_DIR, "Burden_Distribution_Styled.png"))
  } else message("Cases_Burden missing -> skip distribution")
}

## ---------- RUN ----------
res <- NULL
tryCatch({
  res <- main_burden_analysis()
  create_summary_plots(res, top_n = TOP_LABEL_N)
  message("Analysis completed successfully.")
}, error = function(e) {
  message("ERROR: ", e$message)
}, finally = {
  if (!is.null(cl)) {
    try(stopCluster(cl), silent = TRUE)
    message("Cluster stopped.")
  }
  message("Session info:"); print(sessionInfo())
})
