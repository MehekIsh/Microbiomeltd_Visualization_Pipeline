plot_microbiome <- function(
    feature_table,
    taxonomy_file,
    metadata_file,
    tree_file = NULL,
    output_dir = "microbiome_outputs",
    plot_types = c("abundance","rarefaction","alpha","beta"),
    tax_ranks = c("Phylum","Class","Order","Family","Genus","Species"),
    dpi = 300,
    width = 10,
    height = 6,
    group_by = NULL,
    alpha_show_p = TRUE   
) {
  # ---------- packages ----------
  suppressPackageStartupMessages({
    req <- c("phyloseq","tidyverse","vegan","openxlsx","RColorBrewer","ggpubr","viridisLite","rstatix", "stringr")
    missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      stop(sprintf("Missing packages: %s\nInstall them, then retry.", paste(missing, collapse = ", ")))
    }
    lapply(req, function(p) library(p, character.only = TRUE))
  })
  
  # ---------- helpers ----------
  dir_create_all <- function(paths) lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE)
  
  # Prefer a column in metadata (not melted dfs). Falls back to a reasonable categorical.
  pick_group_name_from_metadata <- function(metadata, prefer = NULL) {
    cn <- colnames(metadata)
    if (!is.null(prefer) && prefer %in% cn) return(prefer)
    cand <- grep("group|site|location|body|condition|treatment|subject|sampletype",
                 cn, ignore.case = TRUE, value = TRUE)
    if (length(cand) == 1) return(cand)
    if (length(cand) > 1) {
      n_vals <- sapply(metadata[, cand, drop = FALSE], function(x) length(unique(as.character(x))))
      return(cand[which.max(n_vals)])
    }
    for (col in cn) {
      v <- metadata[[col]]
      if (!is.numeric(v) && length(unique(as.character(v))) <= nrow(metadata) / 2) return(col)
    }
    stop("Could not determine a grouping column. Pass `group_by=` explicitly.")
  }
  
  # Attach a single canonical grouping column to a sample-level df: `.__GROUP__`
  attach_group <- function(df_samples, metadata, group_name) {
    md <- metadata %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, !!group_name) %>%
      dplyr::rename(`.__GROUP__` = !!group_name)
    dplyr::left_join(df_samples, md, by = "Sample")
  }
  
  fmt_pct <- function(x) {
    s <- format(round(x, 2), nsmall = 2, trim = TRUE)
    sub("\\.?0+$", "", s)
  }
  
  col_fix <- function(v) {
    v <- tolower(v)
    v <- gsub("[^a-z0-9]+", "_", v)
    v <- gsub("_+", "_", v)
    v <- gsub("^_|_$", "", v)
    v
  }
  
  # ---------- load inputs ----------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  otu <- read.table(feature_table, sep = "\t", header = TRUE, row.names = 1,
                    skip = 1, comment.char = "", check.names = FALSE)
  taxonomy <- read.table(taxonomy_file, sep = "\t", header = TRUE, row.names = 1,
                         check.names = FALSE)
  metadata <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1,
                         check.names = FALSE)
  
  # taxonomy parsing
  tax <- taxonomy %>%
    dplyr::select(Taxon) %>%
    tidyr::separate(Taxon, into = c("Domain", tax_ranks), sep = "; ", fill = "right") %>%
    dplyr::mutate(across(everything(), ~gsub(".__", "", .))) %>%
    replace(is.na(.), "") %>%
    as.data.frame()
  rownames(tax) <- rownames(taxonomy)
  
  if ("Domain" %in% colnames(tax)) {
    tax$Domain[tax$Domain == ""] <- "Unclassified"
  }
  
  # backfill children from parent as "<parent>_unclassified"
  for (i in seq_along(colnames(tax))) {
    if (i == 1) next
    parent <- colnames(tax)[i - 1]; me <- colnames(tax)[i]
    empty <- tax[[me]] == "" | is.na(tax[[me]])
    tax[[me]][empty] <- paste0(tax[[parent]][empty], "_unclassified")
  }
  
  # collapse any chained suffixes like "..._unclassified_unclassified" -> "..._unclassified"
  collapse_unclassified <- function(x) {
    x <- gsub("(?:_unclassified)+$", "_unclassified", x, perl = TRUE)
    x
  }
  tax[] <- lapply(tax, collapse_unclassified)
  
  if (all(c("Genus","Species") %in% colnames(tax))) {
    
    simple_species_label <- function(gen, sp, order = c("genus_first","species_first")) {
      order <- match.arg(order)
      gen_raw <- trimws(gen)
      sp_raw  <- trimws(sp)
      
      # if species missing/unknown -> "<Genus>_unclassified"
      if (sp_raw == "" || tolower(sp_raw) == "unclassified") {
        return(paste0(gen_raw, "_unclassified"))
      }
      
      # normalize tokens
      norm <- function(x) {
        x <- gsub("[^A-Za-z0-9]+", "_", x)
        x <- gsub("_+", "_", x)
        x <- gsub("^_|_$", "", x)
        x
      }
      gen_n <- norm(gen_raw)
      sp_n  <- norm(sp_raw)
      
      # drop a leading genus from species if present (case-insensitive)
      sp_parts <- unlist(strsplit(sp_n, "_"))
      if (length(sp_parts) >= 1 && tolower(sp_parts[1]) == tolower(sub("_.*$", "", gen_n))) {
        sp_parts <- sp_parts[-1]
      }
      
      # pick only the first remaining token as the epithet (short label)
      if (!length(sp_parts)) return(paste0(gen_n, "_unclassified"))
      epithet <- sp_parts[1]
      
      # build the final short label
      lab <- if (order == "genus_first") paste0(gen_n, "_", epithet) else paste0(epithet, "_", gen_n)
      # collapse any trailing "_unclassified_unclassified"
      lab <- gsub("(_unclassified)+$", "_unclassified", lab, ignore.case = TRUE)
      lab
    }
    
    tax$Species <- mapply(simple_species_label, tax$Genus, tax$Species, USE.NAMES = FALSE)
  }
  
  # phyloseq
  OTU <- phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX <- phyloseq::tax_table(as.matrix(tax))
  SAM <- phyloseq::sample_data(as.data.frame(metadata))
  ps <- if (!is.null(tree_file)) {
    TREE <- phyloseq::read_tree(tree_file)
    phyloseq::phyloseq(OTU, TAX, SAM, TREE)
  } else {
    phyloseq::phyloseq(OTU, TAX, SAM)
  }
  
  # output dirs
  paths <- list(
    abundance_dir   = file.path(output_dir, "abundance_plot"),
    alpha_dir       = file.path(output_dir, "alpha_diversity"),
    beta_dir        = file.path(output_dir, "beta_diversity"),
    rarefaction_dir = file.path(output_dir, "rarefaction_plot")
  )
  dir_create_all(paths)
  
  # ---------- abundance ----------
  do_abundance <- function() {
    ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x) * 100)
    wb_full <- openxlsx::createWorkbook()
    wb_top  <- openxlsx::createWorkbook()
    
    rank_floor <- function(rank) {
      switch(rank,
             "Species" = 2.5,      #
             "Genus"   = 2.5,    # 2.5%
             "Phylum"  = 1,      # keep your lighter floor for higher ranks if you want
             "Class"   = 2.5,
             5)                  # default for others
    }
    
    grp_name <- pick_group_name_from_metadata(metadata, prefer = group_by)
    
    for (rank in tax_ranks) {
      if (!(rank %in% colnames(phyloseq::tax_table(ps_rel))) || rank == "Domain") next
      
      ps_rank <- phyloseq::tax_glom(ps_rel, taxrank = rank)
      df <- phyloseq::psmelt(ps_rank)
      df[[rank]] <- as.character(df[[rank]])
      df <- attach_group(df, metadata, grp_name)  # adds '.__GROUP__'
      
      # workbook: full (no lumping)
      df_full <- df %>%
        dplyr::group_by(Sample, .data[[rank]]) %>%
        dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
      colnames(df_full)[1] <- rank
      openxlsx::addWorksheet(wb_full, rank)
      openxlsx::writeData(wb_full, rank, df_full)
      
      # cutoff selection (auto, floored by rank_floor)
      floor_cut <- rank_floor(rank)
      med_grp <- df %>%
        dplyr::group_by(.data[[".__GROUP__"]], .data[[rank]]) %>%
        dplyr::summarise(med = stats::median(Abundance), .groups = "drop")
      med_glob <- med_grp %>%
        dplyr::group_by(.data[[rank]]) %>%
        dplyr::summarise(center = stats::median(med), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(center))
      auto_cut <- if (nrow(med_glob) > 0) {
        med_glob$cum <- cumsum(med_glob$center) / sum(med_glob$center)
        k <- min(which(med_glob$cum >= 0.95), 20)
        if (k < nrow(med_glob)) med_glob$center[k] else 0
      } else 0
      cutoff <- max(auto_cut, floor_cut)
      other_label <- paste0("< ", fmt_pct(cutoff), "% abundance")
      
      keep_tbl <- df %>%
        dplyr::group_by(.data[[".__GROUP__"]], .data[[rank]]) %>%
        dplyr::summarise(median = stats::median(Abundance), .groups = "drop") %>%
        dplyr::mutate(keep = median > cutoff)
      
      df_lumped <- df %>%
        dplyr::left_join(keep_tbl[, c(".__GROUP__", rank, "keep")],
                         by = c(".__GROUP__", rank)) %>%
        dplyr::mutate(!!rlang::sym(rank) := ifelse(keep, .data[[rank]], other_label))
      
      df_grouped <- df_lumped %>%
        dplyr::group_by(Sample, .data[[".__GROUP__"]], .data[[rank]]) %>%
        dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
        dplyr::group_by(Sample, .data[[".__GROUP__"]]) %>%
        dplyr::mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
        dplyr::ungroup()
      
      # workbook: top (with lumped)
      df_top <- df_grouped %>%
        dplyr::select(Sample, .data[[rank]], Abundance) %>%
        tidyr::pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
      colnames(df_top)[1] <- rank
      openxlsx::addWorksheet(wb_top, rank)
      openxlsx::writeData(wb_top, rank, df_top)
      
      # plot
      levs <- c(sort(setdiff(unique(df_grouped[[rank]]), other_label)), other_label)
      df_grouped[[rank]] <- factor(df_grouped[[rank]], levels = levs)
      n_levels <- length(levs)
      pal <- viridisLite::viridis(n_levels, option = "D")
      
      p <- ggplot2::ggplot(df_grouped, ggplot2::aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::facet_wrap(~ `.__GROUP__`, scales = "free_x", nrow = 1) +
        ggplot2::labs(
          title = paste0("Relative Abundance - ", rank, "  (lumped at ", fmt_pct(cutoff), "%)"),
          y = "Relative Abundance (%)", x = "Sample_ID"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          axis.text.x.bottom = ggplot2::element_text(angle = -90, vjust = 0.5, hjust = 0),
          legend.position = "right"
        ) +
        ggplot2::scale_fill_manual(values = pal[seq_len(n_levels)])
      
      ggplot2::ggsave(
        file.path(paths$abundance_dir, paste0("relative_abundance_barplot_", tolower(rank), ".tiff")),
        plot = p, dpi = dpi, width = width + 2, height = height, compression = "lzw"
      )
    }
    
    openxlsx::saveWorkbook(wb_full, file.path(paths$abundance_dir, "relative_abundance_full.xlsx"), overwrite = TRUE)
    openxlsx::saveWorkbook(wb_top,  file.path(paths$abundance_dir, "relative_abundance_top.xlsx"),  overwrite = TRUE)
    message("✔️ Relative abundance saved to ", paths$abundance_dir)
  }
  
  # ---------- rarefaction ----------
  # Full-page rarefaction with a bottom legend (fixed device size)
  do_rarefaction <- function(width_in = 14, height_in = 10,
                             legend_at = c("bottom", "none"),
                             dots = TRUE) {
    legend_at <- match.arg(legend_at)
    
    otu_mat <- as(phyloseq::otu_table(ps), "matrix")
    if (phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    n_samples <- nrow(otu_mat)
    cols <- viridisLite::viridis(n_samples, option = "D")
    
    # legend layout helpers (bottom legend)
    leg_ncol_for <- function(n) {
      # aim for ~8–10 items per row
      pmax(1L, ceiling(n / 9))
    }
    leg_cex_for <- function(n) {
      if (n <= 24) 0.9 else if (n <= 48) 0.8 else if (n <= 80) 0.75 else 0.7
    }
    
    plot_rares <- function(mat, fname, title) {
      n <- nrow(mat)
      ncol_legend <- leg_ncol_for(n)
      cex_leg <- leg_cex_for(n)
      
      # Fixed device size, no autoscaling
      dev_width  <- width_in
      dev_height <- height_in
      lwd_main   <- 3
      cex_points <- 1.0
      
      grDevices::tiff(file.path(paths$rarefaction_dir, fname),
                      width = dev_width, height = dev_height, units = "in",
                      res = dpi, compression = "lzw")
      on.exit(grDevices::dev.off(), add = TRUE)
      
      if (legend_at == "bottom") {
        layout(matrix(c(1, 2), nrow = 2, byrow = TRUE),
               heights = c(0.82, 0.18))  # adjust for legend space
      } else {
        layout(matrix(1))
      }
      
      ## panel 1: curves (full width)
      par(mar = c(5, 5, 4, 2), xpd = NA)
      out <- vegan::rarecurve(mat, step = 20, label = FALSE, col = cols,
                              lwd = lwd_main,
                              xlab = "Sequencing Depth",
                              ylab = "Expected Species Richness",
                              main = title)
      
      if (isTRUE(dots)) {
        for (i in seq_along(out)) {
          richness <- out[[i]]
          depths <- attr(out[[i]], "Subsample")
          if (length(depths) >= 8) {
            idx <- round(seq(1, length(depths), length.out = 10))
            points(depths[idx], richness[idx], col = cols[i], pch = 16, cex = cex_points)
          }
        }
      }
      
      ## panel 2: legend (bottom, spans full width)
      if (legend_at == "bottom") {
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("center",
               legend = rownames(mat), col = cols,
               lty = 1, lwd = 2, pch = 16,
               ncol = ncol_legend, cex = cex_leg, bty = "n",
               x.intersp = 0.8, y.intersp = 0.9, seg.len = 1.2)
      }
    }
    
    # Unrarefied
    plot_rares(otu_mat, "rarefaction_curve.tiff", "Rarefaction Curves (Unrarefied)")
    
    # Rarefied (if feasible)
    min_depth <- suppressWarnings(min(rowSums(otu_mat)))
    can_rarefy <- is.finite(min_depth) && min_depth >= 50
    if (can_rarefy) {
      rarefied <- try({
        set.seed(123)
        vegan::rrarefy(otu_mat, sample = min_depth)
      }, silent = TRUE)
      
      if (!inherits(rarefied, "try-error")) {
        plot_rares(
          rarefied,
          "rarefaction_curve_rarefied.tiff",
          sprintf("Rarefaction Curves (Rarefied to %d reads)", min_depth)
        )
      } else {
        message("Rarefied plot skipped (rarefaction failed).")
      }
    } else {
      message("Rarefied plot skipped (min depth too low).")
    }
    
    message("✔️ Rarefaction saved to ", paths$rarefaction_dir)
  }
  
  
  # ---------- alpha (plots + single combined workbook; NULL-safe) ----------
  do_alpha <- function() {
    if (phyloseq::nsamples(ps) < 3) {
      message("Fewer than 3 samples. Skipping alpha diversity.")
      return(invisible(NULL))
    }
    
    # ---- Tunables ----
    iqr_k          <- 3              # extreme outliers = 3×IQR per group
    max_brk        <- 12             # cap significance brackets per panel
    base_palette   <- "Set2"         # "Set2" or "viridis"
    chao_group_col <- NULL           # e.g., "group"; NULL uses default Group
    side_widths    <- c(1, 0.02, 0.40)  # [plots, spacer, side note]
    wrap_chars     <- 58
    
    # Fallbacks for globals
    dpi    <- if (exists("dpi", inherits = TRUE))    get("dpi", inherits = TRUE)    else 300
    width  <- if (exists("width", inherits = TRUE))  get("width", inherits = TRUE)  else 8
    height <- if (exists("height", inherits = TRUE)) get("height", inherits = TRUE) else 5
    
    suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
      library(stringr); library(openxlsx); library(rlang)
    })
    
    grp_name <- pick_group_name_from_metadata(metadata, prefer = group_by)
    
    # ---- Outputs (figures; one workbook with everything) ----
    base_shob <- file.path(paths$alpha_dir, "alpha_shannon_observed")
    base_sich <- file.path(paths$alpha_dir, "alpha_simpson_chao1")
    base_chao <- file.path(paths$alpha_dir, "alpha_chao1_overall_box")
    xlsx_combined <- file.path(paths$alpha_dir, "alpha_diversity_combined.xlsx")
    
    measures <- c("Observed","Shannon","Simpson","Chao1")
    
    # ---- Data prep ----
    alpha_df <- phyloseq::estimate_richness(ps, measures = measures)
    alpha_df <- cbind(phyloseq::sample_data(ps), alpha_df) %>%
      tibble::rownames_to_column("Sample")
    
    alpha_long <- alpha_df |>
      pivot_longer(cols = all_of(measures), names_to = "Measure", values_to = "Value")
    alpha_long <- attach_group(alpha_long, metadata, grp_name)
    alpha_long$Group <- factor(alpha_long$`.__GROUP__`)
    
    sample_group <- alpha_long |> distinct(Sample, Group)
    
    # ---- Helpers ----
    trim_extreme <- function(d, value_col = "Value", k = iqr_k) {
      d %>%
        group_by(Group) %>%
        mutate(
          Q1 = stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE),
          Q3 = stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE),
          IQRv = Q3 - Q1, lo = Q1 - k * IQRv, hi = Q3 + k * IQRv
        ) %>%
        ungroup() %>%
        filter(.data[[value_col]] >= lo, .data[[value_col]] <= hi) %>%
        select(-Q1, -Q3, -IQRv, -lo, -hi)
    }
    
    bracket_df_sig <- function(dfm, alpha_level) {
      if (!requireNamespace("rstatix", quietly = TRUE)) return(NULL)
      if (length(unique(dfm$Group)) < 2) return(NULL)
      
      stat <- rstatix::pairwise_wilcox_test(dfm, Value ~ Group, p.adjust.method = "BH")
      p_use <- if ("p.adj" %in% names(stat)) stat$p.adj else stat$p
      stat  <- stat[order(p_use), , drop = FALSE]
      stat  <- stat[p_use <= alpha_level, , drop = FALSE]
      if (!nrow(stat)) return(NULL)
      if (nrow(stat) > max_brk) stat <- stat[seq_len(max_brk), , drop = FALSE]
      
      stat <- rstatix::add_x_position(stat, x = "Group")
      
      rng <- range(dfm$Value, na.rm = TRUE); span <- diff(rng); if (!is.finite(span) || span <= 0) span <- 1
      y0  <- max(dfm$Value, na.rm = TRUE) + 0.06 * span
      step <- 0.07 * span
      
      stat %>%
        arrange(group1, group2, across(any_of(c("p.adj", "p")))) %>%
        mutate(
          y.position = y0 + step * (row_number() - 1),
          p_show = if_else(!is.na(.data[["p.adj"]]), .data[["p.adj"]], .data[["p"]]),
          label  = case_when(p_show < 0.001 ~ "p < 0.001",
                             TRUE ~ sprintf("p = %.3f", p_show))
        )
    }
    
    dynamic_bounds <- function(v, y_extra = numeric(0)) {
      v <- stats::na.omit(as.numeric(v))
      allv <- c(v, stats::na.omit(y_extra))
      if (!length(allv)) return(c(NA_real_, NA_real_))
      lo <- min(allv); hi <- max(allv)
      if (lo == hi) { lo <- lo * 0.98; hi <- hi * 1.02 }
      rng <- max(hi - lo, .Machine$double.eps)
      c(lo - 0.06*rng, hi + 0.08*rng)
    }
    
    scale_cols <- function(n) {
      if (base_palette == "Set2") ggplot2::scale_color_brewer(palette = "Set2")
      else ggplot2::scale_color_manual(values = viridisLite::viridis(n))
    }
    scale_fillz <- function(n) {
      if (base_palette == "Set2") ggplot2::scale_fill_brewer(palette = "Set2")
      else ggplot2::scale_fill_manual(values = viridisLite::viridis(n))
    }
    
    side_note_grob <- function(alpha_level, iqr_k, wrap_chars_local) {
      txt <- paste0(
        "Only significant pairwise differences are shown (Wilcoxon, BH-adjusted, α = ",
        format(alpha_level, nsmall = 3), "). ",
        "Extreme outliers removed per group using ", iqr_k, "×IQR. ",
        "Boxes use translucent fills; non-significant comparisons are omitted."
      )
      ggpubr::text_grob(
        label = stringr::str_wrap(txt, width = wrap_chars_local),
        x = 0, y = 1, hjust = 0, vjust = 1, size = 10
      )
    }
    
    make_panel <- function(df_long, measure_label, alpha_level) {
      dfm <- df_long %>%
        filter(Measure == measure_label) %>%
        select(Sample, Group, Value) %>%
        filter(is.finite(Value))
      dfm$Group <- droplevels(dfm$Group)
      if (nrow(dfm) < 3 || nlevels(dfm$Group) < 2) {
        return(list(
          plot = ggplot() + labs(title = measure_label) + theme_void(),
          trimmed_data = dfm[0,], sig_pairs = NULL
        ))
      }
      
      dfm <- trim_extreme(dfm, "Value", k = iqr_k)
      
      stat <- bracket_df_sig(dfm, alpha_level)
      y_extra <- if (!is.null(stat)) stat$y.position else numeric(0)
      yl <- dynamic_bounds(dfm$Value, y_extra)
      
      n_groups <- nlevels(dfm$Group)
      
      p <- ggplot(dfm, aes(x = Group, y = Value, fill = Group, color = Group)) +
        geom_boxplot(alpha = 0.35, width = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.85) +
        scale_cols(n_groups) + scale_fillz(n_groups) +
        scale_y_continuous(limits = yl, expand = expansion(mult = c(0.02, 0.05))) +
        labs(title = measure_label, x = NULL, y = "Value") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      if (!is.null(stat)) {
        p <- p + ggpubr::stat_pvalue_manual(
          stat,
          label = "label", xmin = "group1", xmax = "group2",
          y.position = "y.position",
          tip.length = 0.01, bracket.size = 0.5, size = 3.2,
          hide.ns = TRUE
        )
      }
      list(plot = p, trimmed_data = dfm, sig_pairs = stat)
    }
    
    # Safe helper to tag sig tables; returns NULL or a data frame
    safe_mutate <- function(df, measure, alpha_val) {
      if (!is.null(df) && nrow(df) > 0) {
        dplyr::mutate(df, Measure = measure, alpha = alpha_val)
      } else NULL
    }
    
    # Rendering figures and collect significant pairs for a given alpha
    save_alpha_figs <- function(alpha_level, suffix_tag) {
      res_shan <- make_panel(alpha_long, "Shannon", alpha_level)
      res_obsv <- make_panel(alpha_long, "Observed", alpha_level)
      res_simp <- make_panel(alpha_long, "Simpson", alpha_level)
      res_chao <- make_panel(alpha_long, "Chao1",   alpha_level)
      
      fig_core1 <- ggpubr::ggarrange(res_shan$plot, res_obsv$plot, ncol = 2, labels = c("A", "B"))
      fig_core2 <- ggpubr::ggarrange(res_simp$plot, res_chao$plot, ncol = 2, labels = c("A", "B"))
      
      fig1 <- ggpubr::ggarrange(
        fig_core1,
        ggpubr::as_ggplot(grid::rectGrob(gp = grid::gpar(col = NA))),
        ggpubr::as_ggplot(side_note_grob(alpha_level, iqr_k, wrap_chars)),
        ncol = 3, widths = side_widths
      )
      fig2 <- ggpubr::ggarrange(
        fig_core2,
        ggpubr::as_ggplot(grid::rectGrob(gp = grid::gpar(col = NA))),
        ggpubr::as_ggplot(side_note_grob(alpha_level, iqr_k, wrap_chars)),
        ncol = 3, widths = side_widths
      )
      
      ggsave(paste0(base_shob, "_", suffix_tag, ".tiff"),
             plot = fig1, dpi = dpi, width = width + 2, height = height,
             device = "tiff", compression = "lzw")
      ggsave(paste0(base_sich, "_", suffix_tag, ".tiff"),
             plot = fig2, dpi = dpi, width = width + 2, height = height,
             device = "tiff", compression = "lzw")
      
    
      if (!is.null(chao_group_col)) {
        mdf <- metadata %>% as.data.frame() %>% tibble::rownames_to_column("Sample")
        if (chao_group_col %in% names(mdf)) {
          chao_df <- alpha_df %>%
            select(Sample, Value = Chao1) %>%
            left_join(mdf[, c("Sample", chao_group_col)], by = "Sample") %>%
            rename(Group = !!sym(chao_group_col))
        } else {
          warning(sprintf("Column '%s' not found in metadata. Falling back to default Group for Chao1.", chao_group_col))
          chao_df <- alpha_long %>% filter(Measure == "Chao1") %>% select(Sample, Group, Value)
        }
      } else {
        chao_df <- alpha_long %>% filter(Measure == "Chao1") %>% select(Sample, Group, Value)
      }
      
      chao_df <- chao_df %>%
        filter(is.finite(Value)) %>%
        mutate(Group = droplevels(factor(Group))) %>%
        trim_extreme("Value", k = iqr_k)
      
      chao_stat <- bracket_df_sig(chao_df, alpha_level)
      chao_y_extra <- if (!is.null(chao_stat)) chao_stat$y.position else numeric(0)
      yl_ch <- dynamic_bounds(chao_df$Value, chao_y_extra)
      n_groups_ch <- nlevels(chao_df$Group)
      
      p_chao_overall <- ggplot(chao_df, aes(x = Group, y = Value, fill = Group, color = Group)) +
        geom_boxplot(alpha = 0.35, width = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.85) +
        scale_cols(n_groups_ch) + scale_fillz(n_groups_ch) +
        scale_y_continuous(limits = yl_ch, expand = expansion(mult = c(0.02, 0.05))) +
        labs(x = NULL, y = "Chao1", title = "Chao1") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      if (!is.null(chao_stat)) {
        p_chao_overall <- p_chao_overall + ggpubr::stat_pvalue_manual(
          chao_stat, label = "label",
          xmin = "group1", xmax = "group2",
          y.position = "y.position",
          tip.length = 0.01, bracket.size = 0.5, size = 3.2,
          hide.ns = TRUE
        )
      }
      
      p_chao_with_note <- ggpubr::ggarrange(
        p_chao_overall,
        ggpubr::as_ggplot(grid::rectGrob(gp = grid::gpar(col = NA))),
        ggpubr::as_ggplot(side_note_grob(alpha_level, iqr_k, wrap_chars)),
        ncol = 3, widths = side_widths
      )
      
      ggsave(paste0(base_chao, "_", suffix_tag, ".tiff"),
             plot = p_chao_with_note, dpi = dpi, width = width + 1, height = height,
             device = "tiff", compression = "lzw")
      
      # Returning tidy significant pairs, NULL-safe
      dplyr::bind_rows(
        safe_mutate(res_shan$sig_pairs, "Shannon",  alpha_level),
        safe_mutate(res_obsv$sig_pairs, "Observed", alpha_level),
        safe_mutate(res_simp$sig_pairs, "Simpson",  alpha_level),
        safe_mutate(res_chao$sig_pairs, "Chao1",    alpha_level)
      )
    }
    
    # ---- Rendering both sets: p < 0.001 and p < 0.05, and collect sig pairs ----
    sig_p001 <- save_alpha_figs(alpha_level = 0.001, suffix_tag = "p001")
    sig_p05  <- save_alpha_figs(alpha_level = 0.050, suffix_tag = "p05")
    
    # ---- Single combined workbook (front sheet + everything) ----
    wb <- openxlsx::createWorkbook()
    
    # 1) Front sheet: full index set (uses microbiome::alpha if available)
    index_vec <- c(
      "observed","chao1",
      "diversity_inverse_simpson","diversity_gini_simpson","diversity_shannon",
      "diversity_fisher","diversity_coverage",
      "evenness_camargo","evenness_pielou","evenness_simpson","evenness_evar","evenness_bulla",
      "dominance_dbp","dominance_dmn","dominance_absolute","dominance_relative",
      "dominance_simpson","dominance_core_abundance","dominance_gini",
      "rarity_log_modulo_skewness","rarity_low_abundance","rarity_rare_abundance"
    )
    
    if (requireNamespace("microbiome", quietly = TRUE)) {
      alpha_front <- microbiome::alpha(ps, index = index_vec) %>%
        tibble::rownames_to_column("Sample")
    } else {
      message("Package 'microbiome' not found. Computing a partial front sheet; some indices will be NA.")
      basics <- phyloseq::estimate_richness(ps, measures = c("Observed","Shannon","Simpson","Chao1")) %>%
        tibble::rownames_to_column("Sample")
      basics <- basics %>%
        mutate(
          diversity_inverse_simpson = 1 / Simpson,
          diversity_gini_simpson   = 1 - Simpson,
          diversity_shannon        = Shannon,
          observed                 = Observed,
          chao1                    = Chao1
        )
      if (requireNamespace("vegan", quietly = TRUE)) {
        otu_mat <- as(phyloseq::otu_table(ps), "matrix")
        if (phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
        fisher <- apply(otu_mat, 1, function(x) tryCatch(vegan::fisher.alpha(x), error = function(e) NA_real_))
        basics$diversity_fisher <- as.numeric(fisher[match(basics$Sample, names(fisher))])
      }
      alpha_front <- basics %>%
        transmute(
          Sample,
          observed,
          chao1,
          diversity_inverse_simpson,
          diversity_gini_simpson,
          diversity_shannon,
          diversity_fisher = dplyr::coalesce(.data$diversity_fisher, NA_real_),
          diversity_coverage = NA_real_,
          evenness_camargo = NA_real_,
          evenness_pielou = NA_real_,
          evenness_simpson = NA_real_,
          evenness_evar = NA_real_,
          evenness_bulla = NA_real_,
          dominance_dbp = NA_real_,
          dominance_dmn = NA_real_,
          dominance_absolute = NA_real_,
          dominance_relative = NA_real_,
          dominance_simpson = NA_real_,
          dominance_core_abundance = NA_real_,
          dominance_gini = NA_real_,
          rarity_log_modulo_skewness = NA_real_,
          rarity_low_abundance = NA_real_,
          rarity_rare_abundance = NA_real_
        )
    }
    
    front_cols <- c(
      "Sample",
      "observed","chao1",
      "diversity_inverse_simpson","diversity_gini_simpson","diversity_shannon",
      "diversity_fisher","diversity_coverage",
      "evenness_camargo","evenness_pielou","evenness_simpson","evenness_evar","evenness_bulla",
      "dominance_dbp","dominance_dmn","dominance_absolute","dominance_relative",
      "dominance_simpson","dominance_core_abundance","dominance_gini",
      "rarity_log_modulo_skewness","rarity_low_abundance","rarity_rare_abundance"
    )
    alpha_front <- alpha_front[, intersect(front_cols, names(alpha_front)), drop = FALSE]
    alpha_front <- alpha_front[, front_cols[front_cols %in% names(alpha_front)], drop = FALSE]
    
    openxlsx::addWorksheet(wb, "Alpha_All_Indices")
    openxlsx::writeData(wb, "Alpha_All_Indices", alpha_front)
    
    # 2) Data sheets
    openxlsx::addWorksheet(wb, "Alpha_Long"); openxlsx::writeData(wb, "Alpha_Long", alpha_long)
    alpha_wide <- alpha_df %>% select(Sample, everything())
    openxlsx::addWorksheet(wb, "Alpha_Wide"); openxlsx::writeData(wb, "Alpha_Wide", alpha_wide)
    
    # 3) BH pairwise p-value matrices per measure (trimmed data, consistent with plots)
    for (measure in measures) {
      dfm <- alpha_df %>%
        select(Sample, Value = all_of(measure)) %>%
        left_join(sample_group, by = "Sample") %>%
        filter(is.finite(Value))
      if (n_distinct(dfm$Group) >= 2 && nrow(dfm) >= 3) {
        dfm_trim <- trim_extreme(dfm, "Value", k = iqr_k)
        res <- try(stats::pairwise.wilcox.test(dfm_trim$Value, dfm_trim$Group, p.adjust.method = "BH"),
                   silent = TRUE)
        if (!inherits(res, "try-error")) {
          pmat <- as.data.frame.matrix(res$p.value)
          sheet_name <- paste0(measure, "_Wilcoxon_BH")
          openxlsx::addWorksheet(wb, sheet_name)
          openxlsx::writeData(wb, sheet_name, pmat, rowNames = TRUE)
        }
      }
    }
    
    # 4) Significant pairs tables for both alphas (tidy, NULL-safe)
    tidy_write <- function(df, sheet) {
      openxlsx::addWorksheet(wb, sheet)
      if (!is.null(df) && nrow(df)) {
        cols <- intersect(c("Measure","alpha","group1","group2","p","p.adj","p_show","label"), names(df))
        openxlsx::writeData(wb, sheet, df[, cols, drop = FALSE])
      } else {
        openxlsx::writeData(wb, sheet, data.frame(Note = "No significant pairs at this alpha."))
      }
    }
    tidy_write(sig_p001, "Significant_Pairs_p001")
    tidy_write(sig_p05,  "Significant_Pairs_p05")
    
    # Save
    openxlsx::saveWorkbook(wb, xlsx_combined, overwrite = TRUE)
    
    message("✔️ Alpha diversity figures (p<0.001 and p<0.05) and combined workbook saved to ", paths$alpha_dir)
  }
  
  #---------------------------beta-------------------------------
  do_beta <- function() {
    if (phyloseq::nsamples(ps) < 3) {
      message("⚠️ Fewer than 3 samples. Skipping beta diversity.")
      return(invisible(NULL))
    }
    has_tree <- !is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))
    dists <- if (has_tree) c("bray","unifrac","wunifrac") else { message("No tree. Using Bray only."); "bray" }
    
    grp_name <- pick_group_name_from_metadata(metadata, prefer = group_by)
    md <- metadata %>% tibble::rownames_to_column("Sample") %>% dplyr::mutate(Sample = as.character(Sample))
    
    for (dist in dists) {
      # NMDS
      nmds <- phyloseq::ordinate(ps, method = "NMDS", distance = dist)
      nmds_df <- phyloseq::plot_ordination(ps, nmds, justDF = TRUE) %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::mutate(Sample = as.character(Sample)) %>%
        dplyr::left_join(md, by = "Sample") %>%
        attach_group(metadata, grp_name)
      
      p1 <- ggplot2::ggplot(nmds_df, ggplot2::aes(x = NMDS1, y = NMDS2, color = `.__GROUP__`)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::stat_ellipse(type = "norm", level = 0.95) +
        ggplot2::theme_classic() +
        ggplot2::labs(title = paste("NMDS with Ellipses -", dist), x = "NMDS1", y = "NMDS2") +
        ggplot2::theme(legend.position = "bottom")
      ggplot2::ggsave(file.path(paths$beta_dir, paste0("nmds_ellipse_", dist, ".tiff")),
                      plot = p1, dpi = dpi, width = width, height = height, compression = "lzw")
      
      p2 <- ggplot2::ggplot(nmds_df, ggplot2::aes(x = NMDS1, y = NMDS2, color = `.__GROUP__`)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.7) +
        ggplot2::facet_wrap(~ `.__GROUP__`) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Faceted NMDS -", dist), x = "NMDS1", y = "NMDS2") +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(file.path(paths$beta_dir, paste0("nmds_facet_", dist, ".tiff")),
                      plot = p2, dpi = dpi, width = width + 2, height = height, compression = "lzw")
      
      # PCoA
      pcoa <- phyloseq::ordinate(ps, method = "PCoA", distance = dist)
      pcoa_df <- phyloseq::plot_ordination(ps, pcoa, justDF = TRUE) %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::mutate(Sample = as.character(Sample)) %>%
        dplyr::left_join(md, by = "Sample") %>%
        attach_group(metadata, grp_name)
      
      p3 <- ggplot2::ggplot(pcoa_df, ggplot2::aes(x = Axis.1, y = Axis.2, color = `.__GROUP__`)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::stat_ellipse(type = "norm", level = 0.95) +
        ggplot2::theme_light() +
        ggplot2::labs(title = paste("PCoA -", dist), x = "PCoA 1", y = "PCoA 2") +
        ggplot2::theme(legend.position = "bottom")
      ggplot2::ggsave(file.path(paths$beta_dir, paste0("pcoa_", dist, ".tiff")),
                      plot = p3, dpi = dpi, width = width, height = height, compression = "lzw")
    }
    
    message("✔️ Beta diversity saved to ", paths$beta_dir)
  }
  
  # ---------- run selected ----------
  if ("abundance"   %in% plot_types) do_abundance()
  if ("rarefaction" %in% plot_types) do_rarefaction()
  if ("alpha"       %in% plot_types) do_alpha()
  if ("beta"        %in% plot_types) do_beta()
  
  message("✔️ All selected plots and files saved to ", output_dir)
}
