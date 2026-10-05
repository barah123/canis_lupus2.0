# =============================================================================
# CanisLupus 2.0 - Microbiome Analysis Dashboard
# Author: Philip Appiah
# =============================================================================

# ── Package Installation (run once) ──────────────────────────────────────────
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "ggtree", "metacoder", "microbiome", "DESeq2"))
# install.packages(c("shiny","plotly","DT","vegan","ape","ggplot2","tidyverse",
#   "networkD3","heatmaply","patchwork","RColorBrewer","viridis","igraph",
#   "cowplot","ggdendro","dendextend","broom","cluster","reshape2","scales"))

# ── Libraries ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(phyloseq)
  library(microbiome)
  library(plotly)
  library(DT)
  library(vegan)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyverse)
  library(metacoder)
  library(networkD3)
  library(heatmaply)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(igraph)
  library(cowplot)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(scales)
})

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Remove singletons (taxa present in only 1 read across all samples)
remove_singletons <- function(ps) {
  prune_taxa(taxa_sums(ps) > 1, ps)
}

#' CLR transformation (robust for compositional data)
clr_transform <- function(ps) {
  microbiome::transform(ps, "clr")
}

#' Bray-Curtis normalised object (compositional)
compositional <- function(ps) {
  microbiome::transform(ps, "compositional")
}

#' Safe tax_glom: skips if rank missing
safe_tax_glom <- function(ps, rank) {
  if (rank %in% rank_names(ps)) tax_glom(ps, rank) else ps
}

#' Locate a bundled demo dataset file.
#'
#' The Data Upload tab states that demo data loads when nothing is uploaded, so
#' this resolves those files from the example dataset shipped with the app.
#' Errors name the missing file rather than letting read.csv fail obscurely.
demo_file <- function(filename, dir = "skin_data", roots = NULL) {
  if (is.null(roots)) roots <- c(".", "..", "../..")
  for (root in roots) {
    path <- file.path(root, dir, filename)
    if (file.exists(path)) return(path)
  }
  stop("Demo data file not found: ", file.path(dir, filename),
       ". Upload your own files, or restore the bundled ", dir, "/ directory.",
       call. = FALSE)
}

#' Distance metrics offered on the Beta Diversity tab.
#'
#' UniFrac is phylogenetic, so it is only meaningful with a real tree. Without
#' one the app holds a placeholder phylogeny, and UniFrac on that yields numbers
#' that look valid but encode nothing, so those options are withheld.
allowed_beta_distances <- function(tree_is_real) {
  all_metrics <- c("bray", "jaccard", "unifrac", "wunifrac", "euclidean")
  if (isTRUE(tree_is_real)) all_metrics
  else setdiff(all_metrics, c("unifrac", "wunifrac"))
}

#' Resolve the distance metric to compute with, or NULL to refuse.
#'
#' Returns NULL when a phylogenetic metric is requested without a real tree, so
#' callers can explain the refusal instead of silently returning a wrong number.
resolve_beta_distance <- function(selected, tree_is_real) {
  if (is.null(selected) || !nzchar(selected)) return("bray")
  if (selected %in% c("unifrac", "wunifrac") && !isTRUE(tree_is_real)) return(NULL)
  selected
}

#' Is this a complete, unreplicated block design?
#'
#' TRUE when every combination of group and block holds exactly one sample,
#' which is what a Friedman test requires. Repeated-measures data that fails
#' this check (unbalanced, or several samples per cell) must fall back to an
#' unblocked test.
is_complete_block_design <- function(group, block) {
  if (is.null(group) || is.null(block)) return(FALSE)
  if (length(group) != length(block)) return(FALSE)
  if (!length(group)) return(FALSE)
  tab <- table(as.factor(group), as.factor(block))
  length(tab) > 0 && all(tab == 1)
}

#' Decide which group-comparison test the alpha-diversity tab should run.
#'
#' Returns "friedman" only when a usable block variable describes a complete
#' block design; otherwise "kruskal". Running Kruskal-Wallis on repeated
#' measures treats correlated samples as independent draws and can materially
#' mis-state significance, so the blocked case is detected explicitly rather
#' than left to the user.
choose_alpha_test <- function(group, block = NULL, block_var = NULL, group_var = NULL) {
  if (is.null(block) || is.null(block_var) || !nzchar(block_var)) return("kruskal")
  if (!is.null(group_var) && identical(block_var, group_var)) return("kruskal")
  if (is_complete_block_design(group, block)) "friedman" else "kruskal"
}

#' Normalise a phyloseq object before beta-diversity distances.
#'
#' Sequencing depth varies between samples and distances such as Bray-Curtis
#' respond to that, so relative abundance is the default; "none" is for counts
#' that are already normalised or rarefied.
normalise_for_beta <- function(ps, mode = c("relative", "none")) {
  mode <- match.arg(mode)
  if (identical(mode, "none")) return(ps)
  transform_sample_counts(ps, function(x) if (sum(x) > 0) x / sum(x) else x)
}

# =============================================================================
# SHOTGUN / TAXONOMIC-PROFILE INPUT
#
# Shotgun metagenomic profilers (MetaPhlAn, Kraken2/Bracken and similar) emit a
# single table of lineage strings by sample, rather than the separate abundance
# and taxonomy tables that amplicon pipelines produce. These helpers turn such a
# table into the same two matrices the rest of the app already works with.
# =============================================================================

TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#' Guess the separator used in lineage strings.
#'
#' MetaPhlAn uses "|", QIIME-style and several Kraken exports use ";".
detect_lineage_separator <- function(lineages) {
  if (!length(lineages)) return("|")
  n_pipe <- sum(grepl("|", lineages, fixed = TRUE))
  n_semi <- sum(grepl(";", lineages, fixed = TRUE))
  if (n_pipe >= n_semi && n_pipe > 0) "|" else if (n_semi > 0) ";" else "|"
}

#' Strip rank prefixes such as "k__" or "s__" and tidy separators.
clean_taxon_name <- function(x) {
  x <- sub("^[dkpcofgst]__", "", trimws(x))
  x <- gsub("_", " ", x)
  x[!nzchar(x)] <- NA_character_
  x
}

#' Keep only the deepest (leaf) lineages in a hierarchical profile.
#'
#' MetaPhlAn-style tables repeat each clade at every rank, so the same reads
#' appear in a kingdom row, a phylum row, and so on. Summing those rows counts
#' the same organism many times over. A row is a leaf when no other row extends
#' it, which selects one non-overlapping set regardless of how deep individual
#' lineages go.
leaf_lineages <- function(lineages, sep = "|") {
  if (!length(lineages)) return(logical(0))
  parts <- strsplit(lineages, sep, fixed = TRUE)
  ancestors <- new.env(hash = TRUE, parent = emptyenv())
  for (pp in parts) {
    if (length(pp) < 2) next
    for (i in seq_len(length(pp) - 1)) {
      assign(paste(pp[seq_len(i)], collapse = sep), TRUE, envir = ancestors)
    }
  }
  !vapply(lineages, function(l) exists(l, envir = ancestors, inherits = FALSE),
          logical(1), USE.NAMES = FALSE)
}

#' Are these values already relative abundances rather than read counts?
#'
#' Profilers such as MetaPhlAn report percentages, so depth-based QC and
#' singleton removal do not apply to them.
looks_like_relative_abundance <- function(mat) {
  if (!length(mat)) return(FALSE)
  if (any(mat < 0, na.rm = TRUE)) return(FALSE)
  totals <- colSums(mat, na.rm = TRUE)
  totals <- totals[totals > 0]
  if (!length(totals)) return(FALSE)

  # Every column landing on 100 (or on 1) is decisive on its own: read counts
  # would not do that across a whole study.
  on_100 <- all(abs(totals - 100) < 0.5)
  on_1   <- all(abs(totals - 1) < 0.01)
  if (on_100 || on_1) return(TRUE)

  # Otherwise the totals alone are ambiguous, because MetaPhlAn leaves an
  # unclassified fraction out and a real profile often sums to the nineties.
  # Fractional values then settle it: read counts are whole numbers.
  fractional <- any(mat %% 1 != 0, na.rm = TRUE)
  if (!fractional) return(FALSE)

  as_proportion <- all(totals > 0.5 & totals <= 1.01)
  as_percentage <- all(totals > 50 & totals <= 100.5)
  as_proportion || as_percentage
}

#' Read a profiler export, skipping any banner lines above the header.
#'
#' MetaPhlAn writes a version banner such as "#mpa_v30_CHOCOPhlAn_201901" on the
#' first line, which would otherwise be read as the header and push the real
#' header down into the data. The header is the first line that actually
#' contains the separator; anything above it is a banner.
read_profile_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("The taxonomic profile is empty.", call. = FALSE)

  n_tab <- vapply(lines[seq_len(min(20, length(lines)))],
                  function(l) lengths(regmatches(l, gregexpr("\t", l))), integer(1))
  n_com <- vapply(lines[seq_len(min(20, length(lines)))],
                  function(l) lengths(regmatches(l, gregexpr(",", l))), integer(1))
  sep <- if (max(n_tab) >= max(n_com)) "\t" else ","
  counts <- if (identical(sep, "\t")) n_tab else n_com

  header_idx <- which(counts > 0)[1]
  if (is.na(header_idx)) {
    stop("Could not find a header row with more than one column in this file.",
         call. = FALSE)
  }

  df <- utils::read.delim(text = paste(lines[seq(header_idx, length(lines))], collapse = "\n"),
                          sep = sep, check.names = FALSE, comment.char = "",
                          stringsAsFactors = FALSE)
  # Some exports mark the header itself with a leading '#'.
  names(df)[1] <- sub("^#\\s*", "", names(df)[1])
  df
}

#' Parse a shotgun taxonomic profile into abundance and taxonomy matrices.
#'
#' Accepts a data frame whose first column (or row names) holds lineage strings
#' and whose remaining columns are samples. Returns the same pieces an amplicon
#' upload provides, so every downstream module is unchanged.
parse_taxonomic_profile <- function(df, lineage_col = NULL) {
  if (is.null(df) || !nrow(df)) stop("The taxonomic profile is empty.", call. = FALSE)

  if (is.null(lineage_col)) {
    first_is_lineage <- !is.numeric(df[[1]])
    lineages <- if (first_is_lineage) as.character(df[[1]]) else rownames(df)
    abund <- if (first_is_lineage) df[, -1, drop = FALSE] else df
  } else {
    lineages <- as.character(df[[lineage_col]])
    abund <- df[, setdiff(names(df), lineage_col), drop = FALSE]
  }

  if (is.null(lineages) || !length(lineages) || all(is.na(lineages))) {
    stop("Could not find lineage strings. Expected a first column such as ",
         "'clade_name' holding entries like 'k__Bacteria|p__Firmicutes'.", call. = FALSE)
  }

  numeric_cols <- vapply(abund, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    stop("No numeric sample columns found in the taxonomic profile.", call. = FALSE)
  }
  # Drop profiler bookkeeping columns such as NCBI taxid, which are numeric but
  # are not samples.
  drop_names <- grepl("tax(onomy)?_?id|ncbi|clade_taxid", names(abund), ignore.case = TRUE)
  abund <- abund[, numeric_cols & !drop_names, drop = FALSE]
  if (!ncol(abund)) stop("No sample columns remain in the taxonomic profile.", call. = FALSE)

  sep <- detect_lineage_separator(lineages)
  keep <- leaf_lineages(lineages, sep = sep)
  lineages <- lineages[keep]
  abund <- abund[keep, , drop = FALSE]
  if (!nrow(abund)) stop("No taxa remained after collapsing the lineage hierarchy.", call. = FALSE)

  parts <- strsplit(lineages, sep, fixed = TRUE)
  depth <- max(lengths(parts))
  n_ranks <- max(length(TAX_RANKS), depth)
  tax <- t(vapply(parts, function(pp) {
    out <- rep(NA_character_, n_ranks)
    pp <- clean_taxon_name(pp)
    if (length(pp)) out[seq_along(pp)] <- pp
    out
  }, character(n_ranks)))
  colnames(tax) <- if (n_ranks <= length(TAX_RANKS)) {
    TAX_RANKS[seq_len(n_ranks)]
  } else {
    c(TAX_RANKS, paste0("Rank", seq_len(n_ranks - length(TAX_RANKS))))
  }

  # Name each taxon by its deepest resolved rank, kept unique for phyloseq.
  leaf_name <- apply(tax, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r)) r[length(r)] else NA_character_
  })
  leaf_name[is.na(leaf_name)] <- "Unclassified"
  ids <- make.unique(as.character(leaf_name), sep = "_")

  otu <- as.matrix(abund)
  mode(otu) <- "numeric"
  otu[is.na(otu)] <- 0
  rownames(otu) <- ids
  rownames(tax) <- ids

  list(otu = otu, tax = tax,
       is_relative = looks_like_relative_abundance(otu),
       n_ranks = n_ranks, separator = sep)
}

#' Which alpha-diversity indices can be computed from this table?
#'
#' Chao1, ACE and Fisher are estimated from how many taxa are seen exactly once
#' or twice, so they need integer read counts. Relative-abundance profiles, as
#' produced by shotgun profilers, carry no such counts; Shannon, Simpson and
#' inverse Simpson are defined on proportions and remain valid.
COUNT_ONLY_MEASURES <- c("Chao1", "ACE", "Fisher")

available_alpha_measures <- function(is_relative) {
  all_m <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
  if (isTRUE(is_relative)) setdiff(all_m, COUNT_ONLY_MEASURES) else all_m
}

#' Alpha diversity for count or proportion data.
#'
#' phyloseq::estimate_richness refuses non-integer input outright, so proportion
#' data is routed through vegan directly for the indices that remain valid.
compute_alpha <- function(ps, measures, is_relative = FALSE) {
  measures <- intersect(measures, available_alpha_measures(is_relative))
  if (!length(measures)) measures <- "Shannon"

  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)   # vegan wants samples as rows

  if (!isTRUE(is_relative)) {
    df <- suppressWarnings(phyloseq::estimate_richness(ps, measures = measures))
    rownames(df) <- phyloseq::sample_names(ps)
    return(df[, intersect(measures, colnames(df)), drop = FALSE])
  }

  out <- list()
  if ("Observed"   %in% measures) out$Observed   <- rowSums(mat > 0)
  if ("Shannon"    %in% measures) out$Shannon    <- vegan::diversity(mat, index = "shannon")
  if ("Simpson"    %in% measures) out$Simpson    <- vegan::diversity(mat, index = "simpson")
  if ("InvSimpson" %in% measures) out$InvSimpson <- vegan::diversity(mat, index = "invsimpson")
  df <- as.data.frame(out)
  rownames(df) <- rownames(mat)
  df
}

#' Resolve a requested taxonomic rank against what the dataset actually holds.
#'
#' Rank availability varies by source: the bundled amplicon taxonomy stops at
#' Genus, while shotgun profiles usually resolve to Species. Returns NULL when
#' the rank is absent so callers can explain, rather than pushing an
#' unresolvable column name into phyloseq.
safe_rank <- function(ps, rank) {
  ranks <- phyloseq::rank_names(ps)
  if (is.null(rank) || !length(rank) || is.na(rank[1]) || !nzchar(rank[1])) {
    return(if (length(ranks)) ranks[length(ranks)] else NULL)
  }
  if (rank[1] %in% ranks) rank[1] else NULL
}

#' Read a numeric input, falling back when it is missing or unset.
#'
#' Shiny inputs are briefly NULL before the client reports them, and a module
#' that divides by one would otherwise error on first paint.
input_num <- function(x, default) {
  if (is.null(x) || !length(x) || !is.finite(suppressWarnings(as.numeric(x)[1]))) default
  else as.numeric(x)[1]
}

#' Read a character input, falling back when it is missing or unset.
input_chr <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1])) default else as.character(x)[1]
}

# =============================================================================
# VISUAL DESIGN SYSTEM
#
# One palette and one plot theme, applied everywhere, so every figure in the app
# reads as part of the same publication-quality set.
# =============================================================================

# Okabe-Ito: the standard colourblind-safe qualitative palette for scientific
# figures, extended with further distinguishable hues for taxa-heavy plots.
canis_categorical <- c(
  "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#0072B2",
  "#D55E00", "#5D3A9B", "#117733", "#882255", "#F0E442",
  "#44AA99", "#DDCC77", "#AA4499", "#88CCEE", "#999933",
  "#661100", "#6699CC", "#332288", "#AA7744", "#BBBBBB"
)

# Colours assigned by function rather than by taste, so each one has exactly
# one job and the interface cannot drift into clutter. Neither pure white nor
# pure black is used: both strain the eye at length.
#
# The accents are Okabe-Ito colours, the standard colourblind-safe scientific
# set, so the chrome is accessible by construction and agrees with the charts.
# Anthocyanin - indigo violet shading to magenta
# Named for the pigment whose colour shifts with pH. Assertive and modern, and the cleanest separation from the chart palette.
canis_ui <- list(
  primary       = "#4B3E9E",
  secondary     = "#6B6480",
  accent        = "#A31E78",
  success       = "#16613F",
  danger        = "#BF2C30",
  ink           = "#1E1B30",
  muted         = "#63607A",
  surface       = "#FDFDFF",
  canvas        = "#F6F5FB",
  border        = "#E2DFEF",
  border_strong = "#938FA6"
)

# Dark-mode counterparts. Surfaces lift slightly off the background so cards
# stay legible as distinct planes without needing borders.
canis_dark <- list(
  ink        = "#E7E4F2",
  muted      = "#9C98B8",
  surface    = "#1A182C",
  canvas     = "#100E1F",
  border     = "#2C2942",
  primary    = "#A99BF5",
  on_primary = "#120E24"
)

#' Hex to an rgba() string, so derived colours track the palette instead of
#' being hardcoded and silently left behind when the theme changes.
hex_rgba <- function(hex, alpha) {
  v <- strtoi(substring(sub("#", "", hex), c(1, 3, 5), c(2, 4, 6)), 16L)
  sprintf("rgba(%d,%d,%d,%s)", v[1], v[2], v[3], alpha)
}

#' Categorical colours, recycled smoothly when a plot needs more than the set.
canis_colors <- function(n) {
  if (is.null(n) || is.na(n) || n < 1) return(canis_categorical[1])
  if (n <= length(canis_categorical)) return(canis_categorical[seq_len(n)])
  grDevices::colorRampPalette(canis_categorical)(n)
}

scale_fill_canis <- function(...) {
  ggplot2::discrete_scale("fill", palette = function(n) canis_colors(n), ...)
}
scale_color_canis <- function(...) {
  ggplot2::discrete_scale("colour", palette = function(n) canis_colors(n), ...)
}

#' Shared plot theme: light, uncluttered, and readable at publication size.
theme_canis <- function(base_size = 13, dark = FALSE) {
  ink     <- if (dark) canis_dark$ink     else canis_ui$ink
  muted   <- if (dark) canis_dark$muted   else canis_ui$muted
  surface <- if (dark) canis_dark$surface else canis_ui$surface
  border  <- if (dark) canis_dark$border  else canis_ui$border

  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text             = ggplot2::element_text(colour = ink),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.1,
                                               colour = ink,
                                               margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(colour = muted, size = base_size * 0.9),
      axis.title       = ggplot2::element_text(colour = muted, size = base_size * 0.9),
      axis.text        = ggplot2::element_text(colour = muted, size = base_size * 0.82),
      panel.grid.major = ggplot2::element_line(colour = border, linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = surface, colour = NA),
      plot.background  = ggplot2::element_rect(fill = surface, colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = ink,
                                               size = base_size * 0.88),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(colour = muted, size = base_size * 0.85),
      legend.text      = ggplot2::element_text(size = base_size * 0.82),
      plot.margin      = ggplot2::margin(12, 12, 12, 12)
    )
}

#' Hide per-sample axis labels once there are too many to read.
sample_axis_theme <- function(n_samples, max_labels = 25) {
  if (n_samples > max_labels) {
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  } else {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
  }
}

# Kept for backwards compatibility with any external code referencing it.
wolf_pal <- canis_categorical

# =============================================================================
# UI
# =============================================================================

ui <- bslib::page_navbar(
  title = div(
    class = "cl-brand",
    img(src = "canis_logo_128.png", class = "cl-brand-mark",
        alt = "CanisLupus: a wolf's head in a gold ring"),
    div(
      span(class = "cl-brand-name", "CanisLupus"),
      span(class = "cl-brand-sub", "Microbiome Analysis Platform")
    )
  ),
  window_title = "CanisLupus 2.0",
  id = "mainTabs",
  # Per the filling-layouts guidance: fill the tabs built around a single
  # dominant output, and let the output-dense tabs fall back to intrinsic
  # height so nothing is squeezed.
  fillable = c("Alpha diversity", "Beta diversity", "Composition",
               "Composition pie", "Phylogenetic tree", "Heat tree",
               "Abundance heatmap", "Sample clustering",
               "Differential abundance", "Transformations"),
  theme = bslib::bs_theme(
    version = 5,
    base_font    = bslib::font_google("Inter", local = FALSE),
    heading_font = bslib::font_google("Inter", local = FALSE),
    code_font    = bslib::font_google("JetBrains Mono", local = FALSE),
    primary   = canis_ui$primary,
    secondary = canis_ui$secondary,
    success   = canis_ui$success,
    info      = canis_ui$primary,
    warning   = canis_ui$accent,
    danger    = canis_ui$danger,
    bg        = canis_ui$canvas,
    fg        = canis_ui$ink,
    "border-radius"     = "0.6rem",
    "card-border-color" = canis_ui$border,
    # Bootstrap 5.3 colour modes: these drive the dark counterpart, so the
    # toggle switches the whole interface rather than just a few surfaces.
    "body-bg-dark"      = canis_dark$canvas,
    "body-color-dark"   = canis_dark$ink,
    "card-bg-dark"      = canis_dark$surface,
    "border-color-dark" = canis_dark$border
  ),

  header = tags$head(
    tags$link(rel = "icon", type = "image/png", href = "canis_logo_64.png"),
    tags$style(HTML(sprintf("
    :root {
      --cl-primary:%s; --cl-secondary:%s; --cl-accent:%s;
      --cl-ink:%s; --cl-muted:%s; --cl-surface:%s;
      --cl-canvas:%s; --cl-border:%s; --cl-danger:%s;
    }
    body { background:var(--cl-canvas); }

    .cl-brand { display:flex; align-items:center; gap:.7rem; margin-right:1.75rem; }
    /* The mark is already a circle with its own gold ring, so it needs no
       frame of its own; a soft ground just lifts it off the dark navbar. */
    .cl-brand-mark {
      width:42px; height:42px; border-radius:50%%;
      background:#fff; object-fit:cover; display:block;
      box-shadow:0 0 0 1px rgba(0,0,0,.10);
    }
    [data-bs-theme='dark'] .cl-brand-mark { box-shadow:0 0 0 1px rgba(255,255,255,.20); }
    /* Inherit the navbar's own text colour so the brand stays legible whether
       the navbar renders light or dark. */
    .cl-brand-name { display:block; font-weight:650; font-size:1rem; line-height:1.1;
                     color:inherit; }
    .cl-brand-sub  { display:block; font-size:.7rem; opacity:.65; font-weight:400;
                     color:inherit; }

    .navbar { box-shadow:0 1px 0 rgba(0,0,0,.06); border-bottom:2px solid var(--cl-primary); }
    .navbar .nav-link { font-size:.9rem; font-weight:500; }

    /* Cards separate by surface contrast and a soft two-layer shadow rather
       than a hard border: a tight contact shadow plus a wide ambient one. */
    .card {
      border:none; border-radius:8px; background:var(--cl-surface);
      box-shadow:0 3.2px 7.2px 0 rgba(0,0,0,.132), 0 .6px 1.8px 0 rgba(0,0,0,.108);
      margin-bottom:1.25rem;
    }
    .card-body { padding:1.5rem; }
    .card-header {
      background:var(--cl-surface); border-bottom:1px solid var(--cl-border);
      border-radius:8px 8px 0 0;
      font-weight:500; font-size:1.05rem; color:var(--cl-ink);
      padding:1rem 1.5rem;
      display:flex; align-items:center; justify-content:space-between;
    }

    .cl-section-title { font-size:.72rem; font-weight:700; letter-spacing:.09em;
      text-transform:uppercase; color:var(--cl-muted); margin:.15rem 0 .5rem; }

    .info-box {
      background:#FBF7F0; border:1px solid #EBE2D4; border-left:3px solid var(--cl-accent);
      border-radius:.45rem; padding:.6rem .7rem; margin-bottom:.85rem;
      font-size:.79rem; line-height:1.45; color:#41525C;
    }
    .info-box.cl-warn { background:#FBF0EF; border-color:#EFD6D4; border-left-color:var(--cl-danger); }

    .cl-hint { font-size:.75rem; color:var(--cl-muted); margin:-.35rem 0 .8rem; }
    .cl-help {
      display:inline-flex; align-items:center; justify-content:center;
      width:18px; height:18px; margin-left:.45rem; border-radius:50%%;
      border:1px solid var(--cl-border); background:#fff;
      color:var(--cl-muted); font-size:.7rem; font-weight:600;
      cursor:pointer; vertical-align:middle;
    }
    .cl-help:hover { border-color:var(--cl-primary); color:var(--cl-primary); }

    /* One control height and one radius throughout, so a dense sidebar still
       reads as a single system. */
    .form-label, label { font-size:.875rem; font-weight:500; color:var(--cl-ink); margin-bottom:.35rem; }
    .form-control, .form-select {
      font-size:.875rem; min-height:40px; border-radius:4px; border:1px solid %s;
    }
    .form-control:focus, .form-select:focus {
      border-color:var(--cl-primary); box-shadow:0 0 0 .18rem %s;
    }
    .btn {
      font-size:.875rem; font-weight:500; border-radius:4px;
      min-height:40px; padding:8px 16px;
      transition:box-shadow 200ms cubic-bezier(.4,1,.75,.9), background-color 150ms ease;
    }
    .btn-primary { background:var(--cl-primary); border-color:var(--cl-primary); }
    .btn-primary:hover { background:#74491F; border-color:#74491F; }
    .sidebar .form-group, .sidebar .shiny-input-container { margin-bottom:1rem; }

    .irs--shiny .irs-bar, .irs--shiny .irs-single { background:var(--cl-primary); border-color:var(--cl-primary); }
    .irs--shiny .irs-handle { border-color:var(--cl-primary); }

    /* verbatimTextOutput renders as <pre>; textOutput renders as <div> and must
       NOT pick up this panel styling, or every value_box gains a grey box. */
    pre, pre.shiny-text-output {
      background:#FBFCFD; border:1px solid var(--cl-border); border-radius:.45rem;
      padding:.75rem .85rem; font-size:.79rem; line-height:1.5; color:#2A3A44;
    }

    .cl-stat-row { display:flex; flex-wrap:wrap; gap:.75rem; margin-bottom:.25rem; }
    .cl-stat {
      flex:1 1 130px; background:var(--cl-surface); border:1px solid var(--cl-border);
      border-radius:.55rem; padding:.7rem .85rem;
    }
    .cl-stat-label { font-size:.68rem; text-transform:uppercase; letter-spacing:.07em;
      color:var(--cl-muted); font-weight:650; }
    .cl-stat-value { font-size:1.45rem; font-weight:680; color:var(--cl-ink); line-height:1.15; }
    .cl-stat-note  { font-size:.72rem; color:var(--cl-muted); }

    .cl-footer { border-top:1px solid var(--cl-border); margin-top:1.75rem;
      padding:1.1rem 0 1.6rem; color:var(--cl-muted); font-size:.79rem; text-align:center; }
    .cl-footer a { color:var(--cl-primary); text-decoration:none; margin:0 .55rem; font-weight:550; }
    .cl-footer a:hover { text-decoration:underline; }

    /* Bootstrap 5.3 stamps data-bs-theme on <html>; the custom properties have
       to follow or the hand-written components stay stuck in light mode. */
    [data-bs-theme='dark'] {
      --cl-ink:%s; --cl-muted:%s; --cl-surface:%s;
      --cl-canvas:%s; --cl-border:%s;
    }
    [data-bs-theme='dark'] .card {
      background:var(--cl-surface);
      box-shadow:0 3.2px 7.2px 0 rgba(0,0,0,.55), 0 .6px 1.8px 0 rgba(0,0,0,.4);
    }
    [data-bs-theme='dark'] .card-header { background:var(--cl-surface); }
    [data-bs-theme='dark'] .info-box {
      background:#1F1B12; border-color:#3A3323; color:#CBD2D8;
    }
    [data-bs-theme='dark'] pre, [data-bs-theme='dark'] pre.shiny-text-output {
      background:#151515; border-color:var(--cl-border); color:#D7DDE2;
    }
    [data-bs-theme='dark'] .form-control, [data-bs-theme='dark'] .form-select {
      background:#181818; border-color:var(--cl-border); color:var(--cl-ink);
    }
    [data-bs-theme='dark'] .cl-help { background:#181818; }
    /* On a dark canvas the primary lightens, so anything sitting ON it needs
       dark text rather than white. */
    [data-bs-theme='dark'] { --cl-primary:%s; }
    [data-bs-theme='dark'] .btn-primary {
      background:var(--cl-primary); border-color:var(--cl-primary); color:%s;
    }
    [data-bs-theme='dark'] .btn-primary:hover { filter:brightness(1.08); }
    [data-bs-theme='dark'] a { color:var(--cl-primary); }
  ",
  canis_ui$primary, canis_ui$secondary, canis_ui$accent, canis_ui$ink,
  canis_ui$muted, canis_ui$surface, canis_ui$canvas, canis_ui$border, canis_ui$danger,
  canis_ui$border_strong,
  hex_rgba(canis_ui$primary, ".15"),
  canis_dark$ink, canis_dark$muted, canis_dark$surface,
  canis_dark$canvas, canis_dark$border,
  canis_dark$primary, canis_dark$on_primary)))
  ),

  # ===========================================================================
  # DATA
  # ===========================================================================
  bslib::nav_panel(
    "Data",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 330, title = "Input & quality control",

        div(class = "cl-section-title", "Data type"),
        radioButtons(
          "dataType", NULL,
          choices = c("Amplicon (ASV / OTU table)" = "amplicon",
                      "Shotgun (taxonomic profile)" = "shotgun"),
          selected = "amplicon"
        ),

        conditionalPanel(
          "input.dataType == 'amplicon'",
          div(class = "info-box",
              "Upload an abundance table with taxa as rows, plus a matching taxonomy table. ",
              "Leave everything empty to explore the bundled example dataset."),
          fileInput("asvFile", "Abundance table (CSV)", accept = ".csv"),
          fileInput("taxFile", "Taxonomy table (CSV)", accept = ".csv")
        ),
        conditionalPanel(
          "input.dataType == 'shotgun'",
          div(class = "info-box",
              "Upload one merged profile from MetaPhlAn, Kraken2/Bracken or similar: ",
              "lineage strings in the first column, one column per sample. ",
              "Nested ranks are collapsed to their deepest level so nothing is counted twice."),
          fileInput("profileFile", "Taxonomic profile (TSV / CSV)",
                    accept = c(".tsv", ".txt", ".csv"))
        ),

        fileInput("metaFile", "Metadata (CSV)", accept = ".csv"),
        fileInput("treeFile", "Phylogenetic tree (optional)",
                  accept = c(".tree", ".tre", ".nwk", ".txt")),
        div(class = "cl-hint",
            "A tree unlocks UniFrac distances. Shotgun profiles usually have none."),

        tags$hr(),
        div(class = "cl-section-title", "Filtering"),
        checkboxInput("removeSingletons", "Remove singleton taxa", value = TRUE),
        numericInput("minReads", "Minimum reads per sample", value = 1000, min = 0, step = 100),
        div(class = "cl-hint",
            "Both apply to read counts. They are skipped automatically for ",
            "relative-abundance profiles."),

        tags$hr(),
        actionButton("update", "Load & process data", class = "btn-primary w-100"),
        div(class = "mt-2"),
        downloadButton("downloadPS", "Download summary", class = "btn-outline-secondary w-100")
      ),

      # value_box() carries its own theme-aware styling and Bootstrap contrast
      # handling, so the hand-rolled stat tiles are gone. fill = FALSE stops the
      # row claiming vertical space it does not need.
      uiOutput("datasetStatsUI"),
      bslib::layout_columns(
        col_widths = c(7, 5),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Sequencing depth per sample"),
          bslib::card_body(min_height = 300, plotlyOutput("readCountPlot"))
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Dataset summary"),
          verbatimTextOutput("dataSummary")
        )
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Metadata"),
        DTOutput("sampleTable")
      )
    )
  ),

  # ===========================================================================
  # DIVERSITY
  # ===========================================================================
  bslib::nav_menu(
    "Diversity",

    bslib::nav_panel(
      "Rarefaction curve",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Rarefaction",
          div(class = "info-box",
              "Shows whether sequencing depth was deep enough to capture the community. ",
              "A curve that flattens has saturated."),
          selectInput("regionRare", "Colour samples by", choices = NULL),
          numericInput("rareStep", "Step size", value = 500, min = 50, step = 50),
          div(class = "cl-hint",
              "Smaller steps give smoother curves but cost far more computation. ",
              "The step is raised automatically if it would be too fine for the ",
              "sequencing depth in this dataset."),
          numericInput("rareMax", "Maximum depth (0 = auto)", value = 0, min = 0),
          checkboxInput("showRareLine", "Mark a target depth", value = TRUE),
          numericInput("rareDepth", "Target depth", value = 10000, min = 100),
          actionButton("updateRare", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Rarefaction curves"),
          bslib::card_body(min_height = 380, plotlyOutput("rarefactionCurve")),
          bslib::card_body(fill = FALSE, class = "border-top",
                           verbatimTextOutput("rareSummary"))
        )
      )
    ),

    bslib::nav_panel(
      "Alpha diversity",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Alpha diversity",
          bslib::accordion(
            multiple = FALSE, open = "Indices",
            bslib::accordion_panel(
              "Indices",
              checkboxGroupInput(
                "alphaMeasure", NULL,
                choices  = c("Observed", "Chao1", "ACE", "Shannon",
                             "Simpson", "InvSimpson", "Fisher"),
                selected = c("Shannon", "Observed")
              )
            ),
            bslib::accordion_panel(
              "Grouping",
              selectInput("regionAlpha", "Group by", choices = NULL),
              selectInput("blockAlpha", "Block / subject variable (optional)", choices = NULL),
              checkboxInput("alphaStats", "Show group comparison test", value = TRUE)
            )
          ),
          actionButton("updateAlpha", "Update", class = "btn-primary w-100")
        ),
        # Three views of one question belong in one card as tabs, rather than
        # stacked cards each competing for vertical space.
        bslib::navset_card_underline(
          title = tagList(
            "Alpha diversity",
            bslib::popover(
              tags$span(class = "cl-help", "?"),
              title = "Alpha diversity",
              tags$p("Diversity within each sample. Indices weight richness and ",
                     "evenness differently, so reporting more than one is usual."),
              tags$p("Set a block/subject variable when the same subject appears at ",
                     "several levels of the grouping variable, such as repeat ",
                     "timepoints. A Friedman test then replaces Kruskal-Wallis, ",
                     "which assumes independent samples."),
              tags$p("Chao1, ACE and Fisher need integer read counts, so they are ",
                     "withheld for relative-abundance profiles.")
            )
          ),
          full_screen = TRUE,
          bslib::nav_panel("Per sample",  plotlyOutput("alphaDiv")),
          bslib::nav_panel("By group",    plotlyOutput("alphaDivBoxplot")),
          bslib::nav_panel("Group comparison", verbatimTextOutput("alphaStats"))
        )
      )
    ),

    bslib::nav_panel(
      "Beta diversity",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Beta diversity",
          # Accordions render flush inside a sidebar and are the documented way
          # to keep a control-dense sidebar short.
          bslib::accordion(
            multiple = FALSE, open = "Ordination",
            bslib::accordion_panel(
              "Ordination",
              selectInput("betaMethod", "Method",
                          choices = c("PCoA", "NMDS", "RDA", "CCA"), selected = "PCoA"),
              selectInput("betaDistance", "Distance",
                          choices = c("bray", "jaccard", "unifrac", "wunifrac", "euclidean"),
                          selected = "bray"),
              selectInput("betaNormalise", "Normalisation",
                          choices = c("Relative abundance" = "relative",
                                      "None (raw counts)" = "none"),
                          selected = "relative")
            ),
            bslib::accordion_panel(
              "Grouping",
              selectInput("regionBeta", "Colour / group by", choices = NULL),
              selectInput("blockBeta", "Block / subject variable (optional)", choices = NULL)
            ),
            bslib::accordion_panel(
              "Statistics",
              checkboxInput("betaEllipse", "Draw 95% ellipses", value = TRUE),
              checkboxInput("betaPermanova", "Run PERMANOVA", value = TRUE)
            )
          ),
          actionButton("updateBeta", "Update", class = "btn-primary w-100")
        ),
        bslib::navset_card_underline(
          title = tagList(
            "Beta diversity",
            bslib::popover(
              tags$span(class = "cl-help", "?"),
              title = "Beta diversity",
              tags$p("Dissimilarity between samples. PERMANOVA tests whether group ",
                     "centroids differ."),
              tags$p("Samples differ in sequencing depth and Bray-Curtis responds to ",
                     "that, so relative abundance is the safer normalisation."),
              tags$p("Setting a subject variable keeps permutations within each subject, ",
                     "which is what a repeated-measures design requires.")
            )
          ),
          full_screen = TRUE,
          bslib::nav_panel("Ordination", plotlyOutput("betaPlot")),
          bslib::nav_panel("NMDS",       plotlyOutput("nmdsPlot")),
          bslib::nav_panel("PERMANOVA",  verbatimTextOutput("permanovaResult"))
        )
      )
    )
  ),

  # ===========================================================================
  # TAXONOMY
  # ===========================================================================
  bslib::nav_menu(
    "Taxonomy",

    bslib::nav_panel(
      "Composition",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Composition",
          div(class = "info-box",
              "Relative abundance of the most abundant taxa, stacked per sample."),
          selectInput("taxLevelBar", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      selected = "Phylum"),
          selectInput("regionBar", "Facet by", choices = NULL),
          sliderInput("topNBar", "Number of taxa shown", min = 3, max = 30, value = 10),
          selectInput("barTransform", "Abundance",
                      choices = c("Relative" = "relative", "Absolute" = "absolute"),
                      selected = "relative"),
          actionButton("updateBar", "Update", class = "btn-primary w-100")
        ),
        bslib::navset_card_underline(
          title = "Community composition",
          full_screen = TRUE,
          bslib::nav_panel("Per sample",     plotlyOutput("taxaBarplot")),
          bslib::nav_panel("Mean abundance", plotlyOutput("relativeAbundancePlot"))
        )
      )
    ),

    bslib::nav_panel(
      "Composition pie",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Proportions",
          div(class = "info-box",
              "Overall share of each taxon, optionally split across groups."),
          selectInput("taxLevelPie", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      selected = "Phylum"),
          selectInput("regionPie", "Group by", choices = NULL),
          actionButton("updatePie", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Taxonomic proportions"),
          bslib::card_body(min_height = 420, plotlyOutput("pieChart"))
        )
      )
    ),

    bslib::nav_panel(
      "Core microbiome",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Core microbiome",
          div(class = "info-box",
              "Taxa found above a detection threshold in a given fraction of samples. ",
              "These are the consistent community members."),
          selectInput("taxLevelCore", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      selected = "Genus"),
          selectInput("regionCore", "Group by", choices = NULL),
          sliderInput("prevalenceCore", "Prevalence threshold (%)", min = 5, max = 100, value = 50),
          sliderInput("detectionCore", "Detection threshold (%)",
                      min = 0.001, max = 10, value = 0.1, step = 0.01),
          actionButton("updateCore", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Core microbiome"),
          bslib::card_body(min_height = 380, plotlyOutput("coreHeatmap")),
          bslib::card_body(fill = FALSE, class = "border-top",
                           verbatimTextOutput("coreSummary"))
        )
      )
    )
  ),

  # ===========================================================================
  # PHYLOGENETICS
  # ===========================================================================
  bslib::nav_menu(
    "Phylogenetics",

    bslib::nav_panel(
      "Phylogenetic tree",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Tree",
          div(class = "info-box", "Evolutionary relationships among the observed taxa."),
          uiOutput("treeProvenanceNote"),
          selectInput("treeLayout", "Layout",
                      choices = c("rectangular", "circular", "fan", "radial"),
                      selected = "rectangular"),
          selectInput("treeColorBy", "Colour tips by", choices = NULL),
          sliderInput("treeTipSize", "Tip label size", min = 0, max = 5, value = 2, step = 0.5),
          numericInput("treePruneN", "Show top N taxa (0 = all)", value = 50, min = 0),
          actionButton("updateTree", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Phylogeny"),
          bslib::card_body(min_height = 460, plotOutput("phylogeneticTree"))
        )
      )
    ),

    bslib::nav_panel(
      "Heat tree",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Heat tree",
          div(class = "info-box",
              "Taxonomic hierarchy drawn as a tree, with node size and colour scaled by abundance."),
          selectInput("taxLevelHeatTree", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                      selected = "Genus"),
          actionButton("updateHeatTree", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Taxonomic heat tree"),
          bslib::card_body(min_height = 440, plotOutput("heatTree"))
        )
      )
    )
  ),

  # ===========================================================================
  # NETWORK & CLUSTERING
  # ===========================================================================
  bslib::nav_menu(
    "Networks",

    bslib::nav_panel(
      "Abundance heatmap",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Heatmap",
          div(class = "info-box",
              "Clustered abundances across samples and taxa. CLR stabilises variance for ",
              "compositional data and is shown on a diverging scale centred at zero."),
          selectInput("taxLevelHeatmap", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                      selected = "Genus"),
          selectInput("regionHeatmap", "Annotate samples by", choices = NULL),
          selectInput("hmTransform", "Transformation",
                      choices = c("CLR" = "clr", "Relative" = "compositional", "Log10p" = "log10p"),
                      selected = "clr"),
          sliderInput("hmTopN", "Number of taxa", min = 5, max = 50, value = 20),
          actionButton("updateHeatmap", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Abundance heatmap"),
          bslib::card_body(min_height = 460, plotlyOutput("interactiveHeatmap"))
        )
      )
    ),

    bslib::nav_panel(
      "Sample clustering",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Clustering",
          div(class = "info-box",
              "Hierarchical clustering of samples by community composition. ",
              "Bray-Curtis with Ward.D2 linkage is the common microbiome choice."),
          selectInput("distanceMethod", "Distance",
                      choices = c("bray", "jaccard", "euclidean", "manhattan", "canberra"),
                      selected = "bray"),
          selectInput("clusterMethod", "Linkage",
                      choices = c("ward.D2", "complete", "average", "single", "mcquitty"),
                      selected = "ward.D2"),
          selectInput("dendroColorBy", "Colour labels by", choices = NULL),
          checkboxInput("dendroShowBar", "Show composition alongside", value = FALSE),
          actionButton("updateDendro", "Update", class = "btn-primary w-100")
        ),
        bslib::navset_card_underline(
          title = "Sample clustering",
          full_screen = TRUE,
          bslib::nav_panel("Dendrogram",  plotlyOutput("dendrogram")),
          bslib::nav_panel("Composition", plotlyOutput("dendroBar"))
        )
      )
    ),

    bslib::nav_panel(
      "Co-occurrence network",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Network",
          div(class = "info-box",
              "Correlations between taxa on CLR-transformed abundances. Edges are filtered by ",
              "correlation strength only, so treat them as exploratory rather than tested."),
          selectInput("taxLevelNetwork", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                      selected = "Genus"),
          sliderInput("corThreshold", "Minimum |correlation|",
                      min = 0.3, max = 0.99, value = 0.6, step = 0.05),
          selectInput("corMethod", "Correlation method",
                      choices = c("spearman", "pearson"), selected = "spearman"),
          sliderInput("netTopN", "Number of taxa", min = 10, max = 100, value = 40),
          checkboxInput("netNegEdges", "Include negative correlations", value = TRUE),
          actionButton("updateNetwork", "Update", class = "btn-primary w-100")
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Co-occurrence network"),
          bslib::card_body(min_height = 420,
                           networkD3::forceNetworkOutput("correlationNetwork")),
          bslib::card_body(fill = FALSE, class = "border-top",
                           verbatimTextOutput("networkSummary"))
        )
      )
    )
  ),

  # ===========================================================================
  # STATISTICS
  # ===========================================================================
  bslib::nav_menu(
    "Statistics",

    bslib::nav_panel(
      "Differential abundance",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Differential abundance",
          div(class = "info-box",
              "Kruskal-Wallis per taxon with Benjamini-Hochberg correction across taxa, ",
              "identifying those that differ between groups."),
          selectInput("taxLevelDA", "Taxonomic level",
                      choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      selected = "Genus"),
          selectInput("regionDA", "Group by", choices = NULL),
          selectInput("blockDA", "Block / subject variable (optional)", choices = NULL),
          div(class = "info-box",
              "As on the diversity tabs, setting a subject variable switches each taxon to a ",
              "Friedman test so repeated measures are not treated as independent samples."),
          sliderInput("daFDR", "FDR threshold", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
          div(class = "cl-hint",
              "Taxa present in fewer than 10% of samples are excluded before testing, ",
              "so they do not spend the FDR budget."),
          actionButton("updateDA", "Update", class = "btn-primary w-100")
        ),
        bslib::navset_card_underline(
          title = "Differential abundance",
          full_screen = TRUE,
          bslib::nav_panel("Volcano", plotlyOutput("daVolcano")),
          bslib::nav_panel("Results", DTOutput("daTable"))
        )
      )
    ),

    bslib::nav_panel(
      "Transformations",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 310, title = "Transformations",
          div(class = "info-box",
              "Compare how each normalisation reshapes the data before downstream analysis."),
          selectInput("transMethod", "Transformation",
                      choices = c("Raw" = "raw", "Relative" = "compositional", "CLR" = "clr",
                                  "Hellinger" = "hellinger", "Log10p" = "log10p"),
                      selected = "clr"),
          selectInput("transColorBy", "Colour by", choices = NULL),
          actionButton("updateTrans", "Update", class = "btn-primary w-100")
        ),
        bslib::navset_card_underline(
          title = "Transformations",
          full_screen = TRUE,
          bslib::nav_panel("Ordination",   plotlyOutput("transPCoA")),
          bslib::nav_panel("Distribution", plotlyOutput("transHist"))
        )
      )
    )
  ),

  bslib::nav_spacer(),
  bslib::nav_item(bslib::input_dark_mode(id = "dark_mode", mode = "light")),
  bslib::nav_item(
    tags$a(
      href = "https://github.com/barah123/canis_lupus2.0", target = "_blank",
      class = "nav-link", "Source"
    )
  ),

  footer = div(
    class = "cl-footer",
    div(
      a(href = "mailto:pyappiah561@gmail.com", "Email"),
      a(href = "https://www.linkedin.com/in/philip-appiah", target = "_blank", "LinkedIn"),
      a(href = "https://github.com/barah123", target = "_blank", "GitHub")
    ),
    div(class = "mt-1", "CanisLupus 2.0 - open-source microbiome analysis - Philip Appiah")
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # ── Reactive: processed phyloseq object ─────────────────────────────────────
  ps <- reactiveVal(NULL)

  # TRUE only when the phylogeny came from a real tree file. A placeholder tree
  # is generated when none is supplied so that tree-based views still render,
  # but UniFrac distances computed on it are meaningless, so phylogenetic
  # distance metrics are withheld unless this is TRUE.
  tree_is_real <- reactiveVal(FALSE)

  # TRUE when the loaded table is already relative abundance (typical of shotgun
  # profilers), which changes what the depth-based controls can mean.
  data_is_relative <- reactiveVal(FALSE)

  # Plots are drawn server-side and cannot see the client's colour mode, so the
  # toggle is read here and fed into the shared ggplot theme.
  is_dark    <- reactive(identical(input$dark_mode, "dark"))
  plot_theme <- function() theme_canis(dark = isolate(is_dark()))
  
  # ── Update all group selectors when ps changes ────────────────────────────
  observe({
    req(ps())
    vars <- colnames(sample_data(ps()))
    for (id in c("regionBar","regionPie","regionRare","regionAlpha","regionBeta",
                 "regionCore","regionHeatmap","dendroColorBy","regionDA",
                 "transColorBy","blockAlpha","blockBeta","blockDA")) {
      updateSelectInput(session, id, choices = c("None" = "", vars), selected = "")
    }
  })

  # ── Offer only the taxonomic ranks this dataset actually has ───────────────
  # Rank availability varies by source: the bundled amplicon data stops at
  # Genus, while shotgun profiles usually resolve to Species. Offering a rank
  # the taxonomy table lacks would send an unresolvable column name into every
  # downstream module.
  observe({
    req(ps())
    ranks <- rank_names(ps())
    ranks <- ranks[!is.na(ranks) & nzchar(ranks)]
    if (!length(ranks)) return()
    # Tree tips are taxa, so they colour by a taxonomic rank. "None" is offered
    # because at deep ranks every tip is its own colour, which tells you nothing.
    updateSelectInput(session, "treeColorBy",
                      choices = c("None" = "", ranks),
                      selected = isolate(if (!is.null(input$treeColorBy) &&
                                             input$treeColorBy %in% ranks)
                                         input$treeColorBy else "Phylum"))
    for (id in c("taxLevelBar", "taxLevelPie", "taxLevelCore",
                 "taxLevelHeatmap", "taxLevelHeatTree", "taxLevelNetwork",
                 "taxLevelDA")) {
      current   <- isolate(input[[id]])
      preferred <- if (id %in% c("taxLevelBar", "taxLevelPie")) "Phylum" else "Genus"
      selected  <- if (!is.null(current) && current %in% ranks) current
                   else if (preferred %in% ranks) preferred
                   else ranks[length(ranks)]
      updateSelectInput(session, id, choices = ranks, selected = selected)
    }
  })

  # ── Only offer diversity indices this data can support ─────────────────────
  observe({
    req(ps())
    avail <- available_alpha_measures(data_is_relative())
    current <- isolate(input$alphaMeasure)
    keep <- intersect(current, avail)
    if (!length(keep)) keep <- intersect(c("Shannon", "Observed"), avail)
    updateCheckboxGroupInput(session, "alphaMeasure", choices = avail, selected = keep)
    if (isTRUE(data_is_relative())) {
      showNotification(
        paste("Chao1, ACE and Fisher need integer read counts and are unavailable for a",
              "relative-abundance profile. Shannon, Simpson and richness remain valid."),
        type = "warning", duration = 9
      )
    }
  })

  # ── Headline numbers for the data tab ──────────────────────────────────────
  # Headline numbers. Posit recommends textOutput() placeholders inside
  # value_box() so the boxes appear before the values resolve, which avoids the
  # layout shifting once data loads.
  output$datasetStatsUI <- renderUI({
    bslib::layout_columns(
      fill = FALSE,
      col_widths = bslib::breakpoints(sm = c(6, 6, 6, 6, 12),
                                      lg = c(2, 2, 2, 3, 3)),
      bslib::value_box(title = "Samples",         value = textOutput("statSamples"),
                       theme = "primary"),
      bslib::value_box(title = "Taxa",            value = textOutput("statTaxa"),
                       theme = "text-primary"),
      bslib::value_box(title = "Metadata fields", value = textOutput("statMeta"),
                       theme = "text-secondary"),
      bslib::value_box(title = uiOutput("statDepthTitle", inline = TRUE),
                       value = textOutput("statDepth"),
                       theme = "text-secondary",
                       textOutput("statDepthNote")),
      bslib::value_box(title = "Phylogeny",       value = textOutput("statTree"),
                       theme = "text-success",
                       textOutput("statTreeNote"))
    )
  })

  output$statSamples <- renderText({ if (is.null(ps())) "-" else nsamples(ps()) })
  output$statTaxa    <- renderText({
    if (is.null(ps())) "-" else format(ntaxa(ps()), big.mark = ",")
  })
  output$statMeta    <- renderText({
    if (is.null(ps())) "-" else ncol(sample_data(ps()))
  })
  output$statDepthTitle <- renderUI({
    if (isTRUE(data_is_relative())) "Median total" else "Median depth"
  })
  output$statDepth   <- renderText({
    if (is.null(ps())) "-"
    else format(round(stats::median(sample_sums(ps()))), big.mark = ",")
  })
  output$statDepthNote <- renderText({
    if (is.null(ps())) ""
    else if (isTRUE(data_is_relative())) "relative abundance"
    else {
      d <- sample_sums(ps())
      paste0("range ", format(round(min(d)), big.mark = ","), " - ",
             format(round(max(d)), big.mark = ","))
    }
  })
  output$statTree    <- renderText({
    if (is.null(ps())) "-" else if (isTRUE(tree_is_real())) "Provided" else "None"
  })
  output$statTreeNote <- renderText({
    if (is.null(ps())) ""
    else if (isTRUE(tree_is_real())) "UniFrac available" else "UniFrac unavailable"
  })

  # ── Withhold phylogenetic distances when there is no real tree ─────────────
  # UniFrac on the placeholder tree yields numbers that look valid but carry no
  # phylogenetic signal, so those options are only offered with a real tree.
  observe({
    req(ps())
    allowed <- allowed_beta_distances(tree_is_real())

    current  <- isolate(input$betaDistance)
    selected <- if (!is.null(current) && current %in% allowed) current else "bray"

    updateSelectInput(session, "betaDistance", choices = allowed, selected = selected)

    if (!isTRUE(tree_is_real())) {
      showNotification(
        paste("No phylogenetic tree was uploaded, so UniFrac distances are unavailable.",
              "Other analyses are unaffected; upload a tree file to enable them."),
        type = "warning", duration = 10
      )
    }
  })
  
  # ── Load Data ────────────────────────────────────────────────────────────────
  observeEvent(input$update, {
    withProgress(message = "Loading and processing data…", value = 0, {
      
      tryCatch({
        is_shotgun <- identical(input$dataType, "shotgun")
        profile_relative <- FALSE

        if (is_shotgun) {
          incProgress(0.15, detail = "Reading taxonomic profile")
          if (is.null(input$profileFile)) {
            stop("Upload a taxonomic profile, or switch to Amplicon to use the ",
                 "bundled example data.", call. = FALSE)
          }
          # Handles both tab- and comma-separated exports, and skips any
          # version banner MetaPhlAn writes above the header.
          prof <- read_profile_file(input$profileFile$datapath)

          incProgress(0.2, detail = "Collapsing lineage hierarchy")
          parsed <- parse_taxonomic_profile(prof)
          asv <- as.data.frame(parsed$otu, check.names = FALSE)
          tax_mat <- parsed$tax
          profile_relative <- parsed$is_relative
        } else {
          incProgress(0.1, detail = "Reading abundance table")
          asv <- if (!is.null(input$asvFile))
            read.csv(input$asvFile$datapath,  row.names=1, check.names=FALSE)
          else
            read.csv(demo_file("ASV.csv"),     row.names=1, check.names=FALSE)

          incProgress(0.2, detail = "Reading taxonomy")
          tax_mat <- if (!is.null(input$taxFile))
            as.matrix(read.csv(input$taxFile$datapath,  row.names=1))
          else
            as.matrix(read.csv(demo_file("taxonomy.csv"), row.names=1))
        }

        incProgress(0.2, detail = "Reading metadata")
        meta <- if (!is.null(input$metaFile)) {
          read.csv(input$metaFile$datapath, row.names=1)
        } else if (is_shotgun) {
          # Without metadata every sample still needs a row, so that grouping
          # controls have something to offer.
          data.frame(SampleID = colnames(asv), row.names = colnames(asv),
                     stringsAsFactors = FALSE)
        } else {
          read.csv(demo_file("metadata.csv"), row.names=1)
        }

        incProgress(0.1, detail = "Building phyloseq object")
        shared <- intersect(colnames(asv), rownames(meta))
        if (!length(shared)) {
          stop("No sample names are shared between the abundance table and the ",
               "metadata. Column names in the abundance table must match the ",
               "metadata row names.", call. = FALSE)
        }
        asv  <- asv[, shared, drop = FALSE]
        meta <- meta[shared, , drop = FALSE]

        otu  <- otu_table(as.matrix(asv), taxa_are_rows=TRUE)
        tax  <- tax_table(tax_mat)
        sam  <- sample_data(meta)
        ps_obj <- phyloseq(otu, tax, sam)
        
        # Load or generate tree
        using_demo <- !is_shotgun && is.null(input$asvFile)
        demo_tree_path <- if (using_demo) {
          tryCatch(demo_file("tree.txt"), error = function(e) NULL)
        } else NULL
        if (!is.null(input$treeFile)) {
          tree <- read.tree(input$treeFile$datapath)
          ps_obj <- merge_phyloseq(ps_obj, tree)
          tree_is_real(TRUE)
        } else if (!is.null(demo_tree_path)) {
          # Demo data ships with a real tree; use it rather than a random one.
          tree <- read.tree(demo_tree_path)
          ps_obj <- merge_phyloseq(ps_obj, tree)
          tree_is_real(TRUE)
        } else {
          # Placeholder tree so tree-based VIEWS still render. It carries no
          # phylogenetic information, so tree_is_real stays FALSE and the
          # UniFrac distance options are withheld downstream.
          set.seed(42)
          rand_tree <- ape::rtree(ntaxa(ps_obj), rooted=TRUE,
                                  tip.label=taxa_names(ps_obj))
          ps_obj <- merge_phyloseq(ps_obj, rand_tree)
          tree_is_real(FALSE)
        }
        
        incProgress(0.2, detail = "Applying filters")
        # Depth and singleton filters are defined on read counts. A profile
        # already expressed as relative abundance has no depth to filter on and
        # no singletons to remove, so both are skipped rather than applied to
        # numbers where they mean nothing.
        n_before <- ntaxa(ps_obj)
        n_samples_before <- nsamples(ps_obj)
        filters_applied <- !profile_relative

        if (filters_applied) {
          min_reads <- input_num(input$minReads, 0)
          # If the cutoff would remove everything, the table is almost certainly
          # not read counts. Skip the filter and say so, rather than failing.
          if (min_reads > 0 && !any(sample_sums(ps_obj) >= min_reads)) {
            showNotification(
              paste0("Every sample totals less than ", min_reads,
                     ", so the minimum-reads filter was skipped. If this table is ",
                     "not raw read counts, that filter does not apply to it."),
              type = "warning", duration = 10)
          } else if (min_reads > 0) {
            ps_obj <- prune_samples(sample_sums(ps_obj) >= min_reads, ps_obj)
          }
          if (isTRUE(input$removeSingletons)) ps_obj <- remove_singletons(ps_obj)
        }
        # Taxa absent everywhere carry no information in either input type.
        ps_obj <- prune_taxa(taxa_sums(ps_obj) > 0, ps_obj)
        n_after <- ntaxa(ps_obj)

        if (!nsamples(ps_obj) || !ntaxa(ps_obj)) {
          stop("Filtering removed every sample or taxon. Lower the minimum reads ",
               "per sample, or turn off singleton removal.", call. = FALSE)
        }

        data_is_relative(profile_relative)
        ps(ps_obj)

        incProgress(0.2, detail = "Done")
        showNotification(
          paste0("Loaded ", nsamples(ps_obj), " samples and ", n_after, " taxa",
                 if (!filters_applied) " (relative-abundance profile; depth filters skipped)"
                 else if (isTRUE(input$removeSingletons))
                   paste0(" (", n_before - n_after, " taxa removed by filtering)")
                 else "",
                 if (nsamples(ps_obj) < n_samples_before)
                   paste0(", ", n_samples_before - nsamples(ps_obj), " samples below depth cutoff")
                 else ""),
          type = "message", duration = 6
        )
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type="error", duration=10)
      })
    })
  })
  
  # ── Download summary ─────────────────────────────────────────────────────
  output$downloadPS <- downloadHandler(
    filename = function() paste0("canis_lupus_summary_", Sys.Date(), ".txt"),
    content  = function(file) {
      req(ps())
      sink(file)
      cat("=== CanisLupus 2.0 - Filtered Dataset Summary ===\n\n")
      print(ps())
      cat("\nSample read counts:\n")
      print(sort(sample_sums(ps())))
      cat("\nTax ranks:\n")
      print(rank_names(ps()))
      sink()
    }
  )
  
  # ── Data Summary ─────────────────────────────────────────────────────────
  output$dataSummary <- renderPrint({
    req(ps())
    cat("=== Microbiome Dataset ===\n\n")
    cat("Samples    :", nsamples(ps()), "\n")
    cat("Taxa       :", ntaxa(ps()),    "\n")
    cat("Tax ranks  :", paste(rank_names(ps()), collapse=" > "), "\n")
    cat("Tree tips  :", ifelse(!is.null(phy_tree(ps(), errorIfNULL=FALSE)),
                               length(phy_tree(ps())$tip.label), "none"), "\n")
    cat("\nRead depth (min/median/max):",
        min(sample_sums(ps())), "/",
        round(median(sample_sums(ps()))), "/",
        max(sample_sums(ps())), "\n")
    cat("\nSample variables:\n")
    print(colnames(sample_data(ps())))
  })
  
  # ── QC Read Count Plot ───────────────────────────────────────────────────
  output$readCountPlot <- renderPlotly({
    req(ps())
    df <- data.frame(
      Sample = sample_names(ps()),
      Reads  = sample_sums(ps())
    ) %>% arrange(Reads)
    df$Sample <- factor(df$Sample, levels=df$Sample)
    
    plot_ly(df, x=~Sample, y=~Reads, type="bar",
            marker=list(color=canis_ui$primary)) %>%
      layout(title  = "Read Counts per Sample",
             xaxis  = list(tickangle=-45, title=""),
             yaxis  = list(title="Total Reads"),
             shapes = list(list(type="line", x0=0, x1=1, xref="paper",
                                y0=input$minReads, y1=input$minReads,
                                line=list(color="red", dash="dash"))))
  })
  
  output$sampleTable <- renderDT({
    req(ps())
    datatable(as(sample_data(ps()), "data.frame"), fillContainer = TRUE,
              options=list(scrollX=TRUE, pageLength=10),
              class="cell-border stripe")
  })
  
  # ============================================================
  # RAREFACTION CURVE (improved: proper per-sample curves)
  # ============================================================
  output$rarefactionCurve <- renderPlotly({
    req(ps())
    input$updateRare
    input$dark_mode
    isolate({
      # Rarefaction subsamples reads, so it is only defined on integer counts.
      # A relative-abundance profile has no reads left to subsample.
      if (isTRUE(data_is_relative())) {
        return(plot_ly() %>%
                 add_annotations(
                   text = paste0("Rarefaction needs raw read counts.<br><br>",
                                 "This dataset is a relative-abundance profile, so there are no<br>",
                                 "reads to subsample. Sequencing depth was already decided<br>",
                                 "upstream by the profiler."),
                   x = 0.5, y = 0.5, showarrow = FALSE, align = "left") %>%
                 layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }
      withProgress(message="Computing rarefaction curves…", {
        otu_mat <- t(as(otu_table(ps()), "matrix"))
        depths   <- rowSums(otu_mat)
        maxd     <- if (input_num(input$rareMax, 0) > 0) input_num(input$rareMax, 0) else max(depths)

        # Total work is (summed depth / step), and it grows without bound as the
        # step shrinks. On a deeply sequenced study a small step will exhaust
        # memory outright, so the step is raised to whatever keeps the curves
        # under a sane number of points, and the user is told.
        MAX_POINTS <- 8000
        asked_step <- max(input_num(input$rareStep, 100), 1)
        min_step   <- max(1, ceiling(sum(depths) / MAX_POINTS))
        step_sz    <- max(asked_step, min_step)
        if (step_sz > asked_step) {
          showNotification(
            sprintf(paste("A step of %s would need roughly %s points across %d samples.",
                          "Using %s instead to keep the curves computable."),
                    format(asked_step, big.mark = ","),
                    format(round(sum(depths) / asked_step), big.mark = ","),
                    nrow(otu_mat), format(step_sz, big.mark = ",")),
            type = "warning", duration = 9)
        }

        # tidy = TRUE is essential, not cosmetic: the default draws a base-R
        # plot as a side effect and only returns the curves invisibly. Inside
        # renderPlotly there is no device sized for that, so it fails with
        # "figure margins too large". The tidy form touches no device at all.
        rc <- vegan::rarecurve(otu_mat, step = step_sz,
                               sample = input_num(input$rareDepth, 10000),
                               tidy = TRUE)
        names(rc)[match(c("Site", "Sample", "Species"), names(rc))] <-
          c("SampleID", "Reads", "OTUs")

        # rarecurve() has no maximum-depth argument, so the cap is applied to
        # the returned curves; passing one through would be silently ignored.
        if (input_num(input$rareMax, 0) > 0) {
          keep <- rc$Reads <= maxd
          if (any(keep)) rc <- rc[keep, , drop = FALSE]
        }

        # as.data.frame() leaves this an S4 sample_data, whose `[` does not drop
        # to a vector: samp_data[ids, var] would come back as a one-column
        # sample_data, and as.character() would deparse it into the single
        # literal string 'c("BB", "BB", ...)'. Every curve then shared one
        # bogus level and the legend printed that string. as() gives a real
        # data frame, and the column is looked up by name to be certain.
        samp_data <- as(sample_data(ps()), "data.frame")
        grp_var   <- if (!is.null(input$regionRare) && nzchar(input$regionRare) &&
                         input$regionRare %in% colnames(samp_data))
          input$regionRare else NULL

        rc$Group <- if (!is.null(grp_var)) {
          g <- as.character(samp_data[[grp_var]])[
                 match(as.character(rc$SampleID), rownames(samp_data))]
          g[is.na(g) | !nzchar(g)] <- "Unassigned"
          g
        } else "All samples"

        p <- ggplot(rc, aes(x = Reads, y = OTUs, group = SampleID, colour = Group)) +
          geom_line(alpha = 0.75, linewidth = 0.6) +
          plot_theme() +
          scale_color_canis() +
          labs(title = "Rarefaction curves", x = "Sequencing depth (reads)",
               y = "Observed taxa",
               colour = if (!is.null(grp_var)) grp_var else NULL)

        if (isTRUE(input$showRareLine)) {
          p <- p + geom_vline(xintercept = input_num(input$rareDepth, 10000),
                              linetype = "dashed", colour = canis_ui$danger,
                              linewidth = 0.6)
        }
        # One line is drawn per sample, so ggplotly emits one trace per sample
        # and each claims its own legend entry. With 75 samples the legend
        # swamps the panel. Collapse it to one entry per group, and let a
        # click on that entry toggle the whole group.
        gp <- ggplotly(p)
        if (!is.null(grp_var)) {
          shown <- character(0)
          for (i in seq_along(gp$x$data)) {
            tr <- gp$x$data[[i]]
            if (!identical(tr$mode, "lines") || is.null(tr$legendgroup)) next
            key <- as.character(tr$legendgroup)[1]
            gp$x$data[[i]]$name       <- key
            gp$x$data[[i]]$showlegend <- !(key %in% shown)
            shown <- c(shown, key)
          }
        }
        gp
      })
    })
  })
  
  output$rareSummary <- renderPrint({
    req(ps())
    input$updateRare
    input$dark_mode
    isolate({
      ss <- sort(sample_sums(ps()))
      if (isTRUE(data_is_relative())) {
        cat("=== Sample totals ===\n\n")
        cat("This dataset is a relative-abundance profile, so these are proportions\n")
        cat("rather than read counts, and rarefaction does not apply.\n\n")
        cat("Total per sample: min", signif(min(ss), 4),
            "| median", signif(stats::median(ss), 4),
            "| max", signif(max(ss), 4), "\n")
        return(invisible(NULL))
      }
      cat("=== Read-depth summary ===\n")
      cat("Min:", min(ss), "\n")
      cat("1st Qu:", quantile(ss, 0.25), "\n")
      cat("Median:", median(ss), "\n")
      cat("Mean:", round(mean(ss)), "\n")
      cat("3rd Qu:", quantile(ss, 0.75), "\n")
      cat("Max:", max(ss), "\n")
      cat("\nSamples below rarefaction depth (", input$rareDepth, "):",
          sum(ss < input$rareDepth), "\n")
    })
  })
  
  # ============================================================
  # ALPHA DIVERSITY
  # ============================================================
  #' Alpha diversity in long form, shared by both plots and the test output.
  alpha_long <- reactive({
    req(ps())
    measures <- input$alphaMeasure
    if (is.null(measures) || !length(measures)) measures <- "Shannon"
    rich <- compute_alpha(ps(), measures, data_is_relative())
    meta <- as(sample_data(ps()), "data.frame")
    long <- utils::stack(rich)
    names(long) <- c("Value", "Index")
    long$Sample <- rep(rownames(rich), times = ncol(rich))
    grp_var <- input$regionAlpha
    long$Group <- if (!is.null(grp_var) && nzchar(grp_var) && grp_var %in% colnames(meta)) {
      as.character(meta[long$Sample, grp_var])
    } else "All samples"
    list(long = long, wide = rich, measures = colnames(rich), group_var = grp_var)
  })

  output$alphaDiv <- renderPlotly({
    req(ps())
    input$updateAlpha
    input$dark_mode
    isolate({
      a <- alpha_long()
      p <- ggplot(a$long, aes(x = Sample, y = Value, colour = Group)) +
        geom_point(size = 2.6, alpha = 0.9) +
        facet_wrap(~Index, scales = "free_y") +
        plot_theme() +
        scale_color_canis() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(title = "Alpha diversity per sample", x = "Sample", y = NULL,
             colour = if (!is.null(a$group_var) && nzchar(a$group_var)) a$group_var else NULL)
      ggplotly(p)
    })
  })

  output$alphaDivBoxplot <- renderPlotly({
    req(ps())
    input$updateAlpha
    input$dark_mode
    isolate({
      a <- alpha_long()
      if (is.null(a$group_var) || !nzchar(a$group_var) ||
          !a$group_var %in% colnames(as(sample_data(ps()), "data.frame"))) {
        return(plot_ly() %>%
                 add_annotations(text = "Select a grouping variable to compare groups.",
                                 x = 0.5, y = 0.5, showarrow = FALSE) %>%
                 layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }
      p <- ggplot(a$long, aes(x = Group, y = Value, fill = Group)) +
        # Both layers key off Group; letting each contribute a legend makes
        # ggplotly emit combined "(level,1)/(level,2)" entries, so only the
        # fill legend is kept.
        geom_boxplot(alpha = 0.65, outlier.shape = NA, colour = "#5A6C77", linewidth = 0.35) +
        geom_jitter(width = 0.16, size = 1.9, alpha = 0.85,
                    colour = canis_ui$ink, show.legend = FALSE) +
        facet_wrap(~Index, scales = "free_y") +
        plot_theme() +
        scale_fill_canis() +
        scale_color_canis() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        labs(title = paste("Alpha diversity by", a$group_var), x = NULL, y = NULL)
      ggplotly(p)
    })
  })
  
  output$alphaStats <- renderPrint({
    req(ps())
    input$updateAlpha
    input$dark_mode
    isolate({
      if (is.null(input$regionAlpha) || nchar(input$regionAlpha)==0 ||
          !input$alphaStats) return(invisible(NULL))
      grp_var <- input$regionAlpha
      if (!grp_var %in% colnames(sample_data(ps()))) return(invisible(NULL))
      measures <- input$alphaMeasure
      if (length(measures)==0) measures <- "Shannon"
      meta_df <- as(sample_data(ps()), "data.frame")
      rich <- compute_alpha(ps(), measures, data_is_relative())
      measures <- colnames(rich)
      rich$Group <- as.factor(meta_df[[grp_var]])

      # Optional blocking/subject variable. See choose_alpha_test(): a Friedman
      # test replaces Kruskal-Wallis only for a complete block design, e.g. the
      # same individuals sampled at every level of the grouping variable.
      block_var <- input$blockAlpha
      has_block <- !is.null(block_var) && nchar(block_var) > 0 &&
        block_var %in% colnames(meta_df)
      if (has_block) rich$Block <- as.factor(meta_df[[block_var]])

      test_choice <- choose_alpha_test(
        group     = rich$Group,
        block     = if (has_block) rich$Block else NULL,
        block_var = if (has_block) block_var else NULL,
        group_var = grp_var
      )
      use_block <- identical(test_choice, "friedman")

      if (has_block && !use_block) {
        if (identical(block_var, grp_var)) {
          cat("Block/subject variable is the same as the grouping variable -- ignoring it.\n\n")
        } else {
          block_tab <- table(rich$Group, rich$Block)
          n_empty <- sum(block_tab == 0)
          n_multi <- sum(block_tab > 1)
          cat("=== Blocking requested but not usable ===\n\n")
          cat("A Friedman test needs exactly one sample per combination of '", grp_var,
              "' and '", block_var, "' (a complete block design), but ", sep="")
          if (n_empty > 0 && n_multi > 0) {
            cat(n_empty, " combination(s) have no sample and ", n_multi, " have more than one.\n", sep="")
          } else if (n_empty > 0) {
            cat(n_empty, " combination(s) have no sample.\n", sep="")
          } else {
            cat(n_multi, " combination(s) have more than one sample.\n", sep="")
          }
          if (n_empty > 0) {
            cat("\nIf you expected a complete design, check whether the QC settings dropped samples: ",
                "the 'Min reads per sample' filter on the Data Upload tab removes low-depth samples, ",
                "which can leave gaps here. Lowering it (or setting it to 0) restores them.\n", sep="")
          }
          cat("\nFalling back to an unblocked Kruskal-Wallis test below. If your samples really are ",
              "repeated measures on the same subjects, interpret that result with caution.\n\n", sep="")
        }
      }

      if (use_block) {
        cat("=== Friedman Tests (blocked by '", block_var, "') ===\n\n", sep="")
        cat("Each level of '", grp_var, "' was sampled exactly once within each level of '",
            block_var, "', so a Friedman test is used instead of Kruskal-Wallis to account for ",
            "this repeated-measures/blocked design.\n\n", sep="")
        for (m in measures) {
          if (m %in% colnames(rich)) {
            ft <- friedman.test(rich[[m]], groups = rich$Group, blocks = rich$Block)
            cat(m, ": chi-sq =", round(unname(ft$statistic), 3),
                ", df =", ft$parameter, ", p =", round(ft$p.value, 4), "\n")
          }
        }
      } else {
        cat("=== Kruskal-Wallis Tests ===\n\n")
        cat("Treats all samples as independent draws. If the same subjects were sampled more than\n",
            "once across '", grp_var, "' (e.g. multiple body sites or timepoints per person), select\n",
            "a block/subject variable on the left to use a paired-appropriate test instead.\n\n", sep="")
        for (m in measures) {
          if (m %in% colnames(rich)) {
            kw <- kruskal.test(rich[[m]] ~ rich$Group)
            cat(m, ": p =", round(kw$p.value, 4), "\n")
          }
        }
      }
    })
  })
  
  # ============================================================
  # BETA DIVERSITY
  # ============================================================

  #' Phyloseq object to use for all beta-diversity distance/ordination work.
  #' Samples usually differ in sequencing depth, and distances such as
  #' Bray-Curtis are sensitive to that, so the default converts to relative
  #' abundance first; otherwise library size, not composition, drives the
  #' apparent dissimilarity.
  beta_ps <- reactive({
    req(ps())
    normalise_for_beta(ps(), if (identical(input$betaNormalise, "none")) "none" else "relative")
  })

  #' Distance metric to actually use, refusing UniFrac without a real tree.
  #' The selector already hides those options; this is the safety net for a
  #' stale input value, so a placeholder tree can never reach a UniFrac result.
  beta_distance <- reactive({
    resolve_beta_distance(input$betaDistance, tree_is_real())
  })

  output$betaPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    input$dark_mode
    isolate({
      ps_beta <- beta_ps()
      dmetric <- beta_distance()
      validate(need(!is.null(dmetric),
                    "UniFrac needs a phylogenetic tree. Upload one on the Data Upload tab, or choose a non-phylogenetic distance."))
      ord <- ordinate(ps_beta, method=input$betaMethod, distance=dmetric)
      p   <- plot_ordination(ps_beta, ord, type="samples") +
        plot_theme() +
        labs(title=paste(input$betaMethod, "(", dmetric, ")"))
      if (!is.null(input$regionBeta) && nchar(input$regionBeta)>0 &&
          input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color=input$regionBeta), size=3)
        if (input$betaEllipse)
          p <- p + stat_ellipse(aes_string(color=input$regionBeta), level=0.95)
      } else {
        p <- p + geom_point(color=canis_ui$primary, size=3)
      }
      ggplotly(p)
    })
  })
  
  output$nmdsPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    input$dark_mode
    isolate({
      ps_beta <- beta_ps()
      dmetric <- beta_distance()
      validate(need(!is.null(dmetric),
                    "UniFrac needs a phylogenetic tree. Upload one on the Data Upload tab, or choose a non-phylogenetic distance."))
      suppressMessages({
        ord <- ordinate(ps_beta, method="NMDS", distance=dmetric)
      })
      p <- plot_ordination(ps_beta, ord, type="samples") +
        plot_theme() +
        labs(title=paste("NMDS (", dmetric, ")"))
      if (!is.null(input$regionBeta) && nchar(input$regionBeta)>0 &&
          input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color=input$regionBeta), size=3)
        if (input$betaEllipse)
          p <- p + stat_ellipse(aes_string(color=input$regionBeta), level=0.95)
      } else {
        p <- p + geom_point(color=canis_ui$secondary, size=3)
      }
      ggplotly(p)
    })
  })
  
  output$permanovaResult <- renderPrint({
    req(ps())
    input$updateBeta
    input$dark_mode
    isolate({
      if (!input$betaPermanova) return(invisible(NULL))
      grp_var <- input$regionBeta
      if (is.null(grp_var) || nchar(grp_var)==0 ||
          !grp_var %in% colnames(sample_data(ps()))) {
        cat("Select a grouping variable to run PERMANOVA.\n")
        return(invisible(NULL))
      }
      dmetric <- beta_distance()
      if (is.null(dmetric)) {
        cat("UniFrac needs a phylogenetic tree.\n\n")
        cat("No tree was uploaded, so the app is using a placeholder phylogeny that carries no\n")
        cat("real evolutionary information. Running UniFrac on it would produce a result that\n")
        cat("looks valid but is not, so it has been withheld. Upload a tree file on the Data\n")
        cat("Upload tab, or choose a non-phylogenetic distance such as Bray-Curtis.\n")
        return(invisible(NULL))
      }
      dist_mat <- phyloseq::distance(beta_ps(), method=dmetric)
      meta_df  <- as(sample_data(ps()), "data.frame")
      meta_df[[grp_var]] <- as.factor(meta_df[[grp_var]])
      set.seed(42)

      cat("=== PERMANOVA (adonis2) ===\n")
      cat("Response variable:", grp_var, "\n")
      cat("Distance metric:  ", dmetric, "\n")
      cat("Normalisation:    ",
          if (identical(input$betaNormalise, "none")) "none (raw counts)" else "relative abundance", "\n")

      # Optional blocking/subject variable: restrict permutations to within
      # each block (strata) and add the block as a term, instead of the
      # default unrestricted permutation, which assumes every sample is an
      # independent draw. Unrestricted permutation on repeated-measures data
      # (e.g. the same subjects sampled at every group level) tests the wrong
      # null distribution -- see manuscript/validation_analysis.R for a
      # worked example on this app's own bundled skin_data/.
      block_var <- input$blockBeta
      use_block <- !is.null(block_var) && nchar(block_var) > 0 &&
        block_var %in% colnames(meta_df)

      if (use_block && block_var == grp_var) {
        cat("Block/subject variable is the same as the grouping variable -- ignoring it.\n")
        use_block <- FALSE
      }

      if (use_block) {
        meta_df[[block_var]] <- as.factor(meta_df[[block_var]])
        cat("Blocked by:       ", block_var,
            "(permutations restricted within each block; no block:group interaction is tested)\n\n", sep=" ")
        form <- stats::as.formula(paste0("dist_mat ~ `", block_var, "` + `", grp_var, "`"))
        perm <- vegan::adonis2(form, data = meta_df, permutations = 999,
                                strata = meta_df[[block_var]], by = "terms")
      } else {
        cat("\n")
        form <- stats::as.formula(paste0("dist_mat ~ `", grp_var, "`"))
        perm <- vegan::adonis2(form, data = meta_df, permutations = 999)
      }
      print(perm)

      # Homogeneity of dispersion
      bd <- vegan::betadisper(dist_mat, meta_df[[grp_var]])
      cat("\n=== Homogeneity of dispersion (betadisper) ===\n")
      print(anova(bd))
    })
  })
  
  # ============================================================
  # STACKED BAR PLOTS
  # ============================================================
  output$taxaBarplot <- renderPlotly({
    req(ps())
    input$updateBar
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelBar)
      validate(need(!is.null(rank), "That taxonomic rank is not present in this dataset."))
      ps_use <- if (identical(input$barTransform, "relative"))
        transform_sample_counts(ps(), function(x) if (sum(x) > 0) x/sum(x) else x) else ps()
      n_bar <- min(input_num(input$topNBar, 10), ntaxa(ps_use))
      top_t <- names(sort(taxa_sums(ps_use), decreasing=TRUE)[seq_len(n_bar)])
      ps_top <- prune_taxa(top_t, ps_use)
      p <- plot_bar(ps_top, fill=rank) +
        plot_theme() +
        labs(title=paste("Top", n_bar, rank), x = "Sample", y = "Abundance") +
        scale_fill_manual(values = canis_colors(n_bar)) +
        sample_axis_theme(nsamples(ps_top))
      if (!is.null(input$regionBar) && nchar(input$regionBar)>0 &&
          input$regionBar %in% colnames(sample_data(ps())))
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales="free_x")
      ggplotly(p)
    })
  })
  
  output$relativeAbundancePlot <- renderPlotly({
    req(ps())
    input$updateBar
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelBar)
      validate(need(!is.null(rank), "That taxonomic rank is not present in this dataset."))
      ps_glom <- safe_tax_glom(ps(), rank)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) if (sum(x) > 0) x/sum(x) else x)
      p <- plot_bar(ps_rel, fill=rank) +
        plot_theme() +
        labs(title=paste(rank, "relative abundance"), x = "Sample", y = "Relative abundance") +
        scale_fill_manual(values = canis_colors(ntaxa(ps_rel))) +
        sample_axis_theme(nsamples(ps_rel))
      if (!is.null(input$regionBar) && nchar(input$regionBar)>0 &&
          input$regionBar %in% colnames(sample_data(ps())))
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales="free_x")
      ggplotly(p)
    })
  })
  
  # ============================================================
  # PIE CHART
  # ============================================================
  output$pieChart <- renderPlotly({
    req(ps())
    input$updatePie
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelPie)
      validate(need(!is.null(rank), "That taxonomic rank is not present in this dataset."))
      ps_glom <- safe_tax_glom(ps(), rank)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      melted  <- psmelt(ps_rel)
      grp_var <- input$regionPie
      lvl_col <- input$taxLevelPie
      
      if (!is.null(grp_var) && nchar(grp_var)>0 && grp_var %in% colnames(melted)) {
        agg <- melted %>%
          group_by(across(all_of(c(lvl_col, grp_var)))) %>%
          summarise(Abundance=mean(Abundance), .groups="drop")
        plot_ly(agg, labels=~get(lvl_col), values=~Abundance,
                color=~get(grp_var), type="pie",
                textposition="inside", textinfo="label+percent") %>%
          layout(title=paste(lvl_col, "by", grp_var))
      } else {
        agg <- melted %>% group_by(across(all_of(lvl_col))) %>%
          summarise(Abundance=mean(Abundance), .groups="drop")
        plot_ly(agg, labels=~get(lvl_col), values=~Abundance, type="pie",
                marker = list(colors = canis_colors(nrow(agg)),
                              line = list(color = "#FFFFFF", width = 1)),
                textposition="inside", textinfo="label+percent") %>%
          layout(title = list(text = paste(lvl_col, "composition"), x = 0.02,
                              font = list(size = 14)))
      }
    })
  })
  
  # ============================================================
  # CORE MICROBIOME
  # ============================================================
  #' Core taxa, shared by the heatmap and the summary.
  #'
  #' Both thresholds are expressed as percentages in the interface, so detection
  #' must be applied to relative abundance. Applying it to raw counts would make
  #' every setting behave as "at least one read", since the slider maximum of
  #' 10% is still below a count of 1. Prevalence is inclusive so that the 100%
  #' setting can return the taxa present in every sample.
  core_taxa <- reactive({
    req(ps())
    rank <- safe_rank(ps(), input$taxLevelCore)
    if (is.null(rank)) return(NULL)
    ps_glom <- safe_tax_glom(ps(), rank)
    ps_rel  <- transform_sample_counts(ps_glom, function(x) if (sum(x) > 0) x / sum(x) else x)

    rel <- as(otu_table(ps_rel), "matrix")
    if (!taxa_are_rows(ps_rel)) rel <- t(rel)

    detection  <- input_num(input$detectionCore, 0.1) / 100
    prevalence <- input_num(input$prevalenceCore, 50) / 100
    prev <- rowMeans(rel > detection)
    core_ids <- names(which(prev >= prevalence))

    list(ids = core_ids, ps_glom = ps_glom, ps_rel = ps_rel, rank = rank,
         prevalence = prev)
  })

  output$coreHeatmap <- renderPlotly({
    req(ps())
    input$updateCore
    input$dark_mode
    isolate({
      core <- core_taxa()
      validate(need(!is.null(core), "That taxonomic rank is not present in this dataset."))
      if (!length(core$ids)) {
        return(plot_ly() %>%
                 add_annotations(text = "No taxa meet both thresholds. Try lowering them.",
                                 x = 0.5, y = 0.5, showarrow = FALSE) %>%
                 layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }
      # Relative abundance is computed across the whole community before
      # pruning, so the scale reads as share of the community rather than
      # share of the core.
      mat <- as(otu_table(prune_taxa(core$ids, core$ps_rel)), "matrix")
      rownames(mat) <- as.character(tax_table(core$ps_glom)[core$ids, core$rank])
      heatmaply(t(mat), colors = viridis::viridis(256),
                main = paste("Core microbiome -", core$rank),
                xlab = "Taxa", ylab = "Samples")
    })
  })

  output$coreSummary <- renderPrint({
    req(ps())
    input$updateCore
    input$dark_mode
    isolate({
      core <- core_taxa()
      if (is.null(core)) { cat("That taxonomic rank is not present in this dataset.\n"); return(invisible(NULL)) }
      cat("=== Core microbiome ===\n")
      cat("Taxonomic level:     ", core$rank, "\n")
      cat("Detection threshold: ", input_num(input$detectionCore, 0.1), "% relative abundance\n", sep = "")
      cat("Prevalence threshold:", input_num(input$prevalenceCore, 50), "% of samples\n")
      cat("Taxa examined:       ", ntaxa(core$ps_glom), "\n")
      cat("Core taxa found:     ", length(core$ids), "\n\n")
      if (length(core$ids)) {
        nm <- as.character(tax_table(core$ps_glom)[core$ids, core$rank])
        prev_pct <- round(100 * core$prevalence[core$ids])
        ord <- order(prev_pct, decreasing = TRUE)
        cat("Core taxa (prevalence at this detection threshold):\n")
        for (i in ord) cat(sprintf("  %-45s %3d%%\n", nm[i], prev_pct[i]))
      } else {
        cat("No taxon exceeds", input_num(input$detectionCore, 0.1),
            "% relative abundance in at least", input_num(input$prevalenceCore, 50),
            "% of samples.\n")
      }
    })
  })
  
  # ============================================================
  # PHYLOGENETIC TREE (improved with ggtree options)
  # ============================================================

  # Says plainly whether the displayed tree is the user's phylogeny or the
  # placeholder, so no one reads structure into a tree that has none.
  output$treeProvenanceNote <- renderUI({
    req(ps())
    if (isTRUE(tree_is_real())) {
      div(class = "info-box",
          "Using the phylogenetic tree supplied with your data.")
    } else {
      div(class = "info-box",
          style = "border-left:4px solid #B00020; background:#fdf0f0;",
          tags$strong("Placeholder tree: no phylogeny was uploaded."),
          " The branching shown here is arbitrary and carries no evolutionary meaning.",
          " It is drawn only so the layout options remain usable.",
          " Do not interpret clustering in this tree, and upload a tree file to enable",
          " UniFrac distances on the Beta Diversity tab.")
    }
  })

  output$phylogeneticTree <- renderPlot({
    
    req(ps())
    input$updateTree
    input$dark_mode
    
    isolate({
      
      tree_ps <- ps()
      
      # prune taxa if user requests
      if (input$treePruneN > 0 && ntaxa(tree_ps) > input$treePruneN) {
        
        top_t <- names(sort(taxa_sums(tree_ps), decreasing = TRUE))[1:input$treePruneN]
        
        tree_ps <- prune_taxa(top_t, tree_ps)
        
        # IMPORTANT: prune tree to match taxa
        phy_tree(tree_ps) <- ape::keep.tip(
          phy_tree(tree_ps),
          taxa_names(tree_ps)
        )
      }
      
      tree <- phy_tree(tree_ps)

      rank <- safe_rank(ps(), input$treeColorBy)
      colour_by <- !is.null(input$treeColorBy) && nzchar(input$treeColorBy) &&
        !is.null(rank)

      p <- ggtree::ggtree(tree, layout = input$treeLayout)

      if (colour_by) {
        # Attach the taxonomy to the tree so tips can be coloured by rank.
        # ggtree matches on a column named `label`, which holds the tip names.
        lab <- as.character(tax_table(tree_ps)[taxa_names(tree_ps), rank])
        lab[is.na(lab) | !nzchar(lab)] <- "Unclassified"
        tip_data <- data.frame(label = taxa_names(tree_ps),
                               Taxon = lab, stringsAsFactors = FALSE)
        p <- p %<+% tip_data +
          ggtree::geom_tippoint(ggplot2::aes(colour = Taxon), size = 1.8, na.rm = TRUE) +
          ggtree::geom_tiplab(ggplot2::aes(colour = Taxon),
                              size = input_num(input$treeTipSize, 2), na.rm = TRUE) +
          scale_color_canis() +
          ggplot2::labs(colour = rank)
      } else {
        p <- p + ggtree::geom_tiplab(size = input_num(input$treeTipSize, 2))
      }

      p <- p + ggtree::theme_tree2() +
        ggplot2::ggtitle(if (colour_by) paste("Phylogeny, tips coloured by", rank)
                         else "Phylogeny")

      print(p)
      
    })
    
  })
  
  # ============================================================
  # HEAT TREE (metacoder – improved)
  # ============================================================
  output$heatTree <- renderPlot({
    req(ps())
    input$updateHeatTree
    input$dark_mode
    isolate({
      tryCatch({
        obj <- metacoder::parse_phyloseq(ps())
        obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "otu_table")
        obj$data$tax_prop  <- obj$data$tax_abund
        obj$data$tax_prop[,-1] <- obj$data$tax_abund[,-1] /
          colSums(obj$data$tax_abund[,-1], na.rm=TRUE)
        
        metacoder::heat_tree(obj,
                             node_label  = taxon_names,
                             node_size   = n_obs,
                             node_color  = n_obs,
                             node_size_axis_label = "OTU count",
                             node_color_axis_label= "OTU count",
                             initial_layout       = "reingold-tilford",
                             layout               = "davidson-harel",
                             title = paste("Taxonomic Heat Tree –", input$taxLevelHeatTree, "level"))
      }, error = function(e) {
        plot.new()
        title(paste("Heat Tree error:", e$message))
      })
    })
  })
  
  # ============================================================
  # INTERACTIVE HEATMAP (CLR-transformed, improved)
  # ============================================================
  output$interactiveHeatmap <- renderPlotly({
    req(ps())
    input$updateHeatmap
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelHeatmap)
      validate(need(!is.null(rank), "That taxonomic rank is not present in this dataset."))
      ps_glom <- safe_tax_glom(ps(), rank)

      # Select top N taxa
      top_t <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[seq_len(min(input_num(input$hmTopN, 20), ntaxa(ps_glom)))]
      ps_top <- prune_taxa(top_t, ps_glom)
      
      # Apply transformation
      ps_trans <- tryCatch(
        microbiome::transform(ps_top, input$hmTransform),
        error = function(e) transform_sample_counts(ps_top, function(x) x/sum(x))
      )
      
      mat <- as(otu_table(ps_trans), "matrix")
      rownames(mat) <- as.character(tax_table(ps_trans)[, rank])

      # The matrix is plotted transposed, so samples are ROWS and taxa are
      # columns. The grouping variable describes samples, so it belongs on the
      # row side; passing it as a column annotation mismatches the dimensions.
      plot_mat <- t(mat)
      row_ann <- NULL
      grp_var <- input$regionHeatmap
      if (!is.null(grp_var) && nzchar(grp_var) &&
          grp_var %in% colnames(sample_data(ps()))) {
        grps <- as.character(sample_data(ps_trans)[[grp_var]])
        row_ann <- data.frame(Group = grps, row.names = rownames(plot_mat))
      }

      # CLR is centred on zero, so it reads correctly only on a diverging scale;
      # the strictly positive transforms use a sequential one.
      hm_colors <- if (identical(input$hmTransform, "clr")) {
        grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256)
      } else {
        viridis::viridis(256)
      }

      heatmaply(plot_mat,
                row_side_colors = row_ann,
                row_side_palette = function(n) canis_colors(n),
                colors = hm_colors,
                main   = paste(input$hmTransform, "heatmap -", rank),
                xlab   = "Taxa", ylab = "Samples",
                margins = c(80, 80, 40, 40))
    })
  })
  
  # ============================================================
  # DENDROGRAM (Bray-Curtis + Ward, standardised microbiome approach)
  # ============================================================
  output$dendrogram <- renderPlotly({
    req(ps())
    input$updateDendro
    input$dark_mode
    isolate({
      otu_mat <- t(as(otu_table(ps()), "matrix"))
      # Normalise before distance calculation
      otu_norm <- vegan::decostand(otu_mat, "total")
      
      dist_mat <- if (input$distanceMethod %in% c("bray","jaccard"))
        vegan::vegdist(otu_norm, method=input$distanceMethod)
      else
        dist(otu_norm, method=input$distanceMethod)
      
      hc   <- hclust(dist_mat, method=input$clusterMethod)
      dend <- as.dendrogram(hc)
      
      # Color labels by group if selected
      grp_var <- input$dendroColorBy
      if (!is.null(grp_var) && nchar(grp_var)>0 &&
          grp_var %in% colnames(sample_data(ps()))) {
        meta_df  <- as(sample_data(ps()), "data.frame")
        grp_fac  <- as.factor(meta_df[[grp_var]])
        pal      <- setNames(canis_colors(nlevels(grp_fac)), levels(grp_fac))
        label_colors <- pal[grp_fac[hc$order]]
        dend <- dendextend::color_labels(dend, col=label_colors)
      }
      
      # Convert to ggplot via ggdendro
      dend_data <- ggdendro::dendro_data(hc, type="rectangle")
      p <- ggplot() +
        geom_segment(data=dend_data$segments,
                     aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.5) +
        geom_text(data=dend_data$labels,
                  aes(x=x, y=y, label=label), hjust=1.1, size=3, angle=90) +
        scale_y_reverse(expand=c(0.3,0)) +
        theme_minimal() +
        labs(title=paste("Sample Dendrogram |",input$distanceMethod,"distance |",
                         input$clusterMethod,"clustering"),
             x="", y="Height")
      ggplotly(p)
    })
  })
  
  output$dendroBar <- renderPlotly({
    req(ps())
    input$updateDendro
    input$dark_mode
    isolate({
      if (!input$dendroShowBar) return(NULL)
      ps_glom <- safe_tax_glom(ps(), "Phylum")
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      p <- plot_bar(ps_rel, fill="Phylum") +
        plot_theme() +
        labs(title="Phylum Composition (ordered as dendrogram)")
      ggplotly(p)
    })
  })
  
  # ============================================================
  # CORRELATION NETWORK (Spearman, real correlations)
  # ============================================================
  output$correlationNetwork <- renderForceNetwork({
    req(ps())
    input$updateNetwork
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelNetwork)
      validate(need(!is.null(rank), "That taxonomic rank is not present in this dataset."))
      ps_glom <- safe_tax_glom(ps(), rank)

      # Top N taxa
      n_taxa <- min(input_num(input$netTopN, 40), ntaxa(ps_glom))
      validate(need(n_taxa >= 3, "At least three taxa are needed to build a network."))
      top_t  <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[seq_len(n_taxa)]
      ps_top <- prune_taxa(top_t, ps_glom)
      
      # CLR-transform before correlation
      ps_clr <- tryCatch(clr_transform(ps_top),
                         error=function(e) compositional(ps_top))
      
      otu_mat <- as(otu_table(ps_clr), "matrix")   # taxa x samples
      # Spearman/Pearson correlation across taxa
      cor_mat <- cor(t(otu_mat), method=input_chr(input$corMethod, "spearman"))

      # Build edge list. Distinct taxa can share a rank label, so labels are
      # made unique or separate nodes would render identically.
      tax_labels <- as.character(tax_table(ps_top)[, rank])
      tax_labels[is.na(tax_labels) | !nzchar(tax_labels)] <- "Unclassified"
      tax_labels <- make.unique(tax_labels, sep = " ")
      diag(cor_mat) <- 0
      thr <- input_num(input$corThreshold, 0.6)
      idx <- which(abs(cor_mat) >= thr & upper.tri(cor_mat), arr.ind=TRUE)

      validate(need(nrow(idx) > 0, paste0(
        "No taxon pair reaches |r| >= ", thr,
        ". Lower the correlation threshold, or include more taxa.")))
      
      cor_vals <- cor_mat[idx]
      links <- data.frame(
        source = idx[,1] - 1,
        target = idx[,2] - 1,
        value  = abs(cor_vals),
        sign   = ifelse(cor_vals > 0, 1, 2)
      )
      
      if (!isTRUE(input$netNegEdges)) links <- links[links$sign==1,]
      
      # Node abundance as group bin
      abund   <- rowSums(as.matrix(otu_table(ps_top)))
      grp_bin <- as.integer(cut(abund, breaks=5, labels=FALSE))
      nodes   <- data.frame(name=tax_labels, group=grp_bin, stringsAsFactors=FALSE)
      
      forceNetwork(Links=links, Nodes=nodes,
                   Source="source", Target="target",
                   NodeID="name", Group="group",
                   Value="value",
                   opacity=0.85, zoom=TRUE,
                   linkDistance=120, charge=-50,
                   colourScale=JS('d3.scaleOrdinal(d3.schemeCategory10)'))
    })
  })
  
  output$networkSummary <- renderPrint({
    req(ps())
    input$updateNetwork
    input$dark_mode
    isolate({
      rank <- safe_rank(ps(), input$taxLevelNetwork)
      if (is.null(rank)) { cat("That taxonomic rank is not present in this dataset.\n"); return(invisible(NULL)) }
      ps_glom <- safe_tax_glom(ps(), rank)
      n_taxa  <- min(input_num(input$netTopN, 40), ntaxa(ps_glom))
      if (n_taxa < 3) { cat("At least three taxa are needed to build a network.\n"); return(invisible(NULL)) }
      top_t   <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[seq_len(n_taxa)]
      ps_top  <- prune_taxa(top_t, ps_glom)
      ps_clr  <- tryCatch(clr_transform(ps_top), error=function(e) compositional(ps_top))
      otu_mat <- as(otu_table(ps_clr), "matrix")
      cor_mat <- cor(t(otu_mat), method=input_chr(input$corMethod, "spearman"))
      diag(cor_mat) <- 0
      thr   <- input_num(input$corThreshold, 0.6)
      n_pos <- sum(cor_mat >=  thr) / 2
      n_neg <- sum(cor_mat <= -thr) / 2
      show_neg <- isTRUE(input$netNegEdges)
      n_shown  <- n_pos + if (show_neg) n_neg else 0

      # How many of these edges would survive a significance test? Edges are
      # selected on magnitude alone, so this says how much of the network is
      # attributable to correlation strength rather than evidence.
      cm    <- cor_mat[upper.tri(cor_mat)]
      n_obs <- ncol(otu_mat)
      tstat <- suppressWarnings(cm * sqrt((n_obs - 2) / pmax(1 - cm^2, .Machine$double.eps)))
      q     <- stats::p.adjust(2 * stats::pt(abs(tstat), df = n_obs - 2, lower.tail = FALSE), "BH")
      sel   <- if (show_neg) abs(cm) >= thr else cm >= thr
      n_sig <- sum(sel & q < 0.05, na.rm = TRUE)

      cat("=== Co-occurrence network ===\n")
      cat("Taxa analysed      :", n_taxa, "\n")
      cat("Correlation method :", input_chr(input$corMethod, "spearman"), "\n")
      cat("Threshold          : |r| >=", thr, "\n")
      cat("Positive edges     :", n_pos, "\n")
      cat("Negative edges     :", n_neg, if (!show_neg) "(hidden)" else "", "\n")
      cat("Edges shown        :", n_shown, "\n\n")
      cat("Edges are selected by correlation strength alone, with no significance test.\n")
      cat("For reference, of the", n_shown, "edges shown,", n_sig,
          "would also pass a\nBH-corrected test at q < 0.05 across all",
          length(cm), "taxon pairs.\n")
      if (n_shown > n_sig) {
        cat("Treat the remaining", n_shown - n_sig, "as exploratory rather than established.\n")
      }
    })
  })
  
  # ============================================================
  # DIFFERENTIAL ABUNDANCE (Kruskal-Wallis + BH)
  # ============================================================
  #' Differential abundance results, computed once and shared by the plot and
  #' the table so the two can never disagree.
  #'
  #' Works directly on the agglomerated abundance matrix rather than a melted
  #' data frame. psmelt() renames any sample variable that collides with a
  #' taxonomic rank (a metadata column called "Class" becomes "sample_Class"),
  #' which would otherwise silently test the taxonomy column instead of the
  #' grouping variable. Operating on the matrix also keeps one test per
  #' agglomerated taxon: tax_glom separates taxa by full lineage, so grouping by
  #' the rank label alone would re-pool distinct lineages that happen to share a
  #' name and test them as a single inflated sample.
  da_results <- reactive({
    req(ps())
    grp_var <- input$regionDA
    if (is.null(grp_var) || !nzchar(grp_var) ||
        !grp_var %in% colnames(sample_data(ps()))) return(NULL)

    rank <- safe_rank(ps(), input$taxLevelDA)
    if (is.null(rank)) return(list(error = "That taxonomic rank is not present in this dataset."))
    ps_glom <- safe_tax_glom(ps(), rank)
    ps_rel  <- transform_sample_counts(ps_glom, function(x) if (sum(x) > 0) x / sum(x) else x)

    ab <- as(otu_table(ps_rel), "matrix")
    if (!taxa_are_rows(ps_rel)) ab <- t(ab)
    meta <- as(sample_data(ps_rel), "data.frame")
    grp  <- factor(meta[[grp_var]])

    if (nlevels(grp) < 2) return(list(error = "The grouping variable needs at least two levels."))

    # Near-absent taxa cannot yield a finding but still consume FDR budget.
    min_prev <- max(3, ceiling(0.10 * ncol(ab)))
    keep <- rowSums(ab > 0) >= min_prev
    n_dropped <- sum(!keep)
    ab <- ab[keep, , drop = FALSE]
    if (!nrow(ab)) {
      return(list(error = paste0("No taxon is present in at least ", min_prev,
                                 " samples at this rank, so nothing can be tested.")))
    }

    # Repeated measures: reuse the same complete-block check the diversity tabs
    # use, so a paired design is not tested as if it were independent.
    block_var <- input$blockDA
    has_block <- !is.null(block_var) && nzchar(block_var) && block_var %in% colnames(meta)
    blk <- if (has_block) factor(meta[[block_var]]) else NULL
    test_choice <- choose_alpha_test(grp, blk, if (has_block) block_var else NULL, grp_var)
    use_block <- identical(test_choice, "friedman")

    pvals <- apply(ab, 1, function(v) {
      tryCatch({
        if (use_block) friedman.test(v, groups = grp, blocks = blk)$p.value
        else kruskal.test(v, grp)$p.value
      }, error = function(e) NA_real_)
    })

    group_means <- t(apply(ab, 1, function(v) tapply(v, grp, mean)))
    colnames(group_means) <- paste0("mean_", levels(grp))

    # Direction: which group carries the taxon, and how large the gap is.
    pseudo <- min(ab[ab > 0], na.rm = TRUE) / 2
    hi <- max.col(group_means, ties.method = "first")
    lo <- max.col(-group_means, ties.method = "first")
    log2fc <- log2((group_means[cbind(seq_len(nrow(group_means)), hi)] + pseudo) /
                     (group_means[cbind(seq_len(nrow(group_means)), lo)] + pseudo))

    fdr <- p.adjust(pvals, method = "BH")   # adjust the raw p-values, never rounded ones

    res <- data.frame(
      Taxon       = as.character(tax_table(ps_glom)[rownames(ab), rank]),
      TaxonID     = rownames(ab),
      as.data.frame(group_means, check.names = FALSE),
      Enriched_in = levels(grp)[hi],
      log2FC      = log2fc,
      Prevalence  = rowSums(ab > 0) / ncol(ab),
      p_value     = pvals,
      FDR         = fdr,
      Significant = !is.na(fdr) & fdr < input_num(input$daFDR, 0.05),
      row.names   = NULL, check.names = FALSE, stringsAsFactors = FALSE
    )
    res <- res[order(res$FDR, na.last = TRUE), , drop = FALSE]

    list(results = res, rank = rank, group_var = grp_var,
         test = if (use_block) "Friedman (blocked)" else "Kruskal-Wallis",
         blocked = use_block, block_var = if (use_block) block_var else NULL,
         n_tested = nrow(ab), n_dropped = n_dropped, min_prev = min_prev,
         n_groups = nlevels(grp))
  })

  output$daVolcano <- renderPlotly({
    req(ps())
    input$updateDA
    input$dark_mode
    isolate({
      da <- da_results()
      msg <- if (is.null(da)) "Select a grouping variable."
             else if (!is.null(da$error)) da$error else NULL
      if (!is.null(msg)) {
        return(plot_ly() %>%
                 add_annotations(text = msg, x = 0.5, y = 0.5, showarrow = FALSE) %>%
                 layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      res <- da$results
      res$neg_log_fdr <- -log10(pmax(res$FDR, 1e-300))
      res$Status <- ifelse(res$Significant, "Significant", "Not significant")
      fc_label <- if (da$n_groups > 2) "log2 fold change (highest vs lowest group)"
                  else "log2 fold change between groups"

      plot_ly(res,
              x = ~log2FC, y = ~neg_log_fdr,
              color = ~Status,
              colors = stats::setNames(c(canis_ui$primary, "#C3CDD4"),
                                       c("Significant", "Not significant")),
              text = ~paste0("<b>", Taxon, "</b>",
                             "<br>enriched in: ", Enriched_in,
                             "<br>log2FC: ", signif(log2FC, 3),
                             "<br>FDR: ", signif(FDR, 3),
                             "<br>prevalence: ", round(100 * Prevalence), "%"),
              hoverinfo = "text",
              type = "scatter", mode = "markers",
              marker = list(size = 9, opacity = 0.85,
                            line = list(width = 0.5, color = "white"))) %>%
        layout(
          title = list(text = paste0(da$test, " - ", da$rank,
                                     "  (", da$n_tested, " taxa tested)"),
                       x = 0.02, font = list(size = 14)),
          xaxis = list(title = fc_label, zeroline = TRUE, zerolinecolor = "#DFE5EA"),
          yaxis = list(title = "-log10(FDR)"),
          legend = list(orientation = "h", y = -0.18),
          shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                             y0 = -log10(input_num(input$daFDR, 0.05)),
                             y1 = -log10(input_num(input$daFDR, 0.05)),
                             line = list(color = canis_ui$danger, dash = "dash", width = 1))),
          margin = list(t = 50)
        )
    })
  })

  output$daTable <- renderDT({
    req(ps())
    input$updateDA
    input$dark_mode
    isolate({
      da <- da_results()
      if (is.null(da) || !is.null(da$error)) return(NULL)

      res <- da$results
      disp <- res
      for (nm in names(disp)[vapply(disp, is.numeric, logical(1))]) {
        disp[[nm]] <- signif(disp[[nm]], 3)
      }

      caption <- paste0(
        da$test, " on relative abundance at ", da$rank, " level, grouped by '", da$group_var, "'. ",
        da$n_tested, " taxa tested",
        if (da$n_dropped > 0) paste0(" (", da$n_dropped, " excluded: present in fewer than ",
                                     da$min_prev, " samples)") else "",
        ". Benjamini-Hochberg FDR across the tested taxa.",
        if (da$blocked) paste0(" Blocked by '", da$block_var, "'.") else ""
      )

      datatable(disp, caption = caption,
                fillContainer = TRUE,
                options = list(scrollX = TRUE, pageLength = 15),
                rownames = FALSE,
                class = "cell-border stripe") %>%
        formatStyle("Significant",
                    backgroundColor = styleEqual(TRUE, "#F7EFD9"),
                    fontWeight = styleEqual(TRUE, "bold"))
    })
  })
  
  # ============================================================
  # TRANSFORMATION EXPLORER
  # ============================================================
  output$transPCoA <- renderPlotly({
    req(ps())
    input$updateTrans
    input$dark_mode
    isolate({
      ps_use <- if (input$transMethod == "raw") ps() else
        tryCatch(microbiome::transform(ps(), input$transMethod),
                 error=function(e) ps())
      suppressMessages({
        ord <- ordinate(ps_use, method="PCoA", distance="euclidean")
      })
      p <- plot_ordination(ps_use, ord, type="samples") +
        plot_theme() +
        labs(title=paste("PCoA –", input$transMethod, "transformation"))
      grp_var <- input$transColorBy
      if (!is.null(grp_var) && nchar(grp_var)>0 &&
          grp_var %in% colnames(sample_data(ps())))
        p <- p + geom_point(aes_string(color=grp_var), size=3) +
        scale_color_canis()
      else
        p <- p + geom_point(color=canis_ui$primary, size=3)
      ggplotly(p)
    })
  })
  
  output$transHist <- renderPlotly({
    req(ps())
    input$updateTrans
    input$dark_mode
    isolate({
      ps_use <- if (input$transMethod == "raw") ps() else
        tryCatch(microbiome::transform(ps(), input$transMethod),
                 error=function(e) ps())
      vals <- as.vector(as(otu_table(ps_use), "matrix"))
      vals <- vals[is.finite(vals)]

      # CLR is centred on zero, so roughly half its values are negative and
      # dropping them would hide most of the distribution this tab exists to
      # show. Zero-bounded transforms keep the zeros visible instead, since the
      # zero-inflation is itself the point.
      n_total <- length(vals)
      n_zero  <- sum(vals == 0)
      subtitle <- if (n_zero > 0) {
        sprintf("%s values, of which %s are zero (%.1f%%)",
                format(n_total, big.mark = ","), format(n_zero, big.mark = ","),
                100 * n_zero / n_total)
      } else sprintf("%s values", format(n_total, big.mark = ","))

      plot_ly(x = vals, type = "histogram",
              marker = list(color = canis_ui$primary, opacity = 0.85,
                            line = list(width = 0.3, color = "white")),
              nbinsx = 60) %>%
        layout(title = list(text = paste0("Value distribution - ", input$transMethod,
                                          "<br><sup>", subtitle, "</sup>"),
                            x = 0.02, font = list(size = 14)),
               xaxis = list(title = "Value"),
               yaxis = list(title = "Frequency"),
               margin = list(t = 55))
    })
  })
  
}

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)
