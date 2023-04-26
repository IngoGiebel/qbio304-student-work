# About this script-------------------------------------------------------------
#
# TODO general description
#
# R environment data and values which are commonly used in this script
# in different steps (with [species] = onivara|osativa):
#
# Step 1:
#
# - studydesign_df
# - studydesign_[species]_df
# - samples_[species]
#
# Step 2:
#
# - dgelist_[species]
#
# - cpm_log2_[species]_df
# - cpm_log2_[species]_piv
#
# - dgelist_fltr_[species]
# - cpm_fltr_log2_[species]_df
# - cpm_fltr_log2_[species]_piv
#
# - dgelist_fltr_norm_[species]
# - cpm_fltr_norm_log2_[species]
# - cpm_fltr_norm_log2_[species]_df
# - cpm_fltr_norm_log2_[species]_piv
#
# Step 3:
#
# - groups_[species]
# - pca_res_[species]
# - pca_var_[species]
# - pca_per_[species]
# - pca_res_[species]_df
#
# Step 4:
#
# - top_genes_[species]_df
# - dgelist_fltr_norm_voom_[species]
#
# Step 5:
#
# - top_genes_gostres_[species]
# - gmt_[species]
# - camera_[species]_df
#
#
# File structure requirements:
#
# |                                    <- the project's base directory
# |
# |- data                              <- all research data files
# |
# |- results                           <- research result files & figures
# |
# |- scripts                           <- all script files
#      |
#      |- dge-analysis-PRJCA004229.R   <- this R script
#      |
#      ...
#
# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------

# Provides static code analysis for R. It checks for adherence to a given style,
# identifying syntax errors and possible semantic issues
library(lintr)
# For setting the current working directory relative to the filepath
# of this script
library(this.path)
# Tools for Data Copy-Pasta. RStudio addins and R functions that make
# copy-pasting vectors and tables to text painless.
library(datapasta)
# Hadley Wickham's collection of R packages for data science
library(tidyverse)
# To combine multiple plots in one figure
library(cowplot)
# For making interactive plots
library(plotly)
# For making interactive tables
library(DT)
# For annotating the abundance data
library(biomaRt)
# For getting Kallisto results into R
library(tximport)
# For differential expression analysis. Here used for creating DGEList objects
# and for data normalization.
library(edgeR)
# R interface to the g:Profiler tools
library(gprofiler2)

# ------------------------------------------------------------------------------
# Function definitions
# ------------------------------------------------------------------------------

#' Read in the studydesign tsv-file "studydesign-PRJCA004229.tsv".
#'
#' Read in the studydesign tsv-file "studydesign-PRJCA004229.tsv" and
#' create a tibble for the overall studydesign and one tibble for each
#' of the  studydesigns restricted to Oryzy nivara and Oryza sativa.
#' Furthermore, samples vectors for the Oryza nivara and the Oryza sativa
#' samples are created.
create_studydesign_data <- function() {

  studydesign_df <<- readr::read_tsv(
    "studydesign-PRJCA004229.tsv",
    col_types = "icccccfffcficfc")

  studydesign_onivara_df <<- dplyr::filter(
    studydesign_df,
    Organism == "Oryza nivara")
  samples_onivara <<- studydesign_onivara_df$`Sample name`

  studydesign_osativa_df <<- dplyr::filter(
    studydesign_df,
    Organism == "Oryza sativa")
  samples_osativa <<- studydesign_osativa_df$`Sample name`
}

#' Import the kallisto "abundance.tsv" files.
#'
#' The "abundance.tsv" files must be in a subdirectory
#' <studydesign_df$`Run accession`> of the current working
#' directory.
#'
#' See https://bioconductor.org/packages/devel/bioc/vignettes/-
#' tximport/inst/doc/tximport.html
#'
#' Prerequisites:
#' - studydesign_onivara_df
#' - studydesign_osativa_df
import_kallisto_transcript_abundance_estimates <- function() {

  # Function definitions

  get_biomart <- function(dataset) {
    biomaRt::useMart(
      biomart = "plants_mart",
      dataset = dataset,
      host = "https://plants.ensembl.org")
  }

  get_tx2gene <- function(dataset) {
    biomaRt::getBM(
      attributes = c("ensembl_transcript_id", "ensembl_gene_id", "description"),
      mart = get_biomart(dataset)) |>
      tibble::as_tibble() |>
      dplyr::rename(
        target_id = "ensembl_transcript_id",
        gene_name = "ensembl_gene_id") |>
      dplyr::select("target_id", "gene_name")
  }

  get_abundance_filepaths <- function(studydsg) {
    file.path(studydsg$`Run accession`, "abundance.tsv")
  }

  get_txigene <- function(dataset, studydsg) {
    tximport::tximport(
      get_abundance_filepaths(studydsg),
      type = "kallisto",
      tx2gene = get_tx2gene(dataset),
      countsFromAbundance = "lengthScaledTPM")
  }

  # Create txi_gene_onivara and txi_gene_osativa

  txi_gene_onivara <<- get_txigene(
    dataset = "onivara_eg_gene",
    studydsg = studydesign_onivara_df)

  txi_gene_osativa <<- get_txigene(
    dataset = "osativa_eg_gene",
    studydsg = studydesign_osativa_df)
}

#' Plot some TPM statistics about the imported kallisto data for Oryza nivara.
plot_txi_gene_stats_onivara <- function() {
  transform(
    txi_gene_onivara$abundance,
    SD = matrixStats::rowSds(txi_gene_onivara$abundance),
    AVG = rowMeans(txi_gene_onivara$abundance),
    MED = matrixStats::rowMedians(txi_gene_onivara$abundance)
  ) |>
    ggplot2::ggplot() +
    ggplot2::aes(x = SD, y = MED) +
    ggplot2::geom_point(shape = 1, size = 0.5) +
    ggplot2::geom_smooth(method = lm) +
    ggplot2::geom_hex(show.legend = TRUE) +
    ggplot2::labs(
      x = "Standard deviation",
      y = "Median",
      title = "Oryza nivara - Transcripts per million (TPM)",
      subtitle = "unfiltered, non-normalized data"
    ) +
    ggplot2::theme_bw()
}

#' Plot some TPM statistics about the imported kallisto data for Oryza sativa.
plot_txi_gene_stats_osativa <- function() {
  transform(
    txi_gene_osativa$abundance,
    SD = matrixStats::rowSds(txi_gene_osativa$abundance),
    AVG = rowMeans(txi_gene_osativa$abundance),
    MED = matrixStats::rowMedians(txi_gene_osativa$abundance)
  ) |>
    ggplot2::ggplot() +
    ggplot2::aes(x = SD, y = MED) +
    ggplot2::geom_point(shape = 1, size = 0.5) +
    ggplot2::geom_smooth(method = lm) +
    ggplot2::geom_hex(show.legend = TRUE) +
    ggplot2::labs(
      x = "Standard deviation",
      y = "Median",
      title = "Oryza sativa - Transcripts per million (TPM)",
      subtitle = "unfiltered, non-normalized data"
    ) +
    ggplot2::theme_bw()
}

#' Create a log2 CPM tibble
cpm_log2_df <- function(data, samples) {
  edgeR::cpm(data, log = TRUE) |>
    tibble::as_tibble(rownames = "geneID") |>
    dplyr::rename_with(~ c("geneID", samples))
}

#' Pivot the (log2) CPM data frame.
pivot_cpm_df <- function(cpm_df, samples) {
  tidyr::pivot_longer(
    cpm_df,
    cols = tidyselect::all_of(samples),
    names_to = "sample",
    values_to = "expression")
}

#' Create DGElist objects and compute counts per million (CPM) and their
#' respective log2 values.
#'
#' A DGEList object holds the dataset (read counts and associated information)
#' to be analyzed by edgeR and the subsequent calculations performed on the
#' dataset.
create_deglists_and_cpms <- function() {

  # Oryza nivara

  # Create a DGElist object
  dgelist_onivara <<- edgeR::DGEList(txi_gene_onivara$counts)
  # Compute counts per million (CPM) and their respective log2 values
  cpm_onivara <<- edgeR::cpm(dgelist_onivara)
  cpm_log2_onivara_df <<- cpm_log2_df(cpm_onivara, samples_onivara)
  cpm_log2_onivara_piv <<- pivot_cpm_df(cpm_log2_onivara_df, samples_onivara)

  # Oryza sativa

  # Create a DGElist object
  dgelist_osativa <<- edgeR::DGEList(txi_gene_osativa$counts)
  # Compute counts per million (CPM) and their respective log2 values
  cpm_osativa <<- edgeR::cpm(dgelist_osativa)
  cpm_log2_osativa_df <<- cpm_log2_df(cpm_osativa, samples_osativa)
  cpm_log2_osativa_piv <<- pivot_cpm_df(cpm_log2_osativa_df, samples_osativa)
}

#' Create a plot of the CPM data.
plot_cpm <- function(cpm_piv, organism, log2, filtered, normalized) {
  ggplot2::ggplot(cpm_piv) +
    ggplot2::aes(x = sample, y = expression, fill = sample) +
    ggplot2::geom_violin(trim = FALSE, show.legend = FALSE) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      shape = 95,
      size = 10,
      color = "black",
      show.legend = FALSE) +
    ggplot2::labs(
      x = "sample",
      y = paste0(if (log2) "log2 " else "", "expression"),
      title = paste0(
        organism, " - ",
        if (log2) "Log2 " else "",
        "Counts per Million (CPM)"),
      subtitle = paste0(
        if (filtered) "filtered" else "unfiltered", ", ",
        if (normalized) "normalized" else "non-normalized")) +
    ggplot2::theme_bw()
}

#' Create filtered DGElist objects and compute log2 CPM of the filtered data.
#'
#' Low reads (less than <cpm_thr> CPM in at least <sample_perc> of the samples)
#' are filtered out.
filter_low_cpm <- function(cpm_thr, sample_perc) {

  # Oryza nivara

  # Create a DGElist object
  dgelist_fltr_onivara <<- dgelist_onivara[
    rowSums(cpm_onivara >= cpm_thr) >= length(samples_onivara) * sample_perc, ]
  # Compute counts per million (CPM) and their respective log2 values
  cpm_fltr_log2_onivara_df <<- cpm_log2_df(dgelist_fltr_onivara, samples_onivara)
  cpm_log2_onivara_piv <<- pivot_cpm_df(cpm_log2_onivara_df, samples_onivara)

  # Oryza sativa

  # Create a DGElist object
  dgelist_fltr_osativa <<- dgelist_osativa[
    rowSums(cpm_osativa >= cpm_thr) >= length(samples_osativa) * sample_perc, ]
  # Compute counts per million (CPM) and their respective log2 values
  cpm_fltr_log2_osativa_df <<- cpm_log2_df(dgelist_fltr_osativa, samples_osativa)
  cpm_log2_osativa_piv <<- pivot_cpm_df(cpm_log2_osativa_df, samples_osativa)
}


# ------------------------------------------------------------------------------
# Step 1: Import and annotate the Kallisto abundance files
# ------------------------------------------------------------------------------

# Set the project's data subdirectory as the current working directory
setwd(paste0(this.path::here(), "/../data"))

# Read in the studydesign file
create_studydesign_data()

# Import the kallisto "abundance.tsv" files
import_kallisto_transcript_abundance_estimates()

# ------------------------------------------------------------------------------
# Step 2: Filter and normalize the data
# ------------------------------------------------------------------------------

# Plot some TPM statistics about the imported kallisto data
plot_txi_gene_stats_onivara()
plot_txi_gene_stats_osativa()

# Create DGElist objects and compute counts per million (CPM) and their
# respective log2 values
create_deglists_and_cpms()

# Plot the log2 CPM data (unfiltered, non-normalized)
plot_cpm_onivara <- plot_cpm(
  cpm_piv = cpm_log2_onivara_piv,
  organism = "Oryza nivara",
  log2 = TRUE,
  filtered = FALSE,
  normalized = FALSE)
plot_cpm_onivara

plot_cpm_osativa <- plot_cpm(
  cpm_piv = cpm_log2_osativa_piv,
  organism = "Oryza sativa",
  log2 = TRUE,
  filtered = FALSE,
  normalized = FALSE)
plot_cpm_osativa

# Check how many genes have low reads

# Determine how many genes have no reads at all (in none of the samples)
table(rowSums(dgelist_onivara$counts == 0) == length(samples_onivara))
table(rowSums(dgelist_osativa$counts == 0) == length(samples_osativa))

# Determine how many genes have CPMs >= 1 in at least 1, 2, 3, ... samples
sapply(
  1L:length(samples_onivara),
  function(n) sum(rowSums(cpm_onivara >= 1) >= n))
sapply(
  1L:length(samples_osativa),
  function(n) sum(rowSums(cpm_osativa >= 1) >= n))

# Filter out genes with low reads (< 1 CPM in at least half of the samples)
filter_low_cpm(cpm_thr = 1, sample_perc = 0.5)
sprintf(
  "Oryza nivara - genes filtered out: %d out of %d",
  nrow(dgelist_onivara) - nrow(dgelist_fltr_onivara),
  nrow(dgelist_onivara))
sprintf(
  "Oryza sativa - genes filtered out: %d out of %d",
  nrow(dgelist_osativa) - nrow(dgelist_fltr_osativa),
  nrow(dgelist_osativa))

# Plot the log2 CPM data (filtered, non-normalized)
plot_cpm_fltr_onivara <- plot_cpm(
  cpm_piv = cpm_fltr_log2_onivara_piv,
  organism = "Oryza nivara",
  log2 = TRUE,
  filtered = TRUE,
  normalized = FALSE)
plot_cpm_fltr_onivara

plot_cpm_fltr_osativa <- plot_cpm(
  cpm_piv = cpm_fltr_log2_osativa_piv,
  organism = "Oryza sativa",
  log2 = TRUE,
  filtered = TRUE,
  normalized = FALSE)
plot_cpm_fltr_osativa


# --- Normalize the data

# Oryza nivara

dgelist_fltr_norm_onivara <- edgeR::calcNormFactors(
  dgelist_fltr_onivara,
  method = "TMM")
# Get log2 'counts per million'
cpm_fltr_norm_log2_onivara <-
  edgeR::cpm(dgelist_fltr_norm_onivara, log = TRUE)
cpm_fltr_norm_log2_onivara_df <-
  cpm_fltr_norm_log2_onivara |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_onivara))
# Pivot the dataframe
cpm_fltr_norm_log2_onivara_piv <-
  cpm_fltr_norm_log2_onivara_df |>
  tidyr::pivot_longer(
    cols = tidyselect::all_of(samples_onivara),
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_onivara_3 <-
  ggplot2::ggplot(cpm_fltr_norm_log2_onivara_piv) +
  ggplot2::aes(x = sample, y = expression, fill = sample) +
  ggplot2::geom_violin(trim = FALSE, show.legend = FALSE) +
  ggplot2::stat_summary(
    fun = "median",
    geom = "point",
    shape = 95,
    size = 10,
    color = "black",
    show.legend = FALSE) +
  ggplot2::labs(
    y = "log2 expression",
    x = "sample",
    title = "Oryza nivara - Log2 Counts per Million (CPM)",
    subtitle = "filtered, TMM normalized") +
  ggplot2::theme_bw()
plot_onivara_3

# Oryza sativa

dgelist_fltr_norm_osativa <- edgeR::calcNormFactors(
  dgelist_fltr_osativa,
  method = "TMM")
# Get log2 'counts per million'
cpm_fltr_norm_log2_osativa <-
  edgeR::cpm(dgelist_fltr_norm_osativa, log = TRUE)
cpm_fltr_norm_log2_osativa_df <-
  cpm_fltr_norm_log2_osativa |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_osativa))
# Pivot the dataframe
cpm_fltr_norm_log2_osativa_piv <-
  cpm_fltr_norm_log2_osativa_df |>
  tidyr::pivot_longer(
    cols = tidyselect::all_of(samples_osativa),
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_osativa_3 <-
  ggplot2::ggplot(cpm_fltr_norm_log2_osativa_piv) +
  ggplot2::aes(x = sample, y = expression, fill = sample) +
  ggplot2::geom_violin(trim = FALSE, show.legend = FALSE) +
  ggplot2::stat_summary(
    fun = "median",
    geom = "point",
    shape = 95,
    size = 10,
    color = "black",
    show.legend = FALSE) +
  ggplot2::labs(
    y = "log2 expression",
    x = "sample",
    title = "Oryza sativa - Log2 Counts per Million (CPM)",
    subtitle = "filtered, TMM normalized") +
  ggplot2::theme_bw()
plot_osativa_3

# Combine all six violin plots in one overview plot
cowplot::plot_grid(
  plot_onivara_1 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_osativa_1 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  plot_onivara_2 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_osativa_2 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  plot_onivara_3 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_osativa_3 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  ncol = 2,
  rel_widths = c(2, 1),
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 8)


# ------------------------------------------------------------------------------
# Step 3: Principal component analysis (PCA)
# ------------------------------------------------------------------------------

# Oryza nivara

groups_onivara <- studydesign_onivara_df$Condition
pca_res_onivara <-
  cpm_fltr_norm_log2_onivara |>
  t() |>
  prcomp(scale = FALSE, retx = TRUE)
# Eigenvalues from the PCA result
pca_var_onivara <- pca_res_onivara$sdev^2
# Percentage variance explained by each PC
pca_per_onivara <- round(pca_var_onivara * 100 / sum(pca_var_onivara), 1)
# Plot PC1 and PC2 against each other
pca_res_onivara_df <- tibble::as_tibble(pca_res_onivara$x)
pca_plot_onivara <-
  ggplot2::ggplot(pca_res_onivara_df) +
  ggplot2::aes(
    x = PC1,
    y = PC2,
    label = samples_onivara,
    color = groups_onivara) +
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse() +
  ggplot2::xlab(paste0("PC1 (", pca_per_onivara[1], "%", ")")) +
  ggplot2::ylab(paste0("PC2 (", pca_per_onivara[2], "%", ")")) +
  ggplot2::labs(colour = "group") +
  ggplot2::ggtitle("Oryza nivara - PCA plot") +
  ggplot2::theme_bw()
plotly::ggplotly(pca_plot_onivara)
# Make an interactive table
cpm_fltr_log2_onivara_df |>
  dplyr::mutate(
    normal.AVG = (BJ278C1 + BJ278C2 + BJ278C3 + BJ89C1 + BJ89C2 + BJ89C3) / 6,
    stress.AVG = (BJ278P1 + BJ278P2 + BJ278P3 + BJ89P1 + BJ89P2 + BJ89P3) / 6,
    LogFC = (stress.AVG - normal.AVG)) |>
  dplyr::select(geneID, normal.AVG, stress.AVG, LogFC) |>
  dplyr::mutate(
    normal.AVG = round(normal.AVG, 2),
    stress.AVG = round(stress.AVG, 2),
    LogFC = round(LogFC, 2)) |>
  dplyr::arrange(dplyr::desc(LogFC)) |>
  DT::datatable(
    extensions = c("KeyTable", "FixedHeader"),
    filter = "top",
    options = list(
      keys = TRUE,
      searchHighlight = TRUE,
      pageLength = 100000,
      lengthMenu = c("10", "25", "50", "100")))

# Oryza sativa

groups_osativa <- studydesign_osativa_df$Condition
pca_res_osativa <-
  cpm_fltr_norm_log2_osativa |>
  t() |>
  prcomp(scale = FALSE, retx = TRUE)
# Eigenvalues from the PCA result
pca_var_osativa <- pca_res_osativa$sdev^2
# Percentage variance explained by each PC
pca_per_osativa <- round(pca_var_osativa * 100 / sum(pca_var_osativa), 1)
# Plot PC1 and PC2 against each other
pca_res_osativa_df <- tibble::as_tibble(pca_res_osativa$x)
pca_plot_osativa <-
  ggplot2::ggplot(pca_res_osativa_df) +
  ggplot2::aes(
    x = PC1,
    y = PC2,
    label = samples_osativa,
    color = groups_osativa) +
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse() +
  ggplot2::xlab(paste0("PC1 (", pca_per_osativa[1], "%", ")")) +
  ggplot2::ylab(paste0("PC2 (", pca_per_osativa[2], "%", ")")) +
  ggplot2::labs(colour = "group") +
  ggplot2::ggtitle("Oryza sativa - PCA plot") +
  ggplot2::theme_bw()
plotly::ggplotly(pca_plot_osativa)
# Make an interactive table
cpm_fltr_log2_osativa_df |>
  dplyr::mutate(
    normal.AVG = (NC1 + NC2 + NC1) / 3,
    stress.AVG = (NP1 + NP2 + NP3) / 3,
    LogFC = (stress.AVG - normal.AVG)) |>
  dplyr::select(geneID, normal.AVG, stress.AVG, LogFC) |>
  dplyr::mutate(
    normal.AVG = round(normal.AVG, 2),
    stress.AVG = round(stress.AVG, 2),
    LogFC = round(LogFC, 2)) |>
  dplyr::arrange(dplyr::desc(LogFC)) |>
  DT::datatable(
    extensions = c("KeyTable", "FixedHeader"),
    filter = "top",
    options = list(
      keys = TRUE,
      searchHighlight = TRUE,
      pageLength = 100000,
      lengthMenu = c("10", "25", "50", "100")))


# ------------------------------------------------------------------------------
# Step 4: Identify differentially expressed genes (DEGs)
# ------------------------------------------------------------------------------

# Set up the design and contrast matrices
model_mtx_onivara <- model.matrix(~0 + groups_onivara)
colnames(model_mtx_onivara) <- make.names(levels(groups_onivara))
contrast_mtx_onivara <- limma::makeContrasts(
  stress = drought.stress.condition - normal.condition,
  levels = model_mtx_onivara)
model_mtx_osativa <- model.matrix(~0 + groups_osativa)
colnames(model_mtx_osativa) <- make.names(levels(groups_osativa))
contrast_mtx_osativa <- limma::makeContrasts(
  stress = drought.stress.condition - normal.condition,
  levels = model_mtx_osativa)

# Model mean-variance trend and fit a linear model to the data
dgelist_fltr_norm_voom_onivara <- limma::voom(
  dgelist_fltr_norm_onivara,
  design = model_mtx_onivara,
  plot = TRUE)
dgelist_fltr_norm_lmfit_onivara <- limma::lmFit(
  dgelist_fltr_norm_voom_onivara,
  design = model_mtx_onivara)
dgelist_fltr_norm_voom_osativa <- limma::voom(
  dgelist_fltr_norm_osativa,
  design = model_mtx_osativa,
  plot = TRUE)
dgelist_fltr_norm_lmfit_osativa <- limma::lmFit(
  dgelist_fltr_norm_voom_osativa,
  design = model_mtx_osativa)

# Get Bayesian stats for the contrasts from the linear model fits
ebayes_onivara <-
  dgelist_fltr_norm_lmfit_onivara |>
  limma::contrasts.fit(contrasts = contrast_mtx_onivara) |>
  limma::eBayes()
ebayes_osativa <-
  dgelist_fltr_norm_lmfit_osativa |>
  limma::contrasts.fit(contrasts = contrast_mtx_osativa) |>
  limma::eBayes()

# --- View DEGs - tables of "top-ranked" genes from the Bayesian stats

# Computed statistics:
# - log FC      : log2-fold-change corresponding to the contrast
# - AveExpr     : average log2-expression
# - t           : moderated t-statistic: ratio of the logFC to the standard
#                 error (where the error has been moderated across all genes
#                 because of Bayesian approach)
# - P.Value     : raw p-value
# - adj.P.Value : adjusted p-value (by BH)
# - B           : log-odds that the gene is differentially expressed

# Oryza nivara

# Determine the top-ranked genes
top_genes_onivara_df <-
  limma::topTable(
    ebayes_onivara,
    adjust.method = "BH",
    coef = 1,
    number = 100000,
    sort.by = "logFC") |>
  tibble::as_tibble(rownames = "geneID")
# Create an interactive vulcano plot
(top_genes_onivara_df |>
  ggplot2::ggplot() +
  ggplot2::aes(x = logFC, y = -log10(adj.P.Val), text = geneID) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::geom_hline(
    yintercept = -log10(0.01),
    linetype = "longdash",
    colour= "grey",
    linewidth = 1) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype ="longdash",
    colour = "coral",
    linewidth = 1) +
  ggplot2::geom_vline(
    xintercept = -1,
    linetype = "longdash",
    colour = "cadetblue",
    linewidth = 1) +
  ggplot2::annotate(
    "rect",
    xmin = 1,
    xmax = 12,
    ymin = -log10(0.01),
    ymax = 7.5,
    alpha = .2,
    fill = "coral") +
  ggplot2::annotate(
    "rect",
    xmin = -1,
    xmax = -12,
    ymin = -log10(0.01),
    ymax = 7.5,
    alpha=.2,
    fill = "cadetblue") +
  ggplot2::ggtitle("Oryza nivara - volcano plot") +
  ggplot2::theme_bw()) |>
  plotly::ggplotly()

# Oryza sativa

# Determine the top-ranked genes
top_genes_osativa_df <-
  limma::topTable(
    ebayes_osativa,
    adjust.method = "BH",
    coef = 1,
    number = 100000,
    sort.by = "logFC") |>
  tibble::as_tibble(rownames = "geneID")
# Create an interactive vulcano plot
(top_genes_osativa_df |>
    ggplot2::ggplot() +
    ggplot2::aes(x = logFC, y = -log10(adj.P.Val), text = geneID) +
    ggplot2::geom_point(size = 0.2) +
    ggplot2::geom_hline(
      yintercept = -log10(0.01),
      linetype = "longdash",
      colour= "grey",
      linewidth = 1) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype ="longdash",
      colour = "coral",
      linewidth = 1) +
    ggplot2::geom_vline(
      xintercept = -1,
      linetype = "longdash",
      colour = "cadetblue",
      linewidth = 1) +
    ggplot2::annotate(
      "rect",
      xmin = 1,
      xmax = 12,
      ymin = -log10(0.01),
      ymax = 7.5,
      alpha = .2,
      fill = "coral") +
    ggplot2::annotate(
      "rect",
      xmin = -1,
      xmax = -12,
      ymin = -log10(0.01),
      ymax = 7.5,
      alpha=.2,
      fill = "cadetblue") +
    ggplot2::ggtitle("Oryza sativa - volcano plot") +
    ggplot2::theme_bw()) |>
  plotly::ggplotly()

# --- Make a Venn diagram of the DEGs

# Oryza nivara
test_results_oryza_nivara <- limma::decideTests(
  ebayes_onivara,
  method = "global",
  adjust.method = "BH",
  p.value = 0.01,
  lfc = 7)
head(test_results_oryza_nivara)
summary(test_results_oryza_nivara)
limma::vennDiagram(test_results_oryza_nivara, include = "both")

# Oryza sativa
test_results_oryza_sativa <- limma::decideTests(
  ebayes_osativa,
  method = "global",
  adjust.method = "BH",
  p.value = 0.01,
  lfc = 7)
head(test_results_oryza_sativa)
summary(test_results_oryza_sativa)
limma::vennDiagram(test_results_oryza_sativa, include = "both")

#--- Create an interactive table and a heatmap of the DEGs

# Oryza nivara

# Extract the expression data of the DEGs
colnames(dgelist_fltr_norm_voom_onivara$E) <- samples_onivara
deg_onivara_df <-
  dgelist_fltr_norm_voom_onivara$E[
    test_results_oryza_nivara[, 1] != 0, ] |>
  as_tibble(rownames = "geneID")
# Save the DEGs as a text file
readr::write_tsv(deg_onivara_df, "../results/degs-oryza-nivara.txt")
# Create an interactive table
deg_onivara_df |>
DT::datatable(
  extensions = c("KeyTable", "FixedHeader"),
  caption = "DEGs in Oryza nivara",
  options = list(
    keys = TRUE,
    searchHighlight = TRUE,
    pageLength = 100000,
    lengthMenu = c("10", "25", "50", "100"))) |>
  DT::formatRound(columns = c(2:13), digits = 2)
# Create a heatmap
deg_onivara_df |>
  dplyr::select(!1) |>
  heatmaply::heatmaply(
    xlab = "Samples",
    ylab = "DEGs",
    main = "DEGs in Oryza nivara",
    scale = "column",
    margins = c(60, 100, 40, 20),
    grid_color = "white",
    grid_width = 0.0000001,
    titleX = TRUE,
    titleY = TRUE,
    hide_colorbar = TRUE,
    branches_lwd = 0.1,
    label_names = c("Gene", "Sample:", "Value"),
    fontsize_row = 5,
    fontsize_col = 5,
    labCol = samples_onivara,
    labRow = deg_onivara_df$geneID,
    heatmap_layers = theme(axis.line = element_blank()))

# Oryza sativa

# Extract the expression data of the DEGs
colnames(dgelist_fltr_norm_voom_osativa$E) <- samples_osativa
deg_osativa_df <-
  dgelist_fltr_norm_voom_osativa$E[
    test_results_oryza_sativa[, 1] != 0, ] |>
  as_tibble(rownames = "geneID")
# Save the DEGs as a text file
readr::write_tsv(deg_osativa_df, "../results/degs-oryza-sativa.txt")
# Create an interactive table
deg_osativa_df |>
  DT::datatable(
    extensions = c("KeyTable", "FixedHeader"),
    caption = "DEGs in Oryza sativa",
    options = list(
      keys = TRUE,
      searchHighlight = TRUE,
      pageLength = 100000,
      lengthMenu = c("10", "25", "50", "100"))) |>
  DT::formatRound(columns = c(2:7), digits = 2)
# Create a heatmap
deg_osativa_df |>
  dplyr::select(!1) |>
  heatmaply::heatmaply(
    xlab = "Samples",
    ylab = "DEGs",
    main = "DEGs in Oryza sativa",
    scale = "column",
    margins = c(60, 100, 40, 20),
    grid_color = "white",
    grid_width = 0.0000001,
    titleX = TRUE,
    titleY = TRUE,
    hide_colorbar = TRUE,
    branches_lwd = 0.1,
    label_names = c("Gene", "Sample:", "Value"),
    fontsize_row = 5,
    fontsize_col = 5,
    labCol = samples_osativa,
    labRow = deg_osativa_df$geneID,
    heatmap_layers = theme(axis.line = element_blank()))


# ------------------------------------------------------------------------------
# Step 5: Gene sets enrichment analysis (GSEA)
# ------------------------------------------------------------------------------

# --- GSEA using g:Profiler

# Oryza nivara

# Functional enrichment analysis of the 100 top-ranked genes
top_genes_gostres_onivara <- gprofiler2::gost(
  top_genes_onivara_df$geneID[1:100],
  organism = "onivara",
  correction_method = "fdr")
# Produce an interactive manhattan plot of the enriched GO terms
gprofiler2::gostplot(
  top_genes_gostres_onivara,
  interactive = TRUE,
  capped = FALSE)
# Produce a static publication quality manhattan plot
# with the first 10 top-ranked GO terms highlighted.
gprofiler2::gostplot(
  top_genes_gostres_onivara,
  interactive = FALSE,
  capped = FALSE) |>
  gprofiler2::publish_gostplot(
    highlight_terms = top_genes_gostres_onivara$result$term_id[1:10])
# Generate a table of the gost results of the first 20 top-ranked GO terms
gprofiler2::publish_gosttable(
  top_genes_gostres_onivara,
  highlight_terms = top_genes_gostres_onivara$result$term_id[1:20],
  show_columns = c("source", "term_name", "term_size", "intersection_size"))

# Oryza sativa

# Functional enrichment analysis of the 100 top-ranked genes
top_genes_gostres_osativa <- gprofiler2::gost(
  top_genes_osativa_df$geneID[1:100],
  organism = "osativa",
  correction_method = "fdr")
# Produce an interactive manhattan plot of the enriched GO terms
gprofiler2::gostplot(
  top_genes_gostres_osativa,
  interactive = TRUE,
  capped = FALSE)
# Produce a static publication quality manhattan plot
# with the first 10 top-ranked GO terms highlighted.
gprofiler2::gostplot(
  top_genes_gostres_osativa,
  interactive = FALSE,
  capped = FALSE) |>
  gprofiler2::publish_gostplot(
    highlight_terms = top_genes_gostres_osativa$result$term_id[1:10])
# Generate a table of the gost results of the first 20 top-ranked GO terms
gprofiler2::publish_gosttable(
  top_genes_gostres_osativa,
  highlight_terms = top_genes_gostres_osativa$result$term_id[1:20],
  show_columns = c("source", "term_name", "term_size", "intersection_size"))
