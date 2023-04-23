# About this script-------------------------------------------------------------
#
# TODO general description
#
# R environment data and values which are commonly used in this script
# in different steps (with [species] = oryza_nivara|oryza_sativa):
#
# Step 1:
#
# - studydesign_df
# - studydesign_[species]_df
# - abundance_filepaths_[species]
# - mart_[species]
# - tx2gene_[species]
# - txi_[species]
#
# Step 2:
#
# - samples_[species]
# - dgelist_[species]
# - dge_cpm_log2_[species]_df
# - dge_cpm_log2_[species]_df_pivot
# - dgelist_filtered_[species]
# - dge_cpm_filtered_log2_[species]_df
# - dge_cpm_filtered_log2_[species]_df_pivot
# - dgelist_filtered_norm_[species]
# - dge_cpm_filtered_norm_log2_[species]
# - dge_cpm_filtered_norm_log2_[species]_df
# - dge_cpm_filtered_norm_log2_[species]_df_pivot
#
# Step 3:
#
# - groups_[species]
# - pca_res_[species]
# - pca_var_[species]
# - pca_per_[species]
# - pca_res_[species]_df

# File structure requirements:
#
# |                                    <- the project's base directory
# |
# |- data                              <- all data files
# |
# |- scripts                           <- all script files
#      |
#      |- dge-analysis-PRJCA004229.R   <- this R script
#      |
#      ...

# ------------------------------------------------------------------------------
# Step 0: Load required packages
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

# ------------------------------------------------------------------------------
# Step 1: Import and annotate the Kallisto abundance files
# ------------------------------------------------------------------------------

# Set the project's data subdirectory as the current working directory
setwd(paste0(this.path::here(), "/../data"))

# --- Read in the studydesign file

studydesign_df <- readr::read_tsv(
  "studydesign-PRJCA004229.tsv",
  col_types = "icccccfffcficfc")

studydesign_oryza_nivara_df <-
  studydesign_df |>
  dplyr::filter(Organism == "Oryza nivara")

studydesign_oryza_sativa_df <-
  studydesign_df |>
  dplyr::filter(Organism == "Oryza sativa")

# --- File paths to the Kallisto abundance.tsv files

# - of the Oryza nivara samples
abundance_filepaths_oryza_nivara <- file.path(
  studydesign_oryza_nivara_df$`Run accession`,
  "abundance.tsv")

# - of the Oryza sativa samples
abundance_filepaths_oryza_sativa <- file.path(
  studydesign_oryza_sativa_df$`Run accession`,
  "abundance.tsv")

# --- Get annotations using BiomaRt

# BioMaRt dataset for Oryza nivara
mart_oryza_nivara <- biomaRt::useMart(
  biomart = "plants_mart",
  dataset = "onivara_eg_gene",
  host = "https://plants.ensembl.org")

tx2gene_oryza_nivara <-
  biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "description"),
    mart = mart_oryza_nivara) |>
  tibble::as_tibble() |>
  dplyr::rename(
    target_id = "ensembl_transcript_id",
    gene_name = "ensembl_gene_id") |>
  dplyr::select("target_id", "gene_name")

# BioMaRt dataset for Oryza sativa
mart_oryza_sativa <- biomaRt::useMart(
  biomart = "plants_mart",
  dataset = "osativa_eg_gene",
  host = "https://plants.ensembl.org")

tx2gene_oryza_sativa <-
  biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "description"),
    mart = mart_oryza_sativa) |>
  tibble::as_tibble() |>
  dplyr::rename(
    target_id = "ensembl_transcript_id",
    gene_name = "ensembl_gene_id") |>
  dplyr::select("target_id", "gene_name")

# --- Import Kallisto transcript counts

txi_gene_oryza_nivara <- tximport::tximport(
  abundance_filepaths_oryza_nivara,
  type = "kallisto",
  tx2gene = tx2gene_oryza_nivara,
  countsFromAbundance = "lengthScaledTPM")

txi_gene_oryza_sativa <- tximport::tximport(
  abundance_filepaths_oryza_sativa,
  type = "kallisto",
  tx2gene = tx2gene_oryza_sativa,
  countsFromAbundance = "lengthScaledTPM")

# ------------------------------------------------------------------------------
# Step 2: Filter and normalize the data
# ------------------------------------------------------------------------------

# --- Some TPM statistics

# Oryza nivara
transform(
  txi_gene_oryza_nivara$abundance,
  SD = matrixStats::rowSds(txi_gene_oryza_nivara$abundance),
  AVG = rowMeans(txi_gene_oryza_nivara$abundance),
  MED = matrixStats::rowMedians(txi_gene_oryza_nivara$abundance)
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

# Oryza sativa
transform(
  txi_gene_oryza_sativa$abundance,
  SD = matrixStats::rowSds(txi_gene_oryza_sativa$abundance),
  AVG = rowMeans(txi_gene_oryza_sativa$abundance),
  MED = matrixStats::rowMedians(txi_gene_oryza_sativa$abundance)
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

# --- DGElist objects

# Oryza nivara

samples_oryza_nivara <- studydesign_oryza_nivara_df$`Sample name`
# Create and save a DGE list for Oryza nivara
dgelist_oryza_nivara <- edgeR::DGEList(txi_gene_oryza_nivara$counts)
save(dgelist_oryza_nivara, file = "dgelist_oryza_nivara")
# Get log2 'counts per million'
dge_cpm_oryza_nivara <- edgeR::cpm(dgelist_oryza_nivara)
dge_cpm_log2_oryza_nivara_df <-
  edgeR::cpm(dge_cpm_oryza_nivara, log = TRUE) |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_nivara))
# Pivot the dataframe
dge_cpm_log2_oryza_nivara_df_pivot <-
  dge_cpm_log2_oryza_nivara_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_nivara,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_nivara_1 <-
  ggplot2::ggplot(dge_cpm_log2_oryza_nivara_df_pivot) +
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
    y = "log2 expression",
    title = "Oryza nivara - Log2 Counts per Million (CPM)",
    subtitle = "unfiltered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_nivara_1

# Oryza sativa

samples_oryza_sativa <- studydesign_oryza_sativa_df$`Sample name`
# Create and save a DGE list for Oryza nivara
dgelist_oryza_sativa <- edgeR::DGEList(txi_gene_oryza_sativa$counts)
save(dgelist_oryza_sativa, file = "dgelist_oryza_sativa")
# Get log2 'counts per million'
dge_cpm_oryza_sativa <- edgeR::cpm(dgelist_oryza_sativa)
dge_cpm_log2_oryza_sativa_df <-
  edgeR::cpm(dgelist_oryza_sativa, log = TRUE) |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_sativa))
# Pivot the dataframe
dge_cpm_log2_oryza_sativa_df_pivot <-
  dge_cpm_log2_oryza_sativa_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_sativa,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_sativa_1 <-
  ggplot2::ggplot(dge_cpm_log2_oryza_sativa_df_pivot) +
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
    y = "log2 expression",
    title = "Oryza sativa - Log2 Counts per Million (CPM)",
    subtitle = "unfiltered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_sativa_1

# --- Check how many genes have low reads

# Determine how many genes have no reads at all
# Oryza nivara:
table(rowSums(dgelist_oryza_nivara$counts == 0) == length(samples_oryza_nivara))
# Oryza sativa:
table(rowSums(dgelist_oryza_sativa$counts == 0) == length(samples_oryza_sativa))

# Check how many genes have CPMs >= 1 in at least 1, 2, 3, ... samples
# Oryza nivara:
sapply(
  1L:length(samples_oryza_nivara),
  function(n) sum(rowSums(dge_cpm_oryza_nivara >= 1) >= n))
# Oryza sativa:
sapply(
  1L:length(samples_oryza_sativa),
  function(n) sum(rowSums(dge_cpm_oryza_sativa >= 1) >= n))

# --- Filter out genes with low reads (< 1 CPM in at least half of the samples)

# Oryza nivara

dgelist_filtered_oryza_nivara <- dgelist_oryza_nivara[
  rowSums(dge_cpm_oryza_nivara >= 1) >= length(samples_oryza_nivara) / 2,
  ]
dim(dgelist_oryza_nivara)
dim(dgelist_filtered_oryza_nivara)
# Get log2 'counts per million'
dge_cpm_filtered_log2_oryza_nivara_df <-
  edgeR::cpm(dgelist_filtered_oryza_nivara, log = TRUE) |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_nivara))
# Pivot the dataframe
dge_cpm_filtered_log2_oryza_nivara_df_pivot <-
  dge_cpm_filtered_log2_oryza_nivara_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_nivara,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_nivara_2 <-
  ggplot2::ggplot(dge_cpm_filtered_log2_oryza_nivara_df_pivot) +
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
    subtitle = "filtered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_nivara_2

# Oryza sativa

dgelist_filtered_oryza_sativa <- dgelist_oryza_sativa[
  rowSums(dge_cpm_oryza_sativa >= 1) >= length(samples_oryza_sativa) / 2,
]
dim(dgelist_oryza_sativa)
dim(dgelist_filtered_oryza_sativa)
# Get log2 'counts per million'
dge_cpm_filtered_log2_oryza_sativa_df <-
  edgeR::cpm(dgelist_filtered_oryza_sativa, log = TRUE) |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_sativa))
# Pivot the dataframe
dge_cpm_filtered_log2_oryza_sativa_df_pivot <-
  dge_cpm_filtered_log2_oryza_sativa_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_sativa,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_sativa_2 <-
  ggplot2::ggplot(dge_cpm_filtered_log2_oryza_sativa_df_pivot) +
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
    subtitle = "filtered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_sativa_2

# --- Normalize the data

# Oryza nivara

dgelist_filtered_norm_oryza_nivara <- edgeR::calcNormFactors(
  dgelist_filtered_oryza_nivara,
  method = "TMM")
# Get log2 'counts per million'
dge_cpm_filtered_norm_log2_oryza_nivara <-
  edgeR::cpm(dgelist_filtered_norm_oryza_nivara, log = TRUE)
dge_cpm_filtered_norm_log2_oryza_nivara_df <-
  dge_cpm_filtered_norm_log2_oryza_nivara |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_nivara))
# Pivot the dataframe
dge_cpm_filtered_norm_log2_oryza_nivara_df_pivot <-
  dge_cpm_filtered_norm_log2_oryza_nivara_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_nivara,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_nivara_3 <-
  ggplot2::ggplot(dge_cpm_filtered_norm_log2_oryza_nivara_df_pivot) +
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
plot_oryza_nivara_3

# Oryza sativa

dgelist_filtered_norm_oryza_sativa <- edgeR::calcNormFactors(
  dgelist_filtered_oryza_sativa,
  method = "TMM")
# Get log2 'counts per million'
dge_cpm_filtered_norm_log2_oryza_sativa <-
  edgeR::cpm(dgelist_filtered_norm_oryza_sativa, log = TRUE)
dge_cpm_filtered_norm_log2_oryza_sativa_df <-
  dge_cpm_filtered_norm_log2_oryza_sativa |>
  tibble::as_tibble(rownames = "geneID") |>
  dplyr::rename_with(~ c("geneID", samples_oryza_sativa))
# Pivot the dataframe
dge_cpm_filtered_norm_log2_oryza_sativa_df_pivot <-
  dge_cpm_filtered_norm_log2_oryza_sativa_df |>
  tidyr::pivot_longer(
    cols = samples_oryza_sativa,
    names_to = "sample",
    values_to = "expression")
# Plot this pivoted data
plot_oryza_sativa_3 <-
  ggplot2::ggplot(dge_cpm_filtered_norm_log2_oryza_sativa_df_pivot) +
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
plot_oryza_sativa_3

# Combine all six violin plots in one overview plot
cowplot::plot_grid(
  plot_oryza_nivara_1 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_oryza_sativa_1 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  plot_oryza_nivara_2 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_oryza_sativa_2 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  plot_oryza_nivara_3 + ggplot2::ggtitle("Oryza nivara - Log2 CPM"),
  plot_oryza_sativa_3 + ggplot2::ggtitle("Oryza sativa - Log2 CPM"),
  ncol = 2,
  rel_widths = c(2, 1),
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 8)

# ------------------------------------------------------------------------------
# Step 3: Principal component analysis (PCA)
# ------------------------------------------------------------------------------

# Oryza nivara

groups_oryza_nivara <- studydesign_oryza_nivara_df$Condition
pca_res_oryza_nivara <-
  dge_cpm_filtered_norm_log2_oryza_nivara |>
  t() |>
  prcomp(scale = FALSE, retx = TRUE)
# Eigenvalues from the PCA result
pca_var_oryza_nivara <- pca_res_oryza_nivara$sdev^2
# Percentage variance explained by each PC
pca_per_oryza_nivara <-
  round(pca_var_oryza_nivara * 100 / sum(pca_var_oryza_nivara), 1)
# Plot PC1 and PC2 against each other
pca_res_oryza_nivara_df <- tibble::as_tibble(pca_res_oryza_nivara$x)
pca_plot_oryza_nivara <-
  ggplot2::ggplot(pca_res_oryza_nivara_df) +
  ggplot2::aes(
    x = PC1,
    y = PC2,
    label = samples_oryza_nivara,
    color = groups_oryza_nivara) +
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse() +
  ggplot2::xlab(paste0("PC1 (", pca_per_oryza_nivara[1], "%", ")")) +
  ggplot2::ylab(paste0("PC2 (", pca_per_oryza_nivara[2], "%", ")")) +
  ggplot2::labs(colour = "group") +
  ggplot2::ggtitle("Oryza nivara - PCA plot") +
  ggplot2::theme_bw()
plotly::ggplotly(pca_plot_oryza_nivara)
# Make an interactive table
dge_cpm_filtered_log2_oryza_nivara_df |>
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

groups_oryza_sativa <- studydesign_oryza_sativa_df$Condition
pca_res_oryza_sativa <-
  dge_cpm_filtered_norm_log2_oryza_sativa |>
  t() |>
  prcomp(scale = FALSE, retx = TRUE)
# Eigenvalues from the PCA result
pca_var_oryza_sativa <- pca_res_oryza_sativa$sdev^2
# Percentage variance explained by each PC
pca_per_oryza_sativa <-
  round(pca_var_oryza_sativa * 100 / sum(pca_var_oryza_sativa), 1)
# Plot PC1 and PC2 against each other
pca_res_oryza_sativa_df <- tibble::as_tibble(pca_res_oryza_sativa$x)
pca_plot_oryza_sativa <-
  ggplot2::ggplot(pca_res_oryza_sativa_df) +
  ggplot2::aes(
    x = PC1,
    y = PC2,
    label = samples_oryza_sativa,
    color = groups_oryza_sativa) +
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse() +
  ggplot2::xlab(paste0("PC1 (", pca_per_oryza_sativa[1], "%", ")")) +
  ggplot2::ylab(paste0("PC2 (", pca_per_oryza_sativa[2], "%", ")")) +
  ggplot2::labs(colour = "group") +
  ggplot2::ggtitle("Oryza sativa - PCA plot") +
  ggplot2::theme_bw()
plotly::ggplotly(pca_plot_oryza_sativa)
# Make an interactive table
dge_cpm_filtered_log2_oryza_sativa_df |>
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
