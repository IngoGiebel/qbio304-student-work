# About this script-------------------------------------------------------------
#
# TODO general description
#
# R environment data and values used in this script which are commonly used
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
# - dge_cpm_[species]
# - dge_cpm_log2_[species]_df
# - dge_cpm_log2_[species]_df_pivot




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
# Get 'counts per million'
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
plot_oryza_nivara_1 <- ggplot2::ggplot(dge_cpm_log2_oryza_nivara_df_pivot) +
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
    subtitle = "unfiltered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_nivara_1

# Oryza sativa

samples_oryza_sativa <- studydesign_oryza_sativa_df$`Sample name`
# Create and save a DGE list for Oryza nivara
dgelist_oryza_sativa <- edgeR::DGEList(txi_gene_oryza_sativa$counts)
save(dgelist_oryza_sativa, file = "dgelist_oryza_sativa")
# Get 'counts per million'
dge_cpm_oryza_sativa <- edgeR::cpm(dgelist_oryza_sativa)
dge_cpm_log2_oryza_sativa_df <-
  edgeR::cpm(dge_cpm_oryza_sativa, log = TRUE) |>
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
plot_oryza_sativa_1 <- ggplot2::ggplot(dge_cpm_log2_oryza_sativa_df_pivot) +
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
    subtitle = "unfiltered, non-normalized") +
  ggplot2::theme_bw()
plot_oryza_sativa_1

