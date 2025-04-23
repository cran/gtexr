## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = (Sys.getenv("RUN_VIGNETTES") != "")
)

## ----setup, message = FALSE---------------------------------------------------
library(gtexr)
library(dplyr)
library(purrr)

## -----------------------------------------------------------------------------
get_eqtl_genes("Whole_Blood")

## -----------------------------------------------------------------------------
# to retrieve the first 3 pages, with default setting of 250 items per page
1:3 |>
  map(\(page) get_eqtl_genes("Whole_Blood", page = page, .verbose = FALSE) |>
        suppressWarnings()) |>
  bind_rows()

## -----------------------------------------------------------------------------
get_variant(snpId = "rs1410858") |>
  tidyr::separate(
    col = b37VariantId,
    into = c(
      "chromosome",
      "position",
      "reference_allele",
      "alternative_allele",
      "genome_build"
    ),
    sep = "_",
    remove = FALSE
  ) |>
  select(snpId:genome_build)

## ----get-genes----------------------------------------------------------------
get_genes("CRP") |>
  select(geneSymbol, gencodeId)

## ----get-variant--------------------------------------------------------------
get_variant(snpId = "rs1410858") |>
  select(snpId, variantId)

## ----get-significant-single-tissue-eqtls--------------------------------------
gene_symbol_of_interest <- "CRP"

gene_gencodeId_of_interest <- get_genes(gene_symbol_of_interest) |>
  pull(gencodeId) |>
  suppressMessages()

gene_gencodeId_of_interest |>
  get_significant_single_tissue_eqtls() |>
  distinct(geneSymbol, gencodeId, tissueSiteDetailId)

## ----calculate-eqtls----------------------------------------------------------
variants_of_interest <- c("rs12119111", "rs6605071", "rs1053870")

variants_of_interest |>
  set_names() |>
  map(
    \(x) calculate_expression_quantitative_trait_loci(
      tissueSiteDetailId = "Liver",
      gencodeId = "ENSG00000237973.1",
      variantId = x
    )
  ) |>
  bind_rows(.id = "rsid") |>
  # optionally, reformat output - first extract genomic coordinates and alleles
  tidyr::separate(
    col = "variantId",
    into = c(
      "chromosome",
      "position",
      "reference_allele",
      "alternative_allele",
      "genome_build"
    ),
    sep = "_"
  ) |>
  # ...then ascertain alternative_allele frequency
  mutate(
    alt_allele_count = (2 * homoAltCount) + hetCount,
    total_allele_count = 2 * (homoAltCount + hetCount + homoRefCount),
    alternative_allele_frequency = alt_allele_count / total_allele_count
  ) |>
  select(
    rsid,
    beta = nes,
    se = error,
    pValue,
    minor_allele_frequency = maf,
    alternative_allele_frequency,
    chromosome:genome_build,
    tissueSiteDetailId
  )

