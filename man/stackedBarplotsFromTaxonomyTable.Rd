% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Microbiomics.R
\name{stackedBarplotsFromTaxonomyTable}
\alias{stackedBarplotsFromTaxonomyTable}
\title{Analyze metagenome data sets taxonomy based on taxonomy result.}
\usage{
stackedBarplotsFromTaxonomyTable(
  mapping_file,
  taxonomy_file,
  outdir,
  order_files = NULL,
  tax_level = "L6",
  selected_genuses = NULL,
  relative_abundance = FALSE,
  keep_most_n = NULL,
  remove_n = NULL,
  side_by_side = FALSE,
  facets = NULL,
  pretty_display = TRUE,
  by_average = FALSE,
  summarize_lineage = TRUE,
  type = "shotgun_metagenomics",
  pdf = FALSE,
  range = NULL,
  prefix = NULL,
  remove_legend = FALSE,
  pretty_display_showx = FALSE,
  exclude_string = NULL,
  defined_width = NULL,
  defined_height = NULL,
  show_samples_on_labels = NULL,
  show_borders = TRUE,
  order_bars_by_taxon = NULL,
  specific_color_list = NULL,
  sample_order = NULL,
  legend_pos = "right",
  legend_ncol = 1,
  verbose = 0
)
}
\value{
An object containing the ggplot object, summarized taxonomy table.
}
\description{
Make sure that data frame does not contain any illegal values.
Takes taxonomy tax table in input (i.e. Legacy Qiime format - text file in .tsv format).
Either provide a selected genus values or will print 20 most abundant genus by default.
}
\examples{
plotObject <- plotTaxonFromTaxonomyTable(...);
}
