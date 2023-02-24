# Test Microbiomics
# Test data comes from (MacPherson al., (2018), Scientific Reports; PMID:30046129)
source("./R/Microbiomics.R")
setwd("./") # for testing purposes, should be the root directory of the repo.
mapping_file = "./data/mapping_file.tsv"
taxonomy_file = "./data/feature_table_final_normalized_L6.txt"
outdir = "./output"


# Let's start by generating a stacked barplot of all samples. The stackedBarplotsFromTaxonomyTable() function will return a ggplot object.
# The following options are mutually exclusive mapping_file=<path_to_file>, and mapping=<data.frame>
# One of the important options of stackedBarplotsFromTaxonomyTable() and the figure generating functions in the Microbiomics package in
# general is the facets=c("variable1") or facets=c("variable1", "variable2") parameter. If we only specify a 1-element vector, the figure
# will be divided in panels along the X-axis. For instance, for the following line, if we specify facets=c("Visit"), we'll obtain the 
# following figure which will plot the read counts of each taxa of each sample:
my_facets=c("Visit")
tax_level = "L6"

p_object_1 = stackedBarplotsFromTaxonomyTable(
    mapping_file=mapping_file, mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=my_facets,
    tax_level=tax_level,       pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
    selected_taxa=NULL,     relative_abundance=FALSE,    keep_most_n=20,             remove_n=NULL,          show_borders=FALSE, 
    order_bars_by_taxon=NULL,  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE,
    prefix=NULL,               remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
    defined_height=NULL,       show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
    legend_ncol=1,             verbose=0,                   range=NULL 
)
print(p_object_1)

# First, we can convert the raw counts to % values to make profiles of each sample easier to compare.
p_object_2 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file, mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=my_facets,
  tax_level=tax_level,       pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=NULL,     relative_abundance=TRUE,     keep_most_n=20,             remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon=NULL,  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE,
  prefix=NULL,               remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,       show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,             verbose=0,                   range=NULL
  #order_bars_by_taxon="k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__NULL;g__NULL;s__NULL",
)
print(p_object_2)

# Which is a little better. We can also simplify the taxonomic lineages a bit with the summarize_lineage=TRUE argument.
p_object_3 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file, mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=my_facets,
  tax_level=tax_level,       pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=NULL,     relative_abundance=TRUE,     keep_most_n=20,             remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon=NULL,  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, 
  prefix=NULL,               remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,       show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,             verbose=0,                   range=NULL
  #order_bars_by_taxon="k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__NULL;g__NULL;s__NULL",
)
print(p_object_3)

# We should also order the barplot by the abundance of a taxon. Say we want to order by the o__Clostridiales;g__Blautia taxon.
p_object_4 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file,                              mapping=NULL,                taxonomy_file,               outdir=outdir,         facets=my_facets,
  tax_level=tax_level,                                    pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=NULL,                                  relative_abundance=TRUE,     keep_most_n=20,             remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon="o__Clostridiales;g__Blautia",      side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, range=NULL, 
  prefix=NULL,                                            remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,                                    show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,                                          verbose=0
  #order_bars_by_taxon="k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__NULL;g__NULL;s__NULL",
)
print(p_object_4)

# Or g__Coprococcus
p_object_5 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file,                              mapping=NULL,                taxonomy_file,               outdir=outdir,         facets=my_facets,
  tax_level=tax_level,                                    pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=NULL,                                  relative_abundance=TRUE,     keep_most_n=20,             remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon="o__Clostridiales;g__Coprococcus",  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, range=NULL, 
  prefix=NULL,                                            remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,                                    show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,                                          verbose=0
  #order_bars_by_taxon="k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__NULL;g__NULL;s__NULL",
)
print(p_object_5)

# We may then want to focus on certain taxa. Say we would want to show only g__Coproccus, o__Enterobacteriales:Others and g__Prevotella
p_object_6 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file,                              mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=my_facets,
  tax_level=tax_level,                                    pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=c("g__Coprococcus", "o__Enterobacteriales", "g__Prevotella"),          relative_abundance=TRUE,    keep_most_n=20,         remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon="o__Clostridiales;g__Coprococcus",  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, 
  prefix=NULL,                                            remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,                                    show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,                                          verbose=0,                   range=NULL
  #order_bars_by_taxon="k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__NULL;g__NULL;s__NULL",
)
print(p_object_6)

# That's already more interesting. We can also split by groups (individuals who received the probiotic and the ones that did not).
# In that case we would include : facets=c("Visit", "Groups") 
p_object_7 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file,                              mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=c("Groups", "Visit"),
  tax_level=tax_level,                                    pretty_display=FALSE,        by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=c("g__Coprococcus", "o__Enterobacteriales", "g__Prevotella"),          relative_abundance=TRUE,    keep_most_n=20,         remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon="o__Clostridiales;g__Coprococcus",  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, 
  prefix=NULL,                                            remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,                                    show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,                                          verbose=0,                   range=NULL
)
print(p_object_7)

# Which is nice, but we should remove the empty space with pretty_display=TRUE
p_object_8 = stackedBarplotsFromTaxonomyTable(
  mapping_file=mapping_file,                              mapping=NULL,                taxonomy_file,              outdir=outdir,          facets=c("Groups", "Visit"),
  tax_level=tax_level,                                    pretty_display=TRUE,         by_average=FALSE,           summarize_lineage=TRUE, order_files=NULL, 
  selected_taxa=c("g__Coprococcus", "o__Enterobacteriales", "g__Prevotella"),          relative_abundance=TRUE,    keep_most_n=20,         remove_n=NULL,          show_borders=FALSE, 
  order_bars_by_taxon="o__Clostridiales;g__Coprococcus",  side_by_side=FALSE,          type="16S_amplicons",       png=FALSE,              pdf=FALSE, 
  prefix=NULL,                                            remove_legend=FALSE,         pretty_display_showx=FALSE, exclude_string=NULL,    defined_width=NULL,
  defined_height=NULL,                                    show_samples_on_labels=NULL, specific_color_list=NULL,   sample_order=NULL,      legend_pos="right",
  legend_ncol=1,                                          verbose=0,                   range=NULL
)
print(p_object_8)

