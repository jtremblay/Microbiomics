# Microbiomics
#
# This file contains functions for the R Microbiomics package.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)

#' Load mapping file
#'
#'
#' @param file_path to mapping file in .tsv format.
#' @return A data frame for which the row.names is the first column of the supplied .tsv file.
#' @examples
#' mapping = loadMappingFile(mapping_file);
#' @export data_frame
loadMappingFile <- function(mapping_file){
  mapping = data.frame(fread(mapping_file, sep="\t", colClasses = "character", header=TRUE), check.names=FALSE)
  row.names(mapping) = mapping[,1]
  colnames(mapping)[1] = "sample_id"

  # implement sanity checks.
  validated = validateDataFrame(mapping)

  return(mapping)
}

#' validateDataFrame
#'
#' @param data_frame
#' @return  Nothing. If data frame is offensive, will call stop(). Tolerated character in the input data.frame are: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.,_-
#' @examples
#' mapping = validateDataFrame(mapping);
#' @export data_frame
validateDataFrame <- function(data_frame){
  tolerated_characters = unlist(strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.,_-", ""))

  #first check colnames
  lapply(colnames(data_frame), function(x){
    y = unlist(strsplit(x, ""))
    res = !y %in% tolerated_characters
    if(isTRUE(any(res))){
      stop("The headers in your mapping file contain other characters than alpha numeric and '.', ',', '_' and '-'\nThe headers in your mapping file contains the following values:", paste(colnames(data_frame), collapse="  "))
    }
  })

  # Then check the data frame data types.
  # TODO
}

#' Add repetition variable to mapping data frame based on a specific variable (i.e. column).
#'
#' @param data_frame input mapping data.frame
#' @param string column name of that exists in input mapping data.frame
#' @return The mapping data frame specified in input, but with replicate values integrated in a new column labeled Rep.
#' @examples
#' mapping_with_reps = addRepsToMapping(mapping, a_variable)
#' @export data_frame
addRepsToMapping <- function(curr_mapping, curr_variable){
  # Generacilly add replicates

  curr_mapping = curr_mapping[order(curr_mapping[[curr_variable]]), , drop=FALSE]
  treatments = curr_mapping[[curr_variable]]
  sampleids = row.names(curr_mapping)
  replicates = c()
  last = ""
  x = 1
  for(j in 1:length(treatments)){
    curr = treatments[j]
    if(curr != last){
      replicates = c(replicates, 1)
      last = curr
      x = 2
    }else{
      replicates = c(replicates, x)
      last = curr
      x = x + 1
    }
  }

  tmp = data.frame(row.names=sampleids, Rep=replicates)
  curr_mapping$sampleid = row.names(curr_mapping)
  tmp2 = merge(curr_mapping, tmp, by="row.names")
  tmp2$Treatment.y=NULL
  tmp2$Row.names = NULL
  row.names(tmp2) = tmp2$sampleid
  tmp2$sampleid = NULL

  return(tmp2)
}

#' Add repetition variable to mapping data frame based  on two variables (i.e. columns).
#'
#'
#' @param data_frame mapping (metadata) data.frame
#' @param facet1 variable to use in x facets
#' @param facet2 variable to use in y facets
#' @return  The mapping data frame specified in input, but with an additional column containing a unique identifier based on the facet1 and facet2 variables.
#' @examples
#' mapping_with_reps = addRepsToMapping(mapping, a_variable, another_variable)
#' @export data_frame
#'
addRepsToMappingTwoVariables <- function(curr_mapping, facet1, facet2){
  # Generically add replicates

  curr_mapping$Rep = NULL
  curr_mapping = curr_mapping[order(curr_mapping[[facet2]]), , drop=FALSE]

  facets1 = unique(curr_mapping[[ facet1[1] ]])

  final_df = NULL
  k = 1
  for(curr_facet1 in facets1){
    curr_mapping1 = curr_mapping[curr_mapping[[ facet1 ]] == curr_facet1 ,]
    treatments = curr_mapping1[[facet2]]
    sampleids = row.names(curr_mapping1)
    replicates = c()
    last = ""
    x = 1
    for(j in 1:length(treatments)){
      curr = treatments[j]
      if(curr != last){
        replicates = c(replicates, 1)
        last = curr
        x = 2
      }else{
        replicates = c(replicates, x)
        last = curr
        x = x + 1
      }
    }

    tmp = data.frame(row.names=sampleids, Rep=replicates)
    curr_mapping1$sampleid = row.names(curr_mapping1)
    tmp2 = merge(curr_mapping1, tmp, by="row.names")
    tmp2$Treatment.y=NULL
    tmp2$Row.names = NULL
    row.names(tmp2) = tmp2$sampleid
    tmp2$sampleid = NULL

    if(k == 1){
      final_df = tmp2
    }else{
      final_df = rbind(final_df, tmp2)
    }

    k = k + 1
  }

  return(final_df)
}

#' Takes a data.frame and a boolean value (summarize_lineage=<boolean>). If summarize_lineage=TRUE, the taxonomy column in the data.frame will
#' be shorten to include only certain taxa.
#'
#' @param data_frame, boolean
#' @return  data_frame
#' @examples
#' data_frame_with_processed_taxa_lineages = processTaxonomicLineages(data_frame, summarize_lineages=TRUE/FALSE);
#' @export data_frame
#'
# Here artifically populate rows
processTaxonomicLineages <- function(tTaxonomy3, summarize_lineage){
  if(!is.logical(summarize_lineage)){
    stop("summarize_lineage variable needs to be a boolean (TRUE or FALSE)")
  }

  tTaxonomy3$TMP_ID = seq(1,nrow(tTaxonomy3),1)
  for(p in 1:nrow(tTaxonomy3)){
    y = strsplit(tTaxonomy3[[ "Taxon" ]][p], ";", fixed=TRUE)[[1]]

    if(length(y) < 2){
      y[2] = "p__unknown"; y[3] = "c__unknown"; y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";

    }else if(length(y) < 3){
      y[3] = "c__unknown"; y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";

    }else if(length(y) < 4){
      y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";

    }else if(length(y) < 5){
      y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";

    }else if(length(y) < 6){
      y[6] = "g__unknown"; y[7] = "s__unknown";

    }else if(length(y) < 7){
      y[7] = "s__unknown";

    }else if(length(y) < 8){

    }

    if(p == 1){
      tmp_df = data.frame(y)
      tmp_df = t(tmp_df)
      row.names(tmp_df) = tTaxonomy3[p ,c("TMP_ID")]

    }else{
      tmp_df = rbind(tmp_df, y)
      row.names(tmp_df)[p] = tTaxonomy3[p ,c("TMP_ID")]
    }
  }
  y2 = data.frame(tmp_df)

  colnames(y2) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  if(summarize_lineage == TRUE){
    if(tax_level == "L6"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,4], ";", y2[,6]))
    }else if(tax_level == "L5"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,5]))
    }else if(tax_level == "L4"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,4]))
    }else if(tax_level == "L3"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,3]))
    }else if(tax_level == "L2"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2]))
    }else if(tax_level == "L1"){
      tmp_taxons = y2[,1,drop=FALSE]
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,1]))
    }else if(tax_level == "L7"){
      tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,4], ";", y2[,6], ";", y2[,7]))
    }
    tmp_taxons = data.frame(tmp_taxons)
    colnames(tmp_taxons)[2] = "Taxon"
    taxons = unique(tmp_taxons)
    tmp = merge(tTaxonomy3, tmp_taxons, by.x="TMP_ID", by.y="TMP_ID")
    tmp$Taxon.x = NULL
    colnames(tmp)[ncol(tmp)] = "Taxon"
    tmp$TMP_ID = NULL
    tmp_colnames = colnames(tmp)[1:(ncol(tmp)-1)]
    tmp_colnames = c("Taxon", tmp_colnames)
    tmp = tmp[,tmp_colnames]
    tTaxonomy3 = tmp
  }

  return(tTaxonomy3)
}


#' Converts numbers to formatted scientific numbers.
#'
#' @param data_frame
#' @return Converts huge numbers to scientific notation.
#' @examples
#' vectorOFformattedNumbers = fancyScientific(vectorOfNumbers);
#' @export object
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  #l=12500000000
  l <- format(l, scientific=TRUE)
  l <- gsub("^(\\d)e", "\\1.0e", l)
  l <- gsub("^(\\d)\\.(\\d)\\de", "\\1.\\2e", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


#' Given a data frame object, this function will return suggested dimensions to use for axes font sizes and figure dimensions.
#'
#' @param data_frame
#' @return  An object containing the dimensions and font sizes for the plot
#' @examples
#' dimensionObject <- getImageSpecsBarplot(dataframe);
#' @export object
getImageSpecsBarplot <- function(df, xVariable, yVariable=NULL, pretty_display=TRUE, number_of_x_bars=NULL){
  #df = df2
  #xVariable = facets[2]
  #yVariable = facets[1]
  #pretty_display=TRUE

  # Define width
  width = 15
  height = 15
  length = length(unique(df[[xVariable]]))
  print(paste0("Length xVariable: ", length))
  width = 19;   fontSizeX = 3;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0; fontSizeY = 4;  stripFontSizeY = 8; hjustY = NULL; angleY =90;
  if(length > 51 ){width = 19;   fontSizeX = 3;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 50 ){width = 17;   fontSizeX = 3;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 45 ){width = 15;   fontSizeX = 3;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 40 ){width = 13;   fontSizeX = 3;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 35 ){width = 12;   fontSizeX = 4;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 30 ){width = 11;   fontSizeX = 4;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 25 ){width = 10;   fontSizeX = 4;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 20 ){width = 8;    fontSizeX = 4;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 15 ){width = 8;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 10 ){width = 7.5;  fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 9 ) {width = 7.5;  fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 8 ) {width = 7;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 7 ) {width = 6;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 6 ) {width = 6;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 5 ) {width = 5;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 4 ) {width = 4;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length < 3 ) {width = 4;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}
  if(length <= 2 ) {width = 4;    fontSizeX = 5;  stripFontSizeX = 8; angleX = 90; hjustX = NULL; vjustX=0}

  if(!is.null(yVariable)){
    length = length(unique(df[[yVariable]]))
    print(paste0("Length yVariable: ", length))
    if(length > 51 ){width = width; fontSizeY = 4;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 16}
    if(length < 50 ){width = width; fontSizeY = 4.5; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 14}
    if(length < 45 ){width = width; fontSizeY = 5;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
    if(length < 40 ){width = width; fontSizeY = 5;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
    if(length < 35 ){width = width; fontSizeY = 6;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
    if(length < 30 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 11}
    if(length < 28 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 11}
    if(length < 25 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
    if(length < 20 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
    if(length < 15 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
    if(length < 10 ){width = width; fontSizeY = 7;   stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
    if(length < 8 ) {width = width;  fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 7}
    if(length < 6 ) {width = width;  fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 6}
    if(length < 4 ) {width = width;  fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 5}
    if(length < 3 ) {width = width;  fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 4.5}
    if(length < 2 ) {width = width;  fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 4}
  }else{
    height = 4
    fontSizeY = 9
    stripFontSizeY = 9
    angleY = 0
    hjustY = 0
  }

  # Then adjust width according to the number of columns we have on the x axis.
  tmp_width = 1
  print("Number of x bars:")
  print(number_of_x_bars)
  if(!is.null(number_of_x_bars)){
    if(number_of_x_bars <= 1000 ) {tmp_width = width * 4}
    if(number_of_x_bars <= 200 )  {tmp_width = width * 3}
    if(number_of_x_bars <= 100 )  {tmp_width = width * 2}
    if(number_of_x_bars <= 80 )   {tmp_width = width * 2.10}
    if(number_of_x_bars <= 60 )   {tmp_width = width + 1.95}
    if(number_of_x_bars <= 50 )   {tmp_width = width + 1.55}
    if(number_of_x_bars <= 40 )   {tmp_width = width + 2}
    if(number_of_x_bars <= 20 )   {tmp_width = width + 1}
    if(number_of_x_bars <= 10 )   {tmp_width = width + 0.5}
    if(number_of_x_bars <= 2 )    {tmp_width = width + 0}
    width = tmp_width
  }

  legendFontSize = 6
  legendKeySize = 0.3
  number_of_taxons = length(unique(df$Taxon))
  if(number_of_taxons <= 20){ legendFontSize = 8 ; legendKeySize = 0.3}
  if(number_of_taxons <= 16){ legendFontSize = 7 ; legendKeySize = 0.3}
  if(number_of_taxons <= 14){ legendFontSize = 7 ; legendKeySize = 0.3}
  if(number_of_taxons <= 10){ legendFontSize = 7 ; legendKeySize = 0.3}
  if(number_of_taxons <= 8) { legendFontSize = 7 ; legendKeySize = 0.3}
  if(number_of_taxons <= 6) { legendFontSize = 7 ; legendKeySize = 0.3}

  tmp_width = 1
  tax_char_length = max(nchar(as.character(df$Taxon)))
  print(paste0("tax_char_length:", tax_char_length))
  if(tax_char_length <= 240 ) {tmp_width = width + 10}
  if(tax_char_length <= 220 ) {tmp_width = width + 9}
  if(tax_char_length <= 200 ) {tmp_width = width + 8}
  if(tax_char_length <= 180 ) {tmp_width = width + 7.5}
  if(tax_char_length <= 160 ) {tmp_width = width + 7.5}
  if(tax_char_length <= 150 ) {tmp_width = width + 7}
  if(tax_char_length <= 140 ) {tmp_width = width + 7}
  if(tax_char_length <= 130 ) {tmp_width = width + 6}
  if(tax_char_length <= 120 ) {tmp_width = width + 5}
  if(tax_char_length <= 110 ) {tmp_width = width + 4}
  if(tax_char_length <= 100 ) {tmp_width = width + 4}
  if(tax_char_length <= 80  ) {tmp_width = width + 3}
  if(tax_char_length <= 60  ) {tmp_width = width + 3}
  if(tax_char_length <= 40  ) {tmp_width = width + 2}
  if(tax_char_length <= 20  ) {tmp_width = width}
  if(tax_char_length <= 10  ) {tmp_width = width}
  width = tmp_width

  tmp_width = 1
  if(!is.null(yVariable)){
    y_char_length = max(nchar(as.character(df[[yVariable]])))
    #print("y_char_length:")
    #print(y_char_length)
    #print("yVariable")
    #print(yVariable)
    #print(str(df))
    if(y_char_length <= 180 ) {tmp_width = width + 7}
    if(y_char_length <= 160 ) {tmp_width = width + 6}
    if(y_char_length <= 140 ) {tmp_width = width + 5}
    if(y_char_length <= 120 ) {tmp_width = width + 4}
    if(y_char_length <= 100 ) {tmp_width = width + 3}
    if(y_char_length <= 80  ) {tmp_width = width + 2}
    if(y_char_length <= 60  ) {tmp_width = width + 1}
    if(y_char_length <= 40  ) {tmp_width = width + 0}
    if(y_char_length <= 20  ) {tmp_width = width + 0}
    if(y_char_length <= 10  ) {tmp_width = width + 0}
    width = tmp_width
  }

  tmp_height = 1
  if(!is.null(xVariable)){
    x_char_length = max(nchar(as.character(df[[xVariable]])))
    print(paste0("x_char_length:", x_char_length))
    if(x_char_length <= 180 ) {tmp_height = height + 7}
    if(x_char_length <= 160 ) {tmp_height = height + 6}
    if(x_char_length <= 140 ) {tmp_height = height + 5}
    if(x_char_length <= 120 ) {tmp_height = height + 4}
    if(x_char_length <= 100 ) {tmp_height = height + 3.5}
    if(x_char_length <= 80  ) {tmp_height = height + 3}
    if(x_char_length <= 60  ) {tmp_height = height + 2}
    if(x_char_length <= 40  ) {tmp_height = height + 1.5}
    if(x_char_length <= 20  ) {tmp_height = height + 1}
    if(x_char_length <= 10  ) {tmp_height = height + 1}
    height = tmp_height
  }

  if(pretty_display == FALSE){
    width = width * 3
  }

  result = list(width=width, height=height, fontSizeX=fontSizeX, fontSizeY=fontSizeY,
                stripFontSizeX=stripFontSizeX, stipFontSizeY=stripFontSizeY, vjustX=vjustX,
                hjustX=hjustX, hjustY=hjustY, angleX=angleX, angleY=angleY, legendFontSize=legendFontSize, legendKeySize=legendKeySize)

  return(result)
  #if(!is.null(yVariable)){
  #   length = length(unique(df[[yVariable]]))
  #   print(paste0("Length yVariable: ", length))
  #   if(length > 45 ){width = 15; fontSizeY = 5;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 14}
  #   if(length < 45 ){width = 15; fontSizeY = 5;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
  #   if(length < 40 ){width = 13; fontSizeY = 5;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
  #   if(length < 35 ){width = 12; fontSizeY = 6;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 12}
  #   if(length < 30 ){width = 11; fontSizeY = 7;  stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
  #   if(length < 25 ){width = 10; fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
  #   if(length < 20 ){width = 9;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
  #   if(length < 15 ){width = 8;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
  #   if(length < 10 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 8}
  #   if(length < 8 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 7}
  #   if(length < 6 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 6}
  #   if(length < 4 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 5}
  #   if(length < 3 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 3.5}
  #   if(length < 2 ){width = 7;  fontSizeY = 7; stripFontSizeY = 8; angleY = 0; hjustY = 0; height = 3}
  #}
}

#' Get most n abundant taxa from a taxonomy table data frame (or file).
#'
#'
#'
#' @param file taxonomic summary file in tsv format.
#' @param data_frame taxonomic summary in a data.frame format.
#' @param variables if a vector of variables (i.e. names of the mapping or mapping_file column) is specified, will return the most abundant taxa based on these variables.
#' @param mapping mapping (i.e. metadata) data.frame
#' @param mapping_file mapping_file in .tsv format. First rows = sample ids
#' @param first integer value
#' @param last integer value
#' @return  A list. list(taxa=vector_or_most_abundant_taxa, taxa_table=taxonomy_table_data_frame, mapping=mapping_data_frame)
#' @examples
#' plotObject = getMostAbundantTaxa(...);
#' You have to specify at least one of the taxonomy_file=<filepath> or taxonomy_df=<data_frame> arguments.
#' You can only specify either mapping=<data_frame> or mapping_file=<filepath>, but not both.
#' first=<int> has to be smaller than last=<int>
#' if variables=<list> is specified, at least mapping=<data_frame> or mapping_file=<file> must be specified.
#' @export data_frame
getMostAbundantTaxa <- function(taxonomy_file=NULL, taxonomy_df=NULL, variables=NULL, mapping=NULL, mapping_file=NULL, first=1, last=20){

  if(is.null(taxonomy_file) & is.null(taxonomy_df)){
    stop("You have to specify at least one of the taxonomy_file=<filepath> or taxonomy_df=<data_frame> arguments.")
  }else if(!is.null(taxonomy_file) & is.null(taxonomy_df)){
    tTaxonomy = data.frame(fread(taxonomy_file, showProgress=FALSE, sep="\t", header=TRUE), check.names=FALSE)
  }else if(is.null(taxonomy_file) & !is.null(taxonomy_df)){
    tTaxonomy = taxonomy_df
  }

  if(!is.null(mapping_file) & !is.null(mapping)){
    stop("You can only specify either mapping=<data_frame> or mapping_file=<filepath>, but not both.")
  }else if(!is.null(mapping_file)){
    mapping = loadMappingFile(mapping_file)
    validateDataFrame(mapping)
    valid_rownames = colnames(tTaxonomy)[colnames(tTaxonomy) %in% row.names(mapping)]
    tTaxonomy = tTaxonomy[,c("Taxon", valid_rownames),drop=FALSE]
  }# else mapping = mapping...

  if(first > last){
    stop("first=<int> has to be smaller than last=<int>")
  }

  if(!is.null(variables) & is.null(mapping) & is.null(mapping_file)){
    stop("if variables=<list> is specified, at least mapping=<data_frame> or mapping_file=<file> must be specified.")
  }else if(!is.null(variables) & !is.null(mapping)){
    mapping_select = NULL
    for(i in 1:length(variables)){
      colname = names(variables)[i]
      variable_value = variables[[colname]]
      print(paste0("i: ", i))
      print(paste0("colname: ", colname))
      print("variable value: "); print(variable_value);
      # filter mapping based on the variable list.
      if(i == 1){
        mapping_select = mapping[mapping[[colname]] %in% variable_value,]
      }else{
        mapping_select = mapping_select[mapping_select[[colname]] %in% variable_value,]
      }
    }
    mapping_select = mapping_select[!duplicated(mapping_select),]
    tTaxonomy = tTaxonomy[,colnames(tTaxonomy) %in% c("Taxon", row.names(mapping_select))]
  }

  tTaxonomy2 = cbind(tTaxonomy, (rowSums(tTaxonomy[2:(ncol(tTaxonomy))]))/ncol(tTaxonomy) )
  tTaxonomy2 = tTaxonomy2[order(-tTaxonomy2[, ncol(tTaxonomy2)]),,drop=FALSE]
  tTaxonomy2[,ncol(tTaxonomy2)] = NULL

  most_abundant_taxa = tTaxonomy2[(first:last),]$Taxon

  return(list(taxa=most_abundant_taxa, taxa_table=tTaxonomy2, mapping=mapping_select))
}

#' Analyze metagenome data sets taxonomy based on taxonomy result.
#'
#' Make sure that data frame does not contain any illegal values.
#' Takes taxonomy tax table in input (i.e. Legacy Qiime format - text file in .tsv format).
#' Either provide a selected genus values or will print 20 most abundant genus by default.
#'
#' @param mapping_file default=NULL. A .tsv file in which row names (except header) matches with colnames of taxonomy_file
#' @param mapping default=NULL. A data.frame obtained from loadMapping(mapping_file) function. mapping and mapping_file are mutually exclusive.
#' @param taxonomy_file default=NULL. A .tsv taxonomic summary file.
#' @param outdir default=NULL. The output directory where figures will be writte.
#' @param order_files default=NULL. Simple text file containing one value per facet included in the facets=c("variable/column1"," "variable/column2")
#' @param custom_order default=NULL. A vector containing all variablles of a column in the mapping file of each corresponding facets.
#' @param tax_level default="L6". Taxonomic depth to use to generate figures. L1=kingdom;L2=Phylym;L3=Class;L4=Order;L5=Family;L6=Genus;L7=Species
#' @param selected_genuses default=
#' @param relative_abundance default=FALSE
#' @param keep_most_n default=NULL
#' @param remove_n default=NULL
#' @param side_by_side default=FALSE
#' @param facets default=NULL
#' @param pretty_display default=TRUE
#' @param by_average default=FALSE
#' @param summarize_lineage default=TRUE
#' @param type default="metagenomics"
#' @param png default=FALSE
#' @param pdf default=FALSE
#' @param range default=NULL
#' @param prefix default=NULL
#' @param remove_legend default=FALSE
#' @param pretty_display_showx default=FALSE
#' @param exclude_string default=NULL
#' @param defined_width default=NULL
#' @param defined_height default=NULL
#' @param show_samples_on_labels default=NULL
#' @param show_borders default=TRUE
#' @param order_bars_by_taxon default=NULL
#' @param specific_color_list default=NULL
#' @param specific_palette default=NULL
#' @param sample_order default=NULL
#' @param legend_pos default="right"
#' @param legend_ncol default=1
#' @param verbose default=0
#' @param specific_order default=NULL
#' @param angle_strip_labels_x=0 default=90
#' @param angle_strip_labels_y=0 default=90
#' @param scale default="normal
#' @return  An object containing the ggplot object, summarized taxonomy table.
#' @examples
#' plotObject <- plotTaxonFromTaxonomyTable(...);
#' @export data_frame
stackedBarplotsFromTaxonomyTable <- function(
    mapping_file=NULL, mapping=NULL, taxonomy_file=NULL, outdir=NULL, order_files=NULL, custom_order=NULL, tax_level="L6",
    selected_taxa=NULL, relative_abundance=FALSE, keep_most_n=NULL, remove_n=NULL,
    side_by_side=FALSE, facets=NULL, pretty_display=TRUE, by_average=FALSE, summarize_lineage=TRUE,
    type="shotgun_metagenomics", png=FALSE, pdf=FALSE, range=NULL, prefix=NULL, remove_legend=FALSE, pretty_display_showx=FALSE,
    exclude_string=NULL, defined_width=NULL, defined_height=NULL, show_samples_on_labels=NULL, show_borders=TRUE,
    order_bars_by_taxon=NULL, specific_color_list=NULL, specific_palette=NULL, sample_order=NULL, legend_pos="right",
    legend_ncol=1, verbose=0, specific_order=NULL, angleX=90, scale="normal", angle_strip_labels_x=90, angle_strip_labels_y=0){

  # validate options. Exit if mutually exclusive options.
  if(!is.null(angle_strip_labels_x)){
    if(angle_strip_labels_x != 0 & angle_strip_labels_x != 45 & angle_strip_labels_x != 90){
      stop("angle_strip_labels_x=<int> has to be one of the three following values: 0, 45 or 90")
    }
  }

  if(!is.null(angle_strip_labels_y)){
    if(angle_strip_labels_y != 0 & angle_strip_labels_y != 45 & angle_strip_labels_y != 90){
      stop("angle_strip_labels_y=<int> has to be one of the three following values: 0, 45 or 90")
    }
  }

  if(length(facets) > 2 | length(facets) == 0){
    stop("Only 1 or 2 facets can be specified.")
  }

  if(!is.null(mapping) & !is.null(mapping_file)){
    stop("You can only specified a mapping or a mapping_file, but not both.")
  }

  if(is.null(mapping) & is.null(mapping_file)){
    stop("You have to specify at least a mapping data frame or a text (.tsv) mapping_file,")
  }

  if(!is.null(order_files) & !is.null(custom_order)){
    stop("You can only specify a custom_order vector or a text (.tsv) order_file, but not both. You can however put both arguments at NULL.")
  }

  if(!is.null(selected_taxa)){
    if(!is.vector(selected_taxa)){ stop ("selected_taxa has to be NULL or = to a vector of species/genera.")}
  }

  if(pretty_display == TRUE & length(facets) == 1) {
    warning("pretty_display will be set to FALSE, because facets contains only one element.")
    pretty_display = FALSE
  }

  if(!is.null(specific_color_list) & !is.null(specific_palette)){
    stop("Can only specify either specific_color_list or specific_palette. If neither is specified, the default palette will be used.")
  }else if(!is.null(specific_color_list)){
    if(!is.list(specific_color_list)){
      stop("specific_color_list has to be a list data type.")
    }
  }

  if(is.null(facets)){
    stop("Only 1 or 2 facets can be specified.")
  }else if(length(facets) == 1){
    facet1 = facets[1]
  }else if(length(facets) == 2){
    facet1 = facets[1]
    facet2 = facets[2]
    #if(!is.null(order_bars_by_taxon)){
    #  stop("order_bars_by_taxon can only be specified if facets < 2.")
    #}
  }else{
    stop("Can't have more than 2 facets...")
  }

  if(!is.null(specific_order)){
    #TODO make sure all the variable specified in this vector are actually in the mapping
  }

  if(scale != "normal" & scale != "log2" & scale != "log10"){
    stop("scale=<string> argument has to be set to either normal, log2 or log10")
  }

  if(!is.null(range)){
    if(!is.vector(range)){ stop("range=<vector> argument has to be a vector of 2 elements. The two elements have to be greater than 0 and the second element has to be greater than the first." )  }
    if(length(range) != 2){ stop("range=<vector> argument has to be a vector of 2 elements. The two elements have to be greater than 0 and the second element has to be greater than the first." )  }
    if (!isTRUE(all(range == floor(range)))) stop("'range' must only contain integer values")
  }

  # Default colors
  vColors = c(
    "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", "#DAA520",
    "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", "#E6E6FA", "#FF8C00", "#80008B",
    "#8FBC8F", "#00BFFF", "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC",
    "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5",
    "#FFFFFF", "#990000", "#666600", "#006666", "#330066", "#A0A0A0", "#99004C"
  )
  vColors2 = c(
    "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF",
    "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5", "#FFFFFF", "#990000", "#666600", "#006666",
    "#330066", "#A0A0A0", "#99004C"
  )
  vColors3 = unique(c(vColors, vColors2))

  # check if files exist
  if(!is.null(mapping_file)){
    if(!file.exists(mapping_file)){ stop(paste0("Mapping file does not exist: ", mapping_file)) }
  }
  if(!file.exists(taxonomy_file)){ stop(paste0("OTU table file does not exist: ", taxonomy_file)) }
  outdir1 = paste0(outdir, "/taxonomy_", type)
  dir.create(file.path(outdir1), showWarnings=TRUE, recursive=TRUE)

  # Load files
  if(is.null(mapping)){
    mapping = loadMappingFile(mapping_file)
    validateDataFrame(mapping)
  }else{
    if(nrow(mapping) == 0 | ncol(mapping) == 0){
      stop("mapping file looks empty...")
    }
  }

  for(i in 1:length(facets)){
    if(any(colnames(mapping) %in% facets[i]) == FALSE) { # any() returns false if all of the values are false.
      found = FALSE
      stop(paste0(facets[i], " not in mapping file..."))
    }
  }

  # Generically add replicates
  mapping$Rep = NULL
  if(pretty_display == TRUE){
    mapping = addRepsToMappingTwoVariables(mapping, facets[1], facets[2])
  }else{
    mapping = addRepsToMapping(mapping, facet1)
  }
  print(head(mapping))

  tTaxonomy = data.frame(fread(taxonomy_file, showProgress=FALSE, sep="\t", header=TRUE), check.names=FALSE)
  if(is.null(order_files)){
    print("No order_files provided.")
  }else{
    for(order_file in order_files){
      if(!file.exists(order_file)){ stop(paste0("Order file does not exist: ", order_file)) }
    }
  }

  if(!is.null(exclude_string)){
    print("Exclude string:")
    print(exclude_string)
    tTaxonomy = tTaxonomy[!grepl(exclude_string, tTaxonomy$Taxon),]
    prefix = paste0(prefix, "minusExcluded")
  }

  if(!is.null(keep_most_n) && !is.null(range)){
    stop("keep_most_n=int> and range=c(<posint>, <posint>) are mutually exclusive.")
  }

  if(is.null(keep_most_n) && is.null(range)){
    #keep_most_n = nrow(tTaxonomy);
    keep_most_n = 20;
    if(verbose == 1){print(paste0("nrow(tTaxonomy):", nrow(tTaxonomy))); print(dim(tTaxonomy))}
  }else if(is.null(keep_most_n) && !is.null(range)){
    keep_most_n = range[2]
  }
  #print("keep_most_n: ");print(keep_most_n);

  if(verbose == 1){print(head(tTaxonomy)); print(dim(tTaxonomy))}
  valid_rownames = colnames(tTaxonomy)[colnames(tTaxonomy) %in% row.names(mapping)]
  tTaxonomy = tTaxonomy[,c("Taxon", valid_rownames),drop=FALSE]

  if(verbose == 1){print(head(tTaxonomy)); print(dim(tTaxonomy)); print(valid_rownames)}
  tTaxonomy2 = cbind(tTaxonomy, (rowSums(tTaxonomy[2:(ncol(tTaxonomy))]))/ncol(tTaxonomy) )
  tTaxonomy2 = tTaxonomy2[order(-tTaxonomy2[, ncol(tTaxonomy2)]),,drop=FALSE]
  tTaxonomy2[,ncol(tTaxonomy2)] = NULL

  if(relative_abundance == TRUE){
    # Before melt, convert absolute counts into relative counts if necessary.
    tTaxonomy2Perc = prop.table(data.matrix(tTaxonomy2[,2:ncol(tTaxonomy2)]), margin=2)*100
    tTaxonomy2Perc = data.frame(Taxon=tTaxonomy2$Taxon, tTaxonomy2Perc, check.names=FALSE)
    tTaxonomy2 = tTaxonomy2Perc
  }
  if(nrow(mapping) == 1){
    colnames(tTaxonomy2)[2] = row.names(mapping)[1]
  }

  if(is.null(prefix)){
    prefix = ""
  }

  if(keep_most_n > nrow(tTaxonomy2)){
    keep_most_n = nrow(tTaxonomy2)
  }

  if(is.null(selected_taxa)){
    tTaxonomy3 = tTaxonomy2[(1:as.numeric(keep_most_n)),]
  }else{
    print(selected_taxa)
    selected_taxa = paste(selected_taxa, collapse='|')
    tTaxonomy3 = tTaxonomy2[grepl(selected_taxa, tTaxonomy2$Taxon),]

    if(keep_most_n > nrow(tTaxonomy3)){
      keep_most_n = nrow(tTaxonomy3)
    }

    prefix = paste0(prefix, "selectedGenuses")
    tTaxonomy3 = tTaxonomy3[(1:as.numeric(keep_most_n)),]

  }

  if(!is.null(range)){
    #start = as.numeric(do.call(rbind, str_split(range, ":"))[,c(1)])
    #end = as.numeric(do.call(rbind, str_split(range, ":"))[,c(2)])
    start = range[1]
    end = range[2]

    if(start > nrow(tTaxonomy3)){
      stop(paste0("start position:", start, " is higher than number of rows of the tax table ", nrow(tTaxonomy3)))
    }
    if(end > nrow(tTaxonomy3)){
      end = nrow(tTaxonomy3)
    }

    tTaxonomy3 = tTaxonomy3[(start:end),]
  }

  if(is.null(remove_n)){
    print("Do not remove most n abundant.")
  }else{
    tTaxonomy3 = tTaxonomy3[(remove_n+1):nrow(tTaxonomy3),]
    prefix = paste0("less_", remove_n)
  }

  if(relative_abundance == FALSE){prefix = paste0(prefix, "_absolute")}else{prefix = paste0(prefix, "_relative")}

  tTaxonomy3$Taxon = gsub(";\\s+", ";", tTaxonomy3$Taxon)
  tTaxonomy3$Taxon = gsub("\\s+;", ";", tTaxonomy3$Taxon)

  tTaxonomy3 = processTaxonomicLineages(tTaxonomy3, summarize_lineage=summarize_lineage) #summarize_lineage can be TRUE or FALSE for now.

  # Reorder just to make sure
  tTaxonomy3 = cbind(tTaxonomy3, (rowSums(tTaxonomy3[2:(ncol(tTaxonomy3))]))/ncol(tTaxonomy3) )
  tTaxonomy3 = tTaxonomy3[order(-tTaxonomy3[, ncol(tTaxonomy3)]),]
  tTaxonomy3[,ncol(tTaxonomy3)] = NULL
  order_taxa = unique(tTaxonomy3$Taxon)

  df = reshape2::melt(tTaxonomy3)
  df = merge(df, mapping, by.x="variable", by.y="row.names")
  colnames(df)[2] = "Taxon"
  df$Taxon = factor(df$Taxon, levels=rev(order_taxa))

  if(is.null(order_files)){k = 1}else{k = length(order_files)}
  for(i in 1:k){
    curr_order = NULL
    if(!is.null(order_files)){
      curr_order_file = order_files[i]
      order_name = basename(file_path_sans_ext(curr_order_file))
      curr_order = data.frame(fread(curr_order_file, header=FALSE, sep="\t"))$V1
      if(verbose == 1){print("CURR_ORDER:"); print(curr_order);}
    }else{
      order_name = "all_samples"
      if(length(facets) == 1){
        curr_order = sort(unique(mapping[[ facets[1] ]]))
      }else{
        curr_order = sort(unique(mapping[[ facets[2] ]]))
      }
    }

    # Get appropriate font size and image dimensions.
    if(verbose == 1){print(paste0("length facets:", length(facets)))}
    prefix2 = NULL
    if(is.null(facets) | length(facets) == 1 ){
      if(is.null(facets)){
        facets = c("Treatment")
      }
      if(length(facets) == 1){
        facets = facets
      }
      mapping2 = mapping[mapping[[ facet1 ]] %in% curr_order,]
      df2 = df[df[[ facet1 ]] %in% curr_order,]
      df2[[ facet1 ]] = factor(df2[[ facet1 ]], levels=curr_order)

      if(is.null(order_files)){
        curr_order_facet = as.character(unique(df2[[ facets[1] ]]))
        curr_order_facet = curr_order_facet[order(nchar(curr_order_facet), curr_order_facet)]
        df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order_facet)
      }
      colnames(df2)[1] = "variable2"
      prefix2 = facets[1]

      #######################
      # Order bars by taxon #
      #######################
      df_sorted_taxon = df2 %>%
        group_by(Taxon, variable2) %>%
        dplyr::summarise(total=sum(value)) %>%
        as.data.frame()

      if(is.null(order_bars_by_taxon)){
        order_bars_by_taxon = tTaxonomy3$Taxon[1]
      }
      if(!order_bars_by_taxon %in% as.character(df_sorted_taxon$Taxon)){
        stop(paste0(order_bars_by_taxon, " was not found in the input data..."))
      }
      df_sorted_taxon = df_sorted_taxon[df_sorted_taxon$Taxon == order_bars_by_taxon,]
      sorted_variables = df_sorted_taxon[order(-df_sorted_taxon$total),]$variable2
      df2$variable2 = factor(df2$variable2, levels=sorted_variables)

      number_of_x_bars = length(unique(df2$variable2))
      specs = getImageSpecsBarplot(df2, facets[1], number_of_x_bars=number_of_x_bars)

    }else{
      if(pretty_display == TRUE & length(facets) == 2){
        mapping2 = mapping[mapping[[ facets[1] ]] %in% curr_order,]
        df2 = df[df[[ facets[2] ]] %in% curr_order,]
        df2[[ facets[2] ]] = factor(df2[[ facets[2] ]], levels=curr_order)

        # Here the key is to create a variable that contains the x facet and their replicates Rep column.
        df2$variable2 = paste0(df2[[ facets[2] ]], "-", df2$Rep)

        if(is.null(order_files)){
          curr_order_facet = as.character(unique(df2[[ facets[2] ]]))
          curr_order_facet = curr_order_facet[order(nchar(curr_order_facet), curr_order_facet)]
          df2[[ facets[2] ]] = factor(df2[[ facets[2] ]], levels=curr_order_facet)

          curr_order_1 = as.character(unique(df2[[ facets[1] ]]))
          curr_order_1 = curr_order_1[order(nchar(curr_order_1), curr_order_1)]
          df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order_1)

        }else{
          print("curr_order - pretty_display=TRUE, length=2 - with order file provided:")
          print(curr_order)
          df2[[ facets[2] ]] = factor(df2[[ facets[2] ]], levels=curr_order)
          df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order)
        }
        number_of_x_bars = length(unique(df2$variable2))
        specs = getImageSpecsBarplot(df2, facets[2], facets[1], pretty_display=TRUE, number_of_x_bars=number_of_x_bars)

        # Also order samples
        order_variable2 = as.character(unique(df2$variable2))
        print(order_variable2)
        order_variable2 = order_variable2[order(nchar(order_variable2), order_variable2)]
        df2$variable2 = factor(df2$variable2, levels=order_variable2)

        # if order by taxa
        # Here, my first attempt was to order by the first Y facet (facets[1]) of each X facet (facets[2])
        # but it creates problems if there are missing samples in one of the y facets...
        sorted_variables_based_on_taxa_order = c("")
        if(is.null(order_bars_by_taxon)){
          order_bars_by_taxon = tTaxonomy3$Taxon[1]
        }
        facets2_variables = unique(df2[[ facets[2] ]])
        variables2 =  unique(df2$variables2)
        for(n in 1:length(facets2_variables)){
          df2_tmp = df2[df2[[ facets[2] ]] == facets2_variables[n],]
          df_sorted_taxon = df2_tmp %>%
            group_by(Taxon, variable2) %>%
            dplyr::summarise(total=sum(value)) %>%
            as.data.frame()

          if(!order_bars_by_taxon %in% as.character(df_sorted_taxon$Taxon)){
            stop(paste0(order_bars_by_taxon, " was not found in the input data..."))
          }
          # TODO add an argument allowing to chose on which Y (facets[1]) facet to sort by the selected taxa.
          # for the moment only sort based on the first Y facet using variable2.
          df_sorted_taxon = df_sorted_taxon[df_sorted_taxon$Taxon == order_bars_by_taxon,]
          print(head(df_sorted_taxon))
          sorted_variables = as.character(df_sorted_taxon[order(-df_sorted_taxon$total),]$variable2)
          print(sorted_variables)
          sorted_variables_based_on_taxa_order = c(sorted_variables_based_on_taxa_order, sorted_variables)
        }
        df2$variable2 = factor(df2$variable2, levels=sorted_variables_based_on_taxa_order)
        #}

        if(verbose == 1){
          print("checkpoint-1")
          print(head(df2))
          print(df2$Rep)
        }

        ## If show sample names labels:
        if(!is.null(show_samples_on_labels)){
          n = 1
          curr_order_y = curr_order[curr_order %in% unique(mapping[[ facets[2] ]])]
          for(u in 1:length(curr_order_y)){
            df2$pos[df[[ facets[2] ]] %in% curr_order_y[u]] = n*0.5
            n = n + 1
          }
        }

      }else if(length(facets) == 2){
        mapping2 = mapping[mapping[[ facet2 ]] %in% curr_order,]
        df2 = df[df[[ facet2 ]] %in% curr_order,]
        df2[[ facet2 ]] = factor(df2[[ facet2 ]], levels=curr_order)
        specs = getImageSpecsBarplot(df2, facets[2], facets[1], pretty_display=FALSE)
        if(is.null(order_files)){
          df2 = df2[with(df2, order(as.character( df2[[facets[2]]] ))),]
          df2$variable2 = df2$variable
          curr_order_facet = as.character(unique(df2$variable2))
          df2$variable2 = factor(df2$variable2, levels=curr_order_facet)

          curr_order_1 = as.character(unique(df2[[ facets[1] ]]))
          curr_order_1 = curr_order_1[order(nchar(curr_order_1), curr_order_1)]
          df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order_1)

        }else{
          df2[[ facets[2] ]] = factor(df2[[ facets[2] ]], levels=curr_order)
          df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order)
          df2$variable2 = df2$variable
        }

        # try to eliminate gaps.
        new_variable2_order = NULL
        new_taxa_order = NULL
        df2$variable2 = df2$variable
        facets2_variables = unique(df2[[ facets[2] ]])
        facets1_variables =  unique(df2[[ facets[1] ]])
        l = 1
        for(m in 1:length(facets2_variables)){
          for(n in 1:length(facets1_variables)){
            tmp = df2
            tmp$Taxon = NULL; tmp$value = NULL;
            tmp = tmp[!duplicated(tmp),]
            tmp$variable2 = as.character(tmp$variable2)
            tmp = tmp[tmp[[ facets[1] ]] == facets1_variables[n] & tmp[[ facets[2] ]] == facets2_variables[m], ]
            curr_variable2_order = tmp[order(tmp[[ facets[2] ]]),]$variable2
            if(l == 1){
              new_variable2_order = curr_variable2_order
            }else{
              new_variable2_order = c(new_variable2_order, tmp[order(tmp[[ facets[2] ]]),]$variable2)
            }

            #order by taxa
            if(!is.null(order_bars_by_taxon)){
              variable2 = facets[2]
              df2_tmp = df2[df2$variable %in% curr_variable2_order,]
              df_sorted_taxon = df2_tmp %>%
                group_by(Taxon, variable2) %>%
                dplyr::summarise(total=sum(value)) %>%
                as.data.frame()

              if(!order_bars_by_taxon %in% as.character(df_sorted_taxon$Taxon)){
                stop(paste0(order_bars_by_taxon, " was not found in the input data..."))
              }
              df_sorted_taxon = df_sorted_taxon[df_sorted_taxon$Taxon == order_bars_by_taxon,]
              sorted_variables = df_sorted_taxon[order(-df_sorted_taxon$total),]$variable2
              if(l == 1){
                new_taxa_order = sorted_variables
              }else{
                new_taxa_order = c(new_taxa_order, sorted_variables)
              }

            }
            l = l + 1
          }
        }
        df2$variable2 = factor(df2$variable2, levels=new_variable2_order)
        df2$variable2 = factor(df2$variable2, levels=new_taxa_order)
      }
      if(length(facets) == 2){
        prefix2 = paste0(facets[1], "_", facets[2])
      }

    }
    #df2_orig = df2 #back for develepment purposes

    if(by_average == TRUE){
      # First get samples in this category for which we want the average.
      sample_id = data.frame(row.names(mapping))
      avg_values = data.frame(mapping[, facets])
      avg_values2 = data.frame(mapping[, c("Treatment")])
      colnames(avg_values2)[1] = "ORDER"

      avg = cbind(sample_id, avg_values, avg_values2)
      row.names(avg) = avg[,1]
      colnames(avg)[1] = "SampleID"
      # Then here sort the table according to ORDER
      avg = avg[order(avg$ORDER), ] # TODO smart order

      if(length(facets) == 1){
        colnames(avg)[2] = facets[1]
        avg = cbind(avg, gsub(" ", " ", paste0(avg[[facets[1]]], " ")))
      }else{
        avg = cbind(avg, gsub(" ", " ", paste(avg[[facets[1]]], avg[[facets[2]]])))
      }

      colnames(avg)[ncol(avg)] = "combined"
      avg_factors = unique(avg$combined)
      print(paste0("avg_factors: ", avg_factors))

      finalAverageTable = NULL
      finalAverageTable = data.frame(tTaxonomy3$Taxon)
      finalSDTable = NULL
      finalSDTable = data.frame(tTaxonomy3$Taxon)
      for(i in avg_factors){
        curr_samples = avg[avg$combined %in% i, 1]
        tDataAvg = data.frame(tTaxonomy3[, names(tTaxonomy3) %in% curr_samples], check.names=FALSE)

        # Here just to determine if mean is to be computed. if col = 1 , only 1 value, no means.
        if(ncol(tDataAvg) == 1){
          finalAverageTable = cbind(finalAverageTable, tDataAvg[,1])
          finalSDTable = cbind(finalSDTable, rep(0, nrow(tDataAvg)) )
        }else{
          finalAverageTable = cbind(finalAverageTable, rowMeans(tDataAvg))
          colnames(finalAverageTable)[ncol(finalAverageTable)] = i

          finalSDTable = cbind(finalSDTable,  apply(tDataAvg,1,sd) )
          colnames(finalSDTable)[ncol(finalSDTable)] = i
        }
      }

      #Process average table
      colnames(finalAverageTable)[1] = "Taxon"
      avg_factors = c("Taxon", avg_factors)
      colnames(finalAverageTable) = avg_factors
      ordered_factors = unique(avg_factors)
      finalAverageTableTest = finalAverageTable[,ordered_factors]
      finalAverageTable2 = melt(finalAverageTableTest)

      #Process sd table
      colnames(finalSDTable)[1] = "Taxon"
      colnames(finalSDTable) = avg_factors
      ordered_factors = unique(avg_factors)
      finalSDTableTest = finalSDTable[,ordered_factors]
      finalSDTable2 = melt(finalSDTableTest)

      colnames(finalSDTable2)[3] = "sd"

      if(length(facets) == 1){
        uncombined = data.frame(finalAverageTable2$variable)
        colnames(uncombined) = facets[1]
      }else{
        print(head(finalAverageTable2))
        uncombined = data.frame(do.call('rbind', strsplit(as.character(finalAverageTable2$variable),' ',fixed=TRUE)))
        colnames(uncombined) = c(facets[1], facets[2])
      }
      finalAverageTable3 = cbind(finalAverageTable2, uncombined)
      finalAverageTable3$variable = gsub("&&&", " & ", finalAverageTable3$variable)

      if(length(facets) == 1){
        finalAverageTable3[,ncol(finalAverageTable3)] = factor(finalAverageTable3[,ncol(finalAverageTable3)], levels=unique(finalAverageTable3[,ncol(finalAverageTable3)]))
      }else{
        finalAverageTable3[,ncol(finalAverageTable3)] = factor(finalAverageTable3[,ncol(finalAverageTable3)], levels=unique(finalAverageTable3[,ncol(finalAverageTable3)]))
        finalAverageTable3[,ncol(finalAverageTable3)-1] = factor(finalAverageTable3[,ncol(finalAverageTable3)-1], levels=unique(finalAverageTable3[,ncol(finalAverageTable3)-1]))
      }
      df2 = finalAverageTable3
      df2$Taxon = factor(df2$Taxon, levels=rev(order_taxa))

      # For pretty display
      if(pretty_display == TRUE){
        specs = getImageSpecsBarplot(df2, facets[2], facets[1])
      }else{
        specs = getImageSpecsBarplot(df2, facets[1])
      }

      # Here the key is to create a variable that contains the x facet and their replicates Rep column.
      if(pretty_display == TRUE & length(facets) == 2){
        df2$variable2 = paste0(df2[[ facets[2] ]], "-", df2$Rep)
        df2[[ facets[2] ]] = factor(df2[[ facets[2] ]], levels=curr_order)

        if(is.null(order_files)){
          curr_order_1 = as.character(unique(df2[[ facets[1] ]]))
          curr_order_1 = curr_order_1[order(nchar(curr_order_1), curr_order_1)]
          df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order_1)
        }

      }else{
        df2 = cbind(df2, finalSDTable2$sd)
        colnames(df2)[ncol(df2)] = "sd"
        df2$variable2 = df2$variable
      }

      if(length(facets) == 1){
        df2[[ facets[1] ]] = as.character(df2[[ facets[1] ]])
        df2[[ facets[1] ]] = gsub(" ", "", df2[[ facets[1] ]])
        df2[[ facets[1] ]] = factor(df2[[ facets[1] ]], levels=curr_order)
      }
    }

    if(!is.null(specific_color_list)){
      taxa_in_list = as.character(stack(specific_color_list)[,2])
      df2 = df2[df2$Taxon %in% taxa_in_list,]
      print(taxa_in_list)
      print(head(df2))
    }

    if(!is.null(sample_order) & length(facets) == 1){
      df2$variable2 = factor(df2$variable2, levels=sample_order)
    }

    df2$Taxon = factor(df2$Taxon, levels=rev(order_taxa))
    str(df2)

    print(specs)
    print(head(df2, n=5))
    vColors2 = vColors[1:length(order_taxa)]

    if(scale == "log2"){
      df2$value = log2(df2$value + 1)
    }else if(scale == "log10"){
      df2$value = log10(df2$value + 1)
    }
    #print("ordered by 2:"); print(sorted_variables);
    #df2$variable2 = factor(df2$variable2, levels=sorted_variables)
    .e <- environment()
    if(is.null(angle_strip_labels_x)){
      angle_strip_labels_x = specs$angleX
    }
    if(is.null(angle_strip_labels_y)){
      angle_strip_labels_y = specs$angleY
    }
    p <- ggplot(environment=.e, data=df2, aes(x=variable2, y=value, fill=Taxon)) +
      xlab("") +
      theme(
        text=element_text(family="Helvetica"),
        panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=specs$fontSizeY, colour="black"),
        axis.title=element_text(family="Helvetica", size=(specs$fontSizeY+2)),
        plot.title = element_text(lineheight=1.2, face="bold", size=18),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        legend.key.size = unit(specs$legendKeySize, "cm"),
        legend.text = element_text(size=specs$legendFontSize, face="plain"),
        legend.title = element_text(size=(specs$legendFontSize+1), face="bold"),
        legend.spacing = unit(1, "cm"),
        legend.position=legend_pos,
        strip.text.x = element_text(angle=angle_strip_labels_x, hjust=0.5, vjust=0, size=specs$stripFontSizeX, face="bold"),
        strip.text.y = element_text(angle=angle_strip_labels_y, hjust=specs$hjustY, size=specs$stripFontSizeY, face="bold"),
        strip.background =  element_blank()
      )

    if(angle_strip_labels_x != 0 & angle_strip_labels_x != 45){
      p = p + theme(strip.text.x = element_text(angle=angle_strip_labels_x, hjust=0, vjust=0.5, size=specs$stripFontSizeX, face="bold"))
    }else{
      p = p + theme(strip.text.x = element_text(angle=angle_strip_labels_x, hjust=0.5, vjust=0, size=specs$stripFontSizeX, face="bold"))
    }

    if(relative_abundance == TRUE){scale_y_continuous()}else{scale_y_continuous(labels=fancy_scientific)}

    if(is.null(specific_color_list) & is.null(specific_palette)){
      p = p + scale_fill_manual(values=rev(vColors2))
    }else{
      if(!is.null(specific_palette)){
        p = p + scale_fill_manual(values=specific_palette)
      }else if(!is.null(specific_color_list)){
        p = p + scale_fill_manual(values=specific_color_list)
      }
    }

    # Separately specify ylab
    if(relative_abundance == TRUE){ylabel = "Relative abundance"}else{ylabel = "Abundance"}
    p = p + ylab(ylabel)

    # Add facets if needed
    if(length(facets) == 1){
      p = p +  facet_grid(paste0(". ~ ", facets[1]), scales="free_x", space="free_x")
    }else{
      p = p + facet_grid(paste0(facets[1]," ~ ", facets[2]), scales="free_x", space="free_x")
    }

    if(side_by_side == TRUE){
      if(by_average==TRUE){
        p = p + geom_bar(colour="black", size=0.3, stat="identity", show.legend=TRUE, position="dodge")
        p = p + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.25, position=position_dodge(0.9))
      }else{
        p = p + geom_bar(colour="black", size=0.3, stat="identity", show.legend=TRUE, position=position_dodge())
      }
    }else{
      if(show_borders == TRUE){
        p = p + geom_bar(stat="identity") + geom_bar(colour="black", size=0.3, stat="identity", show.legend=FALSE)
      }else{
        p = p + geom_bar(stat="identity") + geom_bar(colour="black", size=0, stat="identity", show.legend=FALSE)
      }
    }

    if(pretty_display == FALSE){
      if(by_average == FALSE){
        p = p + theme(axis.text.x=element_text(size=specs$fontSizeX, colour="black", angle=as.numeric(90), hjust=1))
      }else{
        p = p + theme(axis.text.x=element_text(size=specs$fontSizeX, colour="black", angle=as.numeric(90), hjust=1))
      }
    }else if(pretty_display == TRUE & pretty_display_showx == TRUE){
      p = p + theme(axis.text.x=element_text(size=specs$fontSizeX, colour="black", angle=as.numeric(90), hjust=1))
    }else if(pretty_display == TRUE & pretty_display_showx == FALSE){
      p = p + theme(axis.text.x=element_text(size=0, colour="black", angle=as.numeric(90), hjust=1))
    }

    if(remove_legend == TRUE){
      p = p + theme(legend.position="none")
      prefix2 = paste0(prefix2, "_nolegend")
    }

    if(!is.null(show_samples_on_labels)){
      p = p + geom_text(aes_string(y="pos", label=show_samples_on_labels), vjust=0, angle=90, size=3, color="white", nudge_y=10, nudge_x=0.35)
    }
    p = p + guides(fill=guide_legend(ncol=legend_ncol))

    if(by_average == TRUE){
      prefix2 = paste0(prefix2, "_average")
    }

    curr_width = specs$width
    curr_height = specs$height

    if(!is.null(defined_width)){
      curr_width = defined_width
    }
    if(!is.null(defined_height)){
      curr_height = defined_height
    }

    if(isTRUE(png)){
      png( file=paste0(outdir1, "/taxonomy_", type, "_facets_", order_name, "", prefix, "_", prefix2, "_", tax_level, ".png"), height=curr_height, width=curr_width, units="in", res=300)
      print(p)
      dev.off()
    }
    if(isTRUE(pdf)){
      pdf( file=paste0(outdir1, "/taxonomy_", type, "__facets_", order_name, "", prefix, "_", prefix2, "_", tax_level, ".pdf"), height=curr_height, width=curr_width)
      print(p)
      dev.off()
    }
  }

  return(
    list(
      prefix=paste0(type, "", order_name, "", prefix, "_", prefix2, "_", tax_level),
      data=df2,
      ggplot_object=p
    )
  )
}

#' Using a taxonomic summary and mapping file along with what variable to compare, this function returns a data frame with evalues and BH for each comparison of each taxa.
#'
#' Make sure that data frame does not contain any illegal values.
#' Takes taxonomy tax table in input (i.e. text file in .tsv format).
#'
#' @param data_frame
#' @return  A dataframe object containing evalues and BH for each comparison of each taxa.
#' @examples
#' dataframe <- findDifferentialTaxaByAnova(...);
#' @export data_frame
findDifferentialTaxaByAnova <- function(mapping_file=NULL, mapping=NULL, taxonomy_file=NULL, variables=NULL, convert_to_relative_abundance=FALSE){

  if(!is.null(mapping) & !is.null(mapping_file)){
    stop("You can only specified a mapping or a mapping_file, but not both.")
  }

  # check if files exist
  if(!is.null(mapping_file)){
    if(!file.exists(mapping_file)){ stop(paste0("Mapping file does not exist: ", mapping_file)) }
  }
  if(!file.exists(taxonomy_file)){ stop(paste0("OTU table file does not exist: ", taxonomy_file)) }
  #outdir1 = paste0(outdir, "/taxonomy_", type)
  #dir.create(file.path(outdir1), showWarnings=TRUE, recursive=TRUE)

  # Load files
  if(is.null(mapping)){
    mapping = loadMappingFile(mapping_file)
    validateDataFrame(mapping)
  }else{
    if(nrow(mapping) == 0 | ncol(mapping) == 0){
      stop("mapping file appears empty...")
    }
  }

  anova_df = NULL
  anova_df2 = NULL

  tTaxonomy = data.frame(fread(taxonomy_file, header=TRUE, sep="\t"), check.names=FALSE)
  tTaxonomy = tTaxonomy[,c("Taxon", mapping$sample_id)]
  if(isTRUE(convert_to_relative_abundance)){
    tTaxonomyPerc = prop.table(data.matrix(tTaxonomy[,2:ncol(tTaxonomy)]), margin=2)*100
    tTaxonomyPerc = data.frame(Taxon=tTaxonomy$Taxon, tTaxonomyPerc, check.names=FALSE)
    tTaxonomy = tTaxonomyPerc
  }
  taxa = unique(tTaxonomy$Taxon)

  colname = names(variables)[1]
  variable_value_1 = variables[[colname]][1]
  variable_value_2 = variables[[colname]][2]

  tTaxonomy = melt(tTaxonomy)
  tTaxonomy = merge(tTaxonomy, mapping, by.x="variable", by.y="sample_id")
  print(head(tTaxonomy))

  for(j in 1:length(taxa)){
    print(paste0("   j: ", j))
    tTaxonomy2 = tTaxonomy[tTaxonomy$Taxon == taxa[j],]
    if(nrow(tTaxonomy2) == 0){next;}
    res = aov(as.formula(paste0("value ~ ", colname)), data=tTaxonomy2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]

    mean_1 = mean(tTaxonomy2[tTaxonomy2[[ colname ]] == variable_value_1,]$value)
    mean_2 = mean(tTaxonomy2[tTaxonomy2[[ colname ]] == variable_value_2,]$value)

    if(j == 1){
      anova_df = data.frame(Taxa=taxa[j], pvalue=pvalue, mean_1=mean_1, mean_2=mean_2)
    }else{
      anova_df = rbind(anova_df, data.frame(Taxa=taxa[j], pvalue=pvalue, mean_1=mean_1, mean_2=mean_2))
    }
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "BH")
  return(anova_df)
}


#' Using KO and COG aggregated summary matrix files along with what variable to compare, this function returns a data frame with
#' evalues and BH for each comparison of each KO.
#'
#' Make sure that data frame does not contain any illegal values.
#' Takes taxonomy tax table in input (i.e. text file in .tsv format).
#'
#' @param data_frame
#' @return  A dataframe object containing evalues and BH for each comparison of each taxa.
#' @examples
#' dataframe <- findDifferentialKOByAnova(...);
#' @export data_frame
findDifferentialKOByAnova <- function(mapping_file=NULL, mapping=NULL, COG_file=NULL, KO_file=NULL, variables=NULL){

  if(!is.null(mapping) & !is.null(mapping_file)){
    stop("You can only specified a mapping or a mapping_file, but not both.")
  }

  # check if files exist
  if(!is.null(mapping_file)){
    if(!file.exists(mapping_file)){ stop(paste0("Mapping file does not exist: ", mapping_file)) }
  }
  if(!file.exists(KO_file)){ stop(paste0("KO table file does not exist: ", KO_file)) }

  if(!is.null(COG_file)){
    if(!file.exists(COG_file)){ stop(paste0("COG file does not exist: ", COG_file)) }
  }else{
    stop(paste0("Please provide a COG file."))
  }

  # Load files
  if(is.null(mapping)){
    mapping = loadMappingFile(mapping_file)
    validateDataFrame(mapping)
  }else{
    if(nrow(mapping) == 0 | ncol(mapping) == 0){
      stop("mapping file appears empty...")
    }
  }

  tKO = data.frame(fread(KO_file, header=TRUE, sep="\t"), check.names=FALSE)
  tKO = tKO[tKO$KO != "NULL",]
  valid_colnames = colnames(tKO)[colnames(tKO) %in% mapping$sample_id]
  tKO = tKO[,c("KO", valid_colnames)]
  KOs = unique(tKO$KO)
  row.names(tKO) = tKO$KO; tKO$KO = NULL;

  # read COG file to normalize by recA (COG0468)
  tCOG = data.frame(fread(COG_file, header=TRUE, sep="\t"), check.names=FALSE)
  tCOG = tCOG[tCOG$C == "COG0468", ]
  row.names(tCOG) = tCOG$COG; tCOG$COG = NULL;
  tCOG = tCOG[,colnames(tKO)]
  if(!identical(colnames(tCOG), colnames(tKO))){stop("colnames(tCOG) not identical to colnames(tKO)...")}

  tKO = t(t(tKO) / as.numeric(tCOG[nrow(tCOG),]))
  tKO = data.frame(tKO)

  colname = names(variables)[1]
  variable_value_1 = variables[[colname]][1]
  variable_value_2 = variables[[colname]][2]

  tKO$KO = row.names(tKO);
  tKO = melt(tKO)
  tKO = merge(tKO, mapping, by.x="variable", by.y="sample_id")
  print(head(tKO))

  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    tKO2 = tKO[tKO$KO == KOs[j],]
    if(nrow(tKO2) == 0){next;}
    res = aov(as.formula(paste0("value ~ ", colname)), data=tKO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]

    mean_1 = mean(tKO2[tKO2[[ colname ]] == variable_value_1,]$value)
    mean_2 = mean(tKO2[tKO2[[ colname ]] == variable_value_2,]$value)

    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, mean_1=mean_1, mean_2=mean_2)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, mean_1=mean_1, mean_2=mean_2))
    }
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "BH")
  return(anova_df)
}



#' Returns a filtered feature table metting specified cutoffs.
#' Make sure that data frame does not contain any illegal values.
#' Takes taxonomy tax table in input (i.e. text file in .tsv format).
#'
#' @param data_frame
#' @return  A dataframe object containing evalues and BH for each comparison of each taxa.
#' @examples
#' dataframe <- filterFeatureTableNbyN(...);
#' @export data_frame
filterFeatureTableNbyN <- function(feature_table_file, X, Y, write_outfile=FALSE, outdir=outdir, prefix=NULL){

  #otu_table_file = "~/Projects/Genorem_amplicons/export/otu_tables/otu_table_filtered.tsv"
  data = data.frame(fread(otu_table_file, sep="\t", header=TRUE, skip="#FEATURE_ID"), check.names=FALSE)

  tax = data[,c(1,ncol(data))]

  # Keep OTUs : each OTU (row) must have at least 5 otus in at least 5 samples (5x5)

  row.names(data) = data[,1]
  data[,1] = NULL
  data$taxonomy = NULL


  data$filt<-apply(data, 1, function(x) sum(x>X)) #you can change the greater than to less than if you want to invert the count.
  df = data[data$filt>=Y,] #the 2 is made up by me for the case of wanting 2 or more columns that are .005 or greater.  Change the 2 for your needs
  data$filt = NULL #deleting dummy columns
  df$filt = NULL

  head(df)
  df[["#FEATURE_ID"]] = row.names(df)
  df2 = merge(df, tax, by.x="#FEATURE_ID", by.y="#FEATURE_ID")
  sample_names = colnames(df2)[2:(ncol(df2)-1)]
  df3 = df2[,c("#FEATURE_ID", sample_names, "taxonomy")]
  otus = df3$`#FEATURE_ID`

  if(write_outfile == TRUE){
    base = file_path_sans_ext(basename(otu_table_file))
    if(is.null(prefix)){
      prefix=""
    }else{
      prefix = paste0("_", prefix)
    }
    outfile = paste0(outdir, "/", base, "_", X, "by", Y, prefix, ".tsv")
    write.table(df3, outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  }

  return(otus)
}
