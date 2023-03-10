% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Microbiomics.R
\name{loadMappingFile}
\alias{loadMappingFile}
\title{Load mapping file}
\usage{
loadMappingFile(mapping_file)
}
\arguments{
\item{file_path}{to mapping file in .tsv format.}
}
\value{
A data frame for which the row.names is the #SampleID column.
}
\description{
Takes a text .tsv format mapping file and loads in into a data.frame.
}
\examples{
mapping <- loadMappingFile(mapping_file);
}
