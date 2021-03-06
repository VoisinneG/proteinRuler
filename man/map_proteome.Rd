% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteinRuler.R
\name{map_proteome}
\alias{map_proteome}
\title{Map protein abundance using protein IDs or gene names}
\usage{
map_proteome(
  df,
  col_ID = "Protein.IDs",
  col_names = "names",
  sep_primary = NULL,
  sep_secondary = NULL,
  proteome_dataset,
  pdata_sep = ";",
  pdata_col_ID = "Protein.IDs",
  pdata_col_gene_name = "Gene.names",
  pdata_col_copy_number = "Copy.Number",
  map_gene_name = FALSE,
  return_indexes = FALSE,
  updateProgress = NULL
)
}
\arguments{
\item{df}{a data.frame or a list}

\item{col_ID}{column of \code{df} containing protein IDs}

\item{col_names}{column of \code{df} containing gene names. 
Only used if \code{map_gene_name = TRUE}.}

\item{sep_primary}{Separator between different proteins in \code{df}}

\item{sep_secondary}{Set of separators used sequentially (from right to left) to 
identify protein IDs for each protein in \code{df}}

\item{proteome_dataset}{Dataset containing protein abundances.}

\item{pdata_sep}{Separator between different proteins in \code{proteome_dataset}}

\item{pdata_col_ID}{column of \code{proteome_dataset} containing protein IDs}

\item{pdata_col_gene_name}{column of \code{proteome_dataset} containing gene names}

\item{pdata_col_copy_number}{column of \code{proteome_dataset} containing 
protein abundances}

\item{map_gene_name}{logical, map protein using gene names rather than protein IDs}

\item{return_indexes}{logical, return vector of matching indexes}

\item{updateProgress}{used to display progress in shiny apps
protein abundances (in log10)}
}
\description{
Add protein abundance to an \code{InteRactome}. For multiple identifiers, 
the abundance of the first match in the proteome dataset is returned.
}
\examples{
names <- LETTERS
names[3] <- "Z;A"
copies <- runif(length(LETTERS), 0, 100)
proteome <- data.frame(names, copies, stringsAsFactors = FALSE)
names(proteome) <- c("Gene.names", "Copy.Number")
df <- data.frame("names" = c("C", "M", "Z", "A", "KL Z;A"))
res <- map_proteome(df, col_names = "names", proteome_dataset = proteome, map_gene_name = TRUE)
}
