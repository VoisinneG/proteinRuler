% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteinRuler.R
\name{proteinRuler}
\alias{proteinRuler}
\title{Compute protein abundance using the protein ruler methodology}
\usage{
proteinRuler(
  df,
  col_protein_id = "Protein IDs",
  sep_id = ";",
  pattern_intensity = "^Intensity.",
  col_intensity = NULL,
  mass_per_cell_in_g = NULL,
  DNA_mass_per_cell = 5.5209e-12,
  filtering = TRUE,
  col_names = "Gene names",
  col_score = "Score",
  Score_threshold = 0,
  idx_histones = NULL,
  col_mass = NULL,
  show_progress = TRUE,
  replace_zero_by_na = TRUE,
  max_keys = 200,
  ...
)
}
\arguments{
\item{df}{A data.frame containing protein intensities. By default, 
protein intensity column names start by "Intensity." 
(use parameter \code{pattern_intensity} to change)}

\item{col_protein_id}{Column with protein IDs. When several protein IDs are found, 
only the first one is conserved.}

\item{sep_id}{character string separating different protein IDs}

\item{pattern_intensity}{Pattern (regular exrpression) used to identfy df's columns 
containing protein intensity values}

\item{col_intensity}{Names of intensity columns. Overrides \code{pattern_intensity}.}

\item{mass_per_cell_in_g}{Compute protein abundance using a constant mass per cell (in g).}

\item{DNA_mass_per_cell}{Mass of DNA per cell (in g). Set to 5.5pg 
(corresponding to murine cells) by default}

\item{filtering}{Filter out contaminants and poorly scored proteins.}

\item{col_names}{Column with gene names.}

\item{col_score}{Column with protein identification score}

\item{Score_threshold}{Threshold on protein identification score}

\item{idx_histones}{Row indexes corresponding to histone proteins. If NULL, 
histones are identified using UniProt annotations corresponding to the 
first ID of the protein group.}

\item{col_mass}{Column with protein mass (in kDa). If NULL, protein mass are 
retireved from UniProt using the first ID of the protein group.}

\item{show_progress}{Show progress bar when querrying annotations from UniProt}

\item{replace_zero_by_na}{Replace zero-valued intensities by NA.}

\item{max_keys}{Maximum query length used to retrieve UniProt data (passed to function
\code{pannot::get_annotations_uniprot()})}

\item{...}{additionnal parameters passed to fucntion \code{pannot::get_annotations_uniprot()}}
}
\value{
a data.frame with protein abundances
}
\description{
Compute protein abundance using the protein ruler methodology
}
\examples{
data("proteinGroups_CD4_Tcells")
res <- proteinRuler(proteinGroups_CD4_Tcells, DNA_mass_per_cell = 5.5209e-12)
str(res)
cond <- "CopyNumber_WT_0"
hist(log10(res$copy_number[[cond]]), main = "",
xlab = paste(cond, "(log10)"), col = rgb(1,0,0,0.25))
}
