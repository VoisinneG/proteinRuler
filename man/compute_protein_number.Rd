% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteinRuler.R
\name{compute_protein_number}
\alias{compute_protein_number}
\title{Compute protein abundance using the protein ruler methodology}
\usage{
compute_protein_number(
  df,
  idx_histones,
  pattern_intensity = "^Intensity.",
  col_intensity = NULL,
  col_mass = "Mass",
  col_ID = names(df)[1],
  mass_per_cell_in_g = NULL,
  DNA_mass_per_cell = 5.5209e-12,
  replace_zero_by_na = TRUE
)
}
\arguments{
\item{df}{A data.frame containing protein intensities. By default, 
protein intensity column names start by "Intensity." 
(use parameter \code{pattern_intensity} to change)}

\item{idx_histones}{Row indexes corresponding to histone proteins}

\item{pattern_intensity}{Pattern (regular exrpression) used to identfy df's columns 
containing protein intensity values}

\item{col_intensity}{Names of intensity columns. Overrides \code{pattern_intensity}.}

\item{col_mass}{Column with protein mass (in kDa)}

\item{col_ID}{Column with IDs}

\item{mass_per_cell_in_g}{Compute protein abundance using a constant mass per cell (in g).}

\item{DNA_mass_per_cell}{Mass of DNA per cell (in g)}

\item{replace_zero_by_na}{Replace zero-valued intensities by NA.}
}
\value{
a list containing the following elements :

\code{copy_number} : a data.frame with protein abundances

\code{summary} : a data.frame with summary variables
}
\description{
Compute protein abundance using the protein ruler methodology
}
\examples{
\dontrun{
data("proteinGroups_CD4_Tcells")
idx_h <- grep("^Histone H", proteinGroups_CD4_Tcells$`Protein names`)
col_mass <- "Mol. weight [kDa]"
res <- compute_protein_number(proteinGroups_CD4_Tcells, 
idx_histones = idx_h, 
col_mass = col_mass)
}
}
