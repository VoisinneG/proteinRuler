#' Compute protein abundance using the protein ruler methodology
#' @description Compute protein abundance using the protein ruler methodology
#' @param df A data.frame containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param col_names Column with gene names.
#' @param col_protein_id Column with protein IDs
#' @param col_score Column with protein identification score
#' @param Score_threshold Threshold on protein identification score
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param organism Either "mouse" or "human"
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @param ... additionnal parameters passed to fucntion \code{compute_protein_number()}
#' @return a data.frame with protein abundances
#' @examples
#' data("proteinGroups_CD4_Tcells")
#' res <- proteinRuler(proteinGroups_CD4_Tcells, organism = "mouse", DNA_mass_per_cell = 5.5209e-12)
#' @import pannot
#' @export
proteinRuler <- function(df,
                         col_names = "Gene names",
                         col_protein_id = "Protein IDs",
                         col_score = "Score",
                         Score_threshold = 0,
                         pattern_intensity = "^Intensity.",
                         col_intensity = NULL,
                         organism = "mouse",
                         mass_per_cell_in_pg = NULL,
                         DNA_mass_per_cell = 5.5209e-12,
                         ...){
  
  
  if(! col_names %in% names(df)){
    stop(paste(col_names, "is not a column name\n", sep = ""))
  }
  if(! col_protein_id %in% names(df)){
    stop(paste(col_protein_id, "is not a column name\n", sep = ""))
  }
  ######################### Filter proteins
  
  idx_filter <- c( grep("KRT",df[[col_names]]), 
                   grep("^CON",df[[col_protein_id]]), 
                   grep("^REV",df[[col_protein_id]]), 
                   grep("^Biognosys",df[[col_protein_id]]) )
  
  if(length(idx_filter)>0){
    df <- df[-idx_filter, ]
  }
  
  if(col_score %in% names(df)){
    df <- df[df[[col_score]]>Score_threshold, ]
  }else{
    warning(paste(col_score, "is not a column name\n", sep = ""))
  }
  
  
  ######################### Import Annotations
  
  df_annot <- pannot::get_annotations(df, name_id = col_protein_id, organism = organism, split_param = ";")
  #convert Mass (in kDa) from factors to numeric values
  df_annot$Mass <- as.numeric( sub( "," , "." , as.character(df_annot$Mass) ) );
  df$Mass <- df_annot$Mass
  
  ########################### Identify Histone proteins
  
  u_prot_fams <- unique( strsplit(paste(as.character(df_annot$Protein.families), collapse=", "), split=", ")[[1]] ) ;
  Histone_fams <- u_prot_fams[grep("^Histone H", u_prot_fams)]
  idx_hist_all<-NULL
  for(i in 1:length(Histone_fams)){
    idx_hist <- grep(paste("(, |^)", Histone_fams[i], "($|, )",sep=""), 
                     as.character(df_annot$Protein.families), fixed=FALSE)
    idx_hist_all <- c(idx_hist_all, idx_hist );
  }
  
  
  res <- compute_protein_number(df,
                                idx_histones = idx_hist_all,
                                pattern_intensity = pattern_intensity,
                                col_intensity = col_intensity,
                                col_ID = col_protein_id,
                                mass_per_cell_in_pg = mass_per_cell_in_pg,
                                DNA_mass_per_cell = DNA_mass_per_cell
  )
  
  return(res)
  
}

#' Compute protein abundance using the protein ruler methodology
#' @description Compute protein abundance using the protein ruler methodology
#' @param df A data.frame containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param idx_histones Row indexes corresponding to histone proteins
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param col_ID Column with protein mass (in kDa)
#' @param col_mass Column with IDs 
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @return a list containing the following elements :
#' @return \code{copy_number} : a data.frame with protein abundances
#' @return \code{summary} : a data.frame with summary variables 
#' @examples
#' data("proteinGroups_CD4_Tcells")
#' idx_histones <- grep("^Histone H", proteinGroups_CD4_Tcells$`Protein names`)
#' col_mass <- "Mol. weight [kDa]"
#' res <- proteinRuler(proteinGroups_CD4_Tcells, idx_histones = idx_histones, col_mass = col_mass)
#' @importFrom stats median
#' @export
compute_protein_number <- function(df,
                                   idx_histones,
                                   pattern_intensity = "^Intensity.",
                                   col_intensity = NULL, 
                                   col_mass = "Mass",
                                   col_ID = names(df)[1],
                                   mass_per_cell_in_pg = NULL,
                                   DNA_mass_per_cell = 5.5209e-12 ){
  
  df_int <- as.data.frame(df)
  if(is.null(col_intensity)){
    col_intensity <- names(df)[grep(pattern_intensity,names(df))]
  }
  
  copy_number <- vector("list", length = length(col_intensity))
  
  I_tot <- rep(0, length(col_intensity))
  I_hist_tot <- rep(0, length(col_intensity))
  median_I_tot <- rep(0, length(col_intensity))
  percentage_mass_hist <- rep(0, length(col_intensity))
  prot_mass_per_cell_pg <- rep(0, length(col_intensity))
  
  for( i in 1:length(col_intensity) ){
    
    median_I_tot[i] <- median(df_int[ , col_intensity[i]], na.rm=TRUE);
    I_hist_tot[i] <- sum(df_int[idx_histones, col_intensity[i]], na.rm=TRUE);
    I_tot[i] <- sum(df_int[, col_intensity[i]], na.rm=TRUE);
    percentage_mass_hist[i] <- I_hist_tot[i]/I_tot[i]*100;
    prot_mass_per_cell_pg[i] <- I_tot[i]/I_hist_tot[i]*DNA_mass_per_cell*1e12; #protein mass in pg
    
    
    if(is.null(mass_per_cell_in_pg)){
      copy_number[[col_intensity[i]]] = 6.022e23 * df_int[, col_intensity[i]] * DNA_mass_per_cell /
        (1e3*df[[col_mass]]*I_hist_tot[i]);
    }else{
      copy_number[[col_intensity[i]]] = 6.022e23 * df_int[, col_intensity[i]] * 40e-12 / 
        (1e3*df[[col_mass]]*I_tot[i])
    }
  }
  
  df_copy_number <- as.data.frame(do.call(cbind, copy_number))
  
  if(is.null(col_intensity)){
    names(df_copy_number) <- gsub(pattern_intensity, "CopyNumber", col_intensity)
  }else{
    names(df_copy_number) <- paste("CopyNumber_", col_intensity, sep= "")
  }
  
  df_copy_number <- cbind( df_int[col_ID], df_copy_number)
  #names(df_copy_number) <- col_ID
  
  summary <- data.frame(column = col_intensity,
                        median_Intensity = median_I_tot,
                        sum_Intensity = I_tot,
                        sum_Intensity_histones = I_hist_tot,
                        percentage_mass_histones = percentage_mass_hist,
                        protein_mass_per_cell_pg = prot_mass_per_cell_pg
  )
  
  return(list("copy_number" = df_copy_number, 
              "summary" = summary
  )
  )
  
}
