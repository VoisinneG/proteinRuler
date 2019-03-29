#' Compute protein abundance using the protein ruler methodology
#' @description Compute protein abundance using the protein ruler methodology
#' @param df A data.frame containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param col_protein_id Column with protein IDs. When several protein IDs are found, only the first one is conserved.
#' @param sep_id character string separating different protein IDs
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @param filtering Filter out contaminants and poorly scored proteins.
#' @param col_names Column with gene names.
#' @param col_score Column with protein identification score
#' @param Score_threshold Threshold on protein identification score
#' @param col_mass Column with protein mass (in kDa)
#' @param idx_histones Row indexes corresponding to histone proteins
#' @param show_progress Show progress bar when querrying annotations from UniProt
#' @param ... additionnal parameters passed to fucntion \code{compute_protein_number()}
#' @return a data.frame with protein abundances
#' @examples
#' data("proteinGroups_CD4_Tcells")
#' res <- proteinRuler(proteinGroups_CD4_Tcells, DNA_mass_per_cell = 5.5209e-12)
#' str(res)
#' cond <- "CopyNumber_WT_0"
#' hist(log10(res$copy_number[[cond]]), main = "",
#' xlab = paste(cond, "(log10)"), col = rgb(1,0,0,0.25))
#' @import queryup
#' @export
proteinRuler <- function(df,
                         col_protein_id = "Protein IDs",
                         sep_id = ";",
                         pattern_intensity = "^Intensity.",
                         col_intensity = NULL,
                         mass_per_cell_in_pg = NULL,
                         DNA_mass_per_cell = 5.5209e-12,
                         filtering = TRUE,
                         col_names = "Gene names",
                         col_score = "Score",
                         Score_threshold = 0,
                         idx_histones = NULL,
                         col_mass = NULL,
                         show_progress = TRUE,
                         ...){
  
  
  
  if(! col_protein_id %in% names(df)){
    stop(paste(col_protein_id, " is not a column name\n", sep = ""))
  }
  
  histone_ids <- df[[col_protein_id]][idx_histones]
  
  ######################### Filter proteins
  
  if(filtering){
    if(! col_names %in% names(df)){
      stop(paste(col_names, " is not a column name\n", sep = ""))
    }
    
    idx_filter <- c( grep("KRT",df[[col_names]]), 
                     grep("^CON",df[[col_protein_id]]), 
                     grep("^REV",df[[col_protein_id]]), 
                     grep("^Biognosys",df[[col_protein_id]]),
                     which(is.na(df[[col_protein_id]])))
    
    if(length(idx_filter)>0){
      df <- df[-idx_filter, ]
    }
    
    if(col_score %in% names(df)){
      df <- df[df[[col_score]]>Score_threshold, ]
    }else{
      warning(paste(col_score, " is not a column name\n", sep = ""))
    }
  }
  
  ######################### Import Annotations
  
  if(is.null(col_mass) | is.null(idx_histones)){
    
    id = sapply(df[[col_protein_id]], function(x){ strsplit(x, split = sep_id)[[1]][1]})
    df_annot <- queryup::get_annotations_uniprot(id = id, 
                                                 columns = c("genes", "families", "mass"), 
                                                 show_progress = show_progress)
    
    if(!is.null(df_annot)){
      #convert Mass from factors to numeric values
      df_annot$Mass <- as.numeric( gsub( "," , "" , as.character(df_annot$Mass) ) )
      df_annot$Mass <- df_annot$Mass/1e3 # Mass in kDa
      
      
      ########################### Identify Histone proteins
      
      u_prot_fams <- unique( strsplit(paste(as.character(df_annot$Protein.families), collapse=", "), split=", ")[[1]] ) ;
      Histone_fams <- u_prot_fams[grep("^Histone H", u_prot_fams)]
      idx_hist_all<-NULL
      for(i in 1:length(Histone_fams)){
        idx_hist <- grep(paste("(, |^)", Histone_fams[i], "($|, )",sep=""), 
                         as.character(df_annot$Protein.families), fixed=FALSE)
        idx_hist_all <- c(idx_hist_all, idx_hist);
      }
      df_annot$is_histone <- FALSE
      df_annot$is_histone[idx_hist_all] <- TRUE
      
    }else{
      stop(paste("Querying UniProt failed. Protein mass and/or histones identity could not be inputed.", "
                 Please try again later or use parameters 'idx_histones' and 'col_mass'instead.", sep = ""))
    }
    
  }
  
  if(is.null(col_mass)){
    col_mass <- "Mass"
    df$Mass <- df_annot$Mass
  }
  
  if(is.null(idx_histones)){
    idx_histones <- idx_hist_all
  }else{
    idx_histones <- match(histone_ids, df[[col_protein_id]])
  }
  
  
  res <- compute_protein_number(df,
                                col_mass = col_mass,
                                idx_histones = idx_histones,
                                pattern_intensity = pattern_intensity,
                                col_intensity = col_intensity,
                                col_ID = col_protein_id,
                                mass_per_cell_in_pg = mass_per_cell_in_pg,
                                DNA_mass_per_cell = DNA_mass_per_cell,
                                ...)
  
  if(!is.null(df_annot)){
    names(df_annot)[which(names(df_annot) == "Mass")] <- "Mass.kDa"
    res[["annotations"]] <- df_annot
  }
  
  return(res)
  
}

#' Compute protein abundance using the protein ruler methodology
#' @description Compute protein abundance using the protein ruler methodology
#' @param df A data.frame containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param idx_histones Row indexes corresponding to histone proteins
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param col_ID Column with IDs 
#' @param col_mass Column with protein mass (in kDa)
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @return a list containing the following elements :
#' @return \code{copy_number} : a data.frame with protein abundances
#' @return \code{summary} : a data.frame with summary variables 
#' @examples
#' data("proteinGroups_CD4_Tcells")
#' idx_h <- grep("^Histone H", proteinGroups_CD4_Tcells$`Protein names`)
#' col_mass <- "Mol. weight [kDa]"
#' res <- compute_protein_number(proteinGroups_CD4_Tcells, idx_histones = idx_h, col_mass = col_mass)
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
    col_int <- names(df)[grep(pattern_intensity,names(df))]
  }else{
    col_int <- col_intensity
  }
  
  copy_number <- vector("list", length = length(col_int))
  
  I_tot <- rep(0, length(col_int))
  I_hist_tot <- rep(0, length(col_int))
  median_I_tot <- rep(0, length(col_int))
  percentage_mass_hist <- rep(0, length(col_int))
  prot_mass_per_cell_pg <- rep(0, length(col_int))
  
  for( i in 1:length(col_int) ){
    
    median_I_tot[i] <- median(df_int[ , col_int[i]], na.rm=TRUE);
    I_hist_tot[i] <- sum(df_int[idx_histones, col_int[i]], na.rm=TRUE);
    I_tot[i] <- sum(df_int[, col_int[i]], na.rm=TRUE);
    percentage_mass_hist[i] <- I_hist_tot[i]/I_tot[i]*100;
    prot_mass_per_cell_pg[i] <- I_tot[i]/I_hist_tot[i]*DNA_mass_per_cell*1e12; #protein mass in pg
    
    
    if(is.null(mass_per_cell_in_pg)){
      copy_number[[col_int[i]]] = 6.022e23 * df_int[, col_int[i]] * DNA_mass_per_cell /
        (1e3*df[[col_mass]]*I_hist_tot[i]);
    }else{
      copy_number[[col_int[i]]] = 6.022e23 * df_int[, col_int[i]] * 40e-12 / 
        (1e3*df[[col_mass]]*I_tot[i])
    }
  }
  
  df_copy_number <- as.data.frame(do.call(cbind, copy_number))
  
  if(is.null(col_intensity)){
    names(df_copy_number) <- gsub(pattern_intensity, "CopyNumber_", col_int)
  }else{
    names(df_copy_number) <- paste("CopyNumber_", col_int, sep= "")
  }
  
  df_copy_number <- cbind( df_int[col_ID], df_copy_number)
  #names(df_copy_number) <- col_ID
  
  summary <- data.frame(column = col_int,
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
