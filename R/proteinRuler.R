#' Compute protein abundance using the protein ruler methodology
#' @description Compute protein abundance using the protein ruler methodology
#' @param df A data.frame containing protein intensities. By default, 
#' protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param col_protein_id Column with protein IDs. When several protein IDs are found, 
#' only the first one is conserved.
#' @param sep_id character string separating different protein IDs
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns 
#' containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @param filtering Filter out contaminants and poorly scored proteins.
#' @param col_names Column with gene names.
#' @param col_score Column with protein identification score
#' @param Score_threshold Threshold on protein identification score
#' @param col_mass Column with protein mass (in kDa). If NULL, protein mass are 
#' retireved from UniProt using the first ID of the protein group.
#' @param idx_histones Row indexes corresponding to histone proteins. If NULL, 
#' histones are identified using UniProt annotations corresponding to the 
#' first ID of the protein group.
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
#' @import pannot
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
  
  ### Sanity checks #######################################################################
  
  if(!is.data.frame(df)){
    stop("df must be a data.frame")
  }
  
  if(is.null(col_intensity)){
    idx_cols <- grep(pattern_intensity, names(df))
    if(length(idx_cols) == 0){
      stop("Could not find pattern in column names")
    }
  }else{
    if(!all(col_intensity %in% names(df))){
      stop(paste("Some intensity columns could no be found"))
    }
    idx_cols <- match(col_intensity, names(df))
  }
  
  if(! all(sapply(idx_cols, function(i){class(df[[i]])}) %in% "numeric")){
    stop("Some intensity columns are not numeric")
  }
  
  if(! col_protein_id %in% names(df)){
    stop(paste(col_protein_id, " is not a column name\n", sep = ""))
  }
  
  if(!is.null(col_mass)){
    if(! col_mass %in% names(df)){
      stop(paste(col_mass, " is not a column name\n", sep = ""))
    }
    if(!is.numeric(df[[col_mass]])){
      stop(paste(col_mass, " column is not of class numeric\n", sep = ""))
    }
  }
  
  ### Filter proteins #######################################################################
  
  if(filtering){
    if(! col_names %in% names(df)){
      stop(paste(col_names, " is not a column name\n", sep = ""))
    }
    
    idx_filter <- c( grep("^CON", as.character(df[[col_protein_id]])), 
                     grep("^REV", as.character(df[[col_protein_id]])), 
                     grep("^Biognosys", as.character(df[[col_protein_id]])),
                     which(is.na(df[[col_protein_id]])))
    
    if(length(idx_filter)>0){
      df <- df[-idx_filter, ]
    }
    
    if(!is.null(col_score)){
      if(! col_score %in% names(df)){
        stop(paste(col_score, " is not a column name\n", sep = ""))
      }
      if(!is.numeric(df[[col_score]])){
        stop(paste(col_score, " column is not of class numeric\n", sep = ""))
      }
    }
    
    df <- df[df[[col_score]] > Score_threshold, ]
    
  }
  
  ### Retrieve protein annotations from UniProt #############################################
  
  df_annot <- NULL
  
  if(is.null(col_mass) | is.null(idx_histones)){
    
    id = sapply(as.character(df[[col_protein_id]]), 
                function(x){ strsplit(x, split = sep_id)[[1]][1]})
    
    df_annot <- pannot::get_annotations_uniprot(id = id, 
                                                 columns = c("genes", "families", "mass"), 
                                                 show_progress = show_progress)
    
    if(!is.null(df_annot)){
      #convert Mass from factors to numeric values
      df_annot$Mass <- as.numeric( gsub( "," , "" , as.character(df_annot$Mass) ) )
      df_annot$Mass <- df_annot$Mass/1e3 # Mass in kDa
      
      
      #### Identify Histone proteins ########################################################
      
      u_prot_fams <- unique( strsplit(paste(as.character(df_annot$Protein.families), 
                                            collapse=", "), 
                                      split=", ")[[1]] ) ;
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
      stop("Querying UniProt failed. Protein mass and/or histones identity 
      could not be inputed..Please try again later or use parameters 
      'idx_histones' and 'col_mass'instead.")
    }
    
  }
  
  if(is.null(col_mass)){
    col_mass <- "Mass"
    df$Mass <- df_annot$Mass
  }
  
  if(is.null(idx_histones)){
    idx_histones <- idx_hist_all
  }else{
    histone_ids <- df[[col_protein_id]][idx_histones]
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
#' @param df A data.frame containing protein intensities. By default, 
#' protein intensity column names start by "Intensity." 
#' (use parameter \code{pattern_intensity} to change)
#' @param idx_histones Row indexes corresponding to histone proteins
#' @param pattern_intensity Pattern (regular exrpression) used to identfy df's columns 
#' containing protein intensity values
#' @param col_intensity Names of intensity columns. Overrides \code{pattern_intensity}.
#' @param col_ID Column with IDs 
#' @param col_mass Column with protein mass (in kDa)
#' @param mass_per_cell_in_pg Compute protein abundance using a constant mass per cell.
#' @param DNA_mass_per_cell Mass of DNA per cell (in g)
#' @param replace_zero_by_na Replace zero-valued intensities by NA.
#' @return a list containing the following elements :
#' @return \code{copy_number} : a data.frame with protein abundances
#' @return \code{summary} : a data.frame with summary variables 
#' @examples
#' \dontrun{
#' data("proteinGroups_CD4_Tcells")
#' idx_h <- grep("^Histone H", proteinGroups_CD4_Tcells$`Protein names`)
#' col_mass <- "Mol. weight [kDa]"
#' res <- compute_protein_number(proteinGroups_CD4_Tcells, 
#' idx_histones = idx_h, 
#' col_mass = col_mass)
#' }
#' @importFrom stats median
#' @export
compute_protein_number <- function(df,
                                   idx_histones,
                                   pattern_intensity = "^Intensity.",
                                   col_intensity = NULL, 
                                   col_mass = "Mass",
                                   col_ID = names(df)[1],
                                   mass_per_cell_in_pg = NULL,
                                   DNA_mass_per_cell = 5.5209e-12,
                                   replace_zero_by_na = TRUE){
  
  ### Sanity checks #######################################################################
  
  if(!is.data.frame(df)){
    stop("df must be a data.frame")
  }
  
  if(is.null(col_intensity)){
    idx_cols <- grep(pattern_intensity, names(df))
    if(length(idx_cols) == 0){
      stop("Could not find pattern in column names")
    }
  }else{
    if(!all(col_intensity %in% names(df))){
      stop(paste("Some intensity columns could no be found"))
    }
    idx_cols <- match(col_intensity, names(df))
  }
  
  if(! all(sapply(idx_cols, function(i){class(df[[i]])}) %in% "numeric")){
    stop("Some intensity columns are not numeric")
  }
  
  if(! col_ID %in% names(df)){
    stop(paste(col_ID, " is not a column name\n", sep = ""))
  }
  
  if(! col_mass %in% names(df)){
    stop(paste(col_mass, " is not a column name\n", sep = ""))
  }
  
  if(!is.numeric(df[[col_mass]])){
    stop(paste(col_mass, " column is not of class numeric\n", sep = ""))
  }
    

  
  ### Compute protein copy number ###########################################################
  
  df_int <- as.data.frame(df)
  if(is.null(col_intensity)){
    col_int <- names(df)[grep(pattern_intensity, names(df))]
  }else{
    col_int <- col_intensity
  }
  
  if(replace_zero_by_na){
    for( i in 1:length(col_int) ){
      idx <- which( df_int[ , col_int[i] ] == 0 )
      df_int[idx , col_int[i]] <- NA
    }
  }
  
  copy_number <- vector("list", length = length(col_int))
  names(copy_number) <- col_int
    
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
    #protein mass in pg
    prot_mass_per_cell_pg[i] <- I_tot[i]/I_hist_tot[i]*DNA_mass_per_cell*1e12; 
    
    if(is.null(mass_per_cell_in_pg)){
      copy_number[[col_int[i]]] = 6.022e23 * df_int[, col_int[i]] * DNA_mass_per_cell / 
        (1e3*df[[col_mass]]*I_hist_tot[i])
    }else{
      copy_number[[col_int[i]]] = 6.022e23 * df_int[, col_int[i]] * mass_per_cell_in_pg / 
        (1e3*df[[col_mass]]*I_tot[i])
    }
  }
  
  df_copy_number <- as.data.frame(do.call(cbind, copy_number))
  
  if(is.null(col_intensity)){
    names(df_copy_number) <- gsub(pattern_intensity, "CopyNumber_", col_int)
  }else{
    names(df_copy_number) <- paste("CopyNumber_", col_int, sep= "")
  }
  
  df_copy_number <- data.frame(df[[col_ID]], df_copy_number)
  names(df_copy_number)[1] <- col_ID
  
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


#' Map protein abundance using protein IDs or gene names
#' @description Add protein abundance to an \code{InteRactome}. For multiple identifiers, 
#' the abundance of the first match in the proteome dataset is returned.
#' @param df a data.frame or a list
#' @param col_ID column of \code{df} containing protein IDs
#' @param col_names column of \code{df} containing gene names. 
#' Only used if \code{map_gene_name = TRUE}.
#' @param sep_primary Separator between different proteins in \code{df}
#' @param sep_secondary Set of separators used sequentially (from right to left) to 
#' identify protein IDs for each protein in \code{df}
#' @param pdata_sep Separator between different proteins in \code{proteome_dataset}
#' @param proteome_dataset Dataset containing protein abundances.
#' @param pdata_col_ID column of \code{proteome_dataset} containing protein IDs
#' @param pdata_col_gene_name column of \code{proteome_dataset} containing gene names
#' @param pdata_col_copy_number column of \code{proteome_dataset} containing 
#' protein abundances
#' @param map_gene_name logical, map protein using gene names rather than protein IDs
#' @param return_indexes logical, return vector of matching indexes
#' @param updateProgress used to display progress in shiny apps
#' protein abundances (in log10)
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples 
#' names <- LETTERS
#' names[3] <- "Z;A"
#' copies <- runif(length(LETTERS), 0, 100)
#' proteome <- data.frame(names, copies, stringsAsFactors = FALSE)
#' names(proteome) <- c("Gene.names", "Copy.Number")
#' df <- data.frame("names" = c("C", "M", "Z", "A", "KL Z;A"))
#' res <- map_proteome(df, col_names = "names", proteome_dataset = proteome, map_gene_name = TRUE)$df
#' @export
map_proteome <- function( df, 
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
                            updateProgress = NULL){
  
  ### Sanity checks #######################################################################
  
  df_int <- df
  
  
  if(! map_gene_name & ! col_ID %in% names(df)){
    stop(paste(col_ID, "is not defined within the interactome"))
  }
  
  if(map_gene_name & ! col_names %in% names(df)){
    stop(paste(col_names, "is not defined within the interactome"))
  }
  if(! map_gene_name & ! pdata_col_ID %in% names(proteome_dataset)){
    stop(paste(pdata_col_ID, "is not defined within the proteome"))
  }
  if(! pdata_col_copy_number %in% names(proteome_dataset)){
    stop(paste(pdata_col_copy_number, "is not defined within the proteome"))
  }
  if(map_gene_name & ! pdata_col_gene_name %in% names(proteome_dataset)){
    stop(paste(pdata_col_gene_name, "is not defined within the proteome"))
  }
  
  
  pdata <- proteome_dataset
  
  col_map <- col_ID
  pdata_col_map <- pdata_col_ID
  
  if(map_gene_name){
    col_map <- col_names
    pdata_col_map <- pdata_col_gene_name
    
  }
  
  if(is.null(sep_primary)){
    sep_primary_int <- sep_primary
  }
  if(is.null(sep_secondary)){
    sep_secondary_int <- sep_secondary
  }
  
  if(map_gene_name){
    sep_secondary_int <- " "
    sep_primary_int <- ";"
  }else{
    sep_secondary_int <- c("|", "-")
    sep_primary_int <- ";"
  }
  
  names <- df_int[[col_map]]
  
  ### Map proteins #########################################################################
  
  idx_match_all <- rep(NA, length(names));
  Copy_Number <- rep(NA, length(names));
  
  cat("Get protein abundances...\n")
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  
  for( i in 1:length(names) ){
    
    if (is.function(updateProgress)) {
      text <- paste0( i/length(names)*100)
      updateProgress(value = floor(i/length(names)*100), 
                     detail = paste(format(i/length(names)*100, digits = 0), "%", sep = ""))
    }
    setTxtProgressBar(pb, i/length(names)*100)
    
    prot_ids <- strsplit(as.character(names[i]), split = sep_primary_int, fixed = TRUE)[[1]]
    
    for(j in 1:length(prot_ids)){
      
      prot_id_int <- prot_ids[j]
      if(length(sep_secondary_int)>0){
        for(k in 1:length(sep_secondary_int)){
          prot_id_int <- strsplit(prot_id_int, 
                                  split = sep_secondary_int[k], 
                                  fixed = TRUE)[[1]][1]
        }
      }
      
      idx_match <-match_multi(x = toupper(prot_id_int), 
                              y = toupper(as.character(pdata[[pdata_col_map]])),
                              sep = pdata_sep)
      if(!is.na(idx_match)){
        idx_match_all[i] <- idx_match
        break
      }
      
      # idx_match <-grep( paste("(^|", sep_primary_int, ")", 
      #                         toupper(prot_id_int), "($|", sep_primary_int, ")", sep =""),
      #                   toupper(as.character(pdata[[pdata_col_map]])),
      #                   fixed = FALSE)
      # 
      # 
      # if(length(idx_match)>0){
      #   idx_match_all[i] <- idx_match[1]
      #   break
      # }
      
    }
    
  }
  close(pb)
  
  Copy_Number <- pdata[[pdata_col_copy_number]][idx_match_all]
  df_int$Copy_Number <- Copy_Number
  
  if(return_indexes){
    return(idx_match_all)
  }else{
    return(df_int)
  }
  
}

#' Extension of the base 'match' function
#' @description Extension of the base 'match' function to the case where 
#' elements of y can contain multiple items to be matched
#' @param x a character vector
#' @param y a character vector
#' @param sep character separating different items to be matched in a given y element
#' @return a numeric vector with matching indices
#' @examples
#' x <- c("a", "b")
#' y <- c( "bb", "c;b", "aa", "c;a;b")
#' match_multi(x=x, y=y, sep = ";")
#' @export
match_multi <- function(x, y, sep=";"){
  if(length(x)>1){
    return(sapply(x, function(x0){match_multi(x0, y, sep = sep)} ))
  }else{
    idx <- grep( paste("(^|", sep,")", 
                       x, 
                       "($|",sep,")", 
                       sep = ""), 
                 y)
    if(length(idx)>0){
      return(idx[1])
    }else{
      return(NA)
    }
  }
}
