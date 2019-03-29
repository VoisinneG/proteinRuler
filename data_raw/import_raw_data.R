library(readxl)

df <- read_excel("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_CD4_Kinetics/data/proteome_effector_CD4_T_cells.xlsx", 
                 sheet = 1,
                ) 

#df <- df[sample(1:dim(df)[1], size = 1000), 1:64]
df <- df[, 1:64]
proteinGroups_CD4_Tcells <- df

# save in ./R/sysdata.rda  ---------------------------------------------------------------------------------------

usethis::use_data(
  proteinGroups_CD4_Tcells,
  internal = FALSE,
  overwrite = TRUE)

