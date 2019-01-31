library(readxl)

proteinGroups_CD4_Tcells <- read_excel("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_CD4_Kinetics/data/proteome_effector_CD4_T_cells.xlsx", sheet = 1) 

# save in ./R/sysdata.rda  ---------------------------------------------------------------------------------------

devtools::use_data(
  proteinGroups_CD4_Tcells,
  pkg=".",
  internal = FALSE,
  overwrite = TRUE)

