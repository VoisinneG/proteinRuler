[![Travis-CI Build Status](https://travis-ci.org/VoisinneG/InteRact.svg?branch=master)](https://travis-ci.org/VoisinneG/proteinRuler) 

# R Package : proteinRuler

The proteinRuler package allows to compute protein abundance from protein intensities using the protein ruler methodology

Installation
---
Install the package from the github repository using:
```
devtools::install_github("VoisinneG/proteinRuler")
```

Usage
---
Import a dataset containing protein intensities and protein IDs.

```
library(proteinRuler)
data("proteinGroups_CD4_Tcells")
```

Compute protein abundances.

```
res <- proteinRuler(df = proteinGroups_CD4_Tcells,
                    organism = "mouse",
                    DNA_mass_per_cell = 5.5209e-12)
                    
print(res$summary)
print(res$copy_number)
```


