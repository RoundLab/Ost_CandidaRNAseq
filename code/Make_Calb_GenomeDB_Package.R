library(c("biomaRt", "GenomeInfoDb", "GO.db"))
library(AnnotationForge)

makeOrgPackageFromNCBI(version = "0.1", author = "W. Zac Stephens <zac.stephens@path.utah.edu>", maintainer = "W. Zac Stephens <zac.stephens@path.utah.edu>", outputDir = "../data/", tax_id = "237561", genus = "Candida", species =  "albicans")
