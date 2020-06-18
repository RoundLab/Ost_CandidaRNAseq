library(data.table)
library(tidyverse)
library(fgsea)

# Built a genome database for C. alb with annotation forge. Now, install from source.
install.packages("../data/org.Calbicans.eg.db/", repos = NULL, type = "source")
library(org.Calbicans.eg.db)

# Download mappings of ORF19 IDs to Assembly 22 IDs and common gene names
download.file(url = "http://candidagenome.org/download/chromosomal_feature_files/C_albicans_SC5314/ORF19_Assembly22_mapping.tab", destfile = "../data/ORF19_Assembly22_mapping.tab")
map.Calb.GeneNames <- read_table2(file = "../data/ORF19_Assembly22_mapping.tab")


# Read in table of sleuth results of all genes
sleuth_table_wt <- read_delim("../data/sleuth_table_wt.txt", 
	"\t", escape_double = FALSE, col_types = cols(GeneID = col_character()), 
	trim_ws = TRUE)

# Merge tables to add the ORF19 mappings to assembly 22 IDs
sleuth_table_wt <- merge.data.table(as.data.table(sleuth_table_wt), map.Calb.GeneNames, by.x = "target_id", by.y = "ASSEMBLY22_ID", all.x = TRUE)

# Read in the Witchley, et al dataset of genes upregulated in hyphal growth
in.Witchley.UpHyph <- read.delim(file="../data/Witchley_Up-Hyphae.txt", header = FALSE)
# Add the Assembly 22 IDs to Witchley dataset
tmp_Witchley <- merge.data.table(in.Witchley.UpHyph, sleuth_table_wt, by.x = "V1", by.y = "ORF19_ID")
# Creat a gene set by Gene ID for the Witchley dataset. By making list we can add multiple gene sets, though only one is added here.
PubGeneSets <- list(tmp_Witchley$V1)
names(PubGeneSets) <- c("Witchley_Up_Hyphal")
PubGeneSets.byGeneID <- list(tmp_Witchley$Gene)
names(PubGeneSets.byGeneID) <- c("Witchley_Up_Hyphal")
rm(tmp_Witchley)
rm(in.Witchley.UpHyph)

# Subset results to significantly differentially expressed genes (q < 0.05). Create named, sorted vector for GSEA.
gL.q05.all <- filter(sleuth_table_wt, qval <= 0.05) %>% dplyr::select(b, Gene)
gL.q05.all <- setNames(gL.q05.all$b, gL.q05.all$Gene)
gL.q05.all <- sort(gL.q05.all, decreasing = TRUE)

# Run gsea with Witchley dataset and significantly different genes in Rag1. 
gseaRes.q05.all <- fgsea(PubGeneSets.byGeneID, gL.q05.all, nperm = 9999)
# Plot. Inverting log2FC values to reflect positives up in WT (aesthetic adjustment)
plotEnrichment(PubGeneSets.byGeneID[["Witchley_Up_Hyphal"]], (gL.q05.all * -1)) + labs(title = "Witchley: Hyphal Growth")


