### Data folder listing:
1. **GeneID2Alias_reform.txt** : Mapping of candidagenome.org Assembly22 format IDs to NCBI Entrez Gene IDs. Result of ../code/MapEntrezID_to_Assembly22.sh
2. **gene_result.txt** : Input to ../code/MapEntrezID_to_Assembly22.sh. NCBI gene metadata table for C. albicans SC5314.
3. **ORF19_Assembly22_mapping.tab** : Orf19 to Assembly22 mapping file from candidagenome.org
4. **sleuth_table_wt.txt** : Differential expression results table from Sleuth. Positive fold-changes are up in Rag-/- animals relative to wild-type.
5. **Witchley_Up-Hyphae.txt** : List of genes upregulated in hyphae from Witchley, et al. 
6. **org.Calbicans.eg.db/**: C. albicans genome database created with AnnotationForge package. Result of ../code/Make_Calb_GenomeDB_Package.R
