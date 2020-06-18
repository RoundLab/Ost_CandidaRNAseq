#!/bin/bash

# Purpose: Map Entrez Gene IDs to C. albicans Assembly 22 style names. 
#   NCBI symbols have the same information embedded from assembly22, but in a NCBI-specific format (eg. with "CAALFM" prefix and underscore only used to separate prefix).  
# Download NCBI gene table for "Candida albicans SC5314" -> gene_result.txt
# Note that NCBI gene table should only contain entries for the SC5314 strain.

cut -f 3,6 ../data/gene_result.txt > ../data/GeneID2Alias.txt

cat ../data/GeneID2Alias.txt | while read GeneID Symbol
do if [[ "${Symbol}" == CAALFM_*  ]]
then 
	Allele=$(echo ${Symbol} | sed 's/CAALFM_//g' | tail -c 2)
	ChromLoc=$(echo ${Symbol} | sed 's/CAALFM_//g' | cut -c 1,2)
	FeatureID=$(echo ${Symbol} | sed 's/CAALFM_//g' | sed 's/^..//g' | sed 's/.$//g')
	NewSymbol=${ChromLoc}_${FeatureID}_${Allele}
else
	NewSymbol=${Symbol}
fi
echo -e "${NewSymbol}\t${GeneID}"
done > ../data/GeneID2Alias_reform.txt

rm ../data/GeneID2Alias.txt

# R code to add EntrezIDs with results Sleuth results table that has Assembly22 IDs.
# sleuth_table_wt <- data.table::merge.data.table(x = sleuth_table_wt, y = GENEID2candID, by.x = "Gene", by.y = "Symbol", all.x = TRUE)
