library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

df.rf <- read.table(file = "data/data_original_genenames.tsv", sep ="\t",
                    header =T)

all_proteins <- colnames(data)
# remove HC
data_subset <- df.rf[!grepl("HC", rownames(df.rf)), ]
proteins_subset <- colnames(data_subset)

# Perform enrichment analysis
perform_enrichment <- function(proteins) {
  # Convert protein names to ENTREZ IDs
  entrez_ids <- bitr(proteins_subset, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Perform GO enrichment analysis
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'ENTREZID',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  # Simplify results
  ego_simplified <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  return(as.data.frame(ego_simplified))
}

# Run the enrichment analysis
enrichment_results <- perform_enrichment(all_proteins)

# Filter the data to include only entries with ontology "BP"
data_bp605 <- enrichment_results %>% filter(ONTOLOGY == "BP")

# Further filter descriptions related to immune and inflammatory responses
immune_inflammatory_terms <- c("immune", "inflamm", "cytokine", "macrophage", "lymphocyte", "wound healing")
data_bp605_filtered <- data_bp605 %>%
  filter(grepl(paste(immune_inflammatory_terms, collapse = "|"), Description, ignore.case = TRUE))

# write.table(data_bp605_filtered, file = "data/New_Inflammatory_markers_vomGenSetENrich.tsv", sep="\t", quote = F)






