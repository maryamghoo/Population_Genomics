library(topGO)   # GO enrichment analysis
library(GO.db)   # GO term information
library(ggplot2) # plotting


##### define directories #####

# directory containing files prepared by prepare_go_enrichment.sh
data_dir <- "bash_results"

# directory for R results
output_dir <- "r_results"

dir.create(output_dir, showWarnings = FALSE)


##### read background genes #####

all_genes <- scan(
  file.path(
    data_dir,
    "all_genes.txt"
  ),
  what = "character",
  quiet = TRUE
)


##### read candidate genes #####

candidate_genes <- scan(
  file.path(
    data_dir,
    "candidate_genes.txt"
  ),
  what = "character",
  quiet = TRUE
)


##### read gene-to-GO mapping #####

geneID2GO <- readMappings(
  file = file.path(
    data_dir,
    "gene2go_map.txt"
  )
)


##### create candidate/background gene vector #####

# 1 = candidate gene
# 0 = background gene

geneList <- factor(
  as.integer(
    all_genes %in% candidate_genes
  )
)

names(geneList) <- all_genes


##### create topGO data object #####

# BP = Biological Process
# nodeSize = minimum number of annotated genes required for a GO term

GOdata <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = geneList,
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO,
  nodeSize = 3
)


##### run GO enrichment #####

# weight01 accounts for the hierarchical structure of the GO graph.
# Fisher's exact test evaluates enrichment.
result <- runTest(
  GOdata,
  algorithm = "weight01",
  statistic = "fisher"
)


##### obtain all tested GO terms #####

go_scores <- score(result)
length(go_scores)

# Number with nominal/raw p < 0.05
sum(go_scores < 0.05, na.rm = TRUE)


##### create complete GO result table #####

# Include all GO terms returned by the test.
# FDR must be calculated before filtering for significance.

allRes <- GenTable(
  GOdata,
  weightFisher = result,
  orderBy = "weightFisher",
  topNodes = length(go_scores)
)


##### restore exact numeric p-values #####

# GenTable may format very small p-values as strings.
# Use the original numeric scores from topGO instead.

allRes$weightFisher <- as.numeric(go_scores[allRes$GO.ID])


##### add GO term names #####

allRes$Term <- Term(
  GOTERM[allRes$GO.ID]
)


##### calculate FDR #####

# Benjamini-Hochberg FDR correction across ALL GO terms
# returned by the enrichment test.

allRes$FDR <- p.adjust(allRes$weightFisher, method = "fdr")


##### identify raw-significant GO terms #####

raw_sig_GO <- subset(allRes, weightFisher < 0.05)


##### identify FDR-significant GO terms #####

fdr_sig_GO <- subset(allRes, FDR < 0.05)


##### print results #####

cat(
  "Total GO terms tested:",
  nrow(allRes),
  "\n"
)

cat(
  "Raw p < 0.05:",
  nrow(raw_sig_GO),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(fdr_sig_GO),
  "\n"
)

fdr_sig_GO


##### save complete results #####

write.csv(
  allRes,
  file.path(
    output_dir,
    "GO_enrichment_all_results.csv"
  ),
  row.names = FALSE
)


##### save raw-significant GO terms #####

write.csv(
  raw_sig_GO,
  file.path(
    output_dir,
    "GO_enrichment_raw_p05.csv"
  ),
  row.names = FALSE
)


##### save FDR-significant GO terms #####

write.csv(
  fdr_sig_GO,
  file.path(
    output_dir,
    "GO_enrichment_FDR05.csv"
  ),
  row.names = FALSE
)


##### prepare enriched GO terms for plotting #####

# Plot terms with raw p < 0.05.
# Terms that also pass FDR < 0.05 are automatically distinguished.

plot_data <- raw_sig_GO

plot_data$FDR_significant <- plot_data$FDR < 0.05


##### GO enrichment plot #####

go_plot <- ggplot(
  plot_data,
  aes(
    x = reorder(
      Term,
      -log10(weightFisher)
    ),
    y = -log10(weightFisher),
    fill = FDR_significant
  )
) +
  geom_col() +
  
  coord_flip() +
  
  scale_fill_manual(
    values = c(
      "FALSE" = "steelblue",
      "TRUE" = "red"
    ),
    labels = c(
      "Raw p < 0.05",
      "FDR < 0.05"
    ),
    name = "Significance"
  ) +
  
  labs(
    x = "GO Term",
    y = expression(-log[10](italic(p))),
    title = "Enriched GO Biological Processes"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    plot.title = element_text(
      hjust = 0.5
    )
  )

go_plot


##### save GO enrichment plot #####

# PNG
ggsave(
  file.path(
    output_dir,
    "GO_enrichment.png"
  ),
  plot = go_plot,
  width = 9,
  height = 7,
  dpi = 600
)


# PDF
ggsave(
  file.path(
    output_dir,
    "GO_enrichment.pdf"
  ),
  plot = go_plot,
  width = 9,
  height = 7
)


# SVG
ggsave(
  file.path(
    output_dir,
    "GO_enrichment.svg"
  ),
  plot = go_plot,
  width = 9,
  height = 7
)
