## code to prepare the data necessary to run the vignette of DiscoNet

# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("tidyverse")
# BiocManager::install("igraph")
library(org.Hs.eg.db)
library(tidyverse)
library(igraph)


# Get prostate cancer signatures and known targets
# Wang G, Zhao D, Spring DJ, DePinho RA. Genetics and biology of prostate cancer. Genes Dev. 2018;32(17-18):1105-1140. doi:10.1101/gad.315739.118
Mutated_genes = c("APC", "AR", "ATM", "BRCA1", "BRCA2", "CHD1", "ERF", "ERG", "ETS2", "ETVs", "EZh2", "FOXA1", "IDH1", "KMT2A", "KMT2C", "KMT2D", "MYC", "MYCN", "NCOR1", "NCOR2", "NKX3", "PTEN", "RB1", "SETD2", "SETDB1", "SMAD4", "SMARCA1", "SMARCB1", "SPOP", "TP53")

# Myers JS, von Lersner AK, Robbins CJ, Sang QX. Differentially Expressed Genes and Signature Pathways of Human Prostate Cancer. PLoS One. 2015;10(12):e0145322. Published 2015 Dec 18. doi:10.1371/journal.pone.0145322
Under_expressed_genes = c("WFDC9","DEFB125","EDDM3B","PAEP","SEMG2","PATE4","EDDM3A","CRISP1","PATE1","DEFB127","AQP2","TMEM114","GRXCR1","SPINT3","CLDN2","SULT2A1","SPINK2","POU3F3","LCN1","PATE")
Over_expressed_genes = c("ANKRD30A","FEZF2","C6orf10","FOXG1","GC","VAX1","SSX2","FGB","SLC45A2","SPINK1","HOXC12","SCN1A","LOC284661","TFDP3","B3GNT6","FOXB2","NR2E1","XAGE1E","TBX10")

Signature_for_prostate_cancer = list(Mutated_genes = Mutated_genes,
                                     DEG = c(Under_expressed_genes, Over_expressed_genes))
Signature_for_prostate_cancer = lapply(Signature_for_prostate_cancer, function(x) paste0("Gene_", x))

# Yap, Timothy A et al. “Drug discovery in advanced prostate cancer: translating biology into therapy.” Nature reviews. Drug discovery vol. 15,10 (2016): 699-718. doi:10.1038/nrd.2016.120
Known_Target = c("PARP1","PARP2","FLT1","FLT4","EZH2","AR","ATM","SPINK1","CHD1","SPOP","AURKA","RB1","MYC")
Known_Target_prot = mapIds(org.Hs.eg.db, Known_Target, column = "UNIPROT", keytype = "SYMBOL")
Known_Target_prot = setNames(paste0("Prot_", Known_Target_prot), names(Known_Target_prot))


# Create Small example network
# This code creates the small network used as example in the r package DiscoNet
# The network is made up of three different layers:
#   - PPI layer (from StrinDB: https://string-db.org/cgi/download?sessionId=bjv6RHolIynV&species_text=Homo+sapiens, file = 9606.protein.links.v12.0.txt.gz)
# - Gene layer (genes whose protein products are present in tthe PPI network)
# - Gene Ontology layer (biological information to enhance the predictive power of the network)
#
# These layers interact with each other in the following way:
#   - PPI --> Gene (directed: transcription factors regulates the expression of genes)
# - Gene --> PPI (directed: a gene is translated into a protein)
# - GO -- Gene (undirected: the GO terms associated with a gene)

# Retrieve and process data
# Downloaded from stringdb, to run this code, you must download it on your own
PPI_mat = read.delim("../9606.protein.links.detailed.v12.0.txt", sep = " ")
prot_info = read_tsv("../9606.protein.info.v12.0.txt")

# Keep only validated interactions and convert ID name to gene SYMBOL
PPI_mat_clean = PPI_mat %>%
  dplyr::filter(experimental >= 700) %>%
  dplyr::mutate(protein1_symbol = prot_info$preferred_name[match(protein1, prot_info[[1]])]) %>%
  dplyr::mutate(protein2_symbol = prot_info$preferred_name[match(protein2, prot_info[[1]])]) %>%
  dplyr::select(protein1_symbol, protein2_symbol)
PPI_mat_clean = PPI_mat_clean[complete.cases(PPI_mat_clean), ]

# Find the proteins in proximity (order = 1) of known targets
PPI_graph = PPI_mat_clean %>%
  igraph::graph_from_data_frame()
PPI_to_keep = lapply(igraph::ego(PPI_graph, order = 1, nodes = Known_Target), names) %>% unlist %>% unique

# Sort original PPI to keep only targets and nearest nodes and convert to Uniprot IDs
PPI_mat_clean = PPI_mat_clean %>%
  dplyr::filter(protein1_symbol %in% PPI_to_keep | protein2_symbol %in% PPI_to_keep) %>%
  dplyr::mutate(protein1_uniprot = paste0("Prot_", mapIds(org.Hs.eg.db, protein1_symbol, "UNIPROT", "SYMBOL"))) %>%
  dplyr::mutate(protein2_uniprot = paste0("Prot_", mapIds(org.Hs.eg.db, protein2_symbol, "UNIPROT", "SYMBOL"))) %>%
  dplyr::filter(protein1_uniprot != "Prot_NA") %>%
  dplyr::filter(protein2_uniprot != "Prot_NA") %>%
  dplyr::select(protein1_uniprot, protein2_uniprot, protein1_symbol, protein2_symbol) %>%
  unique

# Sort TF to gene information
# packageVersion("dorothea") # ‘1.16.0’
TF_gene_mat = dorothea::entire_database %>%
  dplyr::filter(confidence == "A") %>%
  dplyr::mutate(tf_uniprot = paste0("Prot_", mapIds(org.Hs.eg.db, tf, "UNIPROT", "SYMBOL"))) %>%
  dplyr::mutate(target = paste0("Gene_", target)) %>%
  dplyr::filter(tf_uniprot %in% c(PPI_mat_clean$protein1_uniprot, PPI_mat_clean$protein2_uniprot)) %>%
  dplyr::select(tf_uniprot, target)

# Sort Gene to GO information
# packageVersion("GOfuncR") # ‘1.24.0’
GO_gene_mat = GOfuncR::get_anno_categories(str_remove(unique(c(PPI_mat_clean$protein1_symbol, PPI_mat_clean$protein2_symbol, TF_gene_mat$target)), "Gene_"))
GO_gene_mat = GO_gene_mat %>%
  dplyr::group_by(go_id) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::filter((count < 20) & (count > 5)) %>% # Keep GO of size between 5 and 20
  dplyr::mutate(gene = paste0("Gene_", gene)) %>%
  dplyr::mutate(go_id = str_replace(go_id, ":", "_")) %>%
  dplyr::select(1,2)

# Remove NAs
PPI_mat_clean = PPI_mat_clean[complete.cases(PPI_mat_clean), ] %>% unique
TF_gene_mat = TF_gene_mat[complete.cases(TF_gene_mat), ] %>% unique
GO_gene_mat = GO_gene_mat[complete.cases(GO_gene_mat), ] %>% unique

# Sort gene to protein information: directed interactions from gene to protein
list_of_proteins = unique(c(PPI_mat_clean$protein1_uniprot, PPI_mat_clean$protein2_uniprot, TF_gene_mat$tf_uniprot))
list_of_proteins_converted_to_symbols = mapIds(org.Hs.eg.db, str_remove(list_of_proteins, "Prot_"), "SYMBOL", "UNIPROT")
PGi_mat = utils::stack(list_of_proteins_converted_to_symbols) %>%
  dplyr::mutate(values = paste0("Gene_", values)) %>%
  dplyr::mutate(ind = paste0("Prot_", ind)) %>%
  dplyr::filter(values != "Gene_NA") %>%
  dplyr::filter(ind != "Prot_NA")

# Create the network
PPI_net = PPI_mat_clean[, c(1,2)] %>% igraph::graph_from_data_frame(directed = FALSE)
TFG_net = TF_gene_mat  %>% igraph::graph_from_data_frame(directed = TRUE)
GOG_net = GO_gene_mat %>% igraph::graph_from_data_frame(directed = FALSE)
PGi_net = PGi_mat %>% igraph::graph_from_data_frame(directed = TRUE)

clean_net = function (graph) {
  getcompo = components(graph)
  getverti = names(getcompo$membership[which(getcompo$membership == which.max(getcompo$csize))])
  graph = induced_subgraph(graph, getverti)
  return(simplify(graph))
}
A3LNetwork = as.directed(PPI_net, mode = "mutual")+TFG_net+as.directed(GOG_net, mode = "mutual")+PGi_net
A3LNetwork = clean_net(A3LNetwork)

A3LNetwork = igraph::set_vertex_attr(graph = A3LNetwork, name = "type", index = which(str_detect(V(A3LNetwork)$name, "Prot_")), value = "protein")
A3LNetwork = igraph::set_vertex_attr(graph = A3LNetwork, name = "type", index = which(str_detect(V(A3LNetwork)$name, "Gene_")), value = "gene")
A3LNetwork = igraph::set_vertex_attr(graph = A3LNetwork, name = "type", index = which(str_detect(V(A3LNetwork)$name, "GO_"))  , value = "go")

usethis::use_data(A3LNetwork, overwrite = TRUE)

# Create clusters from the PPI subnetwork.
all_proteins = vnames(A3LNetwork, pattern = "Prot_")
A3LNetwork_prot =igraph::induced_subgraph(A3LNetwork, vids = all_proteins)

# Extract clusters from the PPI network
A3LNetwork_prot_clust = extract_cluster(Graph = A3LNetwork_prot, start_nodes = all_proteins, only_cluster = TRUE)
# Only take clusters of node size higher than 10
A3LNetwork_prot_clust = A3LNetwork_prot_clust[nbrs(A3LNetwork_prot_clust) > 10] # Results in 53 modules
# Convert protein clusters to gene clusters.
Cluster_list_gene = lapply(A3LNetwork_prot_clust, function(clust) suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, str_remove(clust, "Prot_"),"SYMBOL","UNIPROT"))))
Cluster_list_gene = lapply(Cluster_list_gene, function(x) paste0("Gene_", x))

data_prostate_cancer = list(Targets = Known_Target_prot, Signatures = Signature_for_prostate_cancer, Clusters = A3LNetwork_prot_clust, Clusters_gene = Cluster_list_gene)
usethis::use_data(data_prostate_cancer, overwrite = TRUE)










