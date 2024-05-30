## DiscoNet R package

***DiscoNet*** is an R package to automatically extract node features from multi-layered biological networks to be exploited by machine learning algorithms for prediction, classification, unsupervised learning, etc.


Four categories of features can be extracted using  ***DiscoNet*** :  

 * **Propagation-based features**: simple one-to-one distances between *start_nodes* and other nodes in the network.  
The distances can be directed, inversely directed, or undirected. (`extract_by_shp`, `extract_by_shp_inv`, `extract_by_shp_ind`, `extract_by_rwr`,`extract_by_rwr_inv`, `extract_by_rwr_ind`).  
 * **Topological metrics and similarities**: common metrics such as degree, betweenness, or centrality are extracted, as well as node similarity based on node neighborhoods (`extract_topolo`).  
 * **Module-based**: Will identify clusters from the network and compute the distance between each *start_nodes* these clusters, also works with a user-wn list of pre-calculated clusters (`extract_cluster`).  
 * **Signature-based**: Based on a user-own list of nodes (signature), different distances between *start_nodes* and these nodes can be calculated 
 (`extract_by_sig`).  

### Installation
#### Latest GitHub Version
Install the devtools package in R, then load it and install the latest stable version of timeOmics from GitHub

```r
## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
## install DiscoNet
devtools::install_github("MilanPicard/DiscoNet")
```

### Maintainer
Milan Picard (milan.picard.1@ulaval.ca)
