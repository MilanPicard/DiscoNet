## DiscoNet R package
![alt text](https://github.com/MilanPicard/DiscoNet/blob/main/dance-floor.png | width=10)
***DiscoNet*** is an R package to automatically extract node features from multi-layered networks.

<img src="https://github.com/MilanPicard/DiscoNet/blob/main/Image1.png?raw=true" width="150" height="150">

![alt text](https://github.com/MilanPicard/DiscoNet/blob/main/Image1.png?raw=true)

Four categories of features can be extracted using  ***DiscoNet*** :  

 * **Propagation-based features**: simple one-to-one distances between nodes. Two types of algorithms are available, shortest paths and random walk with restart. Distances can be directed, inversely directed, or undirected.  
   Functions: `extract_by_shp`, `extract_by_shp_inv`, `extract_by_shp_ind`, `extract_by_rwr`,`extract_by_rwr_inv`, `extract_by_rwr_ind`  
 * **Topological metrics and similarities**: common metrics such as degree, betweenness, or centrality are extracted, as well as node similarity based on node neighborhoods.  
   Functions: `extract_topolo`.  
 * **Module-based**: clusters are identified and distances between nodes and clusters are computed. Also works with a user-wn list of pre-calculated clusters.  
   Function: `extract_cluster`.  
 * **Signature-based**: Based on a user-own list of nodes (signature), different proximities between *start_nodes* and these nodes can be calculated.  
   Functions: `extract_by_sig`.  

***DiscoNet*** also implements different variable selection methods (*Information Gain* & *Adapative LASSO*) in combination with bootstrapings to help in reducing the number of variable extracted from the network. On a medium size network (n_nodes ~= 10k), ***DiscoNet*** can easily produce more than 70k different variables based on its topology.


### Installation
Install the devtools package in R, then load it and install the latest stable version of DiscoNet from GitHub
```r
## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
## install DiscoNet
devtools::install_github("MilanPicard/DiscoNet")
## install DiscoNet AND the vignette
## Will take a good 5 minutes or more depending on your machine.
install_github("MilanPicard/DiscoNet", build_vignettes = TRUE)
```

### Tutorial
A vignette showing a basic use of the package is available : TargetRepositioningProstateCancer.
It can be accessed using:  
```r
vignette(TargetRepositioningProstateCancer)
```

To analyze your own network using DiscoNet, a simple PPI network will do, but larger multi-layer networks can be exploited as well as some functions can exploit multiple cores in parallel. 


### Maintainer
Milan Picard (milan.picard.1@ulaval.ca)
