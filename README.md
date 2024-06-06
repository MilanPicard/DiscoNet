## <img src="https://github.com/MilanPicard/DiscoNet/blob/main/Meta/dance-floor1.png?raw=true" width="40" height="40"> DiscoNet R package
***DiscoNet*** is an R package to automatically extract node features from multi-layered networks.

<img src="https://github.com/MilanPicard/DiscoNet/blob/main/Meta/Image1.png?raw=true" width="850" height="450">

Four categories of features can be extracted using  ***DiscoNet*** :  

 * **Propagation-based features**: simple one-to-one distances between nodes. Two types of algorithms are available, shortest paths and random walk with restart. Distances can be directed, inversely directed, or undirected.  
   Functions: `extract_by_shp`, `extract_by_shp_inv`, `extract_by_shp_ind`, `extract_by_rwr`,`extract_by_rwr_inv`, `extract_by_rwr_ind`  
 * **Topological metrics and similarities**: common metrics such as degree, betweenness, or centrality are extracted, as well as node similarity based on node neighborhoods.  
   Functions: `extract_topolo`.  
 * **Module-based**: clusters are identified and distances between nodes and clusters are computed. Also works with a user-wn list of pre-calculated clusters.  
   Function: `extract_cluster`.  
 * **Signature-based**: Based on a user-own list of nodes (signature), different proximities between *start_nodes* and these nodes can be calculated.  
   Functions: `extract_by_sig`.  

***DiscoNet*** also implements different **variable selection** methods to reduce the number of variable extracted from the network. They include a method based on Information gain (`InformationGain`) and a selection using Adaptative LASSO (`AdaLASSO`). Both can be used with bootstrapings to make the selection more robust and stable (`InformationGain_Bootstrap`, `AdaLASSO_Bootstrap`).  
 On a medium size network (n_nodes ~= 10k), ***DiscoNet*** can easily produce more than 70k different variables based on its topology.


### Installation
Install the devtools package in R, then load it and install the latest stable version of DiscoNet from GitHub
```r
## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
## install DiscoNet
devtools::install_github("MilanPicard/DiscoNet")
```

### Simple Use
```r
## load library
library(DiscoNet)

## Load your (own) network (igraph object)
A3LNetwork = DiscoNet::A3LNetwork

## Extract common node metrics and topological similarities
SomeFeatures = extract_topolo(Graph = A3LNetwork)

## Extract undirected Random Walks distance between every nodes in the network
RWRFeatures = extract_by_rwr_ind(Graph = A3LNetwork, 
                                 nCores = 1) # Add more cores for parallel computing
```

### Full tutorial
A vignette showing basic usage of the package is available : vignette/TargetRepositioningProstateCancer.Rmd.  
It can be accessed using:  
```r
vignette(TargetRepositioningProstateCancer)

## A the vignette's pdf is also available in vignette/
## you can also install DiscoNet WITH the vignette
## Will take a good 5 minutes or more depending on your machine.
install_github("MilanPicard/DiscoNet", build_vignettes = TRUE)
```

To analyze your own network using DiscoNet, a simple PPI network will do, but larger multi-layer networks can be exploited as well as some functions can exploit multiple cores in parallel. 


### Maintainer
Milan Picard (milan.picard.1@ulaval.ca)
