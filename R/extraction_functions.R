#' Check arguments
#'
#' Check arguments given within called functions
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param Signature_list Optional. Your list of signatures.
#' @param Cluster_list Optional. Your list of clusters previously calculated or not.
#' @param dist_matrix Optional. Your distance matrix previously calculated.
#' @return Stop the called function if some errors are detected, modify some arguments if necessary.
check_for_validity = function(Graph, start_nodes, Signature_list = NULL, Cluster_list = NULL, dist_matrix = NULL){
  # Basic verification
  if(!methods::is(Graph, "igraph")){
    stop("The Graph is not an igrah object.")
  }
  if(!is.character(start_nodes)){
    stop("start_nodes is not a character vector.")
  }
  if(!all(start_nodes %in% vnames(Graph))){
    stop("Not all nodes in start_nodes are present in the network!")
  }
  if(nbr(start_nodes) != nbrunique(start_nodes)){
    stop("start_nodes have multiple identical values.")
  }
  # dist_matrix mparameter
  if(!is.null(dist_matrix)){
    if(!is.data.frame(dist_matrix)){
      stop("dist_matrix is not a dataframe!")
    }
    if(names(dist_matrix)[1] != "Target"){
      stop("The first column of dist_matrix is not 'Target'. Are you sure it is the right matrix")
    }
    if(nbr(dist_matrix$Target) != nbrunique(dist_matrix$Target)){
      stop("Some target in dist_matrix appear multiple times!")
    }
    if((!all(names(dist_matrix)[-1] %in% vnames(Graph))) | !all(dist_matrix$Target %in% vnames(Graph))){
      stop("dist_matrix contains targets or nodes that are not in the network")
    }
    if(!all(dist_matrix$Target %in% start_nodes)){
      stop("dist_matrix and start_nodes do not overlap. Is it the right dist_matrix?")
    }
  }
  # Check Signature_list
  if(!is.null(Signature_list)){
    if(is.character(Signature_list)){
      Signature_list = list(Signature_list)
      names(Signature_list) = "Signature_1"
    } else if(is.list(Signature_list)){
      if(is.null(names(Signature_list)) | any(nchar(names(Signature_list)) == 0)){
        names(Signature_list) = paste0("Signature_", 1:nbr(Signature_list))
      }
    } else {
      stop("Signature_list is not a named list of nodes")
    }
    if(!all(unlist(Signature_list) %in% vnames(Graph))){
      stop("Some nodes in Signature_list are not part of the network")
    }
    if(any(nbrs(Signature_list) < 2)){
      stop("Please remove clusters of size < 2")
    }
  }
  # Check Cluster_list
  if(!is.null(Cluster_list)){
    if(is.character(Cluster_list)){
      Cluster_list = list(Cluster_list)
      names(Cluster_list) = "Cluster_1"
    } else if(is.list(Cluster_list)){
      if(is.null(names(Cluster_list)) | any(nchar(names(Cluster_list)) == 0)){
        names(Cluster_list) = paste0("Cluster_", 1:nbr(Signature_list))
      }
    } else {
      stop("Cluster_list is not a named list of nodes")
    }
    if(!all(unlist(Cluster_list) %in% vnames(Graph))){
      stop("Some nodes in Cluster_list are not part of the network")
    }
    if(any(nbrs(Cluster_list) < 2)){
      stop("Please remove clusters of size < 2")
    }
  }
  return(list(Signature_list = Signature_list, Cluster_list = Cluster_list))
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param Signature_list A named list of nodes that are in the network.
#' @param Cluster_list Optional. A Named list of clusters previously calculated.
#' @param dist_matrix Optional. A distance matrix previously calculated.
#' @param return_clusters Default: FALSE.  Wether or not to return the clusters utilized
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different features extracted from the network.
#' @examples
#' # Don't run
#' # Signature = list(Signature1 = c("node7", "node9", "node13"))
#' # Dist_signature = extract_by_sig(Network, c("node1", "node1"), Signature)
#' @export
extract_by_sig = function(Graph, start_nodes = vnames(Graph), Signature_list, Cluster_list = NULL, dist_matrix = NULL, return_clusters = FALSE, verbose = TRUE){
  if(verbose){message("extract_by_sig is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes, Cluster_list = Cluster_list, Signature_list = Signature_list, dist_matrix = dist_matrix)
  no_cluster = FALSE
  Cluster_list = checker$Cluster_list
  Signature_list = checker$Signature_list
  # If the signature given is a list of gene nodes, then they must have at least some overlap with the clusters given.
  if(!is.null(Cluster_list)){
    if(!any(unlist(Signature_list) %in% unlist(Cluster_list))){
      warning("The signature given and the list of cluster given have no nodes in common! Only a partial set of signature-based features can be calculated...")
      no_cluster = TRUE
    }
  }
  # Create list of cluster if none is given
  if(is.null(Cluster_list)){
    if(verbose){message("Calculating clusters on full network...")}
    Cluster_list = extract_cluster(Graph = Graph, start_nodes = start_nodes, Cluster_list = NULL, dist_matrix = dist_matrix, only_cluster = TRUE, verbose = FALSE)
    if(!any(unlist(Signature_list) %in% unlist(Cluster_list))){
      warning("The signature given and the list of clusters calculated have no nodes in common! Only a partial set of signature-based features can be calculated...")
      no_cluster = TRUE
    }
  }
  # Crate distance matrix if none if given
  if(is.null(dist_matrix)){
    if(verbose){message("Calculating distance_matrix...")}
    dist_matrix = extract_by_shp_ind(Graph = Graph, start_nodes = start_nodes, verbose = FALSE)
  }
  dist_matrix = dist_matrix %>% dplyr::slice(match(start_nodes, Target))
  if(verbose){message("Calculating signature-based features...")}
  Results_for_each_signatures = list()
  for(i in 1:nbr(Signature_list)){
    Signature = Signature_list[[i]]
    Signature_name = names(Signature_list)[i]
    # Create features
    # Are targets part of a signature
    Part_of_signature = as.numeric(start_nodes %in% Signature)
    ## Are target neighborhood part of a signature
    start_nodes_neigh = lapply(igraph::ego(graph = Graph, order = 2, nodes = start_nodes, mode = "out", mindist = 1), names)
    names(start_nodes_neigh) = start_nodes
    Close_to_signature = round(sapply(start_nodes_neigh, function(x) sum(x %in% Signature)/length(Signature))*100, 3)
    # What is the proximity between a target and a full signature
    Distance_to_full_signature = rowMeans(dist_matrix[, intersect(names(dist_matrix), Signature)], na.rm = TRUE)
    # What is the prorximity between a target and each node in the signature
    Distance_to_each_signature_node = dist_matrix[, intersect(Signature, colnames(dist_matrix))]
    # What is the proximity between a target and a nodes of a same signature within the same cluster.
    Clustered_signature = lapply(Cluster_list, function(clust) intersect(clust, Signature))
    Clustered_signature = Clustered_signature[nbrs(Clustered_signature) >= 1]
    if(no_cluster | nbr(Clustered_signature) == 0){
      Results_for_one_signature = data.frame(Part_of_signature, Close_to_signature, Distance_to_full_signature, Distance_to_each_signature_node)
      names(Results_for_one_signature) = paste(Signature_name,names(Results_for_one_signature), sep = "_")
      add(Results_for_each_signatures, Results_for_one_signature)
      break()
    }
    Clustered_signature_av_distance = lapply(Clustered_signature, function(cl) if(nbr(cl) == 1){dist_matrix[, intersect(cl, colnames(dist_matrix))]} else {rowMeans(dist_matrix[, intersect(cl, colnames(dist_matrix))], na.rm = TRUE)})
    Results_for_one_signature = data.frame(Part_of_signature, Close_to_signature, Distance_to_full_signature, Clustered_signature_av_distance, Distance_to_each_signature_node)
    names(Results_for_one_signature) = paste(Signature_name,names(Results_for_one_signature), sep = "_")
    add(Results_for_each_signatures, Results_for_one_signature)
  }
  Results_feature_matrix = data.frame(Target = dist_matrix$Target, Results_for_each_signatures)
  rownames(Results_feature_matrix) = NULL
  if(verbose){message("Done.")}
  if(return_clusters){
    return(list(Results_feature_matrix = Results_feature_matrix, Cluster_list = Cluster_list))
  }
  return(Results_feature_matrix)
}


#' Extract cluster-based features
#'
#' Extract cluster-based features from a biological network, including distance to start_nodes to each clusters.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param Cluster_list Optional. A Named list of clusters previously calculated. If not, it will calculate its own from the network
#' @param dist_matrix Optional. A distance matrix previously calculated.
#' @param only_cluster Default: FALSE.  If TRUE, the function only calculates and return the modules identified without distance calculations.
#' @param return_clusters Default: FALSE.  Weather or not to return the clusters utilized
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different features extracted from the network.
#' @examples
#' # Don't run
#' # Dist_clusters = extract_cluster(Network, c("node1", "node1"), return_clusters = TRUE)
#' @export
extract_cluster = function(Graph, start_nodes = vnames(Graph), Cluster_list = NULL, dist_matrix = NULL, only_cluster = FALSE, return_clusters = FALSE, verbose = TRUE){
  if(verbose){message("extract_cluster is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes, Cluster_list = Cluster_list, dist_matrix = dist_matrix)
  Cluster_list = checker$Cluster_list

  if(is.null(Cluster_list)){
    if(verbose){message("Calculating clusters on full network...")}
    Graph_undirected = igraph::as.undirected(Graph, mode = "collapse")
    # Graph_prot = igraph::as.undirected(igraph::induced.subgraph(Graph, vnames(Graph, type = "Prot_")))
    Cluster_list_walktrap = igraph::cluster_walktrap(graph = Graph_undirected, steps = 4)[TRUE]
    Cluster_list_louvain = igraph::cluster_louvain(graph = Graph_undirected, resolution = 2)[TRUE]
    Cluster_list = c(Cluster_list_walktrap, Cluster_list_louvain)
    Cluster_list = Cluster_list[nbrs(Cluster_list) > 1]
    names(Cluster_list) = paste0("Cluster_", 1:nbr(Cluster_list))
  }
  if(only_cluster){
    if(verbose){message("Done.")}
    return(Cluster_list)
  }
  # Calculating wether starting nodes are in clusters.
  Within_clusters = lapply(Cluster_list, function(cl) as.numeric(start_nodes %in% cl))
  names(Within_clusters) = paste0("IsIn_",names(Cluster_list))

  if(is.null(dist_matrix)){
    if(verbose){message("Calculating distance_matrix...")}
    dist_matrix = extract_by_shp_ind(Graph = Graph, start_nodes = start_nodes, verbose = FALSE)
  }
  dist_matrix = dist_matrix %>% dplyr::slice(match(start_nodes, Target))
  # Calculate distances of each strart nodes to each clusters given by extract_clusters using the median distance of all the shortest path from the start node to each member of the cluster.
  if(verbose){message("Calculating distances to clusters...")}
  Clusters_Mean = lapply(Cluster_list, function(cl) rowMeans(dist_matrix[, intersect(cl, colnames(dist_matrix))], na.rm = TRUE))
  Clusters_Mean = as.data.frame(do.call(cbind, Clusters_Mean))
  names(Clusters_Mean) = paste0("DistanceTo_",names(Clusters_Mean))
  Clusters_results = data.frame(Target = dist_matrix$Target, Within_clusters, Clusters_Mean)
  if(verbose){message("Done.")}
  if(return_clusters){
    return(list(Clusters_results = Clusters_results, Cluster_list = Cluster_list))
  }
  return(Clusters_results)
}


#' Extract topological metrics-based features
#'
#' Extract topological-based features from a biological network including degree, betweenness, similarity, etc.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different features extracted from the network.
#' @examples
#' # Don't run
#' # Topological_metrics = extract_topolo(Network, c("node1", "node1"))
#' @export
extract_topolo = function(Graph, start_nodes = vnames(Graph), verbose = TRUE){
  if(verbose){message("extract_topolo is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  if(verbose){message("Calculating general topological metrics...")}
  results = list(
    Degree_in = igraph::degree(Graph, v = start_nodes, mode = "in"),
    Degree_out = igraph::degree(Graph, v = start_nodes, mode = "out"),
    Degree_all = igraph::degree(Graph, v = start_nodes, mode = "all"),
    Closeness_out = signif(igraph::closeness(Graph, vids = start_nodes, mode = "out"), 3),
    Closeness_all = signif(igraph::closeness(Graph, vids = start_nodes, mode = "all"), 3),
    Eccentricity_out = igraph::eccentricity(Graph, vids = start_nodes, mode = "out"),
    Eccentricity_all = igraph::eccentricity(Graph, vids = start_nodes, mode = "all"),
    Eigen_centrality = signif(igraph::eigen_centrality(Graph, scale = 1)$vector[start_nodes], 3),
    Betweenness = signif(igraph::betweenness(Graph, v = start_nodes), 3),
    Harmonic_centrality = signif(igraph::harmonic_centrality(graph = Graph, vids = start_nodes, mode = "out"), 4))

  # similarity measures
  if(verbose){message("Calculating topological similarities...")}
  Similarity.invlogweighted = igraph::similarity.invlogweighted(graph = Graph, mode = "out", vids = start_nodes)
  rownames(Similarity.invlogweighted) = start_nodes
  colnames(Similarity.invlogweighted) = igraph::V(Graph)$name
  diag(Similarity.invlogweighted[start_nodes, start_nodes]) = 1
  colnames(Similarity.invlogweighted) = paste0("Similarity_", colnames(Similarity.invlogweighted))
  Similarity.invlogweighted = signif(Similarity.invlogweighted, 3)

  names(results) = paste0("Topology_", names(results))
  results = c(results, list(Similarity.invlogweighted))

  results = as.data.frame(do.call(cbind, results))[start_nodes, ]
  results = cbind(Target = rownames(results), results)
  rownames(results) = NULL
  if(verbose){message("Done.")}
  return(results)
}


#' Internal function: shortest path distance
#'
#' An internal function that calculate the shortest path distance from starting_nodes to any other nodes in the network.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param to Optional. Distance to what nodes.
#' @param mode Optional. The directionality of the distance measures, can be "out", "in", or "all.
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different distance features extracted from the network.
extract_by_shp.intern = function(Graph, start_nodes, to = vnames(Graph), mode = "out", verbose = TRUE){
  # if(verbose){message("Calculate shortest path distances...")}
  node_names = igraph::V(Graph)$name
  dist_matrix = igraph::distances(graph = Graph, v = start_nodes, to = to, mode = mode)
  # when no path exist, sometimes create Inf values
  dist_matrix = signif(1/dist_matrix, 5)
  diag(dist_matrix) = 1
  # data.frames are better
  dist_matrix = as.data.frame(dist_matrix)
  colnames(dist_matrix) = node_names
  dist_matrix = cbind(Target = rownames(dist_matrix), dist_matrix)
  rownames(dist_matrix) = NULL
  # if(verbose){message("Done.")}
  return(dist_matrix)
}


#' Shortest path distance: to
#'
#' Calculate the shortest path distance from starting_nodes to all other nodes in the network.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_shp = extract_by_shp(Network, c("node1", "node1"))
#' @export
extract_by_shp = function(Graph, start_nodes, verbose = TRUE){
  if(verbose){message("extract_by_shp is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  result = extract_by_shp.intern(Graph = Graph, start_nodes = start_nodes, to = vnames(Graph), mode = "out", verbose = FALSE)
  if(verbose){message("Done.")}
  return(result)
}


#' Shortest path distance: from
#'
#' Calculate the shortest path distance from all other nodes to starting_nodes in the network.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_shp_inv = extract_by_shp_inv(Network, c("node1", "node1"))
#' @export
extract_by_shp_inv = function(Graph, start_nodes, verbose = TRUE){
  if(verbose){message("extract_by_shp_inv is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  results = extract_by_shp.intern(Graph = Graph, start_nodes = start_nodes, to = vnames(Graph), mode = "in", verbose = FALSE)
  if(verbose){message("Done.")}
  return(results)

}


#' Shortest path distance
#'
#' Calculate the shortest path distance between starting_nodes and all other nodes in the network in a undirectional manner.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_shp_ind = extract_by_shp_ind(Network, c("node1", "node1"))
#' @export
extract_by_shp_ind = function(Graph, start_nodes, verbose = TRUE){
  if(verbose){message("extract_by_shp_ind is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  Graph_ind = igraph::as.undirected(Graph, mode = "collapse")
  results = extract_by_shp.intern(Graph = Graph_ind, start_nodes = start_nodes, to = vnames(Graph), mode = "all", verbose = FALSE)
  if(verbose){message("Done.")}
  return(results)
}


#' Random Walk distance: to
#'
#' Calculate the random walk distance from  starting_nodes to  all other nodes in the network.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @param nCores Default: 1 The number of cores used. Keep 1 if you want sequential calculations, >1 for parallel.
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_rwr = extract_by_rwr(Network, c("node1", "node1"))
#' @export
extract_by_rwr = function(Graph, start_nodes, nCores = 1, verbose = TRUE){
  if(verbose){message("extract_by_rwr is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  extract_by_rwr.intern(Graph = Graph, start_nodes = start_nodes, nCores = nCores, verbose = verbose)
}


#' Random Walk distance: from
#'
#' Calculate the random walk distance from all other nodes to starting_nodes in the network.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @param nCores Default: 1 The number of cores used. Keep 1 if you want sequential calculations, >1 for parallel.
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_rwr_inv = extract_by_rwr_inv(Network, c("node1", "node1"))
#' @export
extract_by_rwr_inv = function(Graph, start_nodes, nCores = 1, verbose = TRUE){
  if(verbose){message("extract_by_rwr_inv is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  Graph_inv = igraph::reverse_edges(Graph)
  extract_by_rwr.intern(Graph = Graph_inv, start_nodes = start_nodes, nCores = nCores, verbose = verbose)
}


#' Random Walk distance
#'
#' Calculate the random walk distance between starting_nodes and all other nodes to starting_nodes in the network in a undirected manner.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @param nCores Default: 1 The number of cores used. Keep 1 if you want sequential calculations, >1 for parallel.
#' @return Return a dataframe where columns are different distance features extracted from the network.
#' @examples
#' # Don't run
#' # extract_by_rwr_inv = extract_by_rwr_inv(Network, c("node1", "node1"))
#' @export
extract_by_rwr_ind = function(Graph, start_nodes, nCores = 1, verbose = TRUE){
  if(verbose){message("extract_by_rwr_ind is running...")}
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  Graph_ind = igraph::as.undirected(Graph, mode = "collapse")
  extract_by_rwr.intern(Graph = Graph_ind, start_nodes = start_nodes, nCores = nCores, verbose = verbose)
}


#' Internal function: Random Walk distance
#'
#' Calculate the random walk between two sets of nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @param nCores Default: 1 The number of cores used. Keep 1 if you want sequential calculations, >1 for parallel.
#' @return Return a dataframe where columns are different distance features extracted from the network.
extract_by_rwr.intern = function(Graph, start_nodes, nCores, verbose){
  checker = check_for_validity(Graph = Graph, start_nodes = start_nodes)
  restart = .5
  # Prepare adjacency matrix for later computation
  if(verbose){message("Calculating adjacency matric...")}
  ADJ = get.adjacency_matrix(Graph, attr = NULL)
  if(verbose){message("Calculating Random Walks with Restart...")}
  if(nCores == 1){
    result_list = list()
    for(i in 1:nbr(start_nodes)){
      if(verbose){
        if(i%%1 == 0){
          cat('\r', paste0(i,'/',nbr(start_nodes), " seeds"))
        }
      }
      res = RWR(Graph = Graph, Seeds = start_nodes[i], adjacency = ADJ, r = restart)
      add(result_list, res$Score[match(vnames(Graph), res$Node)])
    }
    df = do.call(rbind, result_list)
    df = as.data.frame(df)
  } else {
    df = as.data.frame(matrix(ncol = igraph::vcount(Graph), nrow = nbr(start_nodes)))
    results = parallel::mclapply(1:nbr(start_nodes), function(i) {if(i%%50 == 0 & verbose){cat('\r', paste0(i,'/',nbr(start_nodes), " seeds"))} ; RWR(Graph = Graph, Seeds = start_nodes[i], adjacency = ADJ, r = restart)}, mc.cores = nCores)

    if(verbose){message("\n Agglomerating results...")}
    result_list = list()
    for(i in 1:nbr(results)){
      add(result_list, results[[i]]$Score[match(vnames(Graph), results[[i]]$Node)])
    }
    df = do.call(rbind, result_list)
    df = as.data.frame(df)
  }
  rownames(df) = start_nodes
  # Remove seeds from the columns features
  df[is.na(df)] = restart
  df = as.data.frame(scales::rescale(as.matrix(df), c(0,1e6)))
  df = cbind(Target = rownames(df), df)
  colnames(df)[-1] = vnames(Graph)
  rownames(df) = NULL
  if(verbose){message("\nDone.")}
  return(df)
}


#' Internal function: calculate adjacency matrix
#'
#' Calculate the adjacency matrix of a network
#' @param graph An igraph object
#' @param attr A weight attributes for the edges, not functional yet.
#' @return The length of the unlist between v1 and v2
get.adjacency_matrix = function(graph, attr = NULL){
  adjacency = igraph::get.adjacency(graph, type = "both", attr = attr, sparse = igraph::getIgraphOpt("sparsematrices"))
  rown = colnames(adjacency)
  coln = rownames(adjacency)
  cla = class(adjacency)
  coeff = Matrix::rowSums(adjacency)
  coeff[coeff == 0] = 1
  adjacency = adjacency/coeff
  colnames(adjacency) = coln
  rownames(adjacency) = rown
  return(methods::as(SparseM::t(adjacency), cla))
}


#' Internal function: random Walk function
#'
#' Calculate a random walk on the network
#' @param Graph An igraph object
#' @param Seeds Starting_nodes
#' @param r Restart probability
#' @param attr Edge weight, not functional yet.
#' @param nCores Default: 1 The number of cores used. Keep 1 if you want sequential calculations, >1 for parallel.
#' @param adjacency Optional. A precomputed adjacency matrix
#' @param stop_delta Optional. When to stop the random walks
#' @param verbose Wether to print internal calculations
#' @return The length of the unlist between v1 and v2
RWR = function(Graph, Seeds, r = .5, attr = NULL, nCores = 1, adjacency = NULL,  verbose = FALSE, stop_delta = 1e-10){
  if(!igraph::is.igraph(Graph)){
    stop("Not igraph object")
  }
  if(!is.character(Seeds)){
    stop("Seeds not characters")
  }
  Seeds = intersect(vnames(Graph), Seeds)
  if(length(Seeds) == 0){
    stop("Seeds not present in graph")
  }
  if(r < 0 | r > 1){
    stop("Restart probability not recognized")
  }
  if(is.null(adjacency)){
    get.adjacency_matrix = function(graph, attr = NULL){
      adjacency = igraph::get.adjacency(graph, type = "both", attr = attr, sparse = igraph::getIgraphOpt("sparsematrices"))
      rown = colnames(adjacency)
      coln = rownames(adjacency)
      cla = class(adjacency)
      coeff = Matrix::rowSums(adjacency)
      coeff[coeff == 0] = 1
      adjacency = adjacency/coeff
      colnames(adjacency) = coln
      rownames(adjacency) = rown
      return(methods::as(SparseM::t(adjacency), cla))
    }
    adjacency = get.adjacency_matrix(graph = Graph, attr = attr)
  }
  prox_vector = matrix(0, nrow = ncol(adjacency), ncol = 1)
  prox_vector[which(colnames(adjacency) %in% Seeds)] = 1/length(Seeds)
  delta=1
  restart_vector = prox_vector
  while(delta >= stop_delta){
    old_prox_vector = prox_vector
    prox_vector = (1 - r) * (adjacency %*% prox_vector) + r * restart_vector
    delta = sqrt(sum((prox_vector - old_prox_vector)^2))
  }
  if(isTRUE(verbose)){message("Random Walk with restart finished")}
  Results = data.frame(Node = rownames(prox_vector), Score = prox_vector[, 1])
  Results = Results[order(Results$Score, decreasing = T), ]
  `rownames<-`(Results[-which(Results$Node %in% Seeds), ], NULL)
}

#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param Signature_list A named list of nodes that are in the network.
#' @param FC_vectors A named list of fold changes values for each node in the network. 0 are put where fold changes is missing
#' @param topo_results The results returned by the function extract_topolo
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return a dataframe where columns are different features extracted from the network.
#' @examples
#' # Don't run
#' # Signature = list(Signature1 = c("node7", "node9", "node13"))
#' # Dist_signature = extract_more_si_features(Network, c("node1", "node1"), Signature)
#' @export

extract_more_si_features = function(Graph, start_nodes = vnames(Graph), Signature_list, FC_vectors, topo_results, verbose = TRUE){
  if(verbose){message("extract_more_si_features is running...")}
  Signature_results_list = list()
  for(i in 1:nbr(Signature_list)){
    sign = Signature_list[[i]]
    sign_name = names(Signature_list)[i]
    
    # Calculate he interconnectivity between nodes and the signature
    interconnectivity_result = sapply(start_nodes, function(node) interconnectivity(Graph, node, degs = sign)) %>% 
      `names<-`(start_nodes) %>% stack %>%  dplyr::select(2,1) %>% `colnames<-`(c("Target", "Interconnectivity"))
    
    # Are node in signature in the neighborhood
    DEG_in_neighborhood_result   = lapply(start_nodes, function(node) DEG_in_neighborhood(Graph, node = node, degs = sign)) %>% 
      do.call(rbind, .) %>% as.data.frame %>% mutate(Target = start_nodes) %>% dplyr::select(3,1,2)
    
    # What are the topological properties of the node in subgraph of the signature
    subgraph_topo = extract_topolo_signature_graph(Graph = Graph, start_nodes = start_nodes, Signature = sign)
    
    # Merge results
    Signature_results = Reduce(function(x, y) full_join(x, y, by = "Target"), list(interconnectivity_result, DEG_in_neighborhood_result, subgraph_topo))
    names(Signature_results)[-1] = paste0(sign_name, "_", names(Signature_results)[-1])
    add(Signature_results_list, Signature_results)
  }
  Signature_results_list_merged = Reduce(function(x, y) full_join(x, y, by = "Target"), Signature_results_list)
  
  neighbordhood_scoring_result_list = list()
  for(i in 1:nbr(FC_vectors)){
    fc = FC_vectors[[i]]
    fc_name = names(FC_vectors)[i]
    
    fc = fc %>% stack %>% `colnames<-`(c("fc","node"))
    V(Graph)$fc = left_join(data.frame(name = V(Graph)$name), fc, by = c("name" = "node")) %>% pull(fc) %>% replace_na(0)
    
    # Simply retrieve the fold changes
    foldi = data.frame(Target = start_nodes, fold_change = V(Graph)[start_nodes]$fc)
    
    # Composite features based on topological properties and fold change
    topo_features_with_foldchange = topo_metric_with_foldchange(Graph = Graph, start_nodes = start_nodes, topo_results = topo_results) 
    
    # What are the topological properties of the signature nodes in the neighborhood
    neighbordhood_scoring_result_1 = lapply(start_nodes, function(node) topo_metric_of_neighboring_signature(Graph = Graph, node = node, degs = sign, topo_results = topo_results, order = 1)) %>% 
      do.call(rbind, .) %>% as.data.frame %>% mutate(Target = start_nodes) %>% dplyr::select(Target,everything())
    
    neighbordhood_scoring_result_2 = lapply(start_nodes, function(node) topo_metric_of_neighboring_signature(Graph = Graph, node = node, degs = sign, topo_results = topo_results, order = 2)) %>% 
      do.call(rbind, .) %>% as.data.frame %>% mutate(Target = start_nodes) %>% dplyr::select(Target,everything())
    
    names(foldi)[-1] = paste0(fc_name, "_", names(foldi)[-1])
    names(neighbordhood_scoring_result_1)[-1] = paste0(fc_name, "_", names(neighbordhood_scoring_result_1)[-1])
    names(neighbordhood_scoring_result_2)[-1] = paste0(fc_name, "_", names(neighbordhood_scoring_result_2)[-1])
    names(topo_features_with_foldchange)[-1] = paste0(fc_name, "_", names(topo_features_with_foldchange)[-1])
    
    results_per_foldchanges = Reduce(function(x, y) full_join(x, y, by = "Target"), list(foldi, topo_features_with_foldchange, neighbordhood_scoring_result_1, neighbordhood_scoring_result_2))
    
    add(neighbordhood_scoring_result_list, results_per_foldchanges)
  }
  neighbordhood_scoring_result_list_merged = Reduce(function(x, y) full_join(x, y, by = "Target"), neighbordhood_scoring_result_list)
  results = full_join(Signature_results_list_merged, neighbordhood_scoring_result_list_merged, by = "Target")
  if(verbose){message("Done.")}
  return(results)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param topo_results The results returned by the function extract_topolo
#' @return Return something

topo_metric_with_foldchange = function(Graph, start_nodes, topo_results){
  topo_results_metric = topo_results[str_detect(colnames(topo_results), "Similarity_", T)]
  topo_results_metric = topo_results_metric %>% filter(Target %in% start_nodes) %>% arrange(Target)
  start_nodes_fc = V(Graph)[sort(start_nodes)]$fc
  results = do.call(rbind, lapply(1:nrow(topo_results_metric), function(i) topo_results_metric[i, -1] * start_nodes_fc[i])) %>% 
    mutate(Target = sort(start_nodes), .before = 1)
  return(results)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param topo_results The results returned by the function extract_topolo
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something

topo_metric_of_neighboring_signature = function(Graph, node, degs, topo_results = NULL, order = 1){
  ######## get neighbors of node
  if(is.null(topo_results)){stop("Please provide the data frame from extract_topolo")}  
  neighbors = names(ego(graph = Graph, order = order, mindist = 1, nodes = node)[[1]])
  
  ######## Extract topological metrics about the neighbooring DEGS
  topo_results_metric = topo_results[str_detect(colnames(topo_results), "Similarity_", T)]
  topo_results_metric = topo_results_metric %>% filter(Target  %in% intersect(degs, neighbors))
  mean_metric = apply(topo_results_metric[, -1], 2, function(x) suppressWarnings(mean(x, na.rm = TRUE)))
  max_metric = apply(topo_results_metric[, -1], 2, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  
  mean_metric[!is.finite(mean_metric)] = 0
  max_metric[!is.finite(max_metric)] = 0
  
  names(mean_metric) = paste0(names(mean_metric),"_neighboring_signature_order_",order, "_mean")
  names(max_metric) = paste0(names(max_metric),"_neighboring_signature_order_",order, "_max")
  
  ######## Extract their fold changes
  fc_v = V(Graph)[node]$fc # for the node studied.
  fc_nei = V(Graph)[neighbors]$fc
  
  # Discordant Neighboring Factor
  if(fc_v == 0){
    dnf = 0
  } else {
    sign_of_neibhoring = sign(sum(fc_nei)/fc_v)
    dnf = sign_of_neibhoring*log(1+abs(sum(fc_nei)/fc_v), base = 1.15)
  }
  names(dnf) = paste0("discordant_neighborin_factor_order_", order)
  
  # Plain average of neighbors fc, without taking into account their direction
  fc_nei_avg = mean(fc_nei)
  fc_nei_avg_sum = (fc_v+fc_nei_avg)/2
  fc_nei_avg_prod = fc_v*fc_nei_avg
  names(fc_nei_avg_sum) = paste0("fc_nei_avg_sum_factor_order_", order)
  names(fc_nei_avg_prod) = paste0("fc_nei_avg_prod_factor_order_", order)
  
  # Absolute average of neighbors fc
  fc_nei_avg_abs = mean(abs(fc_nei))
  fc_nei_avg_abs_sum = (abs(fc_v)+fc_nei_avg_abs)/2
  fc_nei_avg_abs_prod = abs(fc_v)*fc_nei_avg_abs
  names(fc_nei_avg_abs_sum) = paste0("fc_nei_avg_abs_sum_factor_order_", order)
  names(fc_nei_avg_abs_prod) = paste0("fc_nei_avg_abs_prod_factor_order_", order)
  
  # Averaging only neighboors of the same direction, removing the other ones    
  if(fc_v >= 0){
    fc_nei_avg_dir = fc_nei[fc_nei >= 0]
    if(nbr(fc_nei_avg_dir) == 0){fc_nei_avg_dir = 0}
    fc_nei_avg_dir_sum = (fc_v+mean(fc_nei_avg_dir))/2
    fc_nei_avg_dir_prod = fc_v*mean(fc_nei_avg_dir)
  } else {
    fc_nei_avg_dir = fc_nei[fc_nei < 0]
    if(nbr(fc_nei_avg_dir) == 0){fc_nei_avg_dir = 0}
    fc_nei_avg_dir_sum = (fc_v+mean(fc_nei_avg_dir))/2
    fc_nei_avg_dir_prod = fc_v*mean(fc_nei_avg_dir)
  }
  names(fc_nei_avg_dir_sum) = paste0("fc_nei_avg_dir_sum_factor_order_", order)
  names(fc_nei_avg_dir_prod) = paste0("fc_nei_avg_dir_prod_factor_order_", order)
  
  ######## Combined topological features of neighbors with fold changes.
  mean_metric_fc_nei_avg_sum = (mean_metric*fc_nei_avg_sum) %>% setNames(paste0(names(mean_metric), "_avg_fc_sum"))
  mean_metric_fc_nei_avg_prod = (mean_metric*fc_nei_avg_prod) %>% setNames(paste0(names(mean_metric), "_avg__fc_rod"))
  mean_metric_fc_nei_avg_abs_sum = (mean_metric*fc_nei_avg_abs_sum) %>% setNames(paste0(names(mean_metric), "_avg_abs_fc_sum"))
  mean_metric_fc_nei_avg_abs_prod = (mean_metric*fc_nei_avg_abs_prod) %>% setNames(paste0(names(mean_metric), "_avg_abs__fc_rod"))
  mean_metric_fc_nei_avg_dir_sum = (mean_metric*fc_nei_avg_dir_sum) %>% setNames(paste0(names(mean_metric), "_avg_dir_fc_sum"))
  mean_metric_fc_nei_avg_dir_prod = (mean_metric*fc_nei_avg_dir_prod) %>% setNames(paste0(names(mean_metric), "_avg_dir__fc_rod"))
  
  max_metric_fc_nei_avg_sum = (max_metric*fc_nei_avg_sum) %>% setNames(paste0(names(max_metric), "_avg_fc_sum"))
  max_metric_fc_nei_avg_prod = (max_metric*fc_nei_avg_prod) %>% setNames(paste0(names(max_metric), "_avg__fc_rod"))
  max_metric_fc_nei_avg_abs_sum = (max_metric*fc_nei_avg_abs_sum) %>% setNames(paste0(names(max_metric), "_avg_abs_fc_sum"))
  max_metric_fc_nei_avg_abs_prod = (max_metric*fc_nei_avg_abs_prod) %>% setNames(paste0(names(max_metric), "_avg_abs__fc_rod"))
  max_metric_fc_nei_avg_dir_sum = (max_metric*fc_nei_avg_dir_sum) %>% setNames(paste0(names(max_metric), "_avg_dir_fc_sum"))
  max_metric_fc_nei_avg_dir_prod = (max_metric*fc_nei_avg_dir_prod) %>% setNames(paste0(names(max_metric), "_avg_dir__fc_rod"))
  
  # Extract similarities with DEGS
  topo_results_simi = topo_results[c("Target", str_subset(colnames(topo_results), "Similarity_"))]
  topo_results_simi = topo_results_simi %>% dplyr::filter(Target  %in% node)
  topo_results_simi = topo_results_simi[str_remove(colnames(topo_results_simi), "Similarity_") %in% intersect(degs, neighbors)]
  similarity_with_DEGS = rowMeans(topo_results_simi)
  if(!is.finite(similarity_with_DEGS)){similarity_with_DEGS = 0}
  similarity_with_DEGS = similarity_with_DEGS %>% `names<-`(paste0("similarity_with_DEGS_order_", order))
  # Return results
  results = c(mean_metric, max_metric,
              discordant_neighborin_factor = dnf,
              fc_nei_avg_sum = fc_nei_avg_sum,
              fc_nei_avg_prod = fc_nei_avg_prod,
              fc_nei_avg_abs_sum = fc_nei_avg_abs_sum,
              fc_nei_avg_abs_prod = fc_nei_avg_abs_prod,
              fc_nei_avg_dir_sum = fc_nei_avg_dir_sum,
              fc_nei_avg_dir_prod = fc_nei_avg_dir_prod,
              mean_metric_fc_nei_avg_sum, mean_metric_fc_nei_avg_prod,
              mean_metric_fc_nei_avg_abs_sum, mean_metric_fc_nei_avg_abs_prod,
              mean_metric_fc_nei_avg_dir_sum, mean_metric_fc_nei_avg_dir_prod,
              max_metric_fc_nei_avg_sum, max_metric_fc_nei_avg_prod,
              max_metric_fc_nei_avg_abs_sum, max_metric_fc_nei_avg_abs_prod,
              max_metric_fc_nei_avg_dir_sum, max_metric_fc_nei_avg_dir_prod,
              similarity_with_DEGS)
  return(results)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param topo_results The results returned by the function extract_topolo
#' @param Signature A node signature (not list)
#' @param order The order of neighbors to consider
#' @param verbose Default: TRUE. Wether or not to print the internal calculations
#' @return Return something
 
extract_topolo_signature_graph = function(Graph, start_nodes = vnames(Graph), Signature, verbose = TRUE){
  if(verbose){message("extract_topolo_signature_graph is running...")}
  start_node_nei = lapply(neighborhood(graph = Graph, order = 1, nodes = start_nodes, mindist = 1), names)
  
  #### Simple features based on already connected signature components
  Disconnected_subgraph = subgraph(graph = Graph, vids = intersect(Signature, vnames(Graph)))
  component_Graph_DEG = components(Disconnected_subgraph)
  compo_list = split(names(component_Graph_DEG$membership), component_Graph_DEG$membership)
  compo_list = compo_list[nbrs(compo_list) != 1]
  
  # Get list of connected subgraph to the node
  feature_simple = list()
  for(node_i in 1:nbr(start_nodes)){
    # Get node neighbors
    compo_list_touch = compo_list[sapply(compo_list, function(compo) any(start_node_nei[[node_i]] %in% compo))]
    compo_list_touch = compo_list_touch[nbrs(compo_list_touch) != 0]
    
    # Get metric
    a1 = nbr(compo_list_touch)
    a2 = sum(nbrs(compo_list_touch))
    a3 = (a1/nbr(start_node_nei[[node_i]]))*100
    a4 = (a2/nbr(start_node_nei[[node_i]]))*100
    
    add(feature_simple, c(number_of_DEG_components_touched = a1, number_of_DEG_touched_in_component = a2, number_of_DEG_components_touched_perc = a3, number_of_DEG_touched_in_component_perc = a4))
  }
  feature_simple_df = do.call(rbind, feature_simple) %>% as.data.frame %>% dplyr::mutate(Target = start_nodes, .before = 1)
  
  #### Topological features based on signature subgraph
  Connected_subgraph = extract_subgraph(graph = Graph, terminals = intersect(Signature, vnames(Graph)))
  
  node_remaining = intersect(vnames(Connected_subgraph), start_nodes)
  topo_results = list(
    Degree_in = igraph::degree(Connected_subgraph, v = node_remaining, mode = "in"),
    Degree_out = igraph::degree(Connected_subgraph, v = node_remaining, mode = "out"),
    Degree_all = igraph::degree(Connected_subgraph, v = node_remaining, mode = "all"),
    Closeness_out = signif(igraph::closeness(Connected_subgraph, vids = node_remaining, mode = "out"), 3),
    Closeness_all = signif(igraph::closeness(Connected_subgraph, vids = node_remaining, mode = "all"), 3),
    Eccentricity_out = igraph::eccentricity(Connected_subgraph, vids = node_remaining, mode = "out"),
    Eccentricity_all = igraph::eccentricity(Connected_subgraph, vids = node_remaining, mode = "all"),
    Eigen_centrality = signif(igraph::eigen_centrality(Connected_subgraph)$vector[node_remaining], 3),
    Betweenness = signif(igraph::betweenness(Connected_subgraph, v = node_remaining), 3),
    Harmonic_centrality = signif(igraph::harmonic_centrality(graph = Connected_subgraph, vids = node_remaining, mode = "out"), 4))
  
  # Aggregate the results
  resi = do.call(cbind, topo_results) %>% as.data.frame()
  colnames(resi) = paste0(colnames(resi), "_subgraph")
  resi = rownames_to_column(resi, var = "Target")
  
  results = full_join(feature_simple_df, resi, by = "Target")
  results[is.na(results)] = 0
  return(results)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something
DEG_in_neighborhood = function(Graph, node, degs, order = 1){
  node_neighbors = names(ego(graph = Graph, order = order, mindist = 1, nodes = node)[[1]])
  deg_amount = nbrintersect(node_neighbors, degs)
  deg_perc = (deg_amount/nbr(node_neighbors))*100
  if(!is.finite(deg_perc)){deg_perc = 0}
  return(c(deg_amount = deg_amount, deg_perc = deg_perc))
}

#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something
interconnectivity = function(Graph, node, degs, order = 1){
  icn = function(Graph, vertex1, vertex2){
    friend = as.numeric(are_adjacent(graph = Graph, v1 = vertex1, v2 = vertex2))
    vertex1_neighbors = names(ego(graph = Graph, order = order, mindist = 1, nodes = vertex1)[[1]])
    vertex2_neighbors = names(ego(graph = Graph, order = order, mindist = 1, nodes = vertex2)[[1]])
    friend*((2+nbrintersect(vertex1_neighbors,vertex2_neighbors))/(sqrt(nbr(vertex1_neighbors)*nbr(vertex2_neighbors))))
  }
  degs = setdiff(degs, node)
  deg = degs[1]
  return(sum(sapply(degs, function(deg) icn(Graph = Graph, vertex1 = node, vertex2 = deg)))/nbr(degs))
}

#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something
#' @export
extract_subgraph <- function(graph, terminals){
  # Convert to undirected if necessary
  graph <- as.undirected(graph)
  the_subgraph = subgraph(graph = graph, terminals)
  if(nbr(components(the_subgraph)$csize) == 1){
    return(the_subgraph)
  }
  
  please_some_name = function(x){
    if(any(class(x) %in% c("matrix"))){
      return(rownames(x))
    } else {
      return(names(x))
    }
  }
  
  ### Calculate closest terminals for each terminal
  # calculate distance matrice
  graph_distance_sp = distances(graph = graph, v = terminals, to = terminals)
  diag(graph_distance_sp) = 100
  closest_term_to_each_term = apply(graph_distance_sp, 2, function(x) names(which(x == min(x, na.rm = TRUE))))
  
  # closest_term_to_each_term = lapply(1:ncol(graph_distance_sp), function(i) please_some_name(graph_distance_sp[graph_distance_sp[, i] == min_dist_per_terminals[i], ]))
  # names(closest_term_to_each_term) = terminals
  
  # Get all combinations of terminal pairs
  # terminals_pair_comb = do.call(rbind, lapply(1:nbr(terminals), function(i) data.frame(terminals[i], closest_term_to_each_term[[i]])))
  # terminals_pair_comb = combn(terminals, 2)
  
  # Get all nodes used for shortest paths across all combinations
  sp_list = list()
  for(i in 1:nbr(terminals)){
    # Get all shortest path between the two
    sp <- lapply(all_shortest_paths(graph, from = terminals[i], to = closest_term_to_each_term[[i]])$vpaths, function(x) setdiff(names(x), terminals))
    sp = sp[nbrs(sp) >= 1]
    add(sp_list, sp)
  }
  
  # Retrieve the number of times a node within a sp is found
  # Sort them by importance (frequency)
  all_unique_nodes_present_in_each_pair_sp = lapply(sp_list, function(x) unique(unlist(x)))
  node_importance_for_creating_sp = table(unlist(all_unique_nodes_present_in_each_pair_sp)) %>% sort(decreasing = TRUE)
  node_importance_for_creating_sp = split(names(node_importance_for_creating_sp), as.vector(node_importance_for_creating_sp))
  node_importance_for_creating_sp = node_importance_for_creating_sp[order(as.numeric(names(node_importance_for_creating_sp)), decreasing = TRUE)]
  
  # for each important node, create a subgraph, if the subgraph is not fully connected, continue.
  node_to_add = c()
  for(imp_node in node_importance_for_creating_sp){
    add(node_to_add, imp_node)
    the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
    if(nbr(components(the_subgraph)$csize) == 1){
      return(the_subgraph)
    }
  }
}

