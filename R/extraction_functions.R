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

