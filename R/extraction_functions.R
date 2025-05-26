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

extract_more_si_features = function(Graph, start_nodes = vnames(Graph), Signature_list, FC_vectors, topo_results, nCores = 1, verbose = TRUE){
  if(verbose){message("extract_more_si_features is running...")}
  get_results1 = function(Graph, start_nodes, Signature_list, FC_vectors, i){
    if(verbose){cat('\r', paste0(i,'/',nbr(Signature_list), " signature"))}
    sign = Signature_list[[i]]
    sign_name = names(Signature_list)[i]

    # Calculate he interconnectivity between nodes and the signature
    interconnectivity_result1 = interconnectivity(Graph = Graph, start_nodes = start_nodes, degs = sign, order = 1)
    print("interconnectivity_result1")
    interconnectivity_result2 = interconnectivity(Graph = Graph, start_nodes = start_nodes, degs = sign, order = 2)
    print("interconnectivity_result2")

    # Are node in signature in the neighborhood
    DEG_in_neighborhood_result1   = DEG_in_neighborhood(Graph = Graph, start_nodes = start_nodes, degs = sign, order = 1)
    print("DEG_in_neighborhood_result1")
    DEG_in_neighborhood_result2   = DEG_in_neighborhood(Graph = Graph, start_nodes = start_nodes, degs = sign, order = 2)
    print("DEG_in_neighborhood_result2")

    # What are the topological properties of the node in subgraph of the signature
    subgraph_topo = extract_topolo_signature_graph(Graph = Graph, start_nodes = start_nodes, Signature = sign)
    print("subgraph_topo")

    # FOLD CHANGES BASED FEATURES
    get_results2 = function(Graph, start_nodes, FC_vectors, sign, j){
      Graph2 = Graph

      fc = FC_vectors[[j]]
      fc_name = names(FC_vectors)[j]

      fc = fc %>% stack %>% `colnames<-`(c("fc","node"))
      V(Graph2)$fc = left_join(data.frame(name = V(Graph2)$name), fc, by = c("name" = "node")) %>% pull(fc) %>% replace_na(0)

      # Simply retrieve the fold changes
      foldi = data.frame(Target = start_nodes, fold_change = V(Graph2)[start_nodes]$fc)

      # Composite features based on topological properties and fold change
      topo_features_with_foldchange = topo_metric_with_foldchange(Graph = Graph2, start_nodes = start_nodes, topo_results = topo_results)

      # What are the topological properties of the signature nodes in the neighborhood
      neighbordhood_scoring_result_1 = topo_metric_of_neighboring_signature(Graph = Graph2, start_nodes = start_nodes, degs = sign, topo_results = topo_results, order = 1)
      print("neighbordhood_scoring_result_1")
      neighbordhood_scoring_result_2 = topo_metric_of_neighboring_signature(Graph = Graph2, start_nodes = start_nodes, degs = sign, topo_results = topo_results, order = 2)
      print("neighbordhood_scoring_result_2")

      names(foldi)[-1] = paste0(fc_name, "_", names(foldi)[-1])
      names(neighbordhood_scoring_result_1)[-1] = paste0(fc_name, "_", names(neighbordhood_scoring_result_1)[-1])
      names(neighbordhood_scoring_result_2)[-1] = paste0(fc_name, "_", names(neighbordhood_scoring_result_2)[-1])
      names(topo_features_with_foldchange)[-1] = paste0(fc_name, "_", names(topo_features_with_foldchange)[-1])

      results_per_foldchanges = Reduce(function(x, y) full_join(x, y, by = "Target"), list(foldi, topo_features_with_foldchange, neighbordhood_scoring_result_1, neighbordhood_scoring_result_2))
      return(results_per_foldchanges)
    }

    neighbordhood_scoring_result_list = lapply(1:nbr(FC_vectors), function(j) get_results2(Graph = Graph, start_nodes = start_nodes, FC_vectors = FC_vectors, sign = sign, j=j))
    neighbordhood_scoring_result_list_merged = Reduce(function(x, y) full_join(x, y, by = "Target"), neighbordhood_scoring_result_list)

    # Merge results
    Signature_results = Reduce(function(x, y) full_join(x, y, by = "Target"), list(interconnectivity_result1, interconnectivity_result2, DEG_in_neighborhood_result1, DEG_in_neighborhood_result2, subgraph_topo))

    results = full_join(Signature_results, neighbordhood_scoring_result_list_merged, by = "Target")
    names(results)[-1] = paste0(sign_name, "_", names(results)[-1])
    return(results)
  }

  Signature_results_list = parallel::mclapply(1:nbr(Signature_list), function(i) get_results1(Graph = Graph, start_nodes = start_nodes, Signature_list = Signature_list, i=i, FC_vectors = FC_vectors), mc.cores = nCores)
  Signature_results_list_merged = Reduce(function(x, y) full_join(x, y, by = "Target"), Signature_results_list)
  if(verbose){message("Done.")}
  return(Signature_results_list_merged)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param topo_results The results returned by the function extract_topolo
#' @return Return something

topo_metric_with_foldchange = function(Graph, start_nodes, topo_results, fc_list, verbose = TRUE){
  if(verbose){message("topo_metric_with_foldchange  is running...")}
  # Get topological metrics of each nodes
  topo_results_metric = topo_results[str_detect(colnames(topo_results), "Similarity_", T)]
  topo_results_metric = topo_results_metric %>% filter(Target %in% start_nodes)
  topo_results_metric = topo_results_metric %>% select(-c(Topology_Degree_in, Topology_Degree_out, Topology_Closeness_out, Topology_Eccentricity_out))

  results_list = list()
  for(i in 1:nbr(fc_list)){
    fc = fc_list[[i]]
    fc_name = names(fc_list)[i]

    # Get fold changes for each start_nodes
    df_results = data.frame(Target = start_nodes)
    df_results$Fold_change = fc[start_nodes]
    df_results$Fold_change[is.na(df_results$Fold_change)] = 0

    # Multiply each topological metrics by its fold changes to create composite score.
    df_results = full_join(df_results, topo_results_metric, by = "Target") %>%
      mutate(across(-c(Target, Fold_change), ~ . *Fold_change))

    # Rename features by the fold change name
    names(df_results)[-1] = paste0(names(df_results)[-1], "_", fc_name)

    add(results_list, df_results)
  }
  Final_results = Reduce(function(x, y) full_join(x, y, by = "Target"), results_list)
  if(verbose){message("Done.")}
  return(Final_results)
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

topo_metric_of_neighboring_signature = function(Graph, start_nodes, degs, fc_list, topo_results, order = 1, verbose = TRUE){
  if(verbose){message("topo_metric_of_neighboring_signature is running...")}
  ######## Prep topological data
  topo_results_metric = topo_results[str_detect(colnames(topo_results), "Similarity_", T)]
  topo_results_metric = topo_results_metric %>% dplyr::select(-c(Topology_Degree_in, Topology_Degree_out, Topology_Closeness_out, Topology_Eccentricity_out))
  topo_results_simi = topo_results[c("Target", str_subset(colnames(topo_results), "Similarity_"))]

  ######## Prep node neighbors data
  start_node_neighbors = lapply(ego(graph = Graph, order = order, mindist = 1, nodes = start_nodes), names)
  names(start_node_neighbors) = start_nodes
  start_node_neighbors = lapply(start_node_neighbors, function(X) intersect(X, degs))
  start_node_neighbors = start_node_neighbors[nbrs(start_node_neighbors) != 0]

  results_list = list()
  for(i in 1:nbr(fc_list)){
    fc = fc_list[[i]]
    fc_name = names(fc_list)[i]

    # Get fold changes for each start_nodes
    df_results = data.frame(Target = start_nodes)
    df_results$Fold_change = fc[start_nodes]
    df_results$Fold_change[is.na(df_results$Fold_change)] = 0

    ######## Retrieve fold change properties of DEGs in node neighborhood
    if(verbose){message(paste0("Calculating fold change properties of neighboring DEGS of ", i, "/", nbr(fc_list), " fold changes."))}
    neideg_fc_list = list()
    for(j in 1:nbr(start_node_neighbors)){
      fc_v = df_results$Fold_change[j]
      fc_nei = df_results %>% dplyr::filter(Target %in% start_node_neighbors[[j]]) %>% dplyr::pull(Fold_change)

      # Extract neighbors fold changes
      fc_nei_avg = mean(fc_nei)
      fc_nei_avg_sum = (fc_v+fc_nei_avg)/2
      fc_nei_avg_prod = fc_v*fc_nei_avg
      names(fc_nei_avg_sum) = paste("fc_nei_avg_sum_order", fc_name, order, sep = "_")
      names(fc_nei_avg_prod) = paste("fc_nei_avg_prod_order", fc_name, order, sep = "_")

      # Absolute average of neighbors fc
      fc_nei_avg_abs = mean(abs(fc_nei))
      fc_nei_avg_abs_sum = (abs(fc_v)+fc_nei_avg_abs)/2
      fc_nei_avg_abs_prod = abs(fc_v)*fc_nei_avg_abs
      names(fc_nei_avg_abs_sum) = paste("fc_nei_avg_abs_sum_order", fc_name, order, sep = "_")
      names(fc_nei_avg_abs_prod) = paste("fc_nei_avg_abs_prod_order", fc_name, order, sep = "_")

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
      names(fc_nei_avg_dir_sum) = paste("fc_nei_avg_dir_sum_order", fc_name, order, sep = "_")
      names(fc_nei_avg_dir_prod) = paste("fc_nei_avg_dir_prod_order", fc_name, order, sep = "_")

      add(neideg_fc_list, c(fc_nei_avg_sum,fc_nei_avg_prod,fc_nei_avg_abs_sum,fc_nei_avg_abs_prod,fc_nei_avg_dir_sum,fc_nei_avg_dir_prod))
    }

    # Assemble the results of previous loop
    neideg_fc_results = do.call(rbind, neideg_fc_list) %>% as.data.frame %>% mutate(Target = names(start_node_neighbors), .before = 1)

    ######## Calculating topological properties with and without fold change of neighboring DEGS
    if(verbose){message(paste0("Calculating topological properties with and without fold change of neighboring DEGS of ", i, "/", nbr(fc_list), " fold changes."))}
    top_neideg_fc_list = list()
    for(j in 1:nbr(start_node_neighbors)){
      # Get DEG neighbors topological properties
      topo_results_metric_node = topo_results_metric %>% filter(Target  %in% start_node_neighbors[[j]])
      mean_metric = colMeans(topo_results_metric_node[, -1])
      names(mean_metric) = paste0("Mean_",str_remove(names(mean_metric), "Topology_"),"_neighboring_signature_order_",order)

      mean_metric_fc = lapply(unlist(neideg_fc_results[j, -1]), function(x) x*mean_metric)
      mean_metric_fc = lapply(seq_along(mean_metric_fc), function(k) setNames(mean_metric_fc[[k]], paste0(names(mean_metric_fc[[k]]), names(mean_metric_fc)[k])))
      mean_metric_fc = mean_metric_fc %>% unlist

      # Extract similarities with DEGS
      topo_results_simi_node = topo_results_simi %>% dplyr::filter(Target  %in% names(start_node_neighbors)[j])
      topo_results_simi_node = topo_results_simi_node[paste0("Similarity_", start_node_neighbors[[j]])]
      similarity_with_DEGS = rowMeans(topo_results_simi_node)
      similarity_with_DEGS = similarity_with_DEGS %>% `names<-`(paste0("similarity_with_DEGS_order_", order))

      add(top_neideg_fc_list, c(mean_metric, mean_metric_fc, similarity_with_DEGS))
    }

    # Assemble the results of previous loop
    top_neideg_fc_results = do.call(rbind, top_neideg_fc_list) %>% as.data.frame %>% mutate(Target = names(start_node_neighbors), .before = 1)

    # Assemble all results
    results = full_join(neideg_fc_results, top_neideg_fc_results, by = "Target")
    add(results_list, results)
  }

  Final_results = Reduce(function(x, y) full_join(x, y, by = "Target"), results_list)
  if(verbose){message("Done.")}
  return(Final_results)
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

extract_subgraph_features = function(Graph, start_nodes = vnames(Graph), Signature_list, verbose = TRUE){
  if(verbose){message("extract_subgraph_features is running...")}
  Signature_list = lapply(Signature_list, function(x) intersect(x, vnames(Graph)))

  Signature_results_list = list()
  for(i in 1:nbr(Signature_list)){
    sign = Signature_list[[i]]
    sign_name = names(Signature_list)[i]

    subgraph1 = extract_subgraph(graph = Graph, terminals = sign, type = 1)
    subgraph2 = extract_subgraph(graph = Graph, terminals = sign, type = 2)
    subgraph3 = extract_subgraph(graph = Graph, terminals = sign, type = 3)

    nodes_in_sub1 = intersect(vnames(subgraph1), start_nodes)
    nodes_in_sub2 = intersect(vnames(subgraph2), start_nodes)
    nodes_in_sub3 = intersect(vnames(subgraph3), start_nodes)

    topo_metric1 = list(Degree = igraph::degree(subgraph1, v = nodes_in_sub1, mode = "all"),
                        Closeness = signif(igraph::closeness(subgraph1, vids = nodes_in_sub1, mode = "all"), 3),
                        Eccentricity = igraph::eccentricity(subgraph1, vids = nodes_in_sub1, mode = "all"),
                        Eigen_centrality = signif(igraph::eigen_centrality(subgraph1)$vector[nodes_in_sub1], 3),
                        Betweenness = signif(igraph::betweenness(subgraph1, v = nodes_in_sub1), 3),
                        Harmonic_centrality = signif(igraph::harmonic_centrality(graph = subgraph1, vids = nodes_in_sub1, mode = "all"), 4))

    topo_metric2 = list(Degree = igraph::degree(subgraph2, v = nodes_in_sub2, mode = "all"),
                        Closeness = signif(igraph::closeness(subgraph2, vids = nodes_in_sub2, mode = "all"), 3),
                        Eccentricity = igraph::eccentricity(subgraph2, vids = nodes_in_sub2, mode = "all"),
                        Eigen_centrality = signif(igraph::eigen_centrality(subgraph2)$vector[nodes_in_sub2], 3),
                        Betweenness = signif(igraph::betweenness(subgraph2, v = nodes_in_sub2), 3),
                        Harmonic_centrality = signif(igraph::harmonic_centrality(graph = subgraph2, vids = nodes_in_sub2, mode = "all"), 4))

    topo_metric3 = list(Degree = igraph::degree(subgraph3, v = nodes_in_sub3, mode = "all"),
                        Closeness = signif(igraph::closeness(subgraph3, vids = nodes_in_sub3, mode = "all"), 3),
                        Eccentricity = igraph::eccentricity(subgraph3, vids = nodes_in_sub3, mode = "all"),
                        Eigen_centrality = signif(igraph::eigen_centrality(subgraph3)$vector[nodes_in_sub3], 3),
                        Betweenness = signif(igraph::betweenness(subgraph3, v = nodes_in_sub3), 3),
                        Harmonic_centrality = signif(igraph::harmonic_centrality(graph = subgraph3, vids = nodes_in_sub3, mode = "all"), 4))

    df_results = data.frame(Target = start_nodes)
    topo_metric_df1 = do.call(cbind, topo_metric1) %>% as.data.frame %>% `colnames<-`(paste0("Subgraph_1_", colnames(.))) %>% rownames_to_column("Target")
    topo_metric_df2 = do.call(cbind, topo_metric2) %>% as.data.frame %>% `colnames<-`(paste0("Subgraph_2_", colnames(.))) %>% rownames_to_column("Target")
    topo_metric_df3 = do.call(cbind, topo_metric3) %>% as.data.frame %>% `colnames<-`(paste0("Subgraph_3_", colnames(.))) %>% rownames_to_column("Target")

    # Merge results
    df_results = Reduce(function(x, y) full_join(x, y, by = "Target"), list(df_results, topo_metric_df1, topo_metric_df2, topo_metric_df3))
    df_results[is.na(df_results)] = 0

    names(df_results)[-1] = paste0(sign_name, "_", names(df_results)[-1])
    add(Signature_results_list, df_results)
  }
  Signature_results_list_merged = Reduce(function(x, y) full_join(x, y, by = "Target"), Signature_results_list)

  if(verbose){message("Done.")}
  return(Signature_results_list_merged)
}



#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something
DEG_in_neighborhood = function(Graph, start_nodes, degs_list, order = 1, verbose = TRUE){
  if(verbose){message("DEG_in_neighborhood is running...")}

  results_per_deg_list = list()
  for(i in 1:nbr(degs_list)){
    degs = degs_list[[i]]
    deg_name = names(degs_list)[i]

    # Get start_nodes neighbors
    node_neighbors = lapply(ego(graph = Graph, order = order, mindist = 1, nodes = start_nodes), names)

    # Get number of degs in proximity
    Degree_to_signature = sapply(node_neighbors, function(x) nbrintersect(x, degs))
    # Get percentage of degs in proximity
    Degree_perc_to_signature = (Degree_to_signature/nbrs(node_neighbors))*100
    Degree_perc_to_signature[is.finite(Degree_perc_to_signature) != TRUE] = 0

    results = data.frame(Target = start_nodes, Degree_to_signature, Degree_perc_to_signature)
    names(results)[-1] = paste(names(results)[-1], deg_name, "order", order, sep = "_")

    add(results_per_deg_list, results)
  }

  Final_results = Reduce(function(x, y) full_join(x, y, by = "Target"), results_per_deg_list)

  if(verbose){message("done")}
  return(Final_results)
}


#' Extract signature-based features
#'
#' Extract signature-based features from a biological network based on a list of signatures nodes.
#' @param Graph A igraph object.
#' @param start_nodes Features are extracted for each of these nodes
#' @param degs A node signature (not list)
#' @param order The order of neighbors to consider
#' @return Return something
interconnectivity = function(Graph, start_nodes, degs_list, order = 1, verbose = TRUE){
  if(verbose){message("interconnectivity is running...")}

  icn = function(nei1, nei2){
    (2+nbrintersect(nei1,nei2))/(sqrt(nbr(nei1)*nbr(nei2)))
  }

  results_per_deg_list = list()
  for(i in 1:nbr(degs_list)){
    degs = degs_list[[i]]
    deg_name = names(degs_list)[i]

    # Calculates for each start_nodes, the average interconnectivity to degs nodes.
    # Get neighbors of all start_nodes and degs
    start_nodes_neighbors = lapply(ego(graph = Graph, order = order, mindist = 1, nodes = start_nodes), names)
    names(start_nodes_neighbors) = start_nodes
    degs_neighbors = lapply(ego(graph = Graph, order = order, mindist = 1, nodes = degs), names)
    names(degs_neighbors) = degs

    # Find start_nodes that are neighbors to degs
    dist_mat = distances(graph = Graph, v = start_nodes, to = degs)
    node_that_are_neighbors = which(dist_mat <= order, arr.ind = TRUE)

    # For each pair, calculate icn.
    results = mapply(function(nei1, nei2) icn(nei1, nei2), start_nodes_neighbors[start_nodes[node_that_are_neighbors[, 1]]], degs_neighbors[degs[node_that_are_neighbors[, 2]]])
    results = stack(results) %>% group_by(ind) %>% summarise(Interconnectivity = mean(values))

    # Produce results
    results = full_join(data.frame(Target = start_nodes), results, by = c("Target" = "ind")) %>%
      mutate(Interconnectivity = replace_na(Interconnectivity, 0))

    names(results)[-1] = paste(names(results)[-1], deg_name, "order", order, sep = "_")

    add(results_per_deg_list, results)
  }

  Final_results = Reduce(function(x, y) full_join(x, y, by = "Target"), results_per_deg_list)

  if(verbose){message("done")}
  return(Final_results)
}




extract_subgraph <- function(graph, terminals, type){
  if(type == 1){
    return(extract_subgraph_type1(graph = graph, terminals = terminals))
  }
  if(type == 2){
    return(extract_subgraph_type2(graph = graph, terminals = terminals))
  }
  if(type == 3){
    return(extract_subgraph_type3(graph = graph, terminals = terminals))
  }
}

extract_subgraph_type1 = function(graph, terminals){
  # Convert to undirected if necessary
  graph <- as_undirected(graph)
  the_subgraph = subgraph(graph = graph, terminals)
  if(nbr(components(the_subgraph)$csize) == 1){
    return(the_subgraph)
  }

  which.min.names = function(x){
    if(any(class(x) %in% c("matrix"))){
      whiches = which(x == min(x), arr.ind = TRUE)
      return(list(rownames(x)[unique(whiches[, 1])], colnames(x)[unique(whiches[, 2])])) # does not work properly
      # Does a cross between multiple match of rows and colunms
      # It will return a set of shortest distances
      # And a set of other distances
      # It need further filtering to keep only the minimum, which is why
      # in identify_imp_node_from_sp I re-select shortest paths among shortest paths.
    } else {
      return(list(which(x == min(x)), 1))
    }
  }

  identify_imp_node_from_sp = function(sps){
    # sps must have start and end nodes
    split_sps = split(sps, sapply(sps, tail, 1)) # Gives information about which ters will be connected
    split_sps_clean = lapply(split_sps, function(end_node) lapply(end_node, function(sp) setdiff(sp, terminals)))
    min_sp = lapply(split_sps_clean, function(x) nbrs(x)) %>% unlist %>% min

    split_sps_clean = lapply(split_sps_clean, function(x) x[nbrs(x) == min_sp])

    # In order to consider the path length, the selection of nodes is done walk by walk
    walks = nbrs(split_sps_clean[[1]]) %>% unique # They should all have the same length, cause they are sp
    if(nbr(walks) != 1){stop(sps)}
    node_imp_list = c()
    for(walk in 1:walks){
      if(walk > 1){
        split_sps_clean2 = lapply(split_sps_clean, function(X) X[sapply(X, function(x) x[1:(walk-1)] %in% node_imp_list)])
      } else {
        split_sps_clean2 = split_sps_clean
      }
      sp_node_presence = lapply(split_sps_clean2, function(X) unique(unlist(lapply(X, function(x) x[walk]))))
      node_imp = table(unlist(sp_node_presence))
      node_imp = split(names(node_imp), as.vector(node_imp))
      node_imp = node_imp[order(as.numeric(names(node_imp)), decreasing = TRUE)]

      # Taking the first one does not guaranti every path is confirmed, just the maximum number of paths if only one is taken
      add(node_imp_list, unname(node_imp[[1]]))
    }
    node_imp_list = unique(node_imp_list)
    return(node_imp_list)
  }

  ### Calculate closest terminals for each terminal
  # calculate distance matrice
  graph_distance_sp = distances(graph = graph, v = terminals, to = terminals)
  diag(graph_distance_sp) = 100
  closest_term_to_each_term = apply(graph_distance_sp, 2, function(x) names(which(x == min(x, na.rm = TRUE))))

  # Get all nodes used for shortest paths across all combinations
  sp_list = list()
  for(i in 1:nbr(terminals)){
    # Get all shortest path between the two
    sp <- lapply(all_shortest_paths(graph, from = terminals[i], to = closest_term_to_each_term[[i]])$vpaths, function(x) names(x))
    sp = sp[nbrs(sp) >= 3]
    names(sp) = rep(terminals[i], nbr(sp))
    add(sp_list, sp)
  }
  names(sp_list) = terminals

  node_to_add = unique(unlist(sp_list))
  the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
  if(nbr(components(the_subgraph)$csize) == 1){
    return(the_subgraph)
  }
  # However, some groups of ters will still be deconnected from other groups of ters.
  # To connect all ters, first lets find the components of ters with the most centrality, and link every components to this one.
  # Find the components
  subgraph_components = components(the_subgraph)$membership %>% split(names(.), .)

  node_to_keep_for_compo = c()
  for(compo_i in 1:nbr(subgraph_components)){
    compo = subgraph_components[[compo_i]]
    # Find the minimum distance between two nodes from central component and another component
    dist_between_compo = graph_distance_sp[intersect(unlist(subgraph_components[-compo_i]), terminals), intersect(compo, terminals)]
    min_dist = which.min.names(dist_between_compo)

    # Calculate shortest paths between central components and another component
    sp_comp <- lapply(all_shortest_paths(graph, from = min_dist[[1]], to = min_dist[[2]])$vpaths, function(x) names(x))
    add(node_to_keep_for_compo, identify_imp_node_from_sp(sp_comp))
  }
  packageVersion("igraph")

  node_to_add = c(node_to_add, unique(unlist(node_to_keep_for_compo)))
  the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
  if(nbr(components(the_subgraph)$csize) == 1){
    return(the_subgraph)
  }
  stop("Didnt converge to a network, ask you're mom for help")
}

extract_subgraph_type2 <- function(graph, terminals){
  # Convert to undirected if necessary
  graph <- as_undirected(graph)
  the_subgraph = subgraph(graph = graph, terminals)
  if(nbr(components(the_subgraph)$csize) == 1){
    return(the_subgraph)
  }

  please_some_name = function(x){
    if(any(class(x) %in% c("matrix"))){c
      return(rownames(x))
    } else {
      return(names(x))
    }
  }

  which.min.names = function(x){
    if(any(class(x) %in% c("matrix"))){
      whiches = which(x == min(x), arr.ind = TRUE)
      return(list(rownames(x)[unique(whiches[, 1])], colnames(x)[unique(whiches[, 2])])) # does not work properly
      # Does a cross between multiple match of rows and colunms
      # It will return a set of shortest distances
      # And a set of other distances
      # It need further filtering to keep only the minimum, which is why
      # in identify_imp_node_from_sp I re-select shortest paths among shortest paths.
    } else {
      return(list(which(x == min(x)), 1))
    }
  }

  identify_imp_node_from_sp = function(sps){
    # sps must have start and end nodes
    split_sps = split(sps, sapply(sps, tail, 1)) # Gives information about which ters will be connected
    split_sps_clean = lapply(split_sps, function(end_node) lapply(end_node, function(sp) setdiff(sp, terminals)))
    min_sp = lapply(split_sps_clean, function(x) nbrs(x)) %>% unlist %>% min

    split_sps_clean = lapply(split_sps_clean, function(x) x[nbrs(x) == min_sp])

    # In order to consider the path length, the selection of nodes is done walk by walk
    walks = nbrs(split_sps_clean[[1]]) %>% unique # They should all have the same length, cause they are sp
    if(nbr(walks) != 1){stop(sps)}
    node_imp_list = c()
    for(walk in 1:walks){
      if(walk > 1){
        split_sps_clean2 = lapply(split_sps_clean, function(X) X[sapply(X, function(x) x[1:(walk-1)] %in% node_imp_list)])
      } else {
        split_sps_clean2 = split_sps_clean
      }
      sp_node_presence = lapply(split_sps_clean2, function(X) unique(unlist(lapply(X, function(x) x[walk]))))
      node_imp = table(unlist(sp_node_presence))
      node_imp = split(names(node_imp), as.vector(node_imp))
      node_imp = node_imp[order(as.numeric(names(node_imp)), decreasing = TRUE)]

      # Taking the first one does not guaranti every path is confirmed, just the maximum number of paths if only one is taken
      add(node_imp_list, unname(node_imp[[1]]))
    }
    node_imp_list = unique(node_imp_list)
    return(node_imp_list)
  }

  # Count nodes in sps and return
  count_in_sp = function(sp_list){
    all_unique_nodes_present_in_each_pair_sp = lapply(sp_list, function(x) unique(unlist(x)))
    node_importance_for_creating_sp = table(unlist(all_unique_nodes_present_in_each_pair_sp)) %>% sort(decreasing = TRUE)
    node_importance_for_creating_sp = split(names(node_importance_for_creating_sp), as.vector(node_importance_for_creating_sp))
    node_importance_for_creating_sp = node_importance_for_creating_sp[order(as.numeric(names(node_importance_for_creating_sp)), decreasing = TRUE)]
    return(node_importance_for_creating_sp)
  }

  # calculate distance matrice
  graph_distance_sp = distances(graph = graph, v = terminals, to = terminals)
  diag(graph_distance_sp) = 100
  closest_term_to_each_term = apply(graph_distance_sp, 2, function(x) names(which(x == min(x, na.rm = TRUE))))

  # Get all nodes used for shortest paths across all combinations
  sp_list = list()
  for(i in 1:nbr(terminals)){
    # Get all shortest path between the two
    sp <- lapply(all_shortest_paths(graph, from = terminals[i], to = closest_term_to_each_term[[i]])$vpaths, function(x) names(x))
    sp = sp[nbrs(sp) >= 3]
    names(sp) = rep(terminals[i], nbr(sp))
    add(sp_list, sp)
  }
  names(sp_list) = terminals

  # Keep only sp to sort
  sp_list2 = sp_list[nbrs(sp_list) != 0]

  sp_list_info = lapply(sp_list2, function(x) identify_imp_node_from_sp(x))

  node_importance_for_creating_sp = count_in_sp(sp_list_info)

  # After this, all ters will be connected to its closest ters
  node_to_add = c()
  for(imp_node in node_importance_for_creating_sp){
    add(node_to_add, imp_node)
    the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
    if(nbr(components(the_subgraph)$csize) == 1){
      return(the_subgraph)
    }
  }

  # However, some groups of ters will still be deconnected from other groups of ters.
  # To connect all ters, first lets find the components of ters with the most centrality, and link every components to this one.
  # Find the components
  subgraph_components = components(the_subgraph)$membership %>% split(names(.), .)

  node_to_keep_for_compo = c()
  for(compo_i in 1:nbr(subgraph_components)){
    compo = subgraph_components[[compo_i]]
    # Find the minimum distance between two nodes from central component and another component
    dist_between_compo = graph_distance_sp[intersect(unlist(subgraph_components[-compo_i]), terminals), intersect(compo, terminals)]
    min_dist = which.min.names(dist_between_compo)

    # Calculate shortest paths between central components and another component
    sp_comp <- lapply(all_shortest_paths(graph, from = min_dist[[1]], to = min_dist[[2]])$vpaths, function(x) names(x))
    add(node_to_keep_for_compo, identify_imp_node_from_sp(sp_comp))
  }
  # Retrieve the number of times a node within a sp is found
  # Sort them by importance (frequency)
  node_importance_for_creating_sp_compo = count_in_sp(node_to_keep_for_compo)

  # After this, all ters will be connected to its closest ters
  for(imp_node in node_importance_for_creating_sp_compo){
    add(node_to_add, imp_node)
    the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
    if(nbr(components(the_subgraph)$csize) == 1){
      return(the_subgraph)
    }
  }
  stop("Didnt converge to a network, ask you're mom for help")
}

extract_subgraph_type3 = function(graph, terminals){
  # Convert to undirected if necessary
  graph <- as_undirected(graph)
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

  which.min.names = function(x){
    if(any(class(x) %in% c("matrix"))){
      whiches = which(x == min(x), arr.ind = TRUE)
      return(list(rownames(x)[unique(whiches[, 1])], colnames(x)[unique(whiches[, 2])])) # does not work properly
      # Does a cross between multiple match of rows and colunms
      # It will return a set of shortest distances
      # And a set of other distances
      # It need further filtering to keep only the minimum, which is why
      # in identify_imp_node_from_sp I re-select shortest paths among shortest paths.
    } else {
      return(list(which(x == min(x)), 1))
    }
  }

  # Function to select and rank important nodes from a set of shortest paths
  identify_imp_node_from_sp = function(sps){
    # sps must have start and end nodes
    split_sps = split(sps, sapply(sps, tail, 1)) # Gives information about which ters will be connected
    split_sps_clean = lapply(split_sps, function(end_node) lapply(end_node, function(sp) setdiff(sp, terminals)))
    min_sp = lapply(split_sps_clean, function(x) nbrs(x)) %>% unlist %>% min

    split_sps_clean = lapply(split_sps_clean, function(x) x[nbrs(x) == min_sp])

    # In order to consider the path length, the selection of nodes is done walk by walk
    walks = nbrs(split_sps_clean[[1]]) %>% unique # They should all have the same length, cause they are sp
    if(nbr(walks) != 1){stop(sps)}
    node_imp_list = c()
    for(walk in 1:walks){
      if(walk > 1){
        split_sps_clean2 = lapply(split_sps_clean, function(X) X[sapply(X, function(x) x[1:(walk-1)] %in% node_imp_list)])
      } else {
        split_sps_clean2 = split_sps_clean
      }
      sp_node_presence = lapply(split_sps_clean2, function(X) unique(unlist(lapply(X, function(x) x[walk]))))
      forced_node = sp_node_presence[nbrs(sp_node_presence) == 1] %>% unlist %>% unname

      sp_node_presence = sp_node_presence[!sapply(sp_node_presence, function(X) any(X %in% forced_node))]
      if(nbr(sp_node_presence) == 0){
        add(node_imp_list, c(forced_node) %>% unname)
        next
      }
      node_imp = table(unlist(sp_node_presence))
      node_imp = node_imp[setdiff(names(node_imp), forced_node)]
      node_imp = split(names(node_imp), as.vector(node_imp))
      node_imp = node_imp[order(as.numeric(names(node_imp)), decreasing = TRUE)]

      # Taking the first one does not guaranti every path is confirmed, just the maximum number of paths if only one is taken
      missing_links = sp_node_presence[!sapply(sp_node_presence, function(X) any(X %in% node_imp[[1]]))]
      add(node_imp_list, c(forced_node, node_imp[[1]], unlist(missing_links)) %>% unname)

    }
    node_imp_list = unique(node_imp_list)
    return(node_imp_list)
  }

  # Count nodes in sps and return
  count_in_sp = function(sp_list){
    all_unique_nodes_present_in_each_pair_sp = lapply(sp_list, function(x) unique(unlist(x)))
    node_importance_for_creating_sp = table(unlist(all_unique_nodes_present_in_each_pair_sp)) %>% sort(decreasing = TRUE)
    node_importance_for_creating_sp = split(names(node_importance_for_creating_sp), as.vector(node_importance_for_creating_sp))
    node_importance_for_creating_sp = node_importance_for_creating_sp[order(as.numeric(names(node_importance_for_creating_sp)), decreasing = TRUE)]
    return(node_importance_for_creating_sp)
  }

  ### Calculate closest terminals for each terminal
  # calculate distance matrice
  graph_distance_sp = distances(graph = graph, v = terminals, to = terminals)
  diag(graph_distance_sp) = 100
  closest_term_to_each_term = apply(graph_distance_sp, 2, function(x) names(which(x == min(x, na.rm = TRUE))))

  # Get all nodes used for shortest paths across all combinations
  sp_list = list()
  for(i in 1:nbr(terminals)){
    # Get all shortest path between the two
    sp <- lapply(all_shortest_paths(graph, from = terminals[i], to = closest_term_to_each_term[[i]])$vpaths, function(x) names(x))
    sp = sp[nbrs(sp) >= 3]
    names(sp) = rep(terminals[i], nbr(sp))
    add(sp_list, sp)
  }
  names(sp_list) = terminals

  # Keep only sp to sort
  sp_list2 = sp_list[nbrs(sp_list) != 0]

  # Find important nodes based on the number of times they are needed for sps
  sp_list_info = lapply(sp_list2, function(x) identify_imp_node_from_sp(x))

  # Retrieve the number of times an important is identified
  # Sort them by importance (frequency)
  node_importance_for_creating_sp = count_in_sp(sp_list_info)

  # After this, all ters will be connected to its closest ters
  node_to_add = c()
  for(imp_node in node_importance_for_creating_sp){
    add(node_to_add, imp_node)
    the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
    if(nbr(components(the_subgraph)$csize) == 1){
      return(the_subgraph)
    }
  }
  # However, some groups of ters will still be deconnected from other groups of ters.
  # To connect all ters, first lets find the components of ters with the most centrality, and link every components to this one.
  # Find the components
  subgraph_components = components(the_subgraph)$membership %>% split(names(.), .)

  node_to_keep_for_compo = c()
  for(compo_i in 1:nbr(subgraph_components)){
    compo = subgraph_components[[compo_i]]
    # Find the minimum distance between two nodes from central component and another component
    dist_between_compo = graph_distance_sp[intersect(unlist(subgraph_components[-compo_i]), terminals), intersect(compo, terminals)]
    min_dist = which.min.names(dist_between_compo)

    # Calculate shortest paths between central components and another component
    sp_comp <- lapply(all_shortest_paths(graph, from = min_dist[[1]], to = min_dist[[2]])$vpaths, function(x) names(x))
    add(node_to_keep_for_compo, identify_imp_node_from_sp(sp_comp))
  }
  # Retrieve the number of times a node within a sp is found
  # Sort them by importance (frequency)
  node_importance_for_creating_sp_compo = count_in_sp(node_to_keep_for_compo)

  # After this, all ters will be connected to its closest ters
  for(imp_node in node_importance_for_creating_sp_compo){
    add(node_to_add, imp_node)
    the_subgraph = subgraph(graph = graph, c(terminals, node_to_add))
    if(nbr(components(the_subgraph)$csize) == 1){
      return(the_subgraph)
    }
  }
  stop("Didnt converge to a network, ask you're mom for help")
}




