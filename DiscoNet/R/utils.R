#' Create chunks
#'
#' Divide a vector into chunks of the same size
#' @param x The vector
#' @param n The number of chunks
#' @return A list of chunks from the original vector
#' @examples
#' list_of_letters = chunks(letters, 5)
#' @export
chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


#' Length
#'
#' Identical to length
#' @param x A R object
#' @return The length of the object x
#' @examples
#' Length_of_letters = nbr(letters)
#' @export
nbr = function(x){return(length(x))}


#' Lengths
#'
#' Identical to lengths
#' @param x A R object
#' @return The lengths of the object x
#' @examples
#' ll = nbrs(list(1:3, 1:10))
#' @export
nbrs = function(x){return(lengths(x))}


#' Length of intersect
#'
#' Identical to length(intersect(x,y))
#' @param v1 A vector
#' @param v2 Another vector
#' @return The length of the intersect between v1 and v2
#' @examples
#' ll = nbrintersect(1:10, 1:5)
#' @export
nbrintersect = function(v1,v2){return(length(intersect(v1,v2)))}


#' Length of setdiff
#'
#' Identical to length(setdiff(x,y))
#' @param v1 A vector
#' @param v2 Another vector
#' @return The length of the setdiff between v1 and v2
#' @examples
#' ll = nbrsetdiff(1:10, 1:5)
#' @export
nbrsetdiff = function(v1,v2){return(length(setdiff(v1,v2)))}


#' Length of unique
#'
#' Identical to length(unique(v1))
#' @param v1 A vector
#' @return The length of the unique between v1 and v2
#' @examples
#' ll = nbrunique(c(1,2,1))
#' @export
nbrunique = function(v1){return(length(unique(v1)))}


#' Length of unlist
#'
#' Identical to length(unlist(v1))
#' @param v1 A list
#' @return The length of the unlist between v1 and v2
#' @examples
#' ll = nbrunlist(list(1:3, 2:4))
#' @export
nbrunlist = function(v1){return(length(unlist(v1)))}


#' Concatenate with replacement
#'
#' Add an element to a vector or a list with replacement
#' @param v1 A vector or a list
#' @param x The element to be added
#' @return The length of the unlist between v1 and v2
#' @examples
#' a = 1:10
#' add(a, 11)
#' @export
add = function(v1, x){
  if(is.atomic(v1)){eval.parent(substitute(v1 <- c(v1,x)))}
  if(is.list(v1))  {eval.parent(substitute(v1 <- append(v1,list(x))))}
}


#' Retrieve node names
#'
#' Retrieve the names of nodes from an igraph object
#' @param igraph_obj A igraph object
#' @param pattern Pattern to look for.
#' @param negate If TRUE, return non-matching elements.
#' @return The length of the unlist between v1 and v2
#' @examples
#' # Don't run
#' # all_nodes = vnames(igraph_obj)
#' # protein_nodes = vnames(igraph_obj, type = "protein_")
#' @export
vnames = function(igraph_obj, pattern = NULL, negate = FALSE){
  node_names = igraph::V(igraph_obj)$name
  if(is.null(pattern)){
    return(node_names)
  }
  return(stringr::str_subset(string = node_names, pattern = pattern, negate = negate))
}



