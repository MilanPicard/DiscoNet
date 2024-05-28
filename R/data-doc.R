#' The A3LNetwork graph
#'
#' A medium size directed biological network comprised of three layers:
#' - proteins
#' - genes
#' - GO terms
#' The node names have prefix allowing their identication:
#' Prot_, Gene_, GO_
#' A `type` attribute is also present on nodes for identification based on user preference
#' type = protein, gene, or go
#'
#' @format ## `A3LNetwork`
#' An igraph object containing 7 538 nodes and 53 919 directed edges.
#' @source <StrinDB: https://string-db.org/cgi/download?sessionId=bjv6RHolIynV&species_text=Homo+sapiens>
"A3LNetwork"


#' Data for target repositioning for prostate cancer
#'
#' A list containing known targets for prostate cancer as well as a small gene signature
#' @format ## `data_prostate_cancer`
#' An list of size two, containing:
#'  - 13 known therapeutic targets for prostate cancer
#'  - named list containing gene signatures for prostate cancer
#' \describe{
#'   \item{Targets}{Known targets for prostate cancer}
#'   \item{Signatures}{A named list of signatures (Mutated_genes, DEG)}
#' }
#' @source
#' <Wang G, Zhao D, Spring DJ, DePinho RA. Genetics and biology of prostate cancer. Genes Dev. 2018;32(17-18):1105-1140. doi:10.1101/gad.315739.118>
#' <Myers JS, von Lersner AK, Robbins CJ, Sang QX. Differentially Expressed Genes and Signature Pathways of Human Prostate Cancer. PLoS One. 2015;10(12):e0145322. Published 2015 Dec 18. doi:10.1371/journal.pone.0145322>
#' <Yap, Timothy A et al. “Drug discovery in advanced prostate cancer: translating biology into therapy.” Nature reviews. Drug discovery vol. 15,10 (2016): 699-718. doi:10.1038/nrd.2016.120>
"data_prostate_cancer"
