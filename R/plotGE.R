#' plotGE
#'
#' Plot gene expression for spatially resolved transcriptomics data 
#' 
#' Function to generate gene expression plots for spatially resolved transcriptomics datasets
#' 
#' This function generates a gene expression plot for spot-based spatially resolved
#' transcriptomics data with several
#' options available to adjust the plot type and style.
#'
#' @param spe (SpatialExperiment) Input data object.
#' @param gene_id (character) The id of one gene. Refer to gene_id Column in \code{rowData}  
#' @param gene_name (character) The name of one gene. Refer to gene_name Column in \code{rowData} 
#' @param log_count (logical) Whether to plot the log transfered counts if exists. Default
#'   = FALSE.
#' @param Visium  (logical) Whether to show histology image as background as in \link{ggspavis}{plotVisium}. Default
#'   = FALSE.
#' @param ... Other arguments for \link{ggspavis}{plotVisium} and \link{ggspavis}{plotSpots}
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#'
#' @export
#'
#'
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotGE(spe, gene_id = "ENSG00000243485", Visium = TRUE) 
#' # Equivalently 
#' plotGE(spe, gene_name = "MIR1302-2HG", Visium = TRUE) 
plotGE <- function(spe, gene_id, gene_name,
                   log_count = FALSE,
                   Visium = FALSE,   # TODO: future argment to use plotVisium
                   ...){
  
  miss_gene_id <- missing(gene_id)
  miss_gene_name <- missing(gene_name)
  
  # TODO: Add some error prevention
  # if there is only one sample
  if(colData(spe)$sample_id |> unique() |> length() != 1)
    stop("plotGE currently only supports SpatialExperiment Object with one sample")
  
  # if either gene_id or gene_name contains more than 1 element, stop for now
  if(!miss_gene_id & length(gene_id)>1)
    stop("plotGE currently only supports plotting one gene.")
  if(!miss_gene_name & length(gene_name)>1)
    stop("plotGE currently only supports plotting one gene.") 
  
  
  
  # TODO: import rowData
  no_gene_id <- ifelse(miss_gene_id, 
                       NA, 
                       !gene_id %in% rowData(spe)$gene_id)
  
  no_gene_name <- ifelse(miss_gene_name,
                         NA,
                         !gene_name %in% rowData(spe)$gene_name)

  
  if(miss_gene_id){
    if(!miss_gene_name){
      if(no_gene_name)
        stop("Can't find gene_name in the dataset. Please check gene_name.")
      else
        gene_id <- rowData(spe)[which(rowData(spe)$gene_name==gene_name), "gene_id"]
    } else
      stop("Please provide either gene_id or gene_name")
  } else {

    if(no_gene_id){
      if(miss_gene_name){
        stop("Can't find gene_id in the dataset. Please check gene_id")
      } else { # gene_name & gene_id both exist
        if(no_gene_name) # Didn't find gene_name
          stop("Please check gene_id and gene_name. Both are not found in the dataset.")
        else{
          warning("gene_id not found, use gene_name.")
          gene_id <- rowData(spe)[which(rowData(spe)$gene_name==gene_name), "gene_id"]
        }
      }
    } else { # gene_id exist
      gene_name_by_id <- rowData(spe)[gene_id, "gene_name"]
      if(!miss_gene_name){
        if(
          gene_name != gene_name_by_id ||
           is.na(gene_name_by_id) # Didn't find gene_name
           ) #
          warning("Inconsistent gene_id, gene_name pair. Gene_id overrides gene_name.")
      }
    }
  }
  
  tmp_spe <- spe[gene_id,]
  
  colData(tmp_spe)$UMI <- counts(spe)[gene_id,]
  
  if(log_count){
    stop("not implemented")
    # TODO: Check if logcounts exits 
  }
  
  
  if(Visium)
    plotVisium(tmp_spe, fill = "UMI", ...)
  else
    plotSpots(tmp_spe, annotate = "UMI", ...)
  
  # if(log_count)
    

  # TODO: change the plot title to indicate which gene is plotted and if logcounts are plotted
}