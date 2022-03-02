
# various utility & helper functions

##########
### writeList_to_JSON
##########

#' Write a list of strings (will cast to string) to a json to have a simple conversion to python.
#'
#' @param list_with_rows The list that will be written
#' @param filename the destination without ending!
#'
#' @return nothing
#'
#' @export
#'
#' @importFrom jsonlite toJSON

writeList_to_JSON = function(list_with_rows,filename){
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE, auto_unbox = TRUE,digits = NA)
  writeLines(jsonfile,con = paste0(filename))
}


##########
### downsample_balanced_iterative
##########

#' Helper function for a balanced downsampling procedure.
#'
#' The goal is not to have an equal amount of all cells but instead to sample to a target number with the larger classes being downsampled more strongly.
#'
#' @param metadata metadata with a Cell_ID column and predictor_var
#' @param target_sample_size the target number of cells after downsampling
#' @param predictor_var name of metadata column that will be predicted (celltype)
#' @param stepsize stepsize for iterative downsampling (smaller is more accurate but slightly slower)
#' @param global_seed seed used for random forests
#'
#' @return downsampled metadata
#'
#' @export
#'
#' @import dplyr

downsample_balanced_iterative = function(metadata,target_sample_size,predictor_var,stepsize =100,global_seed=1234){
  message("Running downsample_balanced_iterative")
  metadata = as.data.frame(metadata)
  # group and count full data
  probs_for_sample = metadata %>% dplyr::ungroup() %>% dplyr::group_by(!!dplyr::sym(predictor_var)) %>% dplyr::add_count(name = "per_group") %>% dplyr::select(Cell_ID,per_group,!!dplyr::sym(predictor_var))
  probs_for_sample = as.data.frame(probs_for_sample)
  # first: remove all NAs
  probs_for_sample = probs_for_sample[!is.na(probs_for_sample[,predictor_var]),]
  # init
  sizes = seq(from=nrow(probs_for_sample), to=target_sample_size,by = -stepsize)
  #print(sizes)
  n_classes = length(unique(probs_for_sample[,predictor_var]))
  tmp_samples=probs_for_sample
  # loop over sizes and remove iteratively
  cc=0
  for(size in sizes[1:(length(sizes)-1)]){
    cc=cc+1
    # find classes that should be left out
    which_to_downsample = table(tmp_samples[,predictor_var])>(size/n_classes)
    # addiotionally add some other classes that are close
    which_to_downsample = c(which_to_downsample,table(tmp_samples[,predictor_var])>min(table(tmp_samples[,predictor_var])[names(which_to_downsample[which_to_downsample])]))
    # sample from these two classes
    set.seed(global_seed)
    leave_out = sample(tmp_samples$Cell_ID[tmp_samples[,predictor_var] %in% names(which_to_downsample[which_to_downsample])],stepsize)
    # subset
    tmp_samples = tmp_samples[!tmp_samples$Cell_ID %in% leave_out,]
  }

  # return samples IDs
  return(tmp_samples)
}

##########
### gene_pct_cluster
##########

#' Get per cluster gene pct
#'
#' @param seurat_object seurat object
#' @param genes which genes
#' @param col_name metadata_column
#' @param min_expression min expressio. defaults to 0
#' @param return_long
#'
#' @return matrix or long dataframe
#'
#' @export
#'
#' @import dplyr

gene_pct_cluster = function(seurat_object,genes,col_name,min_expression=0,return_long=FALSE){
  # check
  if(!col_name %in% colnames(seurat_object@meta.data)){stop("Please provide a valid metadata column for grouping !")}
  # get data
  gene_expr = Seurat::FetchData(seurat_object,vars = genes)
  group_factor = seurat_object@meta.data[,col_name]
  # calc pct
  gene_expr[gene_expr > min_expression] <- 1 # set to 1 for occ
  cluster_length = table(group_factor)
  per_Cluster_occ=apply(gene_expr,2,function(x,group_factor){tapply(x,INDEX=group_factor,FUN=sum)},group_factor = group_factor)
  per_Cluster_pct = per_Cluster_occ / as.numeric(cluster_length)
  # return
  if(return_long){
    per_Cluster_pct_long = as.data.frame(per_Cluster_pct) %>% dplyr::mutate(group = rownames(per_Cluster_pct)) %>% tidyr::gather(-group,key="gene",value="pct")
    return(per_Cluster_pct_long)
  }else{
    return(per_Cluster_pct)
  }

}
