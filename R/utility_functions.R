
# various utility & helper functions

##########
### match_sample_names
##########

#' Match samples between two metadata samples based on a third column (Barcodes)
#'
#' Requires a similarly named column with corresponding barcodes in both tables
#'
#' @param table_A a data.frame (metadata from single cell data)
#' @param table_B a second data.frame (metadata from single cell data)
#' @param sample_col_A column name with sample names in A
#' @param sample_col_B column names with sample names in B
#' @param barcode_col column with barcodes (in btoh tables
#' @param min_pct minimum of matching barcodes to set as sample
#' @param silent messages
#'
#' @return a data.frame ith matched samples (NA if no matchin sample in B for A). Will drop samples from B that cannot be matched to A.
#'
#' @export
#'

match_sample_names = function(table_A,table_B,sample_col_A,sample_col_B,barcode_col = "Barcode",min_pct=0.5,silent=FALSE){

  A_barcode_list = split(table_A[,barcode_col],table_A[,sample_col_A])
  B_barcode_list = split(table_B[,barcode_col],table_B[,sample_col_B])

  resvec=vector()
  for(i in 1:length(A_barcode_list)){
    current_bcs = A_barcode_list[[i]]
    overlap_bcs = sapply(B_barcode_list,function(x,el){length(intersect(x,el))},el=current_bcs)
    likely_sample=names(which(overlap_bcs==max(overlap_bcs) & max(overlap_bcs) > min_pct*length(current_bcs)))
    if(length(likely_sample)>0){
      if(!silent){message("Associating ",names(A_barcode_list)[i]," with ",likely_sample, " based on an overlap of ",max(overlap_bcs)," out of ",length(current_bcs)," barcodes.")}
      resvec[names(A_barcode_list)[i]] =  likely_sample
    }else{
      resvec[names(A_barcode_list)[i]] =  NA
    }
  }
  resdf = data.frame(Samples_A = names(resvec),Samples_B = resvec)
  return(resdf)
}

##########
### infer_sex
##########

#' Infer sex per sample based on Xist expression in cells
#'
#' The goal is not to have an equal amount of all cells but instead to sample to a target number with the larger classes being downsampled more strongly.
#'
#' @param seurat_object a seurat object
#' @param sample_column column with sample annotation
#' @param min_xist_female min pct of cells in sample that are xist positive to declare female
#' @param max_xist_male max pct of cells in sample that are xist positive to declare male
#'
#' @return a vector with either F,M or U per cell (but based on the samples)
#'
#' @export
#'
#' @import dplyr Seurat

## function to infer sex per sample:
infer_sex = function(seurat_object,sample_column,min_xist_female = 0.7,max_xist_male = 0.1){
  xist_expr = Seurat::FetchData(seurat_object,"Xist")[,1]
  xist_expr[xist_expr>0]=1
  tmp_df = data.frame(Cell_ID = rownames(seurat_object@meta.data), Sample_ID =  seurat_object@meta.data[,sample_column],xist_expr=xist_expr) %>% dplyr::group_by(Sample_ID) %>%
    dplyr::add_count(name="n_cells")  %>% dplyr::mutate(n_xist = sum(xist_expr)) %>% dplyr::mutate(xist_pct = n_xist / n_cells) %>%
    dplyr::mutate(inferred_sex = dplyr::case_when(xist_pct > min_xist_female ~ "F", xist_pct < max_xist_male ~ "M", TRUE ~ "U"))
  return(tmp_df$inferred_sex)
}

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
### find_children
##########

#' Recursively walk through edgelist and return children of specified node(s).
#'
#' @param nodes vector of node(s) present in edges for which all children will be determined
#' @param edges edgelist with two columns: 'to' and 'from'
#'
#' @return vector with children ids
#'
#' @export
#'

find_children = function(nodes,edges){
  current_children = edges$to[edges$from %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_children)>0){
    all_children = c(current_children,find_children(current_children,edges))
  }else{
    all_children = current_children
  }
  return(all_children)
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
  if(!"Cell_ID" %in% colnames(metadata)){
    warning("Error: metadata requires a column 'Cell_ID'")
    return(NULL)
  }
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
#' @param return_long return data frame in lonf format
#'
#' @return matrix or long dataframe
#'
#' @export
#'
#' @import dplyr
#' @importFrom tidyr gather

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
    colnames(per_Cluster_pct) = make.unique(colnames(per_Cluster_pct))
    per_Cluster_pct_long = as.data.frame(per_Cluster_pct) %>% dplyr::mutate(group = rownames(per_Cluster_pct)) %>% tidyr::gather(-group,key="gene",value="pct")
    return(per_Cluster_pct_long)
  }else{
    return(as.data.frame(per_Cluster_pct))
  }

}

##########
### CalculateMultScore
##########

#' Calculate a gene set score for each cell in a Seurat object based on teh product (if one gene is missing the score will be zero)
#'
#' Useful for highlighting all cells that express two or three genes together.
#' Scores for longer gene sets will likely be zero, but the allowed_zeros(pct) argument can be used to allow a certain number of zeros (and calculation of the product only on the other non-zero genes).
#'
#' @param seurat_object seurat object
#' @param features which genes to calculate score on
#' @param cells subset to tehses cells (defaults to NULL which emans all)
#' @param allowed_zeros how many zeros (integer) are allowed (will be ignored when calculating product)
#' @param allowed_zeros_pct how many zeros (pct of number of features) are allowed. Number between 0 and 1. Rounded with floor. If larger than 0 this will overwrite allowed_zeros
#'
#' @return a numeric vector with scores per cell in seurat_object
#'
#' @export
#'
#' @importFrom Seurat FetchData

## with normalization:
CalculateMultScore = function(seurat_object,features,cells=NULL,allowed_zeros = 0,allowed_zero_pct = 0){
  gene_data = Seurat::FetchData(seurat_object,vars = features,cells = cells)
  if(allowed_zero_pct > 0){
    allowed_zeros = floor(ncol(gene_data) * max(1,allowed_zero_pct))
  }
  if(allowed_zeros == 0){
    multScore = apply(gene_data,1,prod)
    multScore = multScore^(1/ncol(gene_data)) # normalize
  }else{
    multScore = apply(gene_data,1,function(x,subtr){
      zero_idx = which(x == 0) # get all zeros
      # decide which ones to drop (maybe not all)
      if(length(zero_idx)>subtr){
        drop_idx = zero_idx[(length(zero_idx)-subtr+1):length(zero_idx)]
      }else{
        drop_idx = zero_idx
      }
      y = x[-drop_idx] # exclude subtr zeros from vector
      z = prod(y)^(1/length(y)) # normalize
    },subtr=allowed_zeros)
  }
  return(multScore)
}
