# functions for multi-subject differential gene expression based on nebula:
# https://github.com/lhe17/nebula
# https://www.nature.com/articles/s42003-021-02146-6


##########
### function: run_nebula
##########

#' Apply nebula on input data
#'
#' Wrapper around model.matrix, nebula::group_cell, nebula::nebula
#'
#' requires nebula and Rfast to be installed : https://github.com/lhe17/nebula
#'
#' @param counts count matrix Dgc
#' @param metadata ...
#' @param formula_string ...
#' @param subject_id ...
#' @param cell_id ...
#' @param cell_size ...
#' @param ... further params passed directly to nebula
#'
#' @return nebula output list
#'
#' @export
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix colSums

run_nebula = function(counts,metadata, formula_string = "~ treatment",subject_id = "Subject_ID",cell_id="Cell_ID",cell_size = NULL,...){

  # optional use of packages:
  if (!requireNamespace("nebula", quietly = TRUE)) {
    warning("The nebula package must be installed to use this function")
    return(NULL)
  }

  # handle cell_size
  if(is.null(cell_size)){
    cell_size = Matrix::colSums(counts)
  }
  if(is.character(cell_size) & length(cell_size)){
    if(! cell_size %in% colnames(metadata)){
      stop("when providing a string for cell_size, it must be a column name of metadata.")
    }else{
      cell_size = metadata[,cell_size]
    }
  }
  if(is.numeric(cell_size)){
    if(length(cell_size) != ncol(counts)){
      stop("When providing a numeric vector for cell_size, it must have the same length as ncol(counts).")
    }
  }

  # make designamtix
  design_matrix = stats::model.matrix(stats::as.formula(formula_string), metadata)

  # prepare input
  input_list = list(
    count = counts,
    id = metadata[,subject_id],
    pred = design_matrix,
    offset = cell_size
  )
  # ensure right order:
  group_result = nebula::group_cell(input_list$count, input_list$id, input_list$pred, input_list$offset )
  if(!is.null(group_result)){input_list = group_result}
  # run nebula
  nebula_results = nebula::nebula(input_list$count, input_list$id, input_list$pred, input_list$offset,...)
  # return results
  return( nebula_results)
}

##########
### function: FindDEG_nebula
##########

#' Run nebula on an seurat object with condition(s) and multiple subjects
#'
#' Wrapper around run_nebula to execute it for all clusters.
#'
#' requires nebula and Rfast to be installed : https://github.com/lhe17/nebula
#' Requires doParallel and foreach to be installed
#'
#' @param seurat_object a seurat object with all required metadata
#' @param cluster_variable column name in seurat meta.data with the clusters
#' @param sample_variable column name in seurat meta.data with the subject or sample id
#' @param primary_variable column name in seurat meta.data with the main variable to compare (condition, treatment etc.)
#' @param reference_level character string. which level of the primary_variable should be the reference in the model
#' @param other_variables column names in seurat meta.data with other relevant variables to include in nebula model
#' @param features subset to these features (speeds up!). NULL means use all features. Additionally nebula's default setting for lowly expressed gene removal is applied.
#' @param nCounts column name in seurat meta.data with the scaling factor ("offset" in nebula). Can be NULL if all features are used. Then the scaling factor is the total counts number.
#' @param min_cells number of cells (in currently evaluated cluster) for a sample to be valid
#' @param min_pct_valid_samples how many valid samples have to be in cluster in order to include this cluster. 0.5 means 50% of samples must have at least min_cells in the cluster. Set this to 1 to requires all samples to have min_cells.
#' @param min_pct min_pct of cells in cluster expressing a feature
#' @param assay "RNA" or other assay from seurat object
#' @param padjust_method defaults to "fdr" . directly passed to stats::p.adjust
#' @param return_full_result return only marker table or list with all nebula outputs
#' @param numCores number of cores used by foreach
#' @param verbose whether to print verbose messages
#' @param ... further params passed directly to nebula::nebula
#'
#' @return dataframe with DEGs or list with dataframe and list off full nebula results depending on return_full_result
#'
#' @export
#'
#' @importFrom dplyr left_join bind_rows
#' @importFrom Seurat SplitObject
#' @importFrom stats p.adjust

FindDEG_nebula = function(seurat_object,cluster_variable="seurat_clusters",sample_variable="Sample_ID",primary_variable="Condition",reference_level="control",other_variables=c(),features = NULL,nCounts=NULL,min_cells=5,min_pct_valid_samples = 0.5,min_pct=0.05,assay="RNA",padjust_method="fdr",return_full_result = TRUE,numCores = 1,verbose=TRUE,...){

  # optional use of packages:
  if (!requireNamespace("nebula", quietly = TRUE)) {
    warning("The nebula package must be installed to use this function")
    return(NULL)
  }
  # optional use of packages:
  if (!requireNamespace("foreach", quietly = TRUE)) {
    warning("The foreach package must be installed to use this function")
    return(NULL)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    warning("The foreach package must be installed to use this function")
    return(NULL)
  }

  # check that input and parameters are valid
  if(! cluster_variable %in%  colnames(seurat_object@meta.data)){
    stop("Cannot find cluster_variable: '",cluster_variable,"' in meta.data.")
  }
  if(! sample_variable %in%  colnames(seurat_object@meta.data)){
    stop("Cannot find sample_variable: '",sample_variable,"' in meta.data.")
  }
  if(! primary_variable %in%  colnames(seurat_object@meta.data)){
    stop("Cannot find primary_variable: '",primary_variable,"' in meta.data.")
  }
  if(! reference_level %in%  unique(seurat_object@meta.data[,primary_variable])){
    stop("Cannot find reference_level: '",reference_level,"' primary_variable levels in meta.data.")
  }
  if(! all(other_variables %in%  colnames(seurat_object@meta.data))){
    stop("Cannot find one or more of the other_variables in meta.data.")
  }

  # check nCounts
  if(!is.null(nCounts)){
    if(is.character(nCounts) & length(nCounts)==1){
      if(nCounts %in% colnames(seurat_object@meta.data)){
        # do nothing
        #nCounts = seurat_object@meta.data[,nCounts]
      }else{
        nCounts = NULL
      }
    }
  }
  # check features
  if(!is.null(features)){
    features = features[features %in% rownames(seurat_object@assays[[assay]]@counts)]
    if(is.null(nCounts)){
      if("nCount_RNA" %in% colnames(seurat_object@meta.data)){
        nCounts = "nCount_RNA"#seurat_object@meta.data[,"nCount_RNA"]
      }else{
        stop("When subsetting features, please provide a valid meta.data column with library/cell sizes via nCounts.")
      }
    }
  }else{
    features = rownames(seurat_object@assays[[assay]]@counts)
  }

  # split object
  if(verbose){message("Using ",length(features)," features.")}
  if(verbose){message("Splitting seurat object by ",cluster_variable)}
  seurat_object_list = Seurat::SplitObject(seurat_object, split.by = cluster_variable)
  savelength = length(seurat_object_list)
  # prefilter seurat list
  pct_valid_sample_list = sapply(seurat_object_list,function(x,primary_variable,min_cells,min_pct_valid_samples){
    # how many samples have more than min cell in subset
    pct_valid_samples = length(which(table(x@meta.data[,primary_variable]) >= min_cells)) / length(unique(x@meta.data[,primary_variable]))
    return(pct_valid_samples)
  },primary_variable,min_cells,min_pct_valid_samples)
  seurat_object_list = seurat_object_list[pct_valid_sample_list >= min_pct_valid_samples]
  if(verbose){message("Split seurat into ",length(seurat_object_list)," objects (removed: ",savelength-length(seurat_object_list)," with too few cells)")}

  # NEBULA
  if(verbose){message("Running NEBULA in parallel on object list using ",min(numCores,length(seurat_object_list))," cores.")}

  # register cores
  doParallel::registerDoParallel(numCores)  # use multicore, set to the number of our cores
  `%dopar%` = foreach::`%dopar%`
  # run in parallel: I pass errors to be able to set the names of the result list afterwards
  nebula_result_list <- foreach::foreach(subset_seurat_object = seurat_object_list,.errorhandling = "pass",.verbose=FALSE) %dopar% {
    message(subset_seurat_object@project.name)

    # get fcs similar to Seurat
    cells_reference_level = rownames(subset_seurat_object@meta.data)[subset_seurat_object@meta.data[,primary_variable] == reference_level]
    cells_other_levels = rownames(subset_seurat_object@meta.data)[subset_seurat_object@meta.data[,primary_variable] != reference_level]
    foldchange_df = FoldChange_seurat(subset_seurat_object@assays[[assay]]@data,cells.1 = cells_other_levels,cells.2 = cells_reference_level ,features = features)

    # further subset
    alpha.min <- pmax(foldchange_df$pct.1, foldchange_df$pct.2)
    names(x = alpha.min) <- foldchange_df$gene
    subset_features <- names(x = which(x = alpha.min >= min_pct))
    foldchange_df_subset = foldchange_df[foldchange_df$gene %in% subset_features,]

    # prepare input
    subset_metadata = subset_seurat_object@meta.data[,c(sample_variable,primary_variable,other_variables,nCounts),drop=FALSE]
    subset_metadata$Cell_ID = rownames(subset_metadata)
    subset_counts = subset_seurat_object@assays[[assay]]@counts[subset_features,]
    # reorder
    if(reference_level %in% subset_seurat_object@meta.data[,primary_variable] ){
      level_order = c(reference_level,unique(subset_metadata[,primary_variable])[unique(subset_metadata[,primary_variable]) != reference_level])
      # reorder metadata
      subset_metadata = subset_metadata[order(match(subset_metadata[,primary_variable], level_order)),]
      subset_metadata[,primary_variable] = factor(subset_metadata[,primary_variable],levels = level_order) # this is the most important one because stats::model.matrix uses this !
      # change counts to metadata order
      subset_counts = subset_counts[,match(rownames(subset_metadata),colnames(subset_counts))]
    }
    if(length(other_variables) > 0){
      subset_formula = paste0("~ ",paste0(primary_variable," + ",paste0(other_variables,collapse = " + ")))
    }else{
      subset_formula = paste0("~ ",paste0(primary_variable))
    }
    # run nebula
    nebula_res = run_nebula(counts = subset_counts,
                            metadata = subset_metadata,
                            formula = subset_formula,
                            subject_id = sample_variable,
                            cell_id="Cell_ID",
                            cell_size = nCounts,...)

    nebula_res$foldchange_df = foldchange_df_subset
    # result
    nebula_res

  }
  names(nebula_result_list) = names(seurat_object_list)
  if(verbose){message("Finalized NEBULA. Returning results.")}
  # convert nebula results
  deg_dataframe_list = sapply(names(nebula_result_list),function(cluster_name,nebula_res_list,primary_variable){
    current_res = nebula_res_list[[cluster_name]]
    if(length(current_res) >= 2){
      marker_summary = current_res$summary[,c("gene"),drop=FALSE]
      marker_summary$cluster = cluster_name
      marker_summary = as.data.frame(cbind(marker_summary,current_res$summary[,colnames(current_res$summary)[grepl(primary_variable,colnames(current_res$summary))]]))
      marker_fcs = current_res$foldchange_df
      marker_fcs$avg_log2FC = round(marker_fcs$avg_log2FC,5)
      marker_fcs$pct_diff = marker_fcs$pct.1 - marker_fcs$pct.2
      marker_summary = dplyr::left_join(marker_summary,marker_fcs,by="gene")
      # if(length(colnames(current_res$summary)[grepl(paste0("p_",primary_variable),colnames(current_res$summary))]) == 1){
      #   marker_summary$padj = stats::p.adjust(p = marker_summary[,colnames(current_res$summary)[grepl(paste0("p_",primary_variable),colnames(current_res$summary))]],method = padjust_method)
      # }
    }else{
      marker_summary=NULL
    }
    marker_summary
  },nebula_res_list = nebula_result_list,primary_variable=primary_variable,simplify = FALSE)

  # return
  deg_dataframe <- tryCatch({
    deg_dataframe = as.data.frame(do.call(dplyr::bind_rows,deg_dataframe_list))
    deg_dataframe
  },
  error=function(cond) {
    message("Cannot rbind dataframes. Error:",cond)
    return(NULL)
  })

  if(return_full_result){
    return_list = list(degs = deg_dataframe, nebula_result_list = nebula_result_list)
  }else{
    return_list=degs
  }
  return(return_list)
}


##########
### foldchange from seurat
##########

#' Calculate general FC and pcts between two cell groups
#'
#' https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
#'
#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#'
#' @importFrom Matrix rowSums

FoldChange_seurat <- function(
    object,
    cells.1,
    cells.2,
    features = NULL
) {
  if(is.null(features)){
    features <- rownames(x = object)
  }else{
    features <- features[features %in% rownames(x = object)]
  }
  # mean.fxn -> always use log2
  fc.name = "avg_log2FC"
  mean.fxn <- function(x) {
    base = 2
    pseudocount.use = 1
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
  }

  # Calculate percent expressed
  thresh.min <- 0
  pct.1 <- round(
    x = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
      length(x = cells.2),
    digits = 3
  )
  # Calculate fold change
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- data.frame(gene = features,fc = fc, pct.1 = pct.1, pct.2 = pct.2)
  colnames(fc.results)[2] <- fc.name
  return(fc.results)
}

