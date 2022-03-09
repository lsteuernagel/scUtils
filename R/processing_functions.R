
# processing function for singl cell data (mostly seurat based)

##########
### function: Read10xFormat
##########

#' Read10xFormat
#'
#' Similar to Seurat's ReadMtx but with less error handling and options
#' Sometimes ReadMtx does not work --> then use this
#'
#' @param mtx Name of the mtx file
#' @param cells Name of the cells file
#' @param features Name of the features file
#' @param cell.column Specify which column of cells file to use for cell names; default is 1
#' @param feature.column Specify which column of features files to use for feature/gene names; default is 2
#'
#' @return dgCMatrix with row and col names
#'
#' @export
#'
#' @import Matrix
#'
#' @importFrom data.table fread
#' @importFrom methods as

Read10xFormat= function(mtx,cells,features, cell.column = 1,feature.column = 2){
  matrix = Matrix::readMM(mtx)
  matrix = methods::as(matrix, "dgCMatrix")
  barcodes = data.table::fread(cells,data.table = F,header = FALSE)
  features = data.table::fread(features,data.table = F,header = FALSE)
  colnames(matrix) = barcodes[,cell.column]
  rownames(matrix) = features[,feature.column]
  return(matrix)
}


##########
### function: ReadDGEFormat
##########

#' ReadDGEFormat
#'
#' Read dge file
#'
#' @param dge Name of the dge file
#' @param feature.column Specify which column of the dge matrix will be used as rownames (cell column start from feature.column+1)
#'
#' @return dgCMatrix with row and col names
#'
#' @export
#'
#' @import Matrix
#'
#' @importFrom data.table fread
#' @importFrom methods as

ReadDGEFormat= function(dge, feature.column = 1){
  matrix = data.table::fread(dge,data.table = FALSE)
  features = matrix[,feature.column]
  matrix = matrix[,(feature.column+1):ncol(matrix)]
  rownames(matrix) = features
  matrix = methods::as(as.matrix(matrix), "dgCMatrix")
  return(matrix)
}


##########
### function: seurat_recipe
##########

#' Apply a standard seurat workflow to a seurat object
#'
#' @param seurat_object seurat object
#' @param assay which assay
#' @param nfeatures_vst number of features for HVG detection
#' @param sample_column a column in meatdata with sample or batch ids. defaults to NULL which means feature detection will not be batch aware
#' @param remove_hvgs  boolean. TRUE: use genes_to_remove to remove from HVGS
#' @param genes_to_remove a vector with genes to exclude
#' @param normalize_data normalize data with standard logNorm
#' @param npcs_PCA how many pcs to save
#' @param findClusters run neighbor and cluster detection with seurat standard
#' @param calcUMAP calc standard seurat UMAP
#' @param k.param k neighbors for UMAP and SNN
#' @param clusterRes resolution for clustering
#' @param key name/key for dimred (will name PCA like this and then use for all following steps)
#' @param seed random seed
#'
#' @return processed seurat_object
#'
#' @export
#'
#' @import Seurat

seurat_recipe = function(seurat_object,assay="RNA",nfeatures_vst = 1000,sample_column = NULL,remove_hvgs = FALSE,genes_to_remove = character(0),normalize_data=TRUE,npcs_PCA = 50,calcUMAP=TRUE,k.param=20,findClusters =F,clusterRes = 1,key="pca",seed=123){
  message("Running standard seurat recipe.")
  set.seed(seed)
  if(dim(seurat_object@assays[[assay]]@data)[2]==0){normalize_data=TRUE}
  if(normalize_data){
    seurat_object <- Seurat::NormalizeData(object = seurat_object, verbose = F,assay=assay)
  }else{
    message("Warning: Skipping Normalization")
  }

  if(!is.null(sample_column)){
    message("Find HVGs split by ",sample_column)
    var_features = identify_variable_features(seurat_object,n_hvgs_sizes=nfeatures_vst*2,batch_var=sample_column,assay_name=assay,method="vst",ignore_genes_regex=NULL,returnSeurat=FALSE,seed=seed)
  }else{
    message("Find HVGs")
    seurat_object <- Seurat::FindVariableFeatures(object = seurat_object,assay=assay, selection.method = "vst", nfeatures = nfeatures_vst*2, verbose = F)
    # var features
    var_features = seurat_object@assays[[assay]]@var.features
  }
  # remove features:
  if(remove_hvgs){
    var_features = var_features[!var_features %in% genes_to_remove] # remove problematic genes
  }
  var_features = var_features[1:nfeatures_vst] # subset to actual length

  # optional downstream:
  message("Run Scale and PCA")
  seurat_object <- Seurat::ScaleData(object = seurat_object,assay=assay,features = var_features, verbose = F)
  seurat_object <- Seurat::RunPCA(object = seurat_object,assay=assay, npcs = npcs_PCA,reduction.name = key,features = var_features,reduction.key=key, verbose = F,seed.use = seed)
  if(calcUMAP){
    message("Run Umap...")
    if(key=="pca"){umap_key="umap"}else{umap_key=paste0("umap_",key)}
    seurat_object <- Seurat::RunUMAP(object = seurat_object,assay=assay, reduction = key,reduction.name=umap_key,reduction.key=umap_key, dims = 1:npcs_PCA,n.neighbors=k.param,verbose = F,seed.use = seed)
  }
  if(findClusters){
    message("Run Graph & Clustering...")
    seurat_object <- Seurat::FindNeighbors(seurat_object,reduction = key, dims = 1:npcs_PCA,k.param = k.param,verbose = F)
    seurat_object <- Seurat::FindClusters(seurat_object, resolution = clusterRes,verbose = F,random.seed = seed)

  }
  message("umap name: ",umap_key)
  return(seurat_object)
}

##########
### identify_variable_features
##########

#' Find HVGs with vst.log per Batch (or other Seurat FindVariableFeatures method)
#'
#' @param seurat_object seurat object
#' @param n_hvgs_sizes integer or integer vector of desired feature set sizes
#' @param batch_var var to split by
#' @param assay_name seurat object assay
#' @param method which HVG method : selection.method from Seurat::FindVariableFeatures
#' @param ignore_genes_regex  a regex o filter genes from HVGs
#' @param returnSeurat return seurat
#' @param seed random seed
#'
#' @return seurat_object with variable features in misc$var_features
#'
#' @export
#'
#' @import Seurat

identify_variable_features = function(seurat_object,n_hvgs_sizes=1000,batch_var,assay_name="RNA",method="vst",ignore_genes_regex=NULL,returnSeurat=TRUE,seed=123){
  require(data.table)
  require(Seurat)
  # find features for largest set:
  hvg_genes = max(n_hvgs_sizes)
  #message("Splitting object by: ",batch_var)
  seurat_object@assays[[assay_name]]@var.features = character()
  Idents(seurat_object) <- batch_var
  seurat_object.list <- Seurat::SplitObject(seurat_object, split.by = batch_var)
  #message("Finding shared variable features: ",hvg_genes)
  split.features <- Seurat::SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = hvg_genes,fvf.nfeatures=hvg_genes,selection.method =method,assay=rep(assay_name,length(seurat_object.list)),verbose=FALSE)
  if(!is.null(ignore_genes_regex)){split.features= split.features[!grepl(ignore_genes_regex,split.features)]}
  seurat_object@misc$var_features[[paste0(assay_name,".log.",method,".split_",batch_var,".features.",hvg_genes)]] = split.features

  # run for all set size (see params)
  all_features = seurat_object@misc$var_features[[paste0(assay_name,".log.",method,".split_",batch_var,".features.",max(n_hvgs_sizes))]]
  for(i in 1:length(n_hvgs_sizes)){
    hvg_genes = n_hvgs_sizes[i]
    if(hvg_genes != max(n_hvgs_sizes)){
      message("Setting shared variable features: ",hvg_genes)
      seurat_object@misc$var_features[[paste0(assay_name,".log.",method,".split_",batch_var,".features.",hvg_genes)]] = all_features[1:hvg_genes]
    }
  }

  if(returnSeurat){
    return(seurat_object)
  }else{
    return(seurat_object@misc$var_features[[paste0(assay_name,".log.",method,".split_",batch_var,".features.",hvg_genes)]])
  }
}

##########
### determine_cluster_resolution
##########

#' Find resolution closest to target_cluster_number
#'
#' Wrapper around Seurat::FindClusters
#'
#' @param seurat_object seurat object
#' @param target_cluster_number integer of desired unique clusters
#' @param resolutions resolutions to calculate
#' @param min_cells min cells for a cluster to count (default: 5)
#' @param graph_name graph_name in seurat_object (default: 'RNA_snn')
#' @param cluster_col_name  name for clumn in seurat object that will contain the clusters from the best resolution. (Defaults to 'seurat_clusters')
#' @param return_seurat  return cluster resolution name or seurat object with column cluster_col_name
#' @param seed random seed
#' @param ... other parameters passed on to Seurat::FindClusters
#'
#' @return seurat_object with variable features in misc$var_features
#'
#' @export
#'
#' @import Seurat
#'

determine_cluster_resolution <- function(seurat_object,target_cluster_number,resolutions = c(0.5,0.75,1,1.5,2:10),min_cells=5,graph_name = "RNA_snn",cluster_col_name = "seurat_clusters",return_seurat =TRUE,seed = 1234,...){

  # run clustering (louvain)
  seurat_object = Seurat::FindClusters(seurat_object,resolution = resolutions,graph.name = graph_name,random.seed =seed,verbose=FALSE,...)
  # identify res clostest to target
  n_clusters_per_res = apply(seurat_object@meta.data[,paste0(graph_name,"_res.",resolutions)],2,function(x,min_cells){length(table(x)[table(x)>min_cells])},min_cells=min_cells)
  res_with_target_n = paste0(graph_name,"_res.",resolutions)[which(abs(n_clusters_per_res-target_cluster_number)== min(abs(n_clusters_per_res-target_cluster_number)))[1]]
  # return
  if(return_seurat){
    seurat_object@meta.data[,cluster_col_name] = seurat_object@meta.data[,res_with_target_n]
    return(seurat_object)
  }else{
    return(res_with_target_n)
  }
}


##########
### apply_DoubletFinder
##########

#' Run DoubletFinder on an Seurat object
#'
#' Wrapper around Seurat::FindClusters
#'
#' see also: https://github.com/chris-mcginnis-ucsf/DoubletFinder
#'
#' @param seurat_object seurat object
#' @param npcs_PCA how many PCs to use (DoubletFinder wil run its own PCA!)
#' @param pN_fixed pN to pass to DoubletFidner
#' @param pK_max pK will be estimated, this parameter allows to se a maximum cap
#' @param doublet_formation_rate how many doublets are expected. Deafults to 0.075
#' @param adjust_nExp whether to adjust number of doublets with modelHomotypic
#' @param doublet_cluster_tresh  defaults to NULL. if not NULL: Clusters with more than doublet_cluster_tresh pct of cells will completely be considered Doublets
#' @param cluster_column  metadata column with clusters
#' @param return_seurat whether to return seurat with 'Doublet' column
#'
#' @return seurat with 'Doublet' column or named vector
#'
#' @export
#'
#' @import Seurat dplyr
#'

apply_DoubletFinder <- function(seurat_object,npcs_PCA=50,pN_fixed=0.25,pK_max=NULL,doublet_formation_rate=0.075,adjust_nExp = FALSE, doublet_cluster_tresh =NULL,cluster_column="seurat_clusters",return_seurat=TRUE){

  # optional use of packages:
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    warning("The DoubletFinder package must be installed to use this function")
    return(NULL)
  }
  ## pK Identification (no ground-truth)
  sweep.res.seurat_object <- DoubletFinder::paramSweep_v3(seurat_object, PCs = 1:npcs_PCA, sct = FALSE)
  sweep.stats_seurat_object <- DoubletFinder::summarizeSweep(sweep.res.seurat_object, GT = FALSE)
  bcmvn_seurat_object <- DoubletFinder::find.pK(sweep.stats_seurat_object)
  pK_at_BCmetrci_max = as.numeric(as.character(bcmvn_seurat_object$pK[bcmvn_seurat_object$BCmetric == max(bcmvn_seurat_object$BCmetric)]))
  if(!is.null(pK_max)){if(pK_at_BCmetrci_max > pK_max){pK_at_BCmetrci_max = pK_max}}

  ## Homotypic Doublet Proportion Estimate
  nExp_poi <- round(doublet_formation_rate*nrow(seurat_object@meta.data))
  if(adjust_nExp){
    homotypic.prop <- DoubletFinder::modelHomotypic(seurat_object@meta.data[,cluster_column])
    nExp_poi <- round(nExp_poi*(1-homotypic.prop))
  }

  ## Run DoubletFinder
  seurat_object <- DoubletFinder::doubletFinder_v3(seurat_object, PCs = 1:npcs_PCA, pN = pN_fixed, pK = pK_at_BCmetrci_max, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  reuse.pANN = paste0("pANN_",pN_fixed,"_",pK_at_BCmetrci_max,"_",nExp_poi)
  # final Doublet column:
  seurat_object@meta.data$Doublet = seurat_object@meta.data[,paste0("DF.classifications_",gsub("pANN_","",reuse.pANN))]

  # optionally add cells in clusters with many doublets
  if(!is.null(doublet_cluster_tresh)){
    if(doublet_cluster_tresh>1){doublet_cluster_tresh = 1}
    # calculate pct of doublets per cluster
    doublet_stats_per_cluster = seurat_object@meta.data %>% dplyr::select(cluster = !!rlang::sym(cluster_column),doublet_finder_classification = !!rlang::sym(paste0("DF.classifications_",gsub("pANN_","",reuse.pANN)))) %>%
      dplyr::group_by(cluster) %>% dplyr::add_count(name="cells_per_cluster") %>% dplyr::filter(doublet_finder_classification=="Doublet") %>% dplyr::add_count(name="doublet_per_cluster") %>%
      dplyr::distinct(cluster,cells_per_cluster,doublet_per_cluster) %>% dplyr::mutate(doublet_pct = doublet_per_cluster / cells_per_cluster)
    # also filter out full cluster with certain pct!
    cells_in_doublet_clusters = rownames(seurat_object@meta.data)[seurat_object@meta.data[,cluster_column] %in% doublet_stats_per_cluster$cluster[doublet_stats_per_cluster$doublet_pct >= doublet_cluster_tresh]]
    seurat_object@meta.data$Doublet[rownames(seurat_object@meta.data) %in% cells_in_doublet_clusters] = "Doublet"
  }
  if(return_seurat){
    return(seurat_object)
  }else{
    retvec = as.character(seurat_object@meta.data$Doublet)
    names(retvec) = rownames(seurat_object@meta.data)
    return(retvec)
  }

}



