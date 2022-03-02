
# processing function for singl cell data (mostly seurat based)

##########
### function: seurat_recipe
##########

#' Apply a standard seurat workflow to a seurat object
#' @param seurat_object
#' @param assay
#' @param nfeatures_vst
#' @param clean_hvg
#' @param normalize_data
#' @param npcs_PCA
#' @param findClusters
#' @param calcUMAP
#' @param clusterRes
#' @param key
#' @param seed
#'
#' @return processed seurat_object
#'
#' @export
#'
#' @import Seurat

# run seurat standard pipeline on an onject
seurat_recipe = function(seurat_object,assay="RNA",nfeatures_vst = 2000,clean_hvg=F,normalize_data=F,npcs_PCA = 50,findClusters =F,calcUMAP=TRUE,clusterRes = 1,key="pca",seed=123){
  message("Running standard seurat recipe.")
  set.seed(seed)
  if(dim(seurat_object@assays[[assay]]@data)[2]==0){normalize_data=TRUE}
  if(normalize_data){
    seurat_object <- NormalizeData(object = seurat_object, verbose = F,assay=assay)
  }else{
    message("Skipping Normalization")
  }
  message("Find HVGs")
  seurat_object <- FindVariableFeatures(object = seurat_object,assay=assay, selection.method = "vst", nfeatures = nfeatures_vst, verbose = F)
  # ignore mt and rpl/s genes in hvg
  if(clean_hvg){
    idx = grep("mt-|Rpl|Rps",seurat_object@assays[[assay]]@var.features)
    if(length(idx)>0){
      print(paste0("Ignoring ",length(idx)," rp and mt genes in hvgs.."))
      seurat_object@assays[[assay]]@var.features = seurat_object@assays[[assay]]@var.features[-idx]
      seurat_object@assays[[assay]]@meta.features$vst.variable[idx] = FALSE
    }
  }
  message("Run Scale and PCA")
  seurat_object <- ScaleData(object = seurat_object,assay=assay, verbose = F)
  seurat_object <- RunPCA(object = seurat_object,assay=assay, npcs = npcs_PCA,reduction.name = key,reduction.key=key, verbose = F,seed.use = seed)
  message("Run Umap...")
  if(key=="pca"){umap_key="umap"}else{umap_key=paste0("umap_",key)}
  if(calcUMAP){
    seurat_object <- RunUMAP(object = seurat_object,assay=assay, reduction = key,reduction.name=umap_key,reduction.key=umap_key, dims = 1:npcs_PCA,verbose = F,seed.use = seed)
  }
  if(findClusters){
    message("Run Graph & Clustering...")
    seurat_object <- FindNeighbors(seurat_object,reduction = key, dims = 1:npcs_PCA,k.param = 20,verbose = F)
    seurat_object <- FindClusters(seurat_object, resolution = clusterRes,verbose = F,random.seed = seed)

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
#' @param batch_var where to store the results
#' @param assay_name seurat object assay
#' @param method which HVG method : selection.method from Seurat::FindVariableFeatures
#' @param ignore_genes_regex  a regex o filter genes from HVGs
#' @param seed random seed
#'
#' @return seurat_object with variable features in misc$var_features
#'
#' @export
#'
#' @import Seurat

identify_variable_features = function(seurat_object,n_hvgs_sizes=1000,batch_var,assay_name="RNA",method="vst",ignore_genes_regex=NULL,seed=123){
  require(data.table)
  require(Seurat)
  # find features for largest set:
  hvg_genes = max(n_hvgs_sizes)
  message("Splitting object by: ",batch_var)
  seurat_object@assays[[assay_name]]@var.features = character()
  Idents(seurat_object) <- batch_var
  seurat_object.list <- Seurat::SplitObject(seurat_object, split.by = batch_var)
  message("Finding shared variable features: ",hvg_genes)
  split.features <- Seurat::SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = hvg_genes,fvf.nfeatures=hvg_genes,selection.method =method,assay=rep(assay_name,length(seurat_object.list)))
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
  return(seurat_object)
}
