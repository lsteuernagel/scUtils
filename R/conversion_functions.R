# various utility & helper functions

##########
### get_mouse_human_genes_conversion
##########

#' Match samples between two metadata samples based on a third column (Barcodes)
#'
#' Requires a similarly named column with corresponding barcodes in both tables
#'
#' hsapiens_homolog_perc_id:  %id. target Human gene identical to query gene
#' hsapiens_homolog_perc_id_r1: %id. query gene identical to target Human gene
#'
#' @param host defaults to https://feb2021.archive.ensembl.org Can also use https://nov2020.archive.ensembl.org or "https://www.ensembl.org" for latest
#'
#' @return a data.frame with columns mouse_gene and human_gene
#'
#' @export
#'
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr distinct filter
#' @import magrittr
#'

get_mouse_human_genes_conversion = function(host = "https://feb2021.archive.ensembl.org"){
  # or nov2020.archive.ensembl.org
  mart <- biomaRt::useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host=host)
  mouse_to_human = biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name','hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name','hsapiens_homolog_perc_id','hsapiens_homolog_perc_id_r1'),mart = mart)
  mouse_to_human_filt = mouse_to_human %>%
    dplyr::distinct(external_gene_name, hsapiens_homolog_associated_gene_name,.keep_all = TRUE) %>%
    dplyr::rename(mouse_gene = external_gene_name, human_gene = hsapiens_homolog_associated_gene_name) %>%
    dplyr::filter( human_gene != "") %>% as.data.frame()
  return(mouse_to_human_filt)
}


##########
### convert_species_seurat
##########

#' Converts the gene annotation in a seurat object from mouse to human or vice versa
#'
#' At the moment works only with gene names (symbols). Not Ids.
#'
#' @param seurat_object Seurat Object to convert
#' @param add_as_assay boolean: Whether to add the data as a second assay named 'assay'_'species'
#' @param assay The seurat assay to convert. defaults to RNA
#' @param host the ensembl host . Relevant for the genome version to use. Set to "https://www.ensembl.org" for latest
#' @param verbose boolean whether to print messages
#'
#' @return species converted seurat object (see also add_as_assay)
#'
#' @export
#'
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr distinct filter group_by slice_max
#' @importFrom SeuratObject CreateAssayObject Key CreateSeuratObject
#' @importFrom magrittr %>%
#' @importFrom Seurat NormalizeData
#'

convert_species_seurat = function(seurat_object,add_as_assay = FALSE,assay = "RNA", host = "https://feb2021.archive.ensembl.org",verbose=TRUE){

  # get all features
  feature_names = rownames(seurat_object@assays[[assay]]@counts)

  # check if normalized data available:
  if(sum(seurat_object@assays[[assay]]@data[1,]) %% 1 == 0){
    if(verbose){message("Cannot find normalized data in @data. Running Seurat::NormalizeData.")}
    seurat_object = Seurat::NormalizeData(object = seurat_object)
  }

  # get conversion
  if(verbose){message("Retrieving conversion information via biomart using host ",host)}
  conversion_table = get_mouse_human_genes_conversion(host=host)

  # get species and subset to available genes in seurat
  feature_names = rownames(seurat_object@assays[[assay]]@counts)
  # determine direction of conversion
  if( length(intersect(feature_names,conversion_table$mouse_gene)) > length(intersect(feature_names,conversion_table$human_gene))){
    # subset to features in current seurat
    conversion_table = conversion_table[conversion_table$mouse_gene %in% feature_names,]
    # set variable for later
    from_species = "mouse"
    to_species = "human"
  }else{
    # subset to features in current seurat
    conversion_table = conversion_table[conversion_table$human_gene %in% feature_names,]
    # resolve 1 mouse : many human --> Use the one with highest expression (or first for ties)
    #mean_expression = data.frame(human_gene = rownames(seurat_object@assays$RNA@counts),
    #                         mean_expression= Matrix::rowMeans(seurat_object@assays$RNA@data))
    #conversion_table = conversion_table %>% dplyr::left_join(mean_expression,by="human_gene")
    # set variable for later
    from_species = "human"
    to_species = "mouse"
  }
  if(nrow(conversion_table) < 100){
    from_species = "unknown"
    to_species = "unknown"
  }
  if(verbose){message("Converting from ",from_species," to ",to_species,". Resolving 1 to many relations.")}
  # resolve 1 mouse : many human --> Use the one with highest hsapiens_homolog_perc_id
  conversion_table = conversion_table %>% dplyr::group_by(mouse_gene) %>%
    dplyr::slice_max(order_by = hsapiens_homolog_perc_id,n = 1,with_ties=FALSE)
  # resolve many mouse : 1 human --> Use the one with the highest hsapiens_homolog_perc_id_r1
  conversion_table = conversion_table %>% dplyr::group_by(human_gene) %>%
    dplyr::slice_max(order_by = hsapiens_homolog_perc_id_r1,n = 1,with_ties=FALSE)

  # make an expression matrix with mouse gene names based on mapping
  # make a normalized expression matrix with mouse gene names based on mapping - using the true libray sizes before subsetting!!
  if(verbose){message("Creating converted count matrices.")}
  if(from_species == "human" ){
    count_matrix_converted = seurat_object@assays[[assay]]@counts[conversion_table$human_gene,]
    rownames(count_matrix_converted) = conversion_table$mouse_gene
    data_matrix_converted = seurat_object@assays$RNA@data[conversion_table$human_gene,]
    rownames(data_matrix_converted) = conversion_table$mouse_gene
  }else if(from_species == "mouse" ){
    count_matrix_converted = seurat_object@assays[[assay]]@counts[conversion_table$mouse_gene,]
    rownames(count_matrix_converted) = conversion_table$human_gene
    data_matrix_converted = seurat_object@assays$RNA@data[conversion_table$mouse_gene,]
    rownames(data_matrix_converted) = conversion_table$human_gene
  }else{
    stop("Species must be human or mouse (found less than 100 valid genes).")
  }
  # optionally add as a new array to the seurat object
  if(add_as_assay){
    if(verbose){message("Adding converted matrices as new assay to existing seurat.")}
    # add as new assay
    new_assay = SeuratObject::CreateAssayObject(data = data_matrix_converted,min.cells = 1,min.features = 1,) %>% suppressWarnings()
    new_assay@counts = count_matrix_converted
    SeuratObject::Key(new_assay) = paste0(assay,"_",to_species)

    seurat_object@assays[[paste0(assay,"_",to_species)]] = new_assay

    return(seurat_object)
  }else{
    if(verbose){message("Creating new Seurat object containing only the converted matrices.")}
    converted_seurat_object = SeuratObject::CreateSeuratObject(counts = count_matrix_converted,meta.data = seurat_object@meta.data)
    ## add other data back in
    converted_seurat_object@reductions = seurat_object@reductions
    converted_seurat_object@graphs = seurat_object@graphs
    converted_seurat_object@neighbors = seurat_object@neighbors
    converted_seurat_object@misc = seurat_object@misc

    return(converted_seurat_object)
  }
}

# Tests:
#seurat_object = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq_test/human_hypo_raw.rds")
# seurat_object_MOUSE = convert_species_seurat(seurat_object = seurat_object,add_as_assay = FALSE)
# seurat_object_MOUSE2 = convert_species_seurat(seurat_object = seurat_object,add_as_assay = TRUE)
