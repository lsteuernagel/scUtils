

# pretty_violin = function(seurat_object,gene,cells,cluster=NULL,cluster_col=NULL,group_by=NULL,violin_adjust=0.5,pt_size=1,text_size=10){
#
#   cells = rownames(seurat_object@meta.data)[seurat_object@meta.data[,cluster_col] %in% cluster]
#   data_input = FetchData(seurat_object,slot = "data",vars = gene,cells=cells)
#   data_input[,condition] = seurat_object@meta.data[cells,condition]
#   data_input[,cluster_col] = seurat_object@meta.data[cells,cluster_col]
#   data_input = data_input %>% tidyr::gather(key="gene",value="expression",-!!sym(condition))
#
#   plot = ggplot(data_input[data_input$gene == gene,],aes_string(x=cluster_col,y="expression",fill=condition))+
#     geom_violin(adjust=violin_adjust)+
#     geom_point(size=pt_size,position=position_jitterdodge(jitter.height = 0,jitter.width = 0.2,dodge.width = 0.9),aes_string(group=condition))+
#     theme(axis.title.x = element_blank(),axis.title.y = element_blank(),text = element_text(size = text_size))
#
#   return(plot)
# }
