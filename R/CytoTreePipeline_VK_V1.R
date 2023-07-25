#### GROUPS TABLE ####

SaveTemplateMetadataTable <- function(path.to.fcs = getwd(),
                                      out.file = "template.metadata.table.xlsx"){
  
  fcs.files <- list.files(path = fcs.directory, pattern = ".fcs",full.names = F,include.dirs = F)
  fcs.files.full <- list.files(path = fcs.directory, pattern = ".fcs",full.names = T,include.dirs = F)
  
  for(i in 1:length(fcs.files)){
    if(!stringr::str_detect(fcs.files.full[i],fcs.files[i])) stop("File ordering issue - Please contact me if you see this error message")
  }
  
  #test names
  if(min(stringr::str_count(fcs.files, "_")) < 2) stop("Naming not compatible. Naming should include at least 2 underscores (_) and look like this: export_YourFCSFileName.fcs")
  fcs.ids <- stringr::str_split(fcs.files, "_", simplify = T)
  if(any(fcs.ids[,1] != "export")) stop("Naming not compatible. All names should start with export_ ....")
  
  md.raw <- data.frame("path_fcs" = fcs.files.full,
                       "raw_fcs" = gsub(".fcs", "",fcs.files), 
                       "loading_nr" = 1:length(fcs.files),
                       "sample_id" = fcs.ids[,2],
                       "population" = unlist(lapply(split(fcs.ids[,2:ncol(fcs.ids)], seq(nrow(fcs.ids))), function(x){
                         return(gsub(".fcs", "",paste0(x[x != ""], collapse = "_")))
                       })),
                       "group_id" = "")

  openxlsx::write.xlsx(md.raw, file = out.file)
}

#### Filter and adjust markers ####

FilterAdjustMarkerNames <- function(fcs.data, 
                                    keep.non.fluorescent.channels = FALSE){
  
  if(keep.non.fluorescent.channels == FALSE) fcs.data <- fcs.data[,!stringr::str_detect(colnames(fcs.data), "FSC|SSC")]
  ncolnames <- gsub("<|>", "", stringr::str_extract(colnames(fcs.data), "<.*?>"))
  for(i in 1:length(ncolnames)) if(ncolnames[i] == "NA") ncolnames[i] <- gsub("<NA>", "", colnames(fcs.data)[i])
  colnames(fcs.data) <- ncolnames
  message(Sys.time(), paste0(" Updated fcs.data markers are: ", paste0(colnames(fcs.data), collapse = ", ")))
  return(fcs.data)
}


#### Subsample Matrix % of each sample ####
SubsampleFCSmatrix <- function(fcs.data = fcs.data,  
                               md = md,
                               subsample.by.md.column = "raw_fcs",
                               subsample.perc = 0.5){
  
  #expand metadata
  if(!subsample.by.md.column %in% names(md)) stop("subsample.by.md.column was not found in column names of input metadata")
  
  ns <- gsub("_\\d+$", "", rownames(fcs.data)) #get rownmames without cell id
  if(!all(ns %in% md$raw_fcs) | !all(md$raw_fcs %in% ns)) stop("More than one raw_fcs file from metadata not found in fcs.data matrix or vice versa. Please check fcs.data via rownames(fcs.data)")
  mdl <- md[match(ns, md$raw_fcs),]
  mdl$cell_id <- rownames(fcs.data)
  
  #fcs.file.rowname <- gsub("_\\d+$", "", rownames(fcs.data))
  fcs.data <- as.data.frame(fcs.data)
  fcs.data <- lapply(split(fcs.data, unlist(mdl[,subsample.by.md.column])), function(ff, subsample.perc){
    idx <- sample.int(nrow(ff), nrow(ff)*subsample.perc)
    return(ff[idx,])}, subsample.perc)
  fcs.data <- bind_rows(fcs.data)
  
  #subset for final table
  mdl <- mdl[mdl$cell_id %in% rownames(fcs.data),]

  message("N cells per file after subsampling: ")
  print(table(gsub("_\\d+$", "", rownames(fcs.data))))
  message("N cells per subsampling group after subsampling: ")
  print(table(unlist(mdl[,subsample.by.md.column])))
  return(fcs.data)
}

#### Apply transformation ####

ApplyTransformation <- function(fcs.data = fcs.data,
                                shiny.transformation_parameters = transformation_parameters){
  
   if(class(shiny.transformation_parameters) != "data.frame") stop("please provide a dataframe for shiny.transformation_parameters")                               
                                  
  FCSlist[["trans_exprs_frame"]] <- FCSlist[["exprs_frame"]]

    trans.param <- shiny.transformation_parameters
    
    trans.param$m <- trans.param$markers
    trans.param$manual_cf <- as.numeric(trans.param$cofactors)
    trans.param$manual_lower_bound <- as.numeric(trans.param$lower_bounds)
    trans.param$manual_upper_bound <- as.numeric(trans.param$upper_bounds)
    trans.param$manual_thresh <- 0
    
  gc()
  #transform exprs matrix
  for(i in 1:nrow(trans.param)){
    message(paste0("Transforming ", trans.param$m[i], " with asinh and cofactor/base ", trans.param$manual_cf[i], " and lower and upper bounds of ", trans.param$manual_lower_bound[i], ", ", trans.param$manual_upper_bound[i]))
    #finf quantiles
    lq <- quantile(fcs.data[,trans.param$m[i]], trans.param$manual_lower_bound[i], names = FALSE, na.rm = TRUE)
    uq <- quantile(fcs.data[,trans.param$m[i]], trans.param$manual_upper_bound[i], names = FALSE, na.rm = TRUE)
    #apply to marker
    fcs.data[,trans.param$m[i]] <- asinh(fcs.data[,trans.param$m[i]] / trans.param$manual_cf[i])
    fcs.data[,trans.param$m[i]] <- (fcs.data[,trans.param$m[i]] - lq)/(uq - lq)
  }
  message("Done! Returning transformed fcs.data")
  return(fcs.data)
}

#### Expand metadata for match with fcs.data ####

ExpandMetadata <- function(md = md,
                           fcs.data = fcs.data){
  ns <- gsub("_\\d+$", "", rownames(fcs.data)) #get rownmames without cell id
  if(!all(ns %in% md$raw_fcs) | !all(md$raw_fcs %in% ns)) stop("More than one raw_fcs file from metadata not found in fcs.data matrix or vice versa. Please check fcs.data via rownames(fcs.data)")
  mdl <- md[match(ns, md$raw_fcs),]
  mdl$cell <- rownames(fcs.data)
  mdl$stage <- mdl$group_id
  return(mdl)
}

#### Proportions barplot function from Platypus ####
#slightly simplified here

CytoProportionsBarplot <- function(cyt = cyt, 
                                   source.group = "cluster.id",
                                   target.group = "group.id",
                                   stacked.plot = TRUE,
                                   verbose = FALSE){
                                     
unique_samples <- unique(cyt@meta.data[,source.group])
unique_clusters <- unique(cyt@meta.data[,target.group])
                                     
cells_per_cluster_per_sample <- list()
for(i in 1:length(unique_samples)){
cells_per_cluster_per_sample[[i]] <- list()
for(j in 1:length(unique_clusters)){
                                         
        if(stacked.plot == F){ #if normal barplot: get % of source group
                    cells_per_cluster_per_sample[[i]][[j]] <- length(which(cyt@meta.data[,source.group]==unique_samples[i] & cyt@meta.data[,target.group]==unique_clusters[j]))/length(which(cyt@meta.data[,source.group]==unique_samples[i])) * 100
        } else { #if stacked barplot: get % of target group
                    cells_per_cluster_per_sample[[i]][[j]] <- length(which(cyt@meta.data[,source.group]==unique_samples[i] & cyt@meta.data[,target.group]==unique_clusters[j]))/length(which(cyt@meta.data[,target.group]==unique_clusters[j])) * 100
              }
           }
      }
 melting <- as.data.frame(reshape2::melt(cells_per_cluster_per_sample))
                                     
  melting$source.group <- "Unkown"
  melting$target.group <- "Unkown"
                                     
   melting$source.group <- as.character(unique_samples[melting$L1])
   melting$target.group <- as.character(unique_clusters[melting$L2])
   colnames(melting) <- c("value", "L2", "Sample", "source", "target")
                                     
   if("factor" %in% class(cyt@meta.data[,source.group])){
    if(verbose) message("Ordering based on existing source group factor levels")
        melting$source <- ordered(as.factor(melting$source), levels = levels(cyt@meta.data[,source.group]))
        print(unique(melting$source))} else {
                                         
   if(verbose) message("Reordering source group, as original column did not contain factor levels")
     melting$source <- ordered(as.factor(melting$source), levels = unique(cyt@meta.data[,source.group]))
          }
    if("factor" %in% class(cyt@meta.data[,target.group])){
     if(verbose) message("Ordering based on existing target group factor levels")
       melting$target <- ordered(as.factor(melting$target), levels = levels(cyt@meta.data[,target.group]))} else {
     if(verbose) message("Reordering target group, as original column did not contain factor levels")
         melting$target <- ordered(as.factor(melting$target), levels = unique(cyt@meta.data[,target.group]))}
                                     
      if(stacked.plot == F){
      if(verbose) message("Returning standard barplot with y axis = % of cells of source group")
      output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = source, y=value, x=target)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.margin = ggplot2::margin(5, 0, 0, 0, "mm")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab(paste0("% of cells of ", source.group)) + ggplot2::xlab(paste0(target.group))
      }else{
      if(verbose) message("Returning stacked barplot with y axis = % of cells of target group")
      output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = source, y=value, x=target)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black", position = "stack") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.margin = ggplot2::margin(5, 0, 0, 0, "mm")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab(paste0("% of cells of ", source.group)) + ggplot2::xlab(paste0(target.group))
   }
return(output.plot)
}

#### Cluster analysis wrapper ####
AnalyzeClustersWrapper <- function(cyt = cyt, file.prefix = "CYT_clustering_3-2", wh = 10:10,cell.size = 15){
  
  #Heatmaps
plotClusterHeatmap(object = cyt, cluster_rows = FALSE, 
                     cluster_cols = FALSE, 
                     markers.to.use = "all", 
                     scale = "column", 
                     group.by = "raw_fcs", 
                     cell.size = cell.size,
                     pdf.name = paste0(file.prefix,"_HEATMAP_byRawFcs"), #set to "none" to avoid saving plot
                     wh = wh) #height : width
  
plotClusterHeatmap(object = cyt, cluster_rows = TRUE, 
                     cluster_cols = FALSE, 
                     markers.to.use = "all", 
                     scale = "column", 
                     group.by = "cluster.id", 
                     cell.size = cell.size,
                     pdf.name = paste0(file.prefix,"_HEATMAP_byCluster"), #set to "none" to avoid saving plot
                     wh = wh) #height : width
                     
plotClusterHeatmap(object = cyt, cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   markers.to.use = "all", 
                   scale = "column", 
                   group.by = "population", 
                   cell.size = cell.size,
                   pdf.name = paste0(file.prefix,"_HEATMAP_byPopulation"), #set to "none" to avoid saving plot
                   wh = wh) #height : width
                                        
  #VLN plots 
  pdf(paste0(file.prefix,"_VLN_allMarkers_byCluster.pdf"), width = 7, height = 4)
       for(i in cyt@markers){
         print(plotViolin(object = cyt, color.by = "cluster.id", marker = i, text.angle = 60) + labs(title = i) + theme(legend.position = "none"))
       }
       dev.off()
                                            
#FeaturePlots
pdf(paste0(file.prefix,"_FeatureUMAP_allMarkers.pdf"), width = 7, height = 5)
for(i in cyt@markers){
   print(plot2D(object = cyt, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
          size = 1, color.by = i) + labs(title = i) + scale_color_viridis_c(option = "A"))
}
dev.off()
                                                
#Cluster-specific markers (1 cluster vs. all)
cl_markers <- runDiffFlex(object = cyt, grouping.column = "population", p.adjust.method = "BH")
openxlsx::write.xlsx(cl_markers,paste0(file.prefix,"_ClusterMarkers.xlsx"))
                                                                     
#Pairwise markers (each cluster vs. each other)
cl_markers <- runDiffFlex(object = cyt, grouping.column = "population", p.adjust.method = "BH", pairwise = TRUE)
openxlsx::write.xlsx(cl_markers,paste0(file.prefix,"_ClusterMarkersPairwise.xlsx"))
message("Done")                    
}


#### Template Annotation table ####

SaveTemplateAnnotationTable <- function(cyt = cyt,
                            out.file = "template.celltypeAnnotation.table.xlsx"){
  
h
  
  md.raw <- data.frame("path_fcs" = fcs.files.full,
                       "raw_fcs" = gsub(".fcs", "",fcs.files), 
                       "loading_nr" = 1:length(fcs.files),
                       "sample_id" = fcs.ids[,2],
                       "population" = unlist(lapply(split(fcs.ids[,2:ncol(fcs.ids)], seq(nrow(fcs.ids))), function(x){
                         return(gsub(".fcs", "",paste0(x[x != ""], collapse = "_")))
                       })),
                       "group_id" = "")
  
  openxlsx::write.xlsx(md.raw, file = out.file)
  
  
  
}

                     