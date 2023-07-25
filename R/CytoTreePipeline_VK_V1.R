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


#### QC plots ####

QCplots <- function(FCSlist,
                    condition_colors = rainbow(20)){
  
  ggdf <-  data.frame(FCSlist[["exprs_metadata"]],FCSlist[["exprs_frame"]])
  
  ggdf <- pivot_longer(ggdf, cols = c((ncol(FCSlist[["exprs_metadata"]])+1):ncol(ggdf)), names_to = "antigen", values_to = "expression") 
  
  pl_qc <- ggplot(ggdf, aes(x = expression, color = condition, 
                            group = sample_id)) +
    geom_density() +
    facet_wrap(~ antigen, nrow = 4, scales = "free") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = condition_colors)
  
  return(pl_qc)
}

AllvAllplots <- function(FCSlist, subsample = T, 
                         subsample.factor = 0.1, 
                         pdf.filename = "traf1", single.plot = TRUE, bins = 150){
  
  if(subsample){
    ggdf <-  sample_n(FCSlist[["trans_exprs_frame"]], nrow(FCSlist[["trans_exprs_frame"]])*subsample.factor)
  } else {
    ggdf <- FCSlist[["trans_exprs_frame"]]
  }
  if(single.plot){
    pdf(paste0("AllvsAll_singlePlot_", pdf.filename, ".pdf"), width = ncol(ggdf), height = ncol(ggdf))
    print(ggplot(ggdf) +
            geom_hex(aes(x = .panel_x, y = .panel_y), bins = bins) +
            facet_matrix(vars(names(FCSlist[["trans_exprs_frame"]]))) + 
            geom_text(inherit.aes = F,aes(x = .panel_x, y = .panel_y, label = paste0(.panel_x, " X ", .panel_y))) +
            theme_minimal() + 
            theme(axis.text = element_blank(), strip.text = element_text(face = "bold")) + 
            scale_fill_gradientn(colours = c("cadetblue2","gold2","firebrick2"), trans = "log10"))
    dev.off()
  }else {
    combs <- as.data.frame(t(combn(names(ggdf),2)))
    pdf(paste0("AllvsAll_", pdf.filename, ".pdf"), width = 5, height = 5)
    for(i in 1:nrow(combs)){
      print(ggplot(ggdf[,unlist(combs[i,])], aes_string(x = unlist(combs[i,])[1],y = unlist(combs[i,])[2])) +
              geom_hex(bins = bins) +
              scale_fill_gradientn(colours = c("#00007F", "blue", "#007FFF", "cyan",
                                               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), trans = "sqrt") + #Colors from Can! 
              theme_minimal()) 
    }
    dev.off()
  }
}


#### CATALYST CLUSTER own version ####

cluster <- function (x, 
                     features = "type", 
                     xdim = 10, 
                     ydim = 10, 
                     maxK = 20, 
                     verbose = TRUE, 
                     seed = 1,
                     stop.at.SOM = FALSE) 
{
  stopifnot(is(x, "SingleCellExperiment"))
  stopifnot(is.logical(verbose), length(verbose) == 1, vapply(list(xdim, 
                                                                   ydim, maxK, seed), function(arg) is.numeric(arg) && 
                                                                length(arg) == 1, logical(1)))
  features <- .get_features(x, features)
  if (is.null(marker_classes(x))) {
    rowData(x)$marker_class <- factor(c("state", "type")[as.numeric(rownames(x) %in% 
                                                                      features) + 1], levels = c("type", "state", "none"))
  }
  rowData(x)$used_for_clustering <- rownames(x) %in% features
  if (verbose) 
    message("o running FlowSOM clustering...")
  fsom <- ReadInput(flowFrame(t(assay(x, "exprs"))))
  som <- BuildSOM(fsom, colsToUse = features, silent = TRUE, 
                  xdim = xdim, ydim = ydim)
  if(stop.at.SOM == FALSE){
    if (verbose) 
      message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <- suppressWarnings(suppressMessages(ConsensusClusterPlus(t(som$map$codes), 
                                                                 maxK = maxK, reps = 100, distance = "euclidean", seed = seed, 
                                                                 plot = NULL)))
    dev.off()
    k <- xdim * ydim
    mcs <- seq_len(maxK)[-1]
    codes <- data.frame(seq_len(k), map(mc[-1], "consensusClass"))
    codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
    colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", 
                                                      mcs))
    x$cluster_id <- factor(som$map$mapping[, 1])
    metadata(x)$cluster_codes <- codes
    metadata(x)$SOM_codes <- som$map$codes
    metadata(x)$delta_area <- .plot_delta_area(mc)
    return(x)
  } else {
    return(som)
  }
}

#Solution from Helena Crowell https://github.com/HelenaLC/CATALYST/issues/170
ReclusterSpec <- function(sce, 
                          existing.clustering.name = "meta20", 
                          new.clustering.name = "Newmeta20",
                          clusters.to.recluster = c(1:10),
                          features = "type", #passed to CATALYST::cluster
                          xdim = 10, #passed to CATALYST::cluster
                          ydim = 10, #passed to CATALYST::cluster
                          maxK = 20, #passed to CATALYST::cluster
                          seed = 1234){ #passed to CATALYST::cluster
  
  k1 <- existing.clustering.name
  k2 <- new.clustering.name
  kids <- clusters.to.recluster
  
  sce$cell_id <- seq_len(ncol(sce))
  sub <- filterSCE(sce, k = k1, cluster_id %in% kids)
  sub <- CATALYST::cluster(sub, verbose = TRUE, features = features, xdim = xdim, ydim = ydim, maxK = maxK, seed = seed)
  tmp <- sce
  # construct new cluster IDs
  old <- as.character(tmp$cluster_id[sub$cell_id])
  new <- paste(old, sub$cluster_id, sep = ".")
  
  tmp$cluster_id <- as.character(tmp$cluster_id)
  tmp$cluster_id[sub$cell_id] <- new
  
  # construct new cluster codes
  codes_old <- cluster_codes(tmp)[, c("som100", k1)]
  codes_new <- cluster_codes(sub)[, c("som100", k2)]
  
  ss <- strsplit(unique(new), "\\.")
  ss <- lapply(ss, as.numeric)
  codes <- mapply(
    function(i, j) 
    {
      c1 <- codes_old[i, ]
      c2 <- codes_new[j, ]
      data.frame(
        foo = paste(c1$som100, c2$som100, sep = "."),
        meta = paste(c1[[k1]], c2[[k2]], sep = "."))
    }, 
    i = sapply(ss, .subset, 1),
    j = sapply(ss, .subset, 2),
    SIMPLIFY = FALSE)
  codes <- do.call(rbind, codes)
  
  # merge old & new cluster codes
  names(codes) <- names(codes_old)
  for (i in seq(ncol(codes_old)))
    codes_old[[i]] <- as.character(codes_old[[i]])
  codes <- rbind(codes_old, codes)
  
  # reorder factor levels
  for (i in seq(ncol(codes))) {
    o <- order(as.numeric(unique(codes[[i]])))
    codes[[i]] <- factor(codes[[i]], unique(codes[[i]])[o])
  }
  metadata(tmp)$cluster_codes <- codes
  
  print(plotDR(sce, color_by = k1))
  print(plotDR(tmp, color_by = k1))
  
  return(tmp)
}
