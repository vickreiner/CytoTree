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
  
  if(keep.non.fluorescent.channels) fcs.data <- fcs.data[,!stringr::str_detect(colnames(fcs.data), "FSC|SSC")]
  ncolnames <- gsub("<|>", "", stringr::str_extract(colnames(fcs.data), "<.*?>"))
  for(i in 1:length(ncolnames)) if(ncolnames[i] == "NA") ncolnames[i] <- gsub("<NA>", "", colnames(fcs.data)[i])
  colnames(fcs.data) <- ncolnames
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




#### ADD GROUPS ####
FlowAppendMetadata <- function(flowtidy = flowtidy,
                               metadata.table.path = "template.metadata.table.xlsx"){
  
  factors.table <- openxlsx::read.xlsx(metadata.table.path)
  
  #check that all names from factors.table are equal to column names to match
  if(all(factors.table$fcs_name %in% unique(paste0(flowtidy$input_file, " / ", flowtidy$sample))) & nrow(factors.table) == (length(unique(paste0(flowtidy$input_file, " / ", flowtidy$sample))))){
    factors.table <- factors.table[match(paste0(flowtidy$input_file, " / ", flowtidy$sample),factors.table$fcs_name),]
    
    
    for(i in 2:ncol(factors.table)){
      flowtidy[,names(factors.table)[i]] <- unlist(factors.table[,i])
    }
    
    return(flowtidy)
    
  } else {
    stop("Please verify that fcs_name in factors.table match the fcs names of the flowdat frame")
  }
  
}



#Import function
importFCS <- function(fcs.directory,
                      metadata.path,
                      filter.markers = TRUE,
                      subsample.factor = 1){
  
  #search for FCS
  fcs.files <- list.files(path = fcs.directory, pattern = ".fcs",full.names = F,include.dirs = F)
  if(length(fcs.files) == 0) stop("No .fcs files found in fcs.directory")
  #load metadata
  load_error <- tryCatch({
    md <- as.data.frame(readxl::read_excel(metadata.path))}, error = function(e){
      stop(paste0("Loading metadata file failed: ", e))})
  if(!"file_name" %in% names(md)) stop("Metadata table is missing the column file_name corresponding to fcs file names to load")
  
  #add rownames for later
  rownames(md) <- md$file_name
  
  #overlay available files with metdata
  if(all(md$file_name %in% fcs.files)){
    message("Found all .fcs files indicated in metadata. Loading...")
    load_error <- tryCatch({
      fcs_raw <- flowCore::read.flowSet(paste0(fcs.directory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)}, error = function(e){
        stop(paste0("Loading fcs files failed: ", e))})
    
  } else {
    stop(paste0("Could not find the following .fcs files in fcs.directory: ", paste0(md$file_name[!md$file_name %in% fcs.files], collapse = "; ")))
  }
  
  #Add panel data from Flowfile
  pd <- data.frame(pData(parameters(fcs_raw[[1]])))
  pd$name <- as.character(unlist(pd$name))
  pd$desc <- as.character(unlist(pd$desc))
  pd$lineage <- TRUE
  pd$functional <- FALSE
  
  if(filter.markers){ #remove FSC etc.
    pd <- subset(pd, desc %in% markernames(fcs_raw))
    #subset the flowset
    fcs_raw <- fcs_raw[,colnames(fcs_raw) %in% pd$name]
  } else {
    pd$desc[is.na(pd$desc)] <- pd$name[is.na(pd$desc)]
  }
  colnames(fcs_raw) <- pd$desc
  names(pd)[2] <- "antigen"
  pd$antigen <- gsub("-", "_",pd$antigen)
  
  #subsample if neccessary
  if(subsample.factor < 1){
    set.seed(8897)
    fcs_raw <- fsApply(fcs_raw, function(ff) {
      idx <- sample.int(nrow(ff), nrow(ff)*subsample.factor)
      ff[idx,]  # alt. ff[order(idx),]
    })
  }
  
  #Add this to output list and general-all-encompassing format (thats the plan....)
  fcs.list <- list("flowSet" = fcs_raw, "metadata" = as.data.frame(md), "panel" = as.data.frame(pd))
  
  #export exprs data frame
  fcs.list[["exprs_frame"]] <- as.data.frame(fsApply((fcs_raw[,]), exprs))
  
  #export exprs data frame metadata
  md$n_cells_per_sample <- fsApply(fcs_raw, nrow)
  fcs.list[["exprs_metadata"]] <- as.data.frame(lapply(md, rep, md$n_cells_per_sample))
  
  print(paste0("Markers: '", paste0(pd$antigen, collapse = "', '"), "'"))
  
  print("Done")
  return(fcs.list) 
}

#### Add lineage markers ####

defineLineageMarkers <- function(FCSlist, lineage.markers){
  
  FCSlist[["panel"]]$lineage <- FCSlist[["panel"]]$antigen %in% lineage.markers
  FCSlist[["panel"]]$functional <- !FCSlist[["panel"]]$antigen %in% lineage.markers
  
  print(paste0("Lineage markers: ", paste0(FCSlist[["panel"]]$antigen[FCSlist[["panel"]]$lineage], collapse = ", ")))
  print(paste0("Functional markers: ", paste0(FCSlist[["panel"]]$antigen[FCSlist[["panel"]]$functional], collapse = ", ")))
  return(FCSlist) 
}


#### Test Transformation ###

FindTransformation <- function(FCSlist,
                               transformation.methods.to.test = list("asinh" = seq(1, 10000, 1000)), 
                               subsample = TRUE,
                               run.id,
                               output.path){
  
  if(any(!names(transformation.methods.to.test) %in% c("asinh"))) stop("Available transformation methods are: asinh")
  if(missing(run.id)) run.id <- paste0("run",round(rnorm(1,0,20),1))
  
  for(t in 1:length(transformation.methods.to.test)){
    curr.method <- names(transformation.methods.to.test)[t]
    curr.cofactors <- transformation.methods.to.test[[t]]
    
    if(class(curr.cofactors) == "list"){
      length.co <- length(curr.cofactors[[1]])
    }else {
      length.co <- length(curr.cofactors)
    }
    
    print("Subsampling and transforming...")
    
    #Subsample or use all cells
    b <- as.list(FCSlist[[4]])
    if(subsample){
      bcurr <- lapply(b, function(x) sample(x, 0.01*length(x)))
    } else {
      bcurr <- b
    }
    
    #expand and transform
    bcurr <- flatten(lapply(bcurr, function(x){
      x <- lapply(curr.cofactors, function(y,x){asinh(x/y)},x )
    }))
    
    
    print("Calculating fits and scores...")
    
    #fit and return fit scores
    fit.out <- lapply(bcurr, function(x){
      out <- locmodes(x,mod0=2,display=FALSE)
      out <- c(out$locations, out$cbw$bw) #extract relevant fit parameters
      out[5] <- out[4] / (abs(mean(max(x), min(x)) - out[2]) + 0.001) #append score
      out[6] <- 1
      out[7] <- 0
      names(out) <- c("m1","thresh","m2","cbw","sc", "upper_bound", "lower_bound")
      return(out)
    })
    
    print("Compiling summary...")
    
    #return final transformed dataset with thresholds
    fit.out <- bind_rows(fit.out)
    fit.out$cf <- rep(curr.cofactors, length(b))
    fit.out$m <- rep(names(b), each = length(curr.cofactors))
    res <- fit.out %>% group_by(m) %>% dplyr::slice(which.min(sc))
    res <- res[order(factor(res$m, levels = unique(fit.out$m))),]
    
    print("Plotting...")
    
    #print PDF with best calculated transformations and thresholds
    
    btraf <- b
    #open PDF (one PDF per transformation method)
    pdf(paste0(output.path,"/TransTest_", run.id,"_",curr.method, ".pdf"), width = 5, height = 5)
    for(i in 1:length(btraf)){
      btraf[[i]] <- asinh(btraf[[i]]/res$cf[i])
      hist(btraf[[i]], main = paste0(res$m[i], " /cf ", res$cf[i], "/thresh ", round(res$thresh[i],2)))
      abline(v = res$thresh[i], col="red", lwd=3, lty=2)
    }
    dev.off()
    
    #write template dataframe for transform parameter
    res$method <- curr.method
    res$manual_cf <- "none"
    res$manual_thresh <- "none"
    res$manual_lower_bound <- "none"
    res$manual_upper_bound <- "none"
    write.xlsx(res,file = paste0(output.path,"/TransTemplate_",run.id,".xlsx"))
  }
  
  print("Wrote template parameter file to output directory. Please add desired method/cofactors in Excel and read the file using the returned path")
  return(output.path) 
}

#### Apply transformation ####

ApplyTransformation <- function(FCSlist,
                                transformation.parameters.path = "none",
                                shiny.transformation_parameters = "none",
                                output.path,
                                run.id){
  
  FCSlist[["trans_exprs_frame"]] <- FCSlist[["exprs_frame"]]
  
  if(transformation.parameters.path != "none"){
    trans.param <- read.xlsx(transformation.parameters.path)
    
    #organise manual input
    trans.param$manual_cf[trans.param$manual_cf == "none"] <- trans.param$cf[trans.param$manual_cf == "none"]
    trans.param$manual_thresh[trans.param$manual_thresh == "none"] <- trans.param$thresh[trans.param$manual_thresh == "none"]
    trans.param$manual_lower_bound[trans.param$manual_lower_bound == "none"] <- trans.param$lower_bound[trans.param$manual_lower_bound == "none"]
    trans.param$manual_upper_bound[trans.param$manual_upper_bound == "none"] <- trans.param$upper_bound[trans.param$manual_upper_bound == "none"]
    
    trans.param$manual_cf <- as.numeric(trans.param$manual_cf)
    trans.param$manual_thresh <- as.numeric(trans.param$manual_thresh)
    trans.param$lower_bound <- as.numeric(trans.param$lower_bound)
    trans.param$upper_bound <- as.numeric(trans.param$upper_bound)
    
  } else if(class(shiny.transformation_parameters) != "character"){
    
    trans.param <- shiny.transformation_parameters
    
    trans.param$m <- trans.param$markers
    trans.param$manual_cf <- as.numeric(trans.param$cofactors)
    trans.param$manual_lower_bound <- as.numeric(trans.param$lower_bounds)
    trans.param$manual_upper_bound <- as.numeric(trans.param$upper_bounds)
    trans.param$manual_thresh <- 0
    
  } else {stop("Please provide input as either a path to an automatic parameters file or a shiny output transformation_parameters data.frame")}
  
  
  #transform flowset
  FCSlist[["trans_flowSet"]] <- fsApply(FCSlist[["flowSet"]], function(x, markers = trans.param$m, cofactors = trans.param$manual_cf, lower_bound = trans.param$lower_bound, upper_bound = trans.param$upper_bound){
    expr <- exprs(x)
    for(i in 1:length(markers)){
      ### formula by Can! See his Github Transformation script! 
      expr[,markers[i]] <- asinh(expr[,markers[i]] / cofactors[i])
      lq <- quantile(expr[,markers[i]], lower_bound[i], names = FALSE, na.rm = TRUE)
      uq <- quantile(expr[,markers[i]], upper_bound[i], names = FALSE, na.rm = TRUE)
      expr[,markers[i]] <- (expr[,markers[i]] - lq)/(uq - lq)
    }
    exprs(x) <- expr
    return(x)
  })
  
  gc()
  #transform exprs matrix
  for(i in 1:nrow(trans.param)){
    curr.marker <- trans.param$m[i]
    curr.param <- trans.param$manual_cf[i]
    
    print(paste0("Transforming ", curr.marker, " with asinh and cofactor/base ", curr.param, " and lower and upper bounds of ", trans.param$manual_lower_bound[i], ", ", trans.param$manual_upper_bound[i]))
    
    
    
    FCSlist[["trans_exprs_frame"]][,curr.marker] <- asinh(FCSlist[["trans_exprs_frame"]][,curr.marker] / curr.param)
    FCSlist[["trans_exprs_frame"]][,curr.marker] <- (FCSlist[["trans_exprs_frame"]][,curr.marker] - lq)/(uq - lq)
  }
  
  #scale matrix
  gc()
  FCSlist[["scaled_exprs_frame"]] <- as.matrix(FCSlist[["trans_exprs_frame"]])
  rng <- colQuantiles(FCSlist[["scaled_exprs_frame"]], probs = c(0.01, 0.99))
  FCSlist[["scaled_exprs_frame"]] <- t((t(FCSlist[["scaled_exprs_frame"]]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  FCSlist[["scaled_exprs_frame"]][FCSlist[["scaled_exprs_frame"]] < 0] <- 0
  FCSlist[["scaled_exprs_frame"]][FCSlist[["scaled_exprs_frame"]] > 1] <- 1
  
  #setting thresholds
  FCSlist[["panel"]]$transformation_method <- "asinh"
  if(all(trans.param$m %in% FCSlist[["panel"]]$antigen)){
    FCSlist[["panel"]]$transformation_cofactor <- trans.param$manual_cf[match(trans.param$m, FCSlist[["panel"]]$antigen)]
    FCSlist[["panel"]]$positive_threshold <- trans.param$manual_thresh[match(trans.param$m, FCSlist[["panel"]]$antigen)]
  }
  
  #return pdf with final transformations and thresholds
  #open PDF (one PDF per transformation method)
  btraf <- as.list(FCSlist[["exprs_frame"]])
  pdf(paste0(output.path,"/TransApply_", run.id, ".pdf"), width = 5, height = 5)
  for(i in 1:length(btraf)){
    btraf[[i]] <- asinh(btraf[[i]]/trans.param$manual_cf[i])
    hist(btraf[[i]], main = paste0(trans.param$m[i], " /cf ", trans.param$manual_cf[i], "/thresh ", round(trans.param$manual_thresh[i],2)))
    abline(v = trans.param$manual_thresh[i], col="red", lwd=3, lty=2)
  }
  dev.off()
  
  print("Done. Returned results in slot trans_exprs_frame of FCSlist. Scaled data is in the slot scaled_exprs_frame. A transformed flowSet instances is in slot trans_flowSet. Also returned pdf TransApply ... for reference")
  return(FCSlist)
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
