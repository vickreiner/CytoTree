#### GROUPS TABLE ####

SaveTemplateMetadataTable <- function(path.to.fcs = getwd(),
                                      out.file = "template.metadata.table.xlsx"){
  
  fcs.files <- list.files(path = path.to.fcs, pattern = ".fcs",full.names = F,include.dirs = F)
  fcs.files.full <- list.files(path = path.to.fcs, pattern = ".fcs",full.names = T,include.dirs = F)
  
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

####  Set subsampling column ####
SetSubsamplingColumn <- function(cyt = cyt,  
                               subsample.by.md.column = "raw_fcs",
                               subsample.perc = 0.5){
  
  #expand metadata
  if(!subsample.by.md.column %in% names(cyt@meta.data)) stop("subsample.by.md.column was not found in column names of input metadata")
  
  #fcs.file.rowname <- gsub("_\\d+$", "", rownames(fcs.data))
  cyt@meta.data <- as.data.frame(cyt@meta.data)
  subsampled.md <- lapply(split(cyt@meta.data, unlist(cyt@meta.data[,subsample.by.md.column])), function(ff, subsample.perc){
    idx <- sample.int(nrow(ff), nrow(ff)*subsample.perc)
    return(ff[idx,])}, subsample.perc)
  subsampled.md <- bind_rows(subsampled.md)
  
  #set column
  cyt@meta.data$dowsample <- 0
  cyt@meta.data$dowsample[rownames(cyt@meta.data) %in% rownames(subsampled.md)] <- 1
  
  message("N cells selected in dowsample column: ")
  print(table(cyt@meta.data$dowsample))
  return(cyt)
}

#### Apply transformation ####

ApplyTransformation <- function(fcs.data = fcs.data,
                                shiny.transformation_parameters = transformation_parameters){
  
   if(class(shiny.transformation_parameters) != "data.frame") stop("please provide a dataframe for shiny.transformation_parameters")                               
                                  
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
  
  if(!all(ns %in% md$raw_fcs)) stop(paste0("File names from rownames of fcs.data were not found in md$raw_fcs: ", paste0(unique(ns[!ns %in% md$raw_fcs]),collapse = ", ")))
  if(!all(md$raw_fcs %in% ns)) stop(paste0("File names from md$raw_fcs of were not found in rownames of fcs.data : ", paste0(unique(md$raw_fcs[!md$raw_fcs %in% ns]),collapse = ", ")))
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
AnalyzeClustersWrapper <- function(cyt = cyt, file.prefix = "CYT_clustering_3-2", wh = 10:10,cell.size = 15, downsample = FALSE, downsample.dimred = TRUE){
  
  #Heatmaps
plotClusterHeatmap(object = cyt, cluster_rows = FALSE, 
                     cluster_cols = FALSE, 
                     markers.to.use = "all", 
                     downsample = downsample, 
                     scale = "column", 
                     group.by = "raw_fcs", 
                     cell.size = cell.size,
                     pdf.name = paste0(file.prefix,"_HEATMAP_byRawFcs"), #set to "none" to avoid saving plot
                     wh = wh) #height : width
plotClusterHeatmap(object = cyt, cluster_rows = TRUE, 
                     cluster_cols = FALSE, 
                     markers.to.use = "all",  
                     downsample = downsample, 
                     scale = "column", 
                     group.by = "cluster.id", 
                     cell.size = cell.size,
                     pdf.name = paste0(file.prefix,"_HEATMAP_byCluster"), #set to "none" to avoid saving plot
                     wh = wh) #height : width
plotClusterHeatmap(object = cyt, cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   markers.to.use = "all", 
                   downsample = downsample, 
                   scale = "column", 
                   group.by = "population", 
                   cell.size = cell.size,
                   pdf.name = paste0(file.prefix,"_HEATMAP_byPopulation"), #set to "none" to avoid saving plot
                   wh = wh) #height : width
                                  
  #VLN plots 
  pdf(paste0(file.prefix,"_VLN_allMarkers_byCluster.pdf"), width = 7, height = 4)
       for(i in cyt@markers){
         print(plotViolin(object = cyt, color.by = "cluster.id", marker = i, text.angle = 60, downsample = downsample) + labs(title = i) + theme(legend.position = "none"))
       }
       dev.off()
                                            
#FeaturePlots
pdf(paste0(file.prefix,"_FeatureUMAP_allMarkers.pdf"), width = 7, height = 5)
for(i in cyt@markers){
   print(plot2D(object = cyt, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
          size = 1, color.by = i, downsample = downsample.dimred) + labs(title = i) + scale_color_viridis_c(option = "A"))
}
dev.off()
                                                
#Cluster-specific markers (1 cluster vs. all)
cl_markers <- runDiffFlex(object = cyt, grouping.column = "cluster.id", p.adjust.method = "BH")
openxlsx::write.xlsx(cl_markers,paste0(file.prefix,"_ClusterMarkers.xlsx"))

#Cluster-specific markers (1 cluster vs. all)
cl_markers <- runDiffFlex(object = cyt, grouping.column = "population", p.adjust.method = "BH")
openxlsx::write.xlsx(cl_markers,paste0(file.prefix,"_PopulationMarkers.xlsx"))
                                                                     
#Pairwise markers (each cluster vs. each other)
#cl_markers <- runDiffFlex(object = cyt, grouping.column = "population", p.adjust.method = "BH", pairwise = TRUE)
#openxlsx::write.xlsx(cl_markers,paste0(file.prefix,"_ClusterMarkersPairwise.xlsx"))
message("Done")                    
}


#### Template Annotation table ####

SaveTemplateAnnotationTable <- function(cyt = cyt,
                            out.file = "template.celltypeAnnotation.table.xlsx"){
  md.raw <- data.frame("cluster.id" = sort(unique(cyt@meta.data$cluster.id)),
                       "primary.celltype" = "undefined", 
                       "secondary.celltype" = "undefined")
  openxlsx::write.xlsx(md.raw, file = out.file)
}

#### Add cluster annotations ####

AnnotateClusters <- function(cyt, 
                             AnnotationTable.file = "celltypeAnnotation.table.xlsx"){
  
  adf <- openxlsx::read.xlsx(AnnotationTable.file)
  adf <- adf[match(cyt@meta.data$cluster.id,adf$cluster.id),]
  rownames(adf) <- rownames(cyt@meta.data)
  cyt <- addMetaData(cyt, meta.info = adf[,2], name = "primary.celltype")
  cyt <- addMetaData(cyt, meta.info = adf[,3], name = "secondary.celltype")
  return(cyt)
}

#### Subclustering function ####

SubclusterCyt <- function(cyt = cyt, cluster.to.subcluster = c(1), k = 5,...){
  
  curr.cyt <- subsetCYT(cyt, cells = rownames(cyt@meta.data)[cyt@meta.data$cluster.id %in% cluster.to.subcluster])
  print(nrow(curr.cyt@log.data))
  old.cls <- curr.cyt@meta.data$cluster.id
  curr.cyt <- curr.cyt %>% runCluster(k = k,...)
  new.cls <- paste0(old.cls, "_", curr.cyt@meta.data$cluster.id)
  
  #remap to complete cyt
  cyt@meta.data$cluster.id.old <- cyt@meta.data$cluster.id
  cyt@meta.data$cluster.id[rownames(cyt@meta.data) %in% rownames(curr.cyt@meta.data)] <- new.cls
  curr.cyt <- NULL
  return(cyt)
}


#### Generate flowtidy dataframe ####

GenerateFlowtidyDF <- function(cyt = cyt, 
                               cyt.slot = "raw.data", #can be "log.data" or "raw.data"
                               sample.id.column = "sample_id", 
                               grouping.column = "group_id", #analysis across these
                               iteration.column = "primary.celltypes" #iterating over these 
                               ){ #(frequencies as percentage of total nr. of cells of a given sample.id.column entry)

  
  #Get marker expression table
  grouping.vector <- apply(cyt@meta.data[,c(sample.id.column, grouping.column, iteration.column)] ,1, paste , collapse = "-//-")
  print(unique(grouping.vector))
  if(cyt.slot == "raw.data"){
    rawmeans <- lapply(split(as.data.frame(cyt@raw.data), grouping.vector), function(ff) return(colMeans(ff)))
    rawmeans <- bind_rows(rawmeans)  
    colnames(rawmeans) <- paste0("RawExprs_", colnames(rawmeans))
  } else {
    rawmeans <- lapply(split(as.data.frame(cyt@log.data), grouping.vector), function(ff) return(colMeans(ff)))
    rawmeans <- bind_rows(rawmeans)  
    colnames(rawmeans) <- paste0("LogExprs_", colnames(rawmeans))
  }
  
  rawmeans <- cbind(rawmeans, cyt@meta.data[duplicated(grouping.vector) == FALSE,c(sample.id.column, grouping.column, iteration.column)])  
  
  #get frequencies
    cell.nr <- as.data.frame(cyt@meta.data[,c(sample.id.column, grouping.column, iteration.column)] %>% group_by_all() %>% summarise(Count = n()) %>% ungroup() %>% group_by_at(1) %>% mutate(PercentageOfSampleID = round(Count/sum(Count)*100, 3)))
  
  if(all(rawmeans[sample.id.column] == cell.nr[sample.id.column]) & 
     all(rawmeans[grouping.column] == cell.nr[grouping.column]) & 
     all(rawmeans[iteration.column] == cell.nr[iteration.column])){
    rawmeans <- cbind(rawmeans, cell.nr)
  } else{
    rawmeans$ordercol <- apply(rawmeans[,c(sample.id.column, grouping.column, iteration.column)] ,1, paste , collapse = "-//-")
    cell.nr$ordercol <- apply(cell.nr[,c(sample.id.column, grouping.column, iteration.column)] ,1, paste , collapse = "-//-")
    rawmeans <- rawmeans[match(cell.nr$ordercol,rawmeans$ordercol),]
    
   # return(list(rawmeans, cell.nr))
    
    if(all(rawmeans$ordercol == cell.nr$ordercol)){
    rawmeans <- cbind(cell.nr[,!colnames(cell.nr) == "ordercol"],rawmeans[,!colnames(rawmeans) %in% c(sample.id.column, grouping.column, iteration.column, "ordercol")])
    } else {stop("Merging error")}
  }
    
  #make tidy
  rawmeans <- pivot_longer(rawmeans, cols = 4:ncol(rawmeans), names_to = "metric", values_to = "value")
  names(rawmeans)[1:3] <- c("sample", "group", "population")
  rawmeans$level <- 1 #placeholder
  rawmeans$sample <- factor(rawmeans$sample, levels = unique(rawmeans$sample))
  rawmeans$group <- factor(rawmeans$group, levels = unique(rawmeans$group))
  rawmeans$population <- factor(rawmeans$population, levels = unique(rawmeans$population))
  return(rawmeans)
}

#### Testing wrapper flowtidy V3 ####

FlowTestGroups <- function(flowtidy = flowtidy,
                           groups.column = "group", #x axis of resulting plots
                           populations.column = "population", #plot list grouped by this
                           samples.column = "sample",
                           population.levels.column = "level",
                           metrics.columns = "metric",
                           values.column = "value",
                           metrics = "all", # "all" or an array of entries from metrics.column
                           population.levels = "all", #"all" or one or more of entries in population.levels.column
                           populations = "all", #"all" or one or more of entries in populations.column
                           paired.samples = FALSE, #passed to pairwise.wilcox.test etc. 
                           p.adjust.method = "BH", #Adjusting over all generated p values passed to pairwise.wilcox.test etc.. Can be: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                           p.adjust.method.pairwise.wilcoxon = "BH", #Adjusting over all p values within one AOV test Can be: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                           signif.level = 0.05, #used to determine significance on corrected and uncorrected (aov) p values
                           NA.remove = TRUE, #TRUE to remove NAs, FALSE to replace NAs with 0
                           run.id = "run1"){
  
  #Function organisation -> 
  #2 groups => pairwise wilcoxon (if less than n = 3 -> fold change) -> summary table -> p value correction -> volcano + individual boxplots
  #3 groups -> AOV as pretest -> summary table -> proceed only with significant -> post hoc pairwise wilcoxon -> summary table -> 3d volcano + normal boxplots
  #making sure to adjust names
  ft <- data.frame("g" = unlist(flowtidy[groups.column]),
                   "s" = unlist(flowtidy[samples.column]),
                   "p" = unlist(flowtidy[populations.column]),
                   "l" = unlist(flowtidy[population.levels.column]),
                   "m" = unlist(flowtidy[metrics.columns]),
                   "v" = unlist(flowtidy[values.column]))
  
  #for percentile functions later
  percfnc <- list()
  for(i in 1:length(unique(ft$m))){
    percfnc[[unique(ft$m)[i]]] <- ecdf(ft$v[ft$m == unique(ft$m)[i]])
  }
  
  #filtering
  if(is.numeric(population.levels)) ft <- subset(ft, l %in% population.levels)
  if(populations[1] != "all") ft <- subset(ft, p %in% populations)
  if(metrics[1] != "all") ft <- subset(ft, m %in% metrics)
  if(nrow(ft) == 0) stop("No populations left after filtering")
  
  #Dealing with NAs 
  if(NA.remove){ft <- ft[!is.na(ft$v),]} else {ft$v[is.na(ft$v)] <- 0}
  if(NA.remove) message("Removed NAs")
  #get groups
  tg <- ft %>% group_by(g) %>% summarise(nsample = length(unique(s)))
  if(nrow(tg) == 2){ #two groups
    if(sum(tg$nsample > 2) == 2){#Both groups with at least 3 samples -> Pairwise wilcoxon
      message("2 groups with 3 or more samples each. Running pairwise Wilcoxon tests...")
      resraw <- list()
      i <- 1
      for(i in 1:length(unique(ft$m))){
        currmetric <- unique(ft$m)[i]
        message(paste0("-------",currmetric))
        for(j in 1:length(unique(ft$p))){
          currpop <- unique(ft$p)[j]
          message(paste0(currpop))
          currft <- ft[ft$m == currmetric & ft$p == currpop,]
          if(length(unique(currft$g)) > 1){
            eff <- rstatix::wilcox_effsize(currft, v ~ g, 
                                           paired = paired.samples)[,-c(1,7)]
            pval <- wilcox.test(v ~ g, data = currft, 
                                paired = paired.samples)
            eff$p.val <- pval$p.value
            
            eff$m <- currmetric
            eff$p <- currpop
            effm <- currft %>% group_by(g) %>% summarise(mv = mean(v))
            eff$FC <- effm$mv[1]/effm$mv[2]
            #adding percentile to filter hits which are caused by minute changes in negative populations
            eff$percentileOfAll <- percfnc[[currmetric]](max(effm$mv))
            resraw[[length(resraw)+1]] <- eff
          }
        }
      }
      resraw <- bind_rows(resraw)
      resraw$adj.p.val <- p.adjust(resraw$p.val, method = p.adjust.method)
      
      res <- list("AOVSummaryTable" = "none", "PairwiseSummaryTable" = resraw)
      #End summary table -> joined processing later
      
    } else {#At least one group with less than 3 samples -> FCs only (same output processing as Wilcoxon)
      warning("Less then 3 samples per group, cannot calculate p values. Basing feature analysis on fold changes only...")
      
      resraw <- list()
      i <- 1
      for(i in 1:length(unique(ft$m))){
        currmetric <- unique(ft$m)[i]
        for(j in 1:length(unique(ft$p))){
          currpop <- unique(ft$p)[j]
          currft <- ft[ft$m == currmetric & ft$p == currpop,]
          if(length(unique(currft$g)) > 1){
            eff <- currft %>% group_by(g) %>% summarise(mv = mean(v), n = n())
            eff <- data.frame("group1" = eff$g[1], 
                              "group2" = eff$g[2], 
                              "effsize" = NA,
                              "FC" = eff$mv[1]/eff$mv[2],
                              #adding percentile to filter hits which are caused by minute changes in negative populations
                              "percentileOfAll" = percfnc[[currmetric]](max(eff$mv)),
                              "n1" = eff$n[1],
                              "n2" = eff$n[2],
                              "p.val" = 1,
                              "m" = currmetric,
                              "p" = currpop)
            resraw[[length(resraw)+1]] <- eff
          }
        }
      }
      resraw <- bind_rows(resraw)
      resraw$adj.p.val <- 1
      
      res <- list("AOVSummaryTable" = "none", "PairwiseSummaryTable" = resraw)
      #End summary table -> joined processing later
    }
    
  } else if(nrow(tg) > 2){ #three or more groups
    if(sum(tg$nsample > 2) > 2){#All groups with at least 3 samples -> AOV + post hoc
      message("more than 2 groups with 3 or more samples each. Running AOV and pairwise Wilcoxon tests...")
      avraw <- list()
      resraw <- list()
      i <- 1
      for(i in 1:length(unique(ft$m))){
        currmetric <- unique(ft$m)[i]
        message(paste0("-------",currmetric))
        for(j in 1:length(unique(ft$p))){
          currpop <- unique(ft$p)[j]
          message(paste0(currpop))
          currft <- ft[ft$m == currmetric & ft$p == currpop,]
          if(length(unique(currft$g)) > 1){
            #currft$v[9] <- 10
            av <- aov(v ~ g, currft)
            av <- summary(av)[[1]][-2,]
            av$m <- currmetric
            av$p <- currpop
            avraw[[length(avraw)+1]] <- av
            
            eff <- rstatix::wilcox_effsize(currft, v ~ g, 
                                           paired = paired.samples)[,-c(1,7)]
            pval <- suppressWarnings(as.vector(pairwise.wilcox.test(x = currft$v, g = currft$g, 
                                                                    paired = paired.samples, p.adjust.method = p.adjust.method.pairwise.wilcoxon)$p.value))
            
            if(length(pval[!is.na(pval)]) == nrow(eff)){
              eff$p.val <- pval[!is.na(pval)]
              eff$signif.aov <- eff$p.val < signif.level
            } else {
              eff$p.val <- NA
              eff$signif.aov <- NA
              warning("Assigned data `pval[!is.na(pval)]` must be compatible with existing data.")
            }
            
            eff$m <- currmetric
            eff$p <- currpop
            
            #add fold change 
            effm <- currft %>% group_by(g) %>% summarise(mv = mean(v))
            eff$FC <- effm$mv[match(eff$group1, effm$g)] /  effm$mv[match(eff$group2, effm$g)]
            #adding percentile to filter hits which are caused by minute changes in negative populations
            eff$percentileOfAll <- percfnc[[currmetric]](max(effm$mv))
            resraw[[length(resraw)+1]] <- eff
            
          }
        }
      }
      resraw <- bind_rows(resraw)
      resraw$adj.p.val <- p.adjust(resraw$p.val, method = p.adjust.method)
      
      res <- list("AOVSummaryTable" = bind_rows(avraw), "PairwiseSummaryTable" = resraw)
      #End summary table -> joined processing later
      
    } else {#At least one group with less than 3 samples -> FCs only (same output processing as Wilcoxon)
      warning("Less then 3 samples in at least one group, cannot calculate p values. Basing feature analysis on fold changes only...")
      
      resraw <- list()
      i <- 1
      for(i in 1:length(unique(ft$m))){
        currmetric <- unique(ft$m)[i]
        message(paste0("-------",currmetric))
        for(j in 1:length(unique(ft$p))){
          currpop <- unique(ft$p)[j]
          message(paste0(currpop))
          currft <- ft[ft$m == currmetric & ft$p == currpop,]
          effa <- as.data.frame(currft %>% group_by(g) %>% summarise(mv = mean(v), n = n()))
          if(nrow(effa) > 1){
            effp <- t(combn(effa$g,2))
            effp <- data.frame("group1" = effp[,1], "group2" = effp[,2])
            eff <- data.frame("effsize" = NA,
                              "FC" = effa$mv[match(effp$group1, effa$g)]/effa$mv[match(effp$group2, effa$g)],
                              #adding percentile to filter hits which are caused by minute changes in negative populations
                              "percentileOfAll" = percfnc[[currmetric]](max(effa$mv)),
                              "n1" = effa$n[match(effp$group1, effa$g)],
                              "n2" = effa$n[match(effp$group2, effa$g)],
                              "p.val" = 1,
                              "m" = currmetric,
                              "p" = currpop)
            eff <- cbind(effp,eff)
            resraw[[length(resraw)+1]] <- eff
          }
        }
      }
      resraw <- bind_rows(resraw)
      resraw$adj.p.val <- 1
      
      res <- list("AOVSummaryTable" = "none", "PairwiseSummaryTable" = resraw)
      #End summary table -> joined processing later
      
    }
  }
  
}


#### Boxplot wrapper ####

FlowTestBoxplots <- function(FlowTestGroups.table, #the full, or a subset, of the FlowTestGroups results list [[2]]
                             flowtidy, # the input dataframe to the FlowTestGroups call
                             plot.raw.pvalue = FALSE, #Default: FALSE, set to true to label with raw p.value
                             signif.level = 0.1,#p vals above this are not labelled
                             pdf.name = "none", #name of output PDF. if set to "none", no pdf is generated
                             png.prefix = "none", #prefix of output svgs. if set to "none" no pngs are saved
                             plot.save.width = 7,
                             plot.save.height = 7,
                             colors = "none",
                             signif.anno.height.factor = 1,#change to increase or decrease height of signif annotation
                             signif.anno.cex = 4,
                             axis.label.cex = 4,
                             subtitle.size = 18,
                             boxplot.width = 0.8,
                             shape.points.by.column = "none", #set to a column to shape points by e.g. by experiment
                             legend.position = "none",
                             no.points = FALSE, #plot only boxplots
                             label.significant.with.symbols = TRUE){ #size of signif annotations
  
  
  if(!inherits(FlowTestGroups.table, "data.frame")) stop("Please input either direct output of FlowTestGroups or the list [[2]] output from Flowtest groups. Current input is not a dataframe")
  
  df <- FlowTestGroups.table
  raw <- flowtidy
  
  #identify how many plots have to me made based on the number of groups
  #make unique identifier for metric + population and split dataframe based on that
  df <- split(df, paste0(df$m, df$p)) #paste 
  
  #plot iteration
  i <- 1
  pl.list <- list()
  for(i in 1:length(df)){
    
    curr_raw <- subset(raw, population == df[[i]]$p[1] & metric == df[[i]]$m[1])
    
    curr_raw$shape.column <- "noshaping"
    if(shape.points.by.column != "none"){
      if(shape.points.by.column %in% names(curr_raw)){
        curr_raw$shape.column <- curr_raw[,shape.points.by.column]
      } else {
        warning("Not able to shape points, because indicated column is not present in flowtidy input")
      }
    } 
    
    #deal with significance labelling
    
    df[[i]]$adj.p.val[is.na(df[[i]]$adj.p.val)] <- 1
    df[[i]]$p.val[is.na(df[[i]]$p.val)] <- 1
    
    if(!plot.raw.pvalue) df[[i]]$p.val <- df[[i]]$adj.p.val
    if(any(df[[i]]$p.val < signif.level)){ #at least one comparison was significant
      
      curr_signif <- subset(df[[i]], p.val < signif.level)
      
      #assign ranks for plotting
      curr_signif$rank <- 1:nrow(curr_signif)
      
      if(label.significant.with.symbols){
        
        signif.steps <- signif.level / c(1,10,100)
        names(signif.steps) <- c("*", "**", "***")
        
        curr_signif$p.val.label <- "ns"
        curr_signif$p.val.label[curr_signif$p.val < signif.steps[1]] <- names(signif.steps)[1]
        curr_signif$p.val.label[curr_signif$p.val < signif.steps[2]] <- names(signif.steps)[2]
        curr_signif$p.val.label[curr_signif$p.val < signif.steps[3]] <- names(signif.steps)[3]
        curr_signif$p.val <- curr_signif$p.val.label
        
        increase_label_factor <- 0.3
      } else {
        curr_signif$p.val <- round(curr_signif$p.val, 7)
        
        increase_label_factor <- 0.5
      }
      
      #add x coordinate for line by ordering groups in results the same as in raw data
      if(!is.factor(curr_raw$group)) curr_raw$group <- factor(curr_raw$group, levels = unique(raw$group))
      
      curr_signif$group1 <- factor(curr_signif$group1, levels = levels(curr_raw$group))
      curr_signif$group2 <- factor(curr_signif$group2, levels = levels(curr_raw$group))
      
      curr_signif$x <- as.numeric(curr_signif$group1)
      curr_signif$xend <- as.numeric(curr_signif$group2)
      
      curr_signif$x.text <- (curr_signif$x + curr_signif$xend) / 2
      
      #add y coordinates ->
      plot.norm.scale <- (max(curr_raw$value) - min(curr_raw$value)) / 10
      
      curr_signif$y <- max(curr_raw$value) + (curr_signif$rank * signif.anno.height.factor * plot.norm.scale)
      curr_signif$y.text <- max(curr_raw$value) + (curr_signif$rank * signif.anno.height.factor * plot.norm.scale + increase_label_factor * plot.norm.scale)
      
      #fix lower y limit:if metric is counts or percentages, lower limit should be 0
      if(stringr::str_detect(curr_signif$m[1], "Perc|perc|count|Count")){ 
        lower.y <- 0 } 
      else {
        lower.y <- min(curr_raw$value) - plot.norm.scale / 2
      }
      
      ##colors
      
      if(colors[1] == "none"){
        colors <- pals::stepped2(length(unique(curr_raw$group)))
      } 
      
      
      #plot
      pl.list[[i]] <- ggplot(curr_raw, aes(x = group, y = value, color = group, shape = shape.column)) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = colors) +
        geom_boxplot(fill = "transparent", outlier.alpha = 0, width = boxplot.width, inherit.aes = F, data = curr_raw, aes(x = group, y = value, color = group)) +
        geom_segment(inherit.aes = F, data = curr_signif, 
                     aes(x = x, xend = xend, y = y, yend = y)) +
        geom_text(inherit.aes = F, data = curr_signif, 
                  aes(x = x.text, y = y.text, label = p.val), size = signif.anno.cex) +
        scale_y_continuous(limits = c(lower.y, (max(curr_signif$y.text) + plot.norm.scale/2))) +
        theme(legend.position = legend.position, plot.subtitle=element_text(size=subtitle.size), axis.text = element_text(size = axis.label.cex), axis.title = element_text(size = axis.label.cex), axis.text.x = element_text(size = axis.label.cex, angle = 30, hjust = 1), legend.text = element_text(size = axis.label.cex), legend.title = element_text(size = axis.label.cex)) +
        labs(y = curr_signif$m[1], x = "", subtitle = unique(curr_raw$population)) 
      
      
      if(!no.points) pl.list[[i]] <- pl.list[[i]] + geom_point()
      
      if(png.prefix != "none") ggsave(pl.list[[length(pl.list)]], filename = paste0(png.prefix, curr_signif$p[1], "_", curr_signif$m[1], ".png"), width = plot.save.width, height = plot.save.height)
      
    } else { # no comparison was significant
      pl.list[[i]] <- ggplot(curr_raw, aes(x = group, y = value, color = group, shape = shape.column)) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = colors) +
        geom_boxplot(fill = "transparent", outlier.alpha = 0, width = boxplot.width, inherit.aes = F, data = curr_raw, aes(x = group, y = value, color = group)) +
        theme(legend.position = legend.position, plot.subtitle=element_text(size=subtitle.size), axis.text = element_text(size = axis.label.cex), axis.title = element_text(size = axis.label.cex), axis.text.x = element_text(size = axis.label.cex, angle = 30, hjust = 1), legend.text = element_text(size = axis.label.cex), legend.title = element_text(size = axis.label.cex))+
        labs(y = curr_raw$metric[1], x = "", subtitle = unique(curr_raw$population))
      curr_signif <- curr_raw # for svg
      
      if(!no.points) pl.list[[i]] <- pl.list[[i]] + geom_point()
      
      if(png.prefix != "none") ggsave(pl.list[[length(pl.list)]], filename = paste0(png.prefix, curr_raw$population[1], "_", curr_raw$metric[1], ".png"), width = plot.save.width, height = plot.save.height)
    }
    
  }
  if(pdf.name != "none"){
    pdf(pdf.name, width = plot.save.width, height = plot.save.height)
    print(pl.list)
    dev.off()
  } 
  return(pl.list)
}


##### Subset the cyt object ####

filterCYT <- function(cyt, cells, invert.filter = F){
  #get indices of cells
  if(invert.filter){
    cells.to.remove <- rownames(cyt@meta.data)[rownames(cyt@meta.data) %in% cells]
  } else {
    cells.to.remove <- rownames(cyt@meta.data)[!rownames(cyt@meta.data) %in% cells]
  }
  
  #go through object and remove cells

  cyt@raw.data <- cyt@raw.data[!rownames(cyt@raw.data) %in% cells.to.remove,]
  if(any(class(cyt@log.data) == "matrix"))  cyt@log.data <- cyt@log.data[!rownames(cyt@log.data) %in% cells.to.remove,]
  cyt@meta.data <- cyt@meta.data[!rownames(cyt@meta.data) %in% cells.to.remove,]
  
  if(any(class(cyt@cell.name) == "character")) cyt@cell.name <- cyt@cell.name[!cyt@cell.name %in% cells.to.remove]
  if(any(class(cyt@pca.value) == "matrix")) cyt@pca.value <- cyt@pca.value[!rownames(cyt@pca.value) %in% cells.to.remove,]
  if(any(class(cyt@umap.value) == "matrix")) cyt@umap.value <- cyt@umap.value[!rownames(cyt@umap.value) %in% cells.to.remove,]
  
  return(cyt)
  
}



