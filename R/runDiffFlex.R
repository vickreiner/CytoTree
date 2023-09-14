#'
#' Calculate differential expression markers
#'
#' @name runDiffFlex
#'
#' @description Calculating differentially expressed markers for each element of grouping vector
#'
#' @param object a CYT object
#' @param grouping.column Column name of meta.data. Defaults to cluster.id For markers specific to branch check runDiff. 
#' @param group1 one or more elements of grouping.column to use as group1 
#' @param group2 one element of grouping.column to use as group2 (only provide if group1 is provided) 
#' @param verbose logic. Whether to print calculation progress.
#'
#' @seealso  \code{runDiff}
#'
#' @return A table of DEM (differentially expressed markers)
#'
#' @import limma
#' @importFrom stringr str_replace_all fixed
#' @importFrom stats model.matrix
#'
#' @export
#' @return a data.frame with differential expressed markers
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#'
#' DEG.table <- runDiffFlex(cyt)
#'
#' 
#'
#'
#'
runDiffFlex <- function(object, grouping.column = "cluster.id", group1 = NULL, group2 = NULL, verbose = FALSE, p.adjust.method = "BH", pairwise = FALSE, downsample = FALSE,...) {

  if (verbose) message(Sys.time(), " Calculating differentially expressed markers.")
  if (missing(object)) stop(Sys.time(), " CYT object is missing.")
  if (!grouping.column %in% colnames(object@meta.data)) stop(Sys.time(), "input to grouping.column is missing in meta.data. Please check spelling")

  all.branch.ids <- unique(object@meta.data[,grouping.column])
  
  total.deg.list <- NULL
  branch.contrast <- NULL
  ga <- go <- NULL
  if (length(all.branch.ids) == 1) {
    stop(Sys.time(), " There is only one group.")
  } else {
    
    if(downsample){
      pdata <- object@meta.data[which(object@meta.data$dowsample == 1), c("cell", grouping.column)]
      names(pdata)[2] <- "branch.id"
      edata <- object@log.data[which(object@meta.data$dowsample == 1), object@markers.idx]
    } else {
      pdata <- object@meta.data[, c("cell", grouping.column)]
      names(pdata)[2] <- "branch.id"
      edata <- object@log.data[, object@markers.idx]
    }
    
    if (is.null(group1) & is.null(group2)) {
      if(pairwise == FALSE){
          for (bid in all.branch.ids) {
            pdata$contrast <- "go"
            pdata$contrast[which(pdata$branch.id == bid)] = "ga"
            design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
            colnames(design) <- stringr::str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
            fit <- lmFit(t(edata), design)
            contrast <- makeContrasts(ga_go = ga - go,
                                     levels = design)
            fits <- contrasts.fit(fit, contrast)
            ebFit <- eBayes(fits)
            
            deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = p.adjust.method, number = Inf)
            deg_sig_list$branch.contrast <- paste0(bid, "_vs_other")
            deg_sig_list$Gene <- rownames(deg_sig_list)
            total.deg.list <- rbind(total.deg.list, deg_sig_list)
        }
      } else { #PAIRWISE COMPARISON
        for (bid in all.branch.ids) {
          for(bid2 in all.branch.ids){
            if(bid != bid2){
              pdata$contrast <- "gz"
              pdata$contrast[pdata$branch.id %in% bid] = "ga"
              pdata$contrast[pdata$branch.id %in% bid2] = "go"
              design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
              colnames(design) <- stringr::str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
              fit <- lmFit(t(edata), design)
              contrast <- makeContrasts(ga_go = ga - go,
                                      levels = design)
              fits <- contrasts.fit(fit, contrast)
              ebFit <- eBayes(fits)
            
              deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = p.adjust.method, number = Inf)
              deg_sig_list$branch.contrast <- paste0(bid, "_vs_", bid2)
              deg_sig_list$Gene <- rownames(deg_sig_list)
              total.deg.list <- rbind(total.deg.list, deg_sig_list)
            }
          }
        }
      }
    } else if (is.null(group2)) {
      if (verbose) message(Sys.time(), " Some of branches will be calculated.")
      pdata$contrast <- "go"
      pdata$contrast[pdata$branch.id %in% group1] = "ga"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = p.adjust.method, number = Inf)
      deg_sig_list$branch.contrast <- paste0(paste0(group1, collapse = "-"), "_vs_other")
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    } else {
      if (verbose) message(Sys.time(), " Some of branches will be calculated.")
      pdata$contrast <- "gz"
      pdata$contrast[pdata$branch.id %in% group1] = "ga"
      pdata$contrast[pdata$branch.id %in% group2] = "go"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = p.adjust.method, number = Inf)
      deg_sig_list$branch.contrast <- paste0(paste0(group1, collapse = "-"), "_vs_", paste0(group2, collapse = "-"))
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    }
  if (verbose) message(Sys.time(), " Calculating differentially expressed markers completed")
  return(total.deg.list)
  }
}
