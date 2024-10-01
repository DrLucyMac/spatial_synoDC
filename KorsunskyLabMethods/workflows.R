subcluster_cells<-function(counts_mat, metadata, adj_mat_cluster, cluster_col, clusters_to_subtype, cluster_res){
    
    finetype<-map(clusters_to_subtype, function(cluster){
        
        # print(cluster)
        cell_ids<-metadata %>% dplyr::filter(!!rlang::sym(cluster_col) %in% cluster) %>% with(cellID)
        # print(head(cell_ids))
        colnums_cell_ids<-which(colnames(counts_mat) %in% cell_ids)
        # print(head(colnums_cell_ids))
        # adj_mat_cluster[colnums_cell_ids, colnums_cell_ids] %>% dim %>% print
    
        finetype_clusters<-data.frame(
            cellID = cell_ids,
            finetype = RunModularityClustering(adj_mat_cluster[colnums_cell_ids, colnums_cell_ids], resolution = cluster_res),
            ori_clust = cluster
        )
    })
}
        
   
merge_clusters<-function(metadata, finetype_clusters_list, ori_clust_res, finetype_cluster_merge){
    
    all_finetyped_clusters<-map(finetype_clusters_list, ~ .x$ori_clust[1]) %>% unlist
    # print(all_finetyped_clusters)
    
    map(finetype_clusters_list, function(finetype_clusters){
        original_cluster<-finetype_clusters$ori_clust[1] # there should be only one original cluster
        original_cluster_res<-rlang::sym(ori_clust_res)
        finetype_cluster_name<-rlang::sym(finetype_cluster_merge)
        # print(ori_clust_res)

        metadata<-metadata %>% 
            left_join(
                finetype_clusters %>% dplyr::select(c(cellID, ori_clust, finetype_cluster_merge)),
                by = "cellID"
            ) %>% 
            mutate(!!original_cluster_res := if_else(
                !!original_cluster_res == !!original_cluster, 
                paste(!!original_cluster_res, !!finetype_cluster_name, sep = "_"),
                !!original_cluster_res
                )
            )
    }) %>% 
        rbindlist %>% 
        # delete all the original clusters that were finetyped
        dplyr::filter(!(!!rlang::sym(ori_clust_res) %in% c(all_finetyped_clusters))) %>% 
        distinct # delete the duplicate entries of cells
    }
        

find_glmer_markers_amp<-function (obj_input, metadata_input, formula_input, cluster_col_name, 
    metadata_cols, one_tailed_input = FALSE, ncores_input = 20) 
{
    pb <- presto::collapse_counts(obj_input, metadata_input, 
        metadata_cols, min_cells_per_group = 3)
    # pb$meta_data <- pb$meta_data %>% unite(SampleFOV, SampleID:FOV, 
    #     remove = FALSE)
    system.time({
        suppressWarnings({
            presto_res <- presto.presto(formula_input, pb$meta_data, 
                pb$counts_mat, size_varname = "logUMI", effects_cov = cluster_col_name, 
                ncore = ncores_input, min_sigma = 0.05, family = "poisson", 
                nsim = 1000)
        })
    })
    contrasts_mat <- make_contrast.presto(presto_res, cluster_col_name)
    effects_marginal <- contrasts.presto(presto_res, contrasts_mat, 
        one_tailed = one_tailed_input) %>% dplyr::mutate(cluster = contrast) %>% 
        dplyr::mutate(logFC = sign(beta) * log2(exp(abs(beta))), 
            SD = log2(exp(sigma)), zscore = logFC/SD) %>% arrange(pvalue)
    effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, 
        method = "BH")
    return(effects_marginal)
}
QC_harmony_pipeline<-function(
    data_obj, 
    ngenes_threshold, 
    ncounts_threshold, 
    normval, 
    do_umap_before = TRUE,
    do_umap_after = TRUE,
    do_cluster_after = TRUE,
    resolution_clustering = 1.5, 
    clustering_ncores = 5, 
    return_object = FALSE,
    ...
    ){   
        set.seed(9)
        obj<-list()
        # param
        obj$harmony_normval<-normval
        obj$harmony_ngenes_threshold<-ngenes_threshold
        obj$harmony_ncounts_threshold<-ncounts_threshold
        # QC
        obj$counts<-QC_gcmat(data_obj$counts, ngenes_threshold, ncounts_threshold)
        obj$metadata<-data_obj$metadata %>% filter(cellID %in% colnames(obj$counts))
        # PCA
        obj$logcpx<-singlecellmethods::normalizeData((obj$counts), normval, method = 'log')
        obj$pca_res<-singlecellmethods::weighted_pca((obj$logcpx), weights = rep(1, dim(obj$counts)[2]), npc = 20, do_corr = FALSE)
        
        if (do_umap_before) {
            # UMAP 
            obj$umap<-uwot::umap(obj$pca_res$embeddings, min_dist = 0.3, spread = 1, ret_extra = 'fgraph', fast_sgd = TRUE)
        }
        
        # Harmony
        # system.time({
            if (return_object){
                Hobj<-harmony::HarmonyMatrix(
                    obj$pca_res$embeddings[, 1:20], 
                    obj$metadata, 
                    do_pca = FALSE, 
                    verbose = TRUE,
                    plot_convergence = TRUE,         
                    epsilon.cluster = -Inf,
                    epsilon.harmony = -Inf, 
                    return_object = TRUE,
                    ...
                )

               obj$H<-t(as.matrix(Hobj$Z_corr))
                # rownames and colnames should be transferred from the matrix passed to HarmonyMatrix
                rownames(obj$H)<-rownames(obj$pca_res$embeddings[, 1:20])
                colnames(obj$H)<-colnames(obj$pca_res$embeddings[, 1:20])

                obj$Hobj$kmeans_rounds<-Hobj$kmeans_rounds

                obj$Hobj$objective_kmeans<-Hobj$objective_kmeans

            }else {
                 obj$H<-harmony::HarmonyMatrix(
                    obj$pca_res$embeddings[, 1:20], 
                    obj$metadata, 
                    do_pca = FALSE, 
                    verbose = TRUE,
                    plot_convergence = TRUE,         
                    epsilon.cluster = -Inf,
                    epsilon.harmony = -Inf,   
                    ...
                )
            }
        # })
        
        if (do_umap_after) {
            set.seed(9)
            # # UMAP after Harmony
            obj$Humap<-uwot::umap(obj$H, min_dist = 0.3, spread = 1, ret_extra = 'fgraph', fast_sgd = TRUE)
        }
    
        if (do_cluster_after) {
            # # Clustering
            set.seed(9)
            diag(obj$Humap$fgraph)<-1
                obj$Humap$snn<-obj$Humap$fgraph %*% obj$Humap$fgraph
                obj$Humap$snn[which(obj$Humap$snn@x < quantile(obj$Humap$snn@x, .01))]<-0 #prune
                obj$Humap$snn<-Matrix::drop0(obj$Humap$snn)
                system.time({
                    obj$Humap$clusters<-RunModularityClustering(obj$Humap$snn, resolution = resolution_clustering, n_cores = clustering_ncores, print.output = FALSE)
                            })
        }
        return(obj)

}

QC_gcmat<-function(gcmat, gene_thresh, count_thresh){
    ngenes<-colSums(gcmat > 0)
    ncounts<-colSums(gcmat)
    gcmat<-gcmat[, ngenes >= gene_thresh & ncounts >= count_thresh]
}

coarsegrain_harmony_pipeline<-function(data_obj, ngenes_threshold, ncounts_threshold, normval, clustering_ncores = 5, return_object = FALSE,...){
        
        set.seed(9)
        obj<-list()
        # QC
        obj$counts<-QC_gcmat(data_obj$counts, ngenes_threshold, ncounts_threshold)
        obj$metadata<-data_obj$metadata %>% dplyr::filter(cellID %in% colnames(obj$counts))
        # PCA
        obj$logcpx<-singlecellmethods::normalizeData((obj$counts), normval, method = 'log')
        obj$pca<-singlecellmethods::weighted_pca((obj$logcpx), weights = rep(1, dim(obj$counts)[2]), npc = 20, do_corr = FALSE)

        # UMAP 
        obj$umap<-uwot::umap(obj$pca$embeddings, min_dist = 0.3, spread = 1, ret_extra = 'fgraph')

        # Harmony
        if (return_object){
            obj$Hobj<-harmony::HarmonyMatrix(
                obj$pca$embeddings[, 1:20], 
                obj$metadata, 
                do_pca = FALSE, 
                verbose = TRUE,
                plot_convergence = TRUE,         
                epsilon.cluster = -Inf,
                epsilon.harmony = -Inf, 
                return_object = TRUE,
                ...
            )
            
            obj$H<-as.matrix(obj$Hobj$Z_corr)
            # rownames(obj$H)<-rownames(obj$counts)
            # colnames(obj$H)<-colnames(obj$counts)
          
               
        }else {
             obj$H<-harmony::HarmonyMatrix(
                obj$pca$embeddings[, 1:20], 
                obj$metadata, 
                do_pca = FALSE, 
                verbose = TRUE,
                plot_convergence = TRUE,         
                epsilon.cluster = -Inf,
                epsilon.harmony = -Inf,   
                ...
            )
        }
        # UMAP after Harmony
        obj$Humap<-uwot::umap(obj$H, min_dist = 0.3, spread = 1, ret_extra = 'fgraph')
        # Clustering
        diag(obj$Humap$fgraph)<-1
            obj$Humap$snn<-obj$Humap$fgraph %*% obj$Humap$fgraph
            obj$Humap$snn[which(obj$Humap$snn@x < quantile(obj$Humap$snn@x, .01))]<-0 #prune
            obj$Humap$snn<-Matrix::drop0(obj$Humap$snn)
            obj$Humap$clusters<-RunModularityClustering(obj$Humap$snn, resolution = c(seq(0.1, 2.1, 0.2), 3.5, 5.5), n_cores = clustering_ncores, print.output = FALSE)

        return(obj)

}


QC_harmony_pipeline_normval<-function(
    data_obj, 
    ngenes_threshold, 
    ncounts_threshold, 
    normval_type = "median", 
    do_umap_before = TRUE, 
    do_umap_after = TRUE, 
    do_cluster_after = TRUE, 
    resolution_clustering = 1.5, 
    clustering_ncores = 5, 
    return_object = FALSE, 
    ...) 
{
    set.seed(9)
    obj <- list()
    obj$harmony_normval_type <- normval_type
    obj$harmony_ngenes_threshold <- ngenes_threshold
    obj$harmony_ncounts_threshold <- ncounts_threshold
    obj$counts <- QC_gcmat(data_obj$counts, ngenes_threshold, 
        ncounts_threshold)
    obj$metadata <- data_obj$metadata %>% filter(cellID %in% 
        colnames(obj$counts))
    normval<-median(colSums(obj$counts))
    obj$harmony_normval<-normval
    obj$logcpx <- singlecellmethods::normalizeData((obj$counts), 
        normval, method = "log")
    obj$pca_res <- singlecellmethods::weighted_pca((obj$logcpx), 
        weights = rep(1, dim(obj$counts)[2]), npc = 20, do_corr = FALSE)
    if (do_umap_before) {
        obj$umap <- uwot::umap(obj$pca_res$embeddings, min_dist = 0.3, 
            spread = 1, ret_extra = "fgraph", fast_sgd = TRUE)
    }
    system.time({
        if (return_object) {
            Hobj <- harmony::HarmonyMatrix(obj$pca_res$embeddings[, 
                1:20], obj$metadata, do_pca = FALSE, verbose = TRUE, 
                plot_convergence = TRUE, epsilon.cluster = -Inf, 
                epsilon.harmony = -Inf, return_object = TRUE, 
                ...)
            obj$H <- t(as.matrix(Hobj$Z_corr))
            rownames(obj$H) <- rownames(obj$pca_res$embeddings[, 
                1:20])
            colnames(obj$H) <- colnames(obj$pca_res$embeddings[, 
                1:20])
            obj$Hobj$kmeans_rounds <- Hobj$kmeans_rounds
            obj$Hobj$objective_kmeans <- Hobj$objective_kmeans
        }
        else {
            obj$H <- harmony::HarmonyMatrix(obj$pca_res$embeddings[, 
                1:20], obj$metadata, do_pca = FALSE, verbose = TRUE, 
                plot_convergence = TRUE, epsilon.cluster = -Inf, 
                epsilon.harmony = -Inf, ...)
        }
    })
    if (do_umap_after) {
        print("UMAP after Harmony")
        set.seed(9)
        obj$Humap <- uwot::umap(obj$H, min_dist = 0.3, spread = 1, 
            ret_extra = "fgraph", fast_sgd = TRUE)
    }
    if (do_cluster_after) {
        print("Building SNN")
        set.seed(9)
        diag(obj$Humap$fgraph) <- 1
        obj$Humap$snn <- obj$Humap$fgraph %*% obj$Humap$fgraph
        obj$Humap$snn[which(obj$Humap$snn@x < quantile(obj$Humap$snn@x, 
            0.01))] <- 0
        obj$Humap$snn <- Matrix::drop0(obj$Humap$snn)
        print("Clustering")
        system.time({
            obj$Humap$clusters <- RunModularityClustering(obj$Humap$snn, 
                resolution = resolution_clustering, n_cores = clustering_ncores, 
                print.output = FALSE)
        })
    }
    return(obj)
}


find_glmer_markers<-function(obj_input, metadata_input, formula_input, cluster_col_name, metadata_cols, one_tailed_input = FALSE, ncores_input = 20){
    
    pb<-presto::collapse_counts(obj_input, metadata_input, metadata_cols, min_cells_per_group = 3)

    pb$meta_data<-pb$meta_data %>% unite(SampleFOV, SampleID:FOV, remove = FALSE)
    system.time({
        suppressWarnings({
            presto_res <- presto.presto(
                formula_input,
                # y ~ 1 + (1|Clust5.5) + (1|SampleFOV/Clust5.5) + (1|SampleID/Clust5.5) + offset(logUMI), 
                pb$meta_data, 
                pb$counts_mat,
                size_varname = 'logUMI', 
                effects_cov = cluster_col_name,
                ncore = ncores_input, 
                min_sigma = .05,
                family = 'poisson',
                nsim = 1000
            )    
        })
    })

    contrasts_mat <- make_contrast.presto(presto_res, cluster_col_name)
    effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = one_tailed_input) %>% 
        dplyr::mutate(cluster = contrast) %>% 
        dplyr::mutate(
            ## convert stats to log2 for interpretability 
            logFC = sign(beta) * log2(exp(abs(beta))),
            SD = log2(exp(sigma)),
            zscore = logFC / SD
        ) %>% 
        # dplyr::dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) %>% 
        arrange(pvalue)

    effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = 'BH')
    
    return(effects_marginal)

}

createFormula <- function(resp, fixed, rand) {
f <- reformulate(c(fixed,rand),response=resp)
## use the parent (createModel) environment, not the
## environment of this function (which does not contain 'data')
environment(f) <- parent.frame()
f
}

function (clusterobj, do_harmony = FALSE, ...) 
{
  input_var <- list(...)
  if (is.null(input_var[["sigma"]])) 
    message("You did not specify sigma for Harmony")
  else sigma_harmony <- input_var[["sigma"]]
  if (is.null(input_var[["theta"]])) 
    message("You did not specify theta for Harmony")
  else theta_harmony <- input_var[["theta"]]
  if (is.null(input_var[["vars_use"]])) 
    message("You did not specify batch var for Harmony")
  else batch_harmony <- input_var[["vars_use"]]
  wts <- clusterobj$metadata %>% with(table(as.factor(region))) %>% 
    prop.table %>% as.data.frame %>% rename(region = Var1) %>% 
    mutate(wts = 1/Freq) %>% mutate(region = as.character(region))
  clusterobj$metadata <- clusterobj$metadata %>% left_join(wts %>% 
                                                             dplyr::select(region, wts), by = "region")
  head(clusterobj$metadata)
  if (do_harmony) {
    norm_value <- median(colSums(clusterobj$counts))
    objH <- QC_harmony_pipeline(clusterobj, ngenes_threshold = 1, 
                                ncounts_threshold = 1, norm_value, ...)
    objH$sigma_harmony <- sigma_harmony
    objH$vars_use <- batch_harmony
    objH$theta_harmony <- theta_harmony
    clusterobj <- objH
  }
  else {
    set.seed(9)
    print("Clustering on SNN")
    system.time({
      clusterobj$umap$clusters <- RunModularityClustering(clusterobj$umap$fgraph, 
                                                          resolution = c(0.2, 0.5, 0.7, 1), n_cores = 10, 
                                                          print.output = FALSE)
    })
  }
  return(clusterobj)
}
