sourceCpp('./R/foo.cpp')
get_idx <- function(img) {
    res <- get_idx_cpp(img)
    res <- purrr::map(res, matrix, ncol=2, byrow=TRUE)
    return(res)
}


posgcmat<-function(gcmat){
    genes_pos<-rownames(gcmat) %>% 
        setdiff(grep('NegPrb', rownames(gcmat), value = TRUE))
    gcmat.pos<-gcmat[genes_pos, ]
    return(gcmat.pos)
}

findgenes<-function(genes, genelist) {
    ifelse(is_empty(grep(genes, genelist)) , 
       "Sorry, that gene is not in the panel!",  genelist[grepl(genes, genelist)])
}


cellgeoms_baysor<-function(segfile){
    system.time({ 
        plan(multicore, workers = availableCores(constraints = 'multicore') - 1)
        # delete transcripts that are noise
        segfile<-segfile %>% filter(is_noise == FALSE)
        transcriptspercell<-furrr::future_map_dfr(.x = unique(segfile$cell), 
            .f = ~ data.frame(
                cell = .x, 
                num_transcripts = sum(segfile$cell == .x)
                ), 
            .options = furrr_options(seed = TRUE)
            )
        cellidx <- transcriptspercell$cell[transcriptspercell$num_transcripts > 5]
        segfile.new <- furrr::future_map_dfr(.x = cellidx, function(.x) { 
            res <- st_as_sf(segfile[segfile$cell == .x, c('x', 'y')], coords = c('x', 'y')) %>%
                st_union() %>% #dont remove the union. It is needed here.
                ct_triangulate()
            resdf <- data.frame(cell = .x, geometry = res)
            return(resdf)
            }, .options = furrr_options(seed = TRUE))
       
        cellgeoms_final<-segfile.new$geometry %>% 
            furrr::future_map(purrr::reduce, st_union, .options = furrr_options(seed = TRUE)) %>%
            st_sfc() %>%
            as.data.frame()
        
        cellgeoms_final<-cellgeoms_final %>%
            cbind(transcriptspercell[transcriptspercell$cell %in% cellidx, ])
        
        return(cellgeoms_final)
        
        })
    
    }

collapse_and_cluster<-function(metadata, gcmat, verbose = FALSE){
    
    metadata<-metadata %>% 
        subset(tileID %in% colnames(gcmat)) %>% 
        st_sf
    
    gcmat<-gcmat[, as.character(metadata$tileID)]
    
    voronoi_coords_centroid<-metadata %>% 
        st_sf %>% 
        st_set_geometry("voronoi_centroid") %>% 
        st_coordinates("voronoi_centroid") %>% 
        as.data.frame %>% 
        dplyr::select(X, Y) %>% 
        as.matrix
    
    spatial_neighbors<-getSpatialNeighborsDist(voronoi_coords_centroid, dist_thresh_um = Inf)
    spatial_neighbors<-spatial_neighbors$adjmat_pruned
    if (verbose) {
        message("Unique diagonal elements in the adjacency matrix: ", unique(diag(spatial_neighbors)))
    }
    print(paste0("Dim gcmat: ", dim(gcmat)))
    message("Dim adjmat: ", dim(spatial_neighbors))
    collapsed_counts<-gcmat %*% (spatial_neighbors + Diagonal(nrow(spatial_neighbors)))
    message("Dim of spa. neighbors: ", dim(spatial_neighbors))
    message("Dim of cells: ", gcmat %>% ncol)
    
    colnames(collapsed_counts)<-colnames(gcmat)
    return(list(
        neighbors = spatial_neighbors, 
        counts = collapsed_counts, 
        metadata = metadata
        )
    )
    
}
    

    
getSpatialNeighborsDist<-function (coords, dist_thresh_um = 30, pxsz = 0.18, return_weights = FALSE) 
{
    dist_thresh_px <- dist_thresh_um / pxsz
    
    triplets <- geometry::delaunayn(coords)
    pairs <- spatula:::triplets_to_pairs(triplets)
    pairs <- unique(pairs)
    dists <- sqrt(rowSums((coords[pairs[, 1], ] - coords[pairs[, 
        2], ])^2))
    if (is.infinite(dist_thresh_um)) {
        idx_keep <- seq_len(nrow(pairs))
    }
    else {
        idx_keep <- which(dists < dist_thresh_px)
    }
    if (return_weights == FALSE) {
        dists <- rep(1, nrow(pairs))
    }
    adjmat_ori <- Matrix::sparseMatrix(i = pairs[, 1], j = pairs[, 
        2], x = dists, dims = c(nrow(coords), nrow(coords)))
    
    adjmat_pruned <- Matrix::sparseMatrix(i = pairs[idx_keep, 1], j = pairs[idx_keep, 
        2], x = dists[idx_keep], dims = c(nrow(coords), nrow(coords)))
    
    return(list(adjmat_ori = adjmat_ori, adjmat_pruned = adjmat_pruned))
}
