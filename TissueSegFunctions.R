grid_glass<-function(glass_region_input, txfile_input, cellgeoms_input, by_counts = TRUE, ncounts = 100){ 

    if (!by_counts){
        median_voronoi_size<-cellgeoms_input %>% 
            with(sqrt(median(polygon_area)))  
        glass_gridded<-st_make_grid(
            st_rectangle(min(txfile_input$x), max(txfile_input$x), min(txfile_input$y), max(txfile_input$y)), 
            cellsize = median_voronoi_size, 
            square = TRUE
        ) %>% 
        st_intersection(st_sf(glass_region_input)) %>% 
        st_cast("MULTIPOLYGON") %>% 
        st_cast("POLYGON") %>% 
        dplyr::mutate(
            tileID = seq(1, nrow(.)),
            tileID = paste(paste0("G", tileID), txfile_input$SampleFOV[1], sep = "_")
        )
    } else { #check if there are acceptable #tx in the bg. If there are then proceed with tiling  
        txfile_glass<-txfile_input %>% st_sf %>% 
            st_join(st_sf(glass_region_input), st_intersects) %>% 
            subset(region %in% "glass")
        # if there are no tx in the glass region, don't tile it
        if (nrow(txfile_glass) < ncounts/2) {
            glass_gridded<-NULL
            } else {
            system.time({glass_gridded<-spatula::split_tx(txfile_glass, ncounts, 1e3)})
            glass_gridded<-glass_gridded %>% 
                purrr::map(~ .x$bbox_geom %>% st_sfc %>% as.data.frame) %>%
                rbindlist    
            }     
    }   
    if (!is.null(glass_gridded)){
        glass_gridded<-glass_gridded %>%
            st_sf %>% 
            st_intersection(st_sf(glass_region_input)) %>% 
            st_cast("MULTIPOLYGON") %>% 
            st_cast("POLYGON") %>% 
            dplyr::mutate(
                tileID = seq(1, nrow(.)),
                tileID = paste(paste0("G", tileID), txfile_input$SampleFOV[1], sep = "_")
            )
    }
    return(glass_gridded)

}

pts_to_voronoi = function(x, y, eps=0, bbox=NULL) {
    if (is.null(bbox)) {
        bbox = st_rectangle(min(x) - eps, max(x) + eps, min(y) - eps, max(y) + eps)        
    }
    tiles = cbind(x, y) %>% st_multipoint() %>% st_voronoi(bbox) %>% st_collection_extract() 
    pts = st_as_sf(data.frame(x, y), coords = c('x', 'y'))
    tiles = tiles[unlist(st_intersects(pts, tiles))] ## put them in order: THIS MUST BE OPTIMIZED AT SOME POINT 
    tiles = st_crop(tiles, bbox)
    return(tiles)
}

getSpatialNeighbors_obtuse_fix<-function(coords, fix_obtuse = TRUE, 
    dist_thresh = Inf, return_weights = FALSE) 
{
    triplets <- geometry::delaunayn(coords)
    if (fix_obtuse){
        triplet_idx <- NULL
        angles_all<-NULL
        lengths_all<-NULL
        for (i in 1:nrow(triplets)) {
            triangle <- triplets[i, ]
            vertices_triangle <- coords[triangle, ]
            a <- vertices_triangle[1, ]
            b <- vertices_triangle[2, ]
            c <- vertices_triangle[3, ]
            ab <- b - a
            bc <- c - b
            ca <- a - c
            length_edge1<-sqrt(sum(ab^2))
            length_edge2<-sqrt(sum(bc^2))
            length_edge3<-sqrt(sum(ca^2))
            dot_product <- sum(-ab * bc)
            norm_product <- sqrt(sum(ab^2)) * sqrt(sum(bc^2)) # not doing -ab because this is squared anyway
            angle12 <- acos(pmax(pmin(dot_product/norm_product, 1.0), -1.0))
            angle12 <- angle12 * (180/pi)
            dot_product <- sum(ca * -bc)
            norm_product <- sqrt(sum(ca^2)) * sqrt(sum(bc^2))
            angle23 <- acos(pmax(pmin(dot_product/norm_product, 1.0), -1.0))
            angle23 <- angle23 * (180/pi)
            dot_product <- sum(-ca * ab)
            norm_product <- sqrt(sum(ca^2)) * sqrt(sum(ab^2))
            angle31 <- acos(pmax(pmin(dot_product/norm_product, 1.0), -1.0)) # making sure values are bounded b/w [-1, 1]
            angle31 <- angle31 * (180/pi)
            angles_all <- rbind(angles_all, c(angle12, angle23, angle31, triangle)) # contains angles and the point idx of triangles
            lengths_all <- rbind(lengths_all, c(length_edge1, length_edge2, length_edge3, triangle)) # contains lengths and the point idx of triangles

            if (any(angles_all > 150)) {
                triplet_idx <- c(triplet_idx, i)
            }
        }
        message("Before fixing for angles: ", print(dim(triplets)))
        triplets <- triplets[-triplet_idx, ]
        message("After fixing for angles: ", print(dim(triplets)))
    } else {
        lengths_all<-NULL
        angles_all<-NULL
    }
    pairs <- spatula:::triplets_to_pairs(triplets)
    pairs <- unique(pairs)
    dists <- sqrt(rowSums((coords[pairs[, 1], ] - coords[pairs[, 2], ])^2))
    # apply the distance threshold 
    if (is.infinite(dist_thresh)) {
        idx_keep <- seq_len(nrow(pairs))
    }
    else {
        idx_keep <- which(dists < dist_thresh)
    }
    if (return_weights == FALSE) {
        dists <- rep(1, nrow(pairs))
    }
    adjmat <- Matrix::sparseMatrix(i = pairs[idx_keep, 
        1], j = pairs[idx_keep, 2], x = dists[idx_keep], dims = c(nrow(coords), 
        nrow(coords)))
    
    return(list(adjmat = adjmat, lengths = lengths_all, angles = angles_all))
}
diffuse_adjmat<-function (metadata, gcmat, lambda = 0.1, k = 1, pxsize = 1, 
    verbose = TRUE, cut_connections = FALSE, prune_connections = TRUE, 
    prune_communities = TRUE, dist_threshold_cells = 50, community_size_thresh = 4) 
{
    message("Diffusing and cleaning up")
    metadata<-metadata %>% subset(tileID %in% colnames(gcmat)) %>% 
        st_sf
    gcmat<-gcmat[, as.character(metadata$tileID)]
    voronoi_coords_centroid <- metadata %>% st_sf %>% st_set_geometry("polygon_centroid") %>% 
        st_coordinates("polygon_centroid") %>% as.data.frame %>% 
        dplyr::select(X, Y) %>% as.matrix
    adjoutput<-getSpatialNeighbors_obtuse_fix(voronoi_coords_centroid, dist_thresh = Inf, return_weights = TRUE, fix_obtuse = FALSE)
    spatial_neighbors<-adjoutput$adjmat
    if (cut_connections) {
        cells.bg<-metadata %>% subset(region %in% "glass") %>% 
            with(tileID)
        if (length(cells.bg) > 0) {
            cellids.bg <- which(colnames(gcmat) %in% cells.bg)
            cellids.fg <- which(!colnames(gcmat) %in% cells.bg)
            spatial_neighbors[cellids.fg, cellids.bg] = 0
            spatial_neighbors[cellids.bg, cellids.fg] = 0
        }
    }
    # within the tissue region
    if (prune_connections) {
        cells.fg<-metadata %>% subset(region %in% "tissue") %>% 
            with(tileID)
        cellids.fg<-which(colnames(gcmat) %in% cells.fg)
        adjmat.fg<-spatial_neighbors[cellids.fg, cellids.fg]
        cellids.far<-which(adjmat.fg > dist_threshold_cells)
        adjmat.fg[cellids.far] = 0
        spatial_neighbors[cellids.fg, cellids.fg]<-adjmat.fg   
    }
    message("here")
    # adjoutput$neighbors_pruned<-spatial_neighbors
    if (prune_communities) {
        cells.fg<-metadata %>% subset(region %in% "tissue") %>% 
            with(tileID)
        cellids.fg<-which(colnames(gcmat) %in% cells.fg)
        cells.fg<-colnames(gcmat)[cellids.fg]
        adjmat.fg<-spatial_neighbors[cellids.fg, cellids.fg]
        graph_adj1<-igraph::graph_from_adjacency_matrix(adjmat.fg, mode = "undirected")
        comp_graph<-components(graph_adj1)
        comp.delete<-which(comp_graph$csize < community_size_thresh)
        cellids.delete<-which(comp_graph$membership %in% comp.delete)
        if (!is.empty(cellids.delete)){
            cells.delete<-cells.fg[cellids.delete]
            adjmat.fg[cellids.delete, ] = 0
            adjmat.fg[, cellids.delete] = 0
            spatial_neighbors[cellids.fg, cellids.fg]<-adjmat.fg

            cellids.delete.global<-which(colnames(gcmat) %in% cells.delete)
            spatial_neighbors<-spatial_neighbors[-cellids.delete.global, -cellids.delete.global]
            # print(cells.delete)
            # print(dim(gcmat))
            gcmat<-gcmat[, !colnames(gcmat) %in% cells.delete]
            metadata_small<-metadata %>% 
                subset(!tileID %in% colnames(gcmat))
            metadata<-metadata %>% 
                subset(tileID %in% colnames(gcmat))
            gcmat<-gcmat[, as.character(metadata$tileID)]
        } else {metadata_small<-NULL}
        
        # print(dim(gcmat))
    } else {metadata_small<-NULL}
    

    # spatial_neighbors <- spatial_neighbors %>% as("dgCMatrix")
    if (verbose) {
    }
    if (verbose) {
        message("  Unique diagonal elements in the adjacency matrix: ", 
            unique(diag(spatial_neighbors)))
    }
    A <- spatial_neighbors 
    A[A > 0]<-1 # this is done so that A remains dgC and not lgC
    A <- Diagonal(n = nrow(A)) + (lambda * A)
    # small_niches <- which(rowSums(A > 0) <= 4)
    # A[small_niches, ] = 0
    # A[, small_niches] = 0
    A <- Diagonal(x = 1/rowSums(A)) %*% A
    M <- A
    Mlist <- list()
    Mlist[[1]] <- A
    
    for (i in seq_len(k - 1)) {
        message("Sum of adj mat before: ", sum(M))
        M <- A %*% M
        Mlist[[i + 1]] <- M
        message("Sum of adj mat before: ", sum(M))
    }
    message(nrow(Mlist[[1]]))
    message(ncol(Mlist[[1]]))
    message(nrow(gcmat))
    message(ncol(gcmat))
    collapsed_counts <- gcmat %*% t(Mlist[[k - 1]])
    colnames(collapsed_counts) <- colnames(gcmat)
    return(list(neighbors = spatial_neighbors, neighbors_collapsed = Mlist, 
        k_input = k, counts = collapsed_counts, metadata = metadata, metadata_small = metadata_small, counts_raw = gcmat))
}
alpha_shape_cells<-function(cells_voronoi){
    # this functions takes the OP of Voronoi tessallation and alpha shapes the cells the edge cells
    message(" alpha shape cells")
    cells_voronoi<-cells_voronoi %>% 
        dplyr::mutate(
            polygon_centroid = st_centroid(polygon),
            polygon_area = st_area(polygon),
            cell_voronoi_dist = st_distance(cell_centroid, polygon_centroid, by_element = TRUE)
        ) %>% 
        dplyr::select(tileID, cellID, everything()) %>% st_sf()
    # create a 15µm buffer from the cell centroid of each cell in the whole tissue
    buffer_cells<-cells_voronoi %>% st_set_geometry("cell_centroid") %>% 
        dplyr::select(cell_centroid) %>% 
        st_buffer(15) %>% 
        st_union  
    # every cell is an intersection b/w itself and the buffer - this way, the "edge cells" will become a circular shape 
    cells_voronoi_final<-cells_voronoi %>% st_set_geometry("polygon") %>% 
        st_intersection(buffer_cells)
    # find the "weird cells" - because we want to change the shapes of edge cells and the 
        ##cells whose cell centroid and voronoi centroid are > 7µm apart
    boundary_tissue<-cells_voronoi_final %>%   
        st_set_geometry("polygon") %>% 
        st_union() %>% 
        st_boundary %>% 
        as.data.frame %>% 
        dplyr::mutate(boundary = TRUE) %>% st_sf
    # correct the shapes of "edge cells" and the cells whose cell centroid and voronoi centroid are > 7µm apart
    cells_edge_cases<-st_join(
        cells_voronoi_final, 
        boundary_tissue,
        join = st_touches
    ) %>% 
        subset(boundary == TRUE) %>% 
        bind_rows(cells_voronoi_final %>% subset(cell_voronoi_dist > 7)) %>% # also change cells whose vornoi shapes are very different from the cell shapes    
        dplyr::select(-boundary) %>% 
        distinct
    # Voronoi shapes are frozen. Don't change shapes beyond this.
    cells_voronoi_final<-cells_voronoi_final %>% 
        subset(!cellID %in% cells_edge_cases$cellID) %>% 
        bind_rows(cells_edge_cases) 
    # re-calculate properties of the new Voronoi polygons
    cells_voronoi_final<-cells_voronoi_final %>%
        dplyr::mutate(
            polygon_centroid = st_centroid(polygon),
            polygon_area = st_area(polygon),
            cell_voronoi_dist = st_distance(cell_centroid, polygon_centroid, by_element = TRUE)
        ) %>% 
        dplyr::select(tileID, cellID, everything()) 
    return(cells_voronoi_final)
    }
    
tessallate_and_build_gcmat<-function(tfile, cellgeoms_input){
    # this function tessallates the tissue, alpha shapes cells, and builds a gene-cell matrix
    message("Voronoi tessallation")
    bbox_tx<-st_rectangle(min(tfile$x), max(tfile$x), min(tfile$y), max(tfile$y))
    message(" tessallate")
    cells_voronoi<-pts_to_voronoi(cellgeoms_input$cell_center_x, cellgeoms_input$cell_center_y, bbox = bbox_tx)
    # let us add cell ids to the Voronoi polygons so any changing of the rows of this dataframe is resistant to cellIDs
    cells_voronoi<-cells_voronoi %>% 
        as.data.frame %>% 
        cbind(cellgeoms_input %>% dplyr::select(cellID, cell_centroid, SampleFOV, SampleID)) %>% 
        dplyr::rename(polygon = geometry) %>% 
        dplyr::mutate(tileID = cellID)
    # alpha shape cells
    cells_voronoi_alpha_shaped<-alpha_shape_cells(cells_voronoi)
    # generate gcmat
    tfile_voronoi<-tfile[cell != 0] %>% 
        st_sf %>% 
        st_join(
            cells_voronoi_alpha_shaped %>% st_sf %>% 
                st_set_geometry("polygon") %>% 
                dplyr::select(cellID, tileID, polygon),
            .predicate = st_intersects
        ) 
    # remove tx not belonging to any cells
    # message(" # tx not belonging to any voronoi polygon: ", data.table(tfile_voronoi_tissue)[is.na(tileID), .N])
    tfile_voronoi<-data.table(tfile_voronoi)[!is.na(tileID)]    
    gcmat_voronoi<-spatula::tx_to_counts(
        tfile_voronoi$gene, 
        tfile_voronoi$cellID, 
        remove_bg = TRUE
    )
    voronoi_obj<-list()
    voronoi_obj$metadata<-cells_voronoi_alpha_shaped  %>% as.data.frame
    voronoi_obj$counts<-gcmat_voronoi[, as.character(voronoi_obj$metadata$cellID)]
    return(voronoi_obj)
    message("done")
}

check_glass<-function(obj, polygon_col_name, buffer_size = 30, tfile, bbox = NULL){
    # this function checks for the presence of glass beyond the buffered tissue regions on the slide
   
    tissue_region<-obj$metadata %>% 
        st_sf %>% 
        st_set_geometry(polygon_col_name) %>% 
        st_union %>% 
        st_cast("MULTIPOLYGON") %>% 
        st_buffer(buffer_size)  
   
    if (is.null(bbox))
        fov_region<-st_bbox(c(xmin = min(tfile$x), xmax = max(tfile$x), ymin = min(tfile$y), ymax = max(tfile$y))) %>% 
            st_as_sfc 
    else fov_region<-st_bbox(c(xmin = bbox[1], xmax = bbox[2], ymin = bbox[3], ymax = bbox[4])) %>% st_as_sfc
    
    glass_region<-st_difference(fov_region, tissue_region) %>% 
        as.data.frame %>% 
        dplyr::mutate(region = "glass") 

    res<-!(is.empty(glass_region$geometry))
    obj<-list()
    obj$tissue_region<-tissue_region
    obj$glass_region<-glass_region
    obj$glass_exists<-res

    return(obj)
    
}