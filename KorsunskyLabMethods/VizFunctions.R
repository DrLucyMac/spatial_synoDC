markers_umap<-function (dim_red_embeddings, metadata, cell_id_colname, color_by, 
    label_by, plot_title, dim_red_type, legend_posn = "right", 
    shape_points = ".", size_points = 0.1) 
{
    plt_df <- dim_red_embeddings[, 1:2] %>% as.data.frame %>% 
        purrr::set_names("V1", "V2") %>% rownames_to_column(cell_id_colname) %>% 
        left_join(metadata %>% dplyr::select(-starts_with("V")), 
            by = cell_id_colname, multiple = "all") %>% sample_frac(1)
    range_col <- plt_df %>% with(geneval) %>% 
        range
    p <- plt_df %>% ggplot() + geom_point(aes(V1, V2, color = geneval), 
        shape = shape_points) + ggtitle(plot_title) + scale_color_gradient2_tableau(limits = c(-range_col[2], 
        range_col[2])) + xlab(paste0(dim_red_type, "_1")) + ylab(paste0(dim_red_type, 
        "_2")) + theme(plot.title = element_textbox_simple()) + 
        theme(legend.position = legend_posn) + facet_wrap(~ target) + NULL
}



plot_dim_red_cont_amp<-function (dim_red_embeddings, metadata, cell_id_colname, color_by, 
    label_by, plot_title, dim_red_type, legend_posn = "right", 
    shape_points = ".", size_points = 0.1) 
{
    plt_df <- dim_red_embeddings[, 1:2] %>% as.data.frame %>% 
        purrr::set_names("V1", "V2") %>% rownames_to_column(cell_id_colname) %>% 
        left_join(metadata %>% dplyr::select(-starts_with("V")), 
            by = cell_id_colname, multiple = "all") %>% subset(target %in% 
        color_by) %>% sample_frac(1)
    range_col <- plt_df %>% with(geneval) %>% 
        range
    p <- plt_df %>% ggplot() + geom_point(aes(V1, V2, color = geneval), 
        size = size_points) + ggtitle(plot_title) + scale_color_gradient2_tableau(limits = c(-range_col[2], 
        range_col[2])) + xlab(paste0(dim_red_type, "_1")) + ylab(paste0(dim_red_type, 
        "_2")) + theme(plot.title = element_textbox_simple()) + 
        theme(legend.position = legend_posn) + NULL
}

plot_marker_counts<-function(counts_mat, obj, markers_use, cluster_col_name){
    
    markers_exp<-counts_mat[
        markers_use, 
        colnames(counts_mat) %in% obj$metadata$cellID
    ] 
    markers_exp_df<-markers_exp %>% 
        t() %>% 
        as.matrix %>% 
        as.data.frame %>%
        rownames_to_column("cellID") %>% 
        gather(target, geneval, -cellID) %>% 
        # left_join(obj$metadata %>% dplyr::select(c(cellID, as.symbol(cluster_col_name))), by = "cellID") %>% 
        identity()
    
    fig.size(3, 4)
    map(markers_use, ~
        plot_dim_red_cont_amp(
            dim_red_embeddings = obj$Humap$embedding[obj$metadata$cellID, ], 
            metadata = markers_exp_df, 
            cell_id_colname = "cellID",
            color_by = .x,
            plot_title = .x, 
            dim_red_type = "UMAP"
        ) 
    )    
}

generate_colors_tableau<-function (input_vec, verbose = FALSE) 
{
    col_by_values <- unique(input_vec)
    if (verbose) {
        message("Number of unique values: ", length(col_by_values))
    }
    if (length(col_by_values) < 10) pal_name <- "Tableau 10" else pal_name <- "Tableau 20"
            
    if (verbose) 
        message("Palette chosen: ", pal_name)
    pal <- tableau_color_pal(pal_name)
    max_n <- attr(pal, "max_n")
    if (length(col_by_values) < max_n) {
        cols_plot <- tableau_color_pal(pal_name)(length(col_by_values))
        names(cols_plot) <- sort(col_by_values)
    }
    else {
        cols_plot <- colorRampPalette(tableau_color_pal(pal_name)(max_n))(length(col_by_values))
        names(cols_plot) <- sort(col_by_values)
    }
    return(cols_plot)
}

plot_dimred_harmony<-function(obj_input, clustercolname = NULL){
    
    
    p_before<-plot_dim_red(obj_input$pca_res$embeddings, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "sample", "PCA before Harmony colored by sample", "PCA")
    p_after<-plot_dim_red(obj_input$H, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "sample", "PCA after Harmony colored by sample", "PCA")
    
    fig.size(5, 10)
    print(p_before + theme(legend.position = "none") |
    p_after + theme(legend.position = "none"))
    
    p_before<-plot_dim_red(obj_input$pca_res$embeddings, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "fine_type", "PCA before Harmony", "PCA")
    p_after<-plot_dim_red(obj_input$H, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "fine_type", "PCA after Harmony", "PCA")
    
    fig.size(5, 15)
    print(p_before |
    p_after) 
    
    fig.size(7, 20)
    print(p_after |
        p_after  + facet_wrap(~ color_col)
    )
    
    p_before<-plot_dim_red(obj_input$umap$embedding, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "sample", "UMAP before Harmony colored by sample", "UMAP")
    p_after<-plot_dim_red(obj_input$Humap$embedding, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "sample", "UMAP after Harmony colored by sample", "UMAP")
    
    fig.size(5, 10)
    print(p_before + theme(legend.position = "none") | 
        p_after + theme(legend.position = "none")
      )
    
    p_before<-plot_dim_red(obj_input$umap$embedding, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "fine_type", "UMAP before Harmony", "UMAP")
    p_after<-plot_dim_red(obj_input$Humap$embedding, rep(1, ncol(obj_input$counts)), obj_input$metadata, "cellID", "fine_type", "UMAP after Harmony", "UMAP")
    
    fig.size(5, 15)
    print(p_before  | 
    p_after
    )
    
    fig.size(7, 20)
    print(p_after |
    p_after + facet_wrap(~ color_col)
    )
    
}

plot_genes_counts<-function (gcmat, genes_use = NULL, gene_int = 10, count_int = 15) 
{
    if (is.null(genes_use)) 
        genes_use <- rownames(gcmat)
    colSums(gcmat[genes_use, ] > 0) %>% as.data.frame %>% ggplot() + 
        geom_histogram(aes(. + 1)) + scale_x_log10() + geom_vline(xintercept = gene_int) + 
        ggtitle("Number of genes") | colSums(gcmat[genes_use, 
        ]) %>% as.data.frame %>% ggplot() + geom_histogram(aes(. + 1)) + 
        scale_x_log10() + geom_vline(xintercept = count_int) + 
        ggtitle("Number of counts")
}

plot_dim_red_cont_V2<-function (dim_red_embeddings, metadata, cell_id_colname, color_by, 
    label_by, plot_title, dim_red_type, legend_posn = "right", 
    shape_points = ".", size_points = 0.1) 
{
    plt_df <- dim_red_embeddings[, 1:2] %>%
        as.data.frame %>% 
        purrr::set_names("V1", "V2") %>% 
        rownames_to_column(cell_id_colname) %>% 
        left_join(metadata %>% 
            dplyr::select(-starts_with("V")), 
            by = cell_id_colname,
            multiple = "all") %>% 
        subset(target %in% color_by) %>% 
        sample_frac(1)
    
   
    
    range_col<-plt_df %>% subset(source %in% "cosmx") %>% with(geneval) %>%  range
    p <- plt_df %>% ggplot() + geom_point(aes(V1, V2, color = geneval), 
        size = size_points) + ggtitle(plot_title) + scale_color_gradient2_tableau(limits = c(-range_col[2], range_col[2])) + 
        xlab(paste0(dim_red_type, 
        "_1")) + ylab(paste0(dim_red_type, "_2")) + theme(plot.title = element_textbox_simple()) + 
        theme(legend.position = legend_posn) + NULL
}
plot_dim_red<-function (dim_red_embeddings, clusters = NULL, metadata, cell_id_colname = "cellID", 
    color_by, plot_title = "title", dim_red_type, legend_posn = "right", 
    shape_points = ".", size_points = 0.1, point_type = "size", 
    legend_title_input = NULL, plot_labels = FALSE, shorten_labels = FALSE) 
{
    if (is.null(clusters)) 
        clusters = rep(1, nrow(dim_red_embeddings))
    plt_df <- dim_red_embeddings[, 1:2] %>% as.data.frame %>% 
        purrr::set_names("V1", "V2") %>% rownames_to_column(cell_id_colname) %>% 
        cbind(clusters) %>% left_join(metadata, by = cell_id_colname) %>% 
        rename(color_col = as.symbol(color_by)) %>% group_by(color_col) %>% 
        mutate(x_mid = median(V1), y_mid = median(V2), label_col = color_col) %>% 
        ungroup() %>% sample_frac(1)
    cols_plot <- generate_colors_tableau(plt_df$color_col)
    if (is.null(legend_title_input)) 
        legend_title_input = "color_col"
    if (plot_labels & shorten_labels) {
        plt_df <- plt_df %>% mutate(label_col = gsub("(.*?):.*", 
            "\\1", color_col))
    }
    p <- ggplot(data = plt_df, aes(V1, V2, color = color_col)) + 
        {
            if (point_type == "size") 
                geom_point(size = size_points)
        } + {
        if (point_type == "shape") 
            geom_point(shape = shape_points)
    } + ggtitle(plot_title) + scale_color_manual(values = cols_plot) + 
        {
            if (plot_labels) 
                geom_label_repel(data = plt_df %>% dplyr::select(starts_with(c("col", 
                  "x_mid", "y_mid", "label"))) %>% distinct, 
                  aes(x = x_mid, y = y_mid, label = label_col, 
                    fill = color_col), color = "black", size = 5, 
                  min.segment.length = 0, max.overlaps = Inf) 
        } + 
    {if (plot_labels)
        scale_fill_manual(values = cols_plot)
    } + {if (plot_labels)
                  guides(fill = guide_legend(title = legend_title_input, 
                    override.aes = list(color = NA)))
    } + 
            {
        if (!plot_labels) 
            guides(color = guide_legend(title = legend_title_input, 
                override.aes = list(shape = 16, alpha = 1, size = 5)))
    } + 
        xlab(paste0(dim_red_type, "_1")) + ylab(paste0(dim_red_type, 
        "_2")) + theme(plot.title = element_textbox_simple(margin = margin(10, 
        0, 10, 0)), legend.position = legend_posn, legend.text = element_textbox_simple(width = unit(5, 
        "cm")))
}


plt_func_prop<-function(df, property_col){
    df %>% 
        ggplot(aes(x = !! sym(property_col))) + 
            geom_density() + 
            scale_x_log10() 
}

plot_dim_red_cont<-function(
    dim_red_embeddings, 
    metadata, 
    cell_id_colname, 
    color_by, 
    label_by,
    plot_title, 
    dim_red_type, 
    legend_posn = "right", 
    shape_points = ".", 
    size_points = 0.1
){
    
plt_df<-dim_red_embeddings[, 1:2] %>%
    as.data.frame %>% 
    purrr::set_names("V1", "V2") %>% 
    rownames_to_column(cell_id_colname) %>% 
    left_join(metadata %>% dplyr::select(-starts_with("V")), by = cell_id_colname) %>% 
    sample_frac(1) 
        
    
# print(head(plt_df))

    
cluster_midpoints<-plt_df %>% 
    group_by(!!as.symbol(label_by)) %>% 
    summarise(
        x_mid = median(V1),
        y_mid = median(V2)
    ) 
    
# print(head(cluster_midpoints))
    
col_by_values<-unique(cluster_midpoints[[label_by]])
cols_plot<-colorRampPalette(tableau_color_pal("Tableau 20")(20))(length(col_by_values))
names(cols_plot)<-col_by_values

    # print(col_by_values)
    
p<-plt_df %>% 
    ggplot() + 
        geom_point(aes(V1, V2, color = !!as.symbol(color_by)), size = size_points) + 
        ggtitle(plot_title) + 
        scale_color_gradient2_tableau(trans = "log10") + 
        # geom_label_repel(data = cluster_midpoints, aes(x = x_mid, y = y_mid, label = label_col, fill = label_col), size = 5, min.segment.length = Inf, max.overlaps = Inf) +
        # scale_fill_manual(values = cols_plot) +
        xlab(paste0(dim_red_type, "_1")) +
        ylab(paste0(dim_red_type, "_2")) +
        # guides(color = guide_legend(override.aes = list(shape = 16, alpha = 1, size = 5, color = NA))) + 
        # guides(fill = guide_legend(override.aes = list(color = NA))) + 
        theme(plot.title = element_textbox_simple()) +
        theme(legend.position = legend_posn) + 
        NULL
    
    
}