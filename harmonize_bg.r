command_args<-commandArgs(trailingOnly = TRUE)

obj_k_path<-command_args[[1]]
k<-command_args[[2]]

print(k)
print(obj_k_path)
system("export OPENBLAS_NUM_THREADS=1", intern = TRUE)

source("../R/utils.R")
start_upR(TRUE)

obj_k<-readRDS(obj_k_path)
obj_k$metadata<-obj_k$metadata %>% 
    as.data.frame %>%
    dplyr::mutate(cellID = tileID)
objH<-dimred_and_cluster(
    obj_k, 
    do_harmony = TRUE, 
    wts = "wts",
    vars_use = c("SampleID", "SampleFOV"), 
    resolution_clustering = c(0.1, 0.5, 0.7, 1, 2),
    theta = c(0, 0),
    sigma = 0.2, 
    max.iter.harmony = 12,
    max.iter.cluster = 40,
    do_QC = TRUE,
    do_cluster_after = TRUE,
    do_cluster_before = TRUE,
    return_object = TRUE,
    do_umap_after = TRUE,
    do_umap_before = TRUE,
    clustering_ncores = 8

)
file_path_save<-"/n/data1/bwh/medicine/korsunsky/lab/rom4535/SpatialMapRAPaper/dataFinal/cache/tissueSegmentation/voronoiObj"
file_path_final<-file.path(file_path_save, paste0("diffused_mat_sample_harmonized_k_", k, ".RDS"))
print(file_path_final)

saveRDS(objH, file_path_final)