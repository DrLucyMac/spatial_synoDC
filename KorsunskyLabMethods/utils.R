dimlist<-function(list_name){
    map(list_name, ~ dim(.x))
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

start_upR<-function(clusterfiles = FALSE){
    source("./R/libs.R")
    source("./R/workflows.R")
    source("./R/utils.R")
    source("./R/VizFunctions.R")
    source("./R/SegFunctions.R")

    set.seed(9)
    theme_set(theme_bw(base_size = 15))
    plan(multicore)
if(clusterfiles){
    source("./R/ModularityClustering/R/modularity_clustering.R")
    sourceCpp("./R/ModularityClustering/src/RModularityOptimizer.cpp")

}
}

