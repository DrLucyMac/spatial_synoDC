suppressPackageStartupMessages({
   
    library(purrr)
    library(tidyr)
    library(data.table)
    library(dplyr)
    library(glue)
    library(tibble)
    library(Matrix)
    library(tidytext)
    library(RANN)
    library(igraph)
    
    
    ## plotting
    
    library(ggplot2)   
    library(ggthemes)
    library(patchwork)
    library(ggrepel)    
    library(scales)
    library(ComplexHeatmap)
    library(viridis)  
    library(tiff)
    library(colorspace)
    library(ggtext)
    library(cowplot)
    ## parallel
    
    library(future)
    library(furrr)
    
    ## spatial 
    
    library(sf)
    library(stars)
    # library(ggsn)
    # library(ggspatial)
    # library(geojsonio)
    library(classInt)
    library(spatula)
    
    library(presto)
    library(singlecellmethods)
    library(harmony)
    library(spatstat)
    library(spatstat.core)
    # library(lisi)
    
    ## glmer
    library(lme4)
    library(arm)
})

fig.size <- function(h, w) {
    options(repr.plot.width=w, repr.plot.height=h)
}