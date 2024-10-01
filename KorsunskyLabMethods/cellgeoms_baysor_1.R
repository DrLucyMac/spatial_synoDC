cellgeoms_baysor1<-function(segfile){
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
                ct_triangulate() %>%
                # purrr::reduce(st_union) %>% 
                # st_sfc() %>%
                # as.data.frame() %>%
                identity()
            resdf <- data.frame(cell = .x, geometry = res)
            return(resdf)
            }, .options = furrr_options(seed = TRUE))
       
        cellgeoms_final<-segfile.new$geometry %>% 
            furrr::future_map(purrr::reduce, st_union, .options = furrr_options(seed = TRUE)) %>%
            st_sfc() %>%
            as.data.frame()
        
        cellgeoms_final<-cellgeoms_final %>%
            cbind(transcriptspercell[transcriptspercell$cell %in% cellidx, ])
        
        # return(cellgeoms_final)
        
        saveRDS(cellgeoms_final,'cellgeoms.RDS')
        
        })
   
    
    }