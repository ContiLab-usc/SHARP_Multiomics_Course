

plot_lucid_in_parallel_plotly<- function(lucidus_fit,
                                         sankey_colors,
                                         text_size = 10, 
                                         n_z_ftrs_to_plot = NULL){
  rownames(lucidus_fit$res_Mu[[1]]) = lucidus_fit$var.names$Znames[[1]]
  rownames(lucidus_fit$res_Mu[[2]]) = lucidus_fit$var.names$Znames[[2]]
  rownames(lucidus_fit$res_Mu[[3]]) = lucidus_fit$var.names$Znames[[3]]
  # Get number of clusters, layers, etc.
  K <- lucidus_fit$K
  dimG <- lucidus_fit$res_Beta$Beta[[1]] %>% ncol()-1
  n_layers   <- length(lucidus_fit$res_Beta$Beta)
  
  # Get top omics features based on effect size
  if(!is.null(n_z_ftrs_to_plot)){
    top_ftrs <- vector("list", n_layers)
    for(i in seq_along(top_ftrs)){
      top_ftrs[[i]] <- lucidus_fit$res_Mu[[i]] %>%
        as_tibble(rownames = "name") %>%
        # rownames_to_column("name") %>%
        mutate(effect_size = abs(V1) + abs(V2)) %>%
        arrange(desc(effect_size))
      top_ftr_nms <- top_ftrs[[i]]$name[1:n_z_ftrs_to_plot[i]]
      lucidus_fit$res_Mu[[i]] <- 
        lucidus_fit$res_Mu[[i]][
          rownames(lucidus_fit$res_Mu[[i]]) %in% top_ftr_nms, ]
    }
  }
  
  mu_lst <- purrr::map(lucidus_fit$res_Mu, 
                       ~as_tibble(.x, rownames = "name"))
  names(mu_lst) <- paste0("layer", c(1:n_layers))
  dimZ <- purrr::map(mu_lst, ncol) %>% as.numeric()-1
  n_features <- purrr::map(mu_lst, nrow) %>% as.numeric()
  names(n_features) <- paste0("layer", c(1:n_layers))
  # Names of features and set order of omics features
  names_features <- bind_rows(mu_lst, .id = "color_group") %>% 
    rowwise() %>%
    mutate(sum = sum(abs(V1)+abs(V2)), 
           pos_c2 = if_else(V2>0, "pos", "neg")) %>%
    group_by(color_group, pos_c2) %>% arrange(-sum, .by_group = TRUE) %>% ungroup() %>% 
    mutate(rnum = row_number()) %>%
    group_by(name) %>% slice_head() %>% ungroup() %>%
    arrange(color_group, rnum) %>%
    dplyr::select(name, color_group)
  
  # Values for g --> x association
  valueGtoX <- c(lapply(lucidus_fit$res_Beta$Beta, 
                        function(x)(x[-1])) %>%
                   unlist(), 
                 rep(0, dimG*n_layers))
  
  # For Cluster 2 (which needs effect estimates): 
  valueGtoX_c1 <- do.call(rbind, lucidus_fit$res_Beta$Beta)[,-1] %>%
    as_tibble() %>%
    dplyr::mutate(layer = str_c("(Layer ", row_number(), ")"),
                  cluster = "Cluster 2") 
  
  # For cluster 1 (ref. cluster, effect est = 0):
  valueGtoX_c2 <- valueGtoX_c1 %>%
    mutate(across(where(is.numeric), ~0), 
           cluster = "Cluster 1")
  
  # combine, pivot longer, and create source and target columns
  GtoX <- bind_rows(valueGtoX_c1, valueGtoX_c2) %>%
    mutate(target = str_c(cluster, layer, sep = " ")) %>%
    pivot_longer(cols = setdiff(colnames(valueGtoX_c1), 
                                c("layer", "cluster")), 
                 names_to = "source", values_to = "value") %>%
    mutate(color_group = as.factor(value > 0), 
           value = abs(value)) %>%
    dplyr::select(source, target, value, color_group) %>%
    as.data.frame()
  
  valueXtoZ <- c(lapply(lucidus_fit$res_Mu, 
                        function(x)x[, 1]) %>% 
                   unlist(), 
                 lapply(lucidus_fit$res_Mu, 
                        function(x)x[, 2]) %>% 
                   unlist())
  
  valueXtoY <- c(rep(0, n_layers), 
                 # rep(lucidus_fit$res_Delta$Delta$mu[1] / n_layers, n_layers),
                 lucidus_fit$res_Gamma$fit$coefficients[c(2,3,4)])
  
  # n features in each layer
  XtoZ <- data.frame(source = c(rep("Cluster 1 (Layer 1)", n_features[1]),
                                rep("Cluster 1 (Layer 2)", n_features[2]),
                                rep("Cluster 1 (Layer 3)", n_features[3]),
                                # rep("Cluster 1 (Layer 4)", n_features[4]),
                                rep("Cluster 2 (Layer 1)", n_features[1]),
                                rep("Cluster 2 (Layer 2)", n_features[2]),
                                rep("Cluster 2 (Layer 3)", n_features[3]) 
                                # rep("Cluster 2 (Layer 4)", n_features[4])
  ), 
  target = rep(c(lapply(lucidus_fit$res_Mu,
                        rownames) %>% unlist()),
               K[1]), 
  value = abs(valueXtoZ), 
  color_group = as.factor(valueXtoZ > 0))
  
  # To change the outcome from left to right hand side, flip source and target
  XtoY <- data.frame(target = rep("outcome", 2*n_layers), 
                     source = c("Cluster 1 (Layer 1)", 
                                "Cluster 1 (Layer 2)",
                                "Cluster 1 (Layer 3)",
                                #"Cluster 1 (Layer 4)",
                                "Cluster 2 (Layer 1)",
                                "Cluster 2 (Layer 2)",
                                "Cluster 2 (Layer 3)"
                                #"Cluster 2 (Layer 4)"

                     ), 
                     value = abs(valueXtoY), 
                     color_group = as.factor(valueXtoY > 0))
  
  # create Sankey diagram
  # Create Links ----
  links <- rbind(GtoX, XtoZ, XtoY) %>%
    mutate(
      # Group: one of exposure, clusters, or outcomes 
      # (doesn't include Z.order by desired order)
      source_group = case_when(
        str_detect(source, "Cluster") ~ "2_Cluster", 
        # source == "Outcome" ~ "3_outcome", # removed when moving outcome to right
        TRUE ~ "1_exposure"), 
      # Source Omics Layer: lc1-lc4 (for omics layers), outcome, or other 
      source_layer = case_when(
        str_detect(source, "Layer 1") ~ "lc1", 
        str_detect(source, "Layer 2") ~ "lc2", 
        str_detect(source, "Layer 3") ~ "lc3", 
    
        # source == "Outcome" ~ str_sub(target, start = -3, end = -2),  # removed when moving outcome to right
        TRUE ~ "exposure"), 
      # Source group_ for color (one of: exposure, : lc1-lc4 (for omics layers), outcome, or other 
      color_group_node = if_else(source == "outcomee", 
                                 "outcome", 
                                 source_layer)) %>%
    group_by(source_group) %>%
    arrange(source_layer, .by_group = TRUE) %>%
    ungroup() %>%
    dplyr::select(source, target, value, color_group, color_group_node)
  
  
  # Create Nodes ----
  nodes <- links %>%
    dplyr::select(source, color_group_node) %>%
    mutate(rownum = row_number()) %>%
    rename(name = source, 
           color_group = color_group_node) %>%
    # Add outcome (only if outcome is on right side)
    bind_rows(data.frame(name = "outcome", color_group = "outcome")) %>%
    group_by(name) %>%
    slice_head() %>%
    ungroup() %>%
    arrange(rownum) %>%
    dplyr::select(-rownum) %>%
    # Add feature names
    bind_rows(names_features) %>% 
    mutate(id = row_number()-1) %>%
    left_join(sankey_colors, by = c( "color_group"= "domain"))
  
  # Join links and nodes for color names -----
  links <- links %>%
    left_join(nodes %>% 
                dplyr::select(id, name), 
              by = c("source" = "name")) %>%
    rename(source_id = id) %>% 
    dplyr::select(source_id, everything()) %>%
    left_join(nodes %>% 
                dplyr::select(id, name), 
              by = c("target" = "name")) %>%
    rename(target_id = id) %>% 
    dplyr::select(source_id,source, target_id,everything()) 
  
  
  # Manually change colors ----
  links <- links  %>%
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4", 
        # Outcome
        str_detect(target, "outcome") &  color_group == TRUE  ~  "red",
        # Methylation 
        str_detect(source, "Layer 1") &  color_group == TRUE  ~  "#bf9000",
        str_detect(source, "Layer 1") &  color_group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(source, "Layer 3") &  color_group == TRUE  ~  "#38761d",
        str_detect(source, "Layer 3") &  color_group == FALSE ~  "#b6d7a8",
        # mirna
        str_detect(source, "Layer 2") &  color_group == TRUE  ~  "#a64d79",
        str_detect(source, "Layer 2") &  color_group == FALSE ~  "#ead1dc",
        
        
        links$color_group == FALSE ~ "#d9d2e9", # Negative association
        links$color_group == TRUE ~  "red"))
  
  ## change node names
  nodes <- nodes %>%
    mutate(name = case_when(name == "Cluster 1 (Layer 1)" ~ "<b>Serum metabolites\nProfile 0</b>",
                            name == "Cluster 2 (Layer 1)" ~ "<b>Serum metabolites\nProfile 1</b>",
                            name == "Cluster 1 (Layer 2)" ~ "<b>Urine metabolites\nProfile 0</b>",
                            name == "Cluster 2 (Layer 2)" ~ "<b>Urine metabolites\nProfile 1</b>",
                            name == "Cluster 1 (Layer 3)" ~ "<b>Protein\nProfile 0</b>",
                            name == "Cluster 2 (Layer 3)" ~ "<b>Protein\nProfile 1</b>",
                            name == "outcome" ~ lucidus_fit$var.names$Ynames,
                            TRUE ~ name),
           x = case_when(
             str_detect(name, "Serum metabolites") |
               str_detect(name, "Urine metabolites") | 
               str_detect(name, "Protein") ~ 1/3, 
             str_detect(name, "cg")|
               str_detect(name, "tc")|
               str_detect(name, "miR")| 
               str_detect(name, "outcome") ~ 2/3))
  
  
 
  (fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    domain = list(
      x =  c(0,0.8),
      y =  c(0,1)),
    # arrangement = "snap",
    node = list(
      label = nodes$name,
      color = nodes$range,
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      ),
      x = nodes$x
    ),
    
    link = list(
      source = links$source_id,
      target = links$target_id,
      value =  links$value+.00000000000000000000001,
      # label = links$source,
      color = links$link_color
    )
  )
  )
  
  fig <- fig %>% layout(
    font = list(
      size = text_size
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )
  
  fig
}

# ----------- reorder LUCID in Parallel model by specifying reference cluster ------------
# note: only works for K = 2 in each omic layer
# reference = c(1,1,2)
# lucidus_fit <- fit_reordered
reorder_lucid_parallel <- function(lucidus_fit,
                                   reference = NULL) {
  if(is.null(reference)) {
    warning("no reference specified, return the original model")
    return(lucidus_fit)
  }
  
  n_omic <- length(reference)
  
  # reorder beta
  GtoX <- lucidus_fit$res_Beta$Beta
  lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
    (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
    # if reference = 2, flip the reference and negate the estimates
  })
  # reorder mu
  XtoZ <- lucidus_fit$res_Mu
  lucidus_fit$res_Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  XtoY <- lucidus_fit$res_Gamma$Gamma$mu
  XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
  XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
  lucidus_fit$res_Gamma$Gamma$mu <- XtoY
  lucidus_fit$res_Gamma$fit$coefficients <- XtoY
  
  # return the object using the new reference
  return(lucidus_fit)
}


## ---- set_color_pal ----

# Set Color Palettes 
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")


# Color pallet for sankey diagrams
sankey_colors_parallel <- matrix(c("exposure", col_pal[6],
                          "lc1",      "#b3d8ff",
                          "lc2",      "#b3d8ff",
                          "lc3",      "#b3d8ff",
                          "lc4",      "#b3d8ff",
                          "layer1",   col_pal[2],
                          "layer2",   col_pal[3],
                          "layer3",   col_pal[1],
                          "layer4",   col_pal[4],
                          "outcome",  col_pal[8],
                      
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red", 
                          "neg_clus_to_out", "#e4e5f2"), 
                        byrow = TRUE, nrow = 14)

# Change to dataframe
colnames(sankey_colors_parallel) <- c("domain", "range")
sankey_colors_parallel <- as_tibble(sankey_colors_parallel)



